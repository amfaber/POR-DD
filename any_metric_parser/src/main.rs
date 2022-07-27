// use std::sync::mpsc::Sender;
// use std::alloc::Global;
use std::{path, fs, io, slice, str};
use flate2::read::GzDecoder;
use std::io::BufRead;
use std::io::prelude::*;
use regex::Regex;
use std::collections::{HashMap};

// use std::time;
// use std::thread;
use threadpool::ThreadPool;
use lazy_static::lazy_static;
use std::sync::{mpsc, Mutex, Arc};

lazy_static!{
    static ref RE: Regex = Regex::new(r"(?:AID\d+)(?P<activity>_(?:in)?active)?(?:_(?P<idx>\d+))?-lig_(?P<protein>\w{4})-rec_docked_rescore$").unwrap();
}

struct MolEntry{
    index : usize,
    name: String,
    file: String,
    vina: f32,
    cnnaffinity: f32,
    cnnscore: f32,
    vsscore: f32,
    activity: String,
    protein: String,
    index_from_file: Option<usize>,
}

// Header: {0: "Rank", 1: "Title", 2: "Vina", 3: "Target", 4: "File", 5: "newdefault_CNNaffinity", 6: "newdefault_CNNscore"}



impl MolEntry{
    fn default() -> MolEntry{
        MolEntry{
            index: 0,
            name: String::new(),
            file: String::new(),
            vina: f32::INFINITY,
            cnnaffinity: f32::NEG_INFINITY,
            cnnscore: f32::NEG_INFINITY,
            vsscore: f32::NEG_INFINITY,
            activity: String::new(),
            protein: String::new(),
            index_from_file: None,
        }
    }

    fn new(line: Vec<&str>) -> MolEntry{
        let index = line[0].parse::<usize>().unwrap()-1;
        let name = line[1].to_string();
        let vina = line[2].parse::<f32>().unwrap();
        let file = line[4].to_string();
        let cnnaffinity = line[5].parse::<f32>().unwrap();
        let cnnscore = line[6].parse::<f32>().unwrap();
        let vsscore = cnnaffinity * cnnscore;

        let properties = match RE.captures(&file){
            Some(res) => res,
            None => {
                panic!("{}", file);
            }
        };
        let activity = match properties.name("activity"){
            Some(res) => res.as_str()[1..].to_string(),
            None => "active".to_string(),
        };
        let protein = properties.name("protein").unwrap().as_str().to_string();

        let index_from_file = match properties.name("idx"){
            Some(res) => Some(res.as_str().parse::<usize>().unwrap()),
            None => None,
        };


        MolEntry{
            index,
            name,
            file,
            vina,
            cnnaffinity,
            cnnscore,
            vsscore,
            activity,
            protein,
            index_from_file,
        }
    }

    fn clone(&self) -> MolEntry
    {
        MolEntry{
            index: self.index,
            name: self.name.clone(),
            file: self.file.clone(),
            vina: self.vina,
            cnnaffinity: self.cnnaffinity,
            cnnscore: self.cnnscore,
            vsscore: self.vsscore,
            activity: self.activity.clone(),
            protein: self.protein.clone(),
            index_from_file: self.index_from_file,
        }
    }

    fn update_to_best(&mut self, other: &Self, compare_type: &str) -> (){

        let new_is_better = match compare_type.to_lowercase().as_str(){
            "vina" => {
                other.vina < self.vina
            },
            "cnnaffinity" => {
                other.cnnaffinity > self.cnnaffinity
            },
            "cnnscore" => {
                other.cnnscore > self.cnnscore
            },
            "vsscore" => {
                other.vsscore > self.vsscore
            },
            _ => {
                panic!("Unsupported comparison type: {}", compare_type);
            }
        };
        if new_is_better{
            *self = other.clone();
        }
    }
}

fn load_mols(input_dir: &path::Path, sdf_name: &str, sdf_contents: &mut String, mols: &mut Vec<&str>) -> Result<(), io::Error>{
    mols.clear();
    sdf_contents.clear();
    let sdf_file_path = input_dir.join(sdf_name.replace("_rescore", ".sdf.gz"));
    let sdf_file = fs::File::open(sdf_file_path)?;
    let mut sdf_decoder = GzDecoder::new(sdf_file);
    // sdf_decoder.read_to_string(sdf_contents)?;
    sdf_decoder.read_to_string(sdf_contents)?;

    let iter = sdf_contents.split("$$$$\n");

    unsafe{
        // let mut i = 0;
        for mol in iter {
            // i += 1;
            if mol.len() == 0 {
                continue;
                // println!("{}", i);
            }
            let slc = str::from_utf8_unchecked(slice::from_raw_parts(mol.as_ptr(), mol.len()));
            mols.push(slc);
        }
    }
    
    Ok(())
}

macro_rules! str_append {
    ($first:expr, $($x:expr),*) => {
        {
            $($first.push_str($x);)*
        }
    };
}

fn open_output_file(
    output_dir: &path::Path,
    activity: &str,
    ) -> Result<fs::File, io::Error>{
    // let mut output_file = None;
    let writer_path = output_dir.
                                join(format!("{}.sdf", activity));
    fs::create_dir_all(writer_path.parent().unwrap()).unwrap();
    let output_file = fs::File::create(writer_path).unwrap();

    Ok(output_file)
}

fn handle_end_of_chunk(
    best_mols: & Vec<MolEntry>,
    target: &str, 
    input_dir: &path::PathBuf, 
    header: &Vec<String>, 
    pool: &Arc<Mutex<ThreadPool>>,
    output_file: &mut Option<fs::File>,
    ) -> Result<(), io::Error>{
    let mut from_each_file = HashMap::<String, Vec<MolEntry>>::new();

    for mol in best_mols.iter(){
        from_each_file.entry(mol.file.clone()).or_insert(Vec::new()).push(mol.clone());
    }

    // let mut iterators = HashMap::<String, _>::new();
    let (tx, rx) = mpsc::channel();

    for (file, molentries) in from_each_file{
        let target_name = target.to_string(); 
        let input_dir = input_dir.clone();
        let header = header.clone();
        let tx = tx.clone();
        // println!("About to load {} in target {}", file, target_name);
        let pool = pool.lock().unwrap();
        (*pool).execute(move || {
            println!("Loading {} in target {}", file, target_name);
            let mut sdf_contents = String::new();
            let mut mols = Vec::new();
            load_mols(&input_dir, &file, &mut sdf_contents, &mut mols).unwrap();
            let mut this_files_structures = Vec::new();
            for molentry in molentries{
                this_files_structures.push(mols[molentry.index].to_string());
                let mol = this_files_structures.iter_mut().next_back().unwrap();
                str_append!(mol, "> <", &header[5], ">\n",
                &molentry.cnnaffinity.to_string(), 
                "\n\n> <", &header[6], ">\n",
                &molentry.cnnscore.to_string(),
                "\n\n> <CNNaffinity_times_score>\n",
                &molentry.vsscore.to_string(),
                "\n\n> <protein>\n",
                &molentry.protein,
                "\n\n$$$$\n");
            }
            tx.send((file, this_files_structures.into_iter())).unwrap();
        });
        drop(pool);
    }
    drop(tx);
    let mut best_structures = HashMap::<String, _>::new();
    for received in rx{
        let (file, structures) = received;
        best_structures.insert(file, structures);
    }

    let mut chunk_to_write = String::new();


    for mol in best_mols.iter(){
        chunk_to_write.push_str(&(best_structures.get_mut(&mol.file).unwrap().next().unwrap()));
    }

    match output_file {
        Some(ref mut f) => {
            f.write_all(chunk_to_write.as_bytes())?;
        },
        None => {
            panic!("No output file opened");
        }
    }
    Ok(())
}

fn read_block(rest_of_lines: &mut io::Lines<io::BufReader<GzDecoder<fs::File>>>,
    leftover: Option<MolEntry>,

    ) -> Result<(
        HashMap<String, usize>,
        Vec<MolEntry>,
        Option<MolEntry>,
        Option<String>,
    ), io::Error> {
    let mut name_to_ind = HashMap::<String, usize>::new();
    let mut best_mols  = Vec::<MolEntry>::new();

    let mut prev_mol = match leftover{
        Some(mol) => {
            best_mols.push(mol.clone());
            name_to_ind.insert(mol.file.clone(), 0);
            mol
        },
        None => MolEntry::default(),
    };
    let mut leftover = None;
    let mut need_new_file = None;
    for line in rest_of_lines{
        let line: Vec<&str> = line.as_ref().unwrap().split_whitespace().collect();

        let new_mol = MolEntry::new(line);

        let at_the_end_of_chunk = new_mol.index_from_file != prev_mol.index_from_file;
        if new_mol.activity != prev_mol.activity {
           need_new_file = Some(new_mol.activity.clone());
        };

        if at_the_end_of_chunk{
            leftover = Some(new_mol);
            // prev_mol = new_mol;
            break;
        }
        
        match name_to_ind.get(&new_mol.name){
            Some(index) => {(&mut best_mols[*index]).update_to_best(&new_mol, "vsscore");}
            None => {
                name_to_ind.insert(new_mol.name.clone(), best_mols.len());
                best_mols.push(new_mol.clone());
            }
        };
        
        prev_mol = new_mol;
    }
    Ok((name_to_ind, best_mols, leftover, need_new_file))
}

fn parse_target(target: &str, input_dir: path::PathBuf, output_dir: path::PathBuf, pool: Arc<Mutex<ThreadPool>>) -> Result<(), io::Error> {
    println!("Starting parse of {}", target);
    let input_dir = input_dir.join(target);
    let output_dir = output_dir.join(target);
    
    // let mut already_opened_files = HashSet::new();

    let summary_file = fs::File::open(input_dir.join("newdefault.summary.gz"))?;
    let summary_decoder = GzDecoder::new(summary_file);
    let mut summary_bufreader = io::BufReader::with_capacity(8*1000, summary_decoder);
    
    let mut header = String::new();
    summary_bufreader.read_line(&mut header)?;
    let header = header.split_whitespace().map(|s| s.to_string()).collect::<Vec<String>>();

    let mut prev_mol = MolEntry::default();
    let mut all_lines = summary_bufreader.lines();

    let leftover = None;
    let (name_to_ind, best_mols, leftover, need_new_file) = read_block(&mut all_lines, leftover)?;
    let mut prev_leftover = leftover;
    let mut output_file = open_output_file(&output_dir, &need_new_file.unwrap())?;

    loop{
        let (name_to_ind, best_mols, leftover, need_new_file) = read_block(&mut all_lines, prev_leftover)?;
        if leftover.is_none(){
            break;
        }
        else {
            prev_leftover = leftover;
        }
    }
    // for line in all_lines{
    //     let line: Vec<&str> = line.as_ref().unwrap().split_whitespace().collect();

    //     let new_mol = MolEntry::new(line);
        
    //     if output_file.is_none() {
    //         output_file = open_output_file(&output_dir, &new_mol.activity)?;
    //     }

    //     let at_the_end_of_chunk = new_mol.index_from_file != prev_mol.index_from_file;

    //     if at_the_end_of_chunk {
    //         handle_end_of_chunk(&mut best_mols, &mut name_to_ind, target, &input_dir, &header, &pool, &mut output_file)?;
    //     }

    //     if new_mol.activity != prev_mol.activity {
    //         output_file = open_output_file(&output_dir, &new_mol.activity)?;
    //     }

        
    //     match name_to_ind.get(&new_mol.name){
    //         Some(index) => {(&mut best_mols[*index]).update_to_best(&new_mol, "vsscore");}
    //         None => {
    //             name_to_ind.insert(new_mol.name.clone(), best_mols.len());
    //             best_mols.push(new_mol.clone());
    //         }
    //     };
        
    //     prev_mol = new_mol;
        
    // }
    // handle_end_of_chunk(&mut best_mols, &mut name_to_ind, target, &input_dir, &header, &pool, &mut output_file)?;
    Ok(())
}

fn parse_all_targets_in_dir(input_dir: path::PathBuf, output_dir: path::PathBuf, threads_inner: usize, threads_outer: usize) -> Result<(), io::Error> {
    let pool_outer = ThreadPool::new(threads_outer);
    let pool_inner = Arc::new(Mutex::new(ThreadPool::new(threads_inner)));
    let targets = fs::read_dir(&input_dir)?;
    // let (tx, rx) = mpsc::channel();
    for target in targets{
        let target = target?;
        if target.file_type()?.is_dir(){
            let target = target.file_name().to_str().unwrap().to_string();
            let input_dir = input_dir.clone();
            let output_dir = output_dir.clone();
            // let tx = tx.clone();
            let pool_inner = pool_inner.clone();
            pool_outer.execute(move || {
                parse_target(&target, input_dir, output_dir, pool_inner).unwrap();
            });
            // parse_target("ADRB2", input_dir, output_dir).unwrap();
            // break;

        }
    }
    pool_outer.join();
    Ok(())
}
fn main() -> Result<(), io::Error> {

    let input_dir = path::PathBuf::from("/home/qzj517/POR-DD/data/gnina_LIT-PCBA_VS_data/lit-pcba");
    let output_dir = path::PathBuf::from("/home/qzj517/POR-DD/data/gnina_LIT-PCBA_VS_data/best_scoring_vsscore");
    parse_all_targets_in_dir(input_dir, output_dir, 7, 2)?;
    Ok(())
}
