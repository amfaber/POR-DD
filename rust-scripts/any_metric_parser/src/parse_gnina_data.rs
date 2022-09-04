use glob::glob;
use std::{path, fs, io, slice, str};
use flate2::read::GzDecoder;
use std::io::BufRead;
use std::io::prelude::*;
use regex::Regex;
use std::collections::{HashMap};


use std::time;
use threadpool::ThreadPool;
use std::sync::{mpsc, Mutex, Arc};


pub fn parse_filename(filename: &str) -> (String, Option<i32>, String){
    let mut iter = filename.split("-");
    let mut iter2 = iter.next().unwrap().split("_");
    iter2.next();
    let activity;
    let index_from_file;
    match iter2.next(){
        Some(res) => {
            match res.parse::<i32>() {
                Ok(index) => {index_from_file = Some(index);
                    activity = "active".to_string();
                },
                Err(_) => {
                    activity = res.to_string();
                    index_from_file = iter2.next().map(|ele| ele.parse::<i32>().unwrap());
                    }
                }
        },
        None => {
            activity = "active".to_string();
            index_from_file = None;
        }

    }
    // drop(iter2);
    let protein = iter.next().unwrap().split("_").nth(1).unwrap().to_string();
    (activity, index_from_file, protein)
}

pub struct MolEntry{
    pub index : usize,
    pub name: String,
    pub file: String,
    pub vina: f32,
    pub cnnaffinity: f32,
    pub cnnscore: f32,
    pub vsscore: f32,
    pub activity: String,
    pub protein: String,
    pub index_from_file: Option<i32>,
}

// Header: {0: "Rank", 1: "Title", 2: "Vina", 3: "Target", 4: "File", 5: "newdefault_CNNaffinity", 6: "newdefault_CNNscore"}



impl MolEntry{
    pub fn default() -> MolEntry{
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
            index_from_file: Some(-1),
        }
    }

    pub fn new(line: Vec<&str>) -> MolEntry{
        let index = line[0].parse::<usize>().unwrap()-1;
        let name = line[1].to_string();
        let mut vina = line[2].parse::<f32>().unwrap();
        let file = line[4].to_string();
        let cnnaffinity: f32;
        let cnnscore: f32;
        let vsscore: f32;
        match line.get(5){
                    Some(s) => {
                        cnnaffinity = s.parse::<f32>().unwrap();
                        cnnscore = line[6].parse::<f32>().unwrap();
                        vsscore = cnnaffinity * cnnscore;
                    },
                    None => {
                        cnnaffinity = f32::NEG_INFINITY;
                        cnnscore = f32::NEG_INFINITY;
                        vsscore = f32::NEG_INFINITY;
                        vina = f32::NEG_INFINITY;
                    },
        };
        

        // let cnnscore = line[6].parse::<f32>().unwrap();
        // let vsscore = cnnaffinity * cnnscore;

        let (activity, index_from_file, protein) = parse_filename(&file);

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

    pub fn clone(&self) -> MolEntry
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

    pub fn update_to_best(&mut self, other: &Self, compare_type: &str, equibind_best: &Option<HashMap<String, String>>) -> (){

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

        let accept_change = match equibind_best.as_ref(){
            Some(hashmap) => {
                match hashmap.get(&self.name){
                    Some(equibind_protein) => {
                        equibind_protein == &other.protein && (equibind_protein != &self.protein || new_is_better)
                    },
                    None => {
                        new_is_better
                    }
                }
            }
            None => {
                new_is_better
            }
        };

        if accept_change{
            *self = other.clone();
        }
    }
}

pub fn load_mols<'a>(input_dir: &path::Path, sdf_name: &str) -> Result<(String, Vec<&'a str>), io::Error>{
    // mols.clear();
    // sdf_contents.clear();
    let mut mols = Vec::new();
    let mut sdf_contents = String::new();
    let sdf_file_path = input_dir.join(sdf_name.replace("_rescore", ".sdf.gz"));
    let sdf_file;
    match fs::File::open(&sdf_file_path){
        Ok(file) => sdf_file = file,
        Err(e) => panic!("Paniced with error {} on filename {}", e, sdf_file_path.to_str().unwrap())
    }
    let mut sdf_decoder = GzDecoder::new(sdf_file);
    // sdf_decoder.read_to_string(sdf_contents)?;
    sdf_decoder.read_to_string(&mut sdf_contents).unwrap();

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

    Ok((sdf_contents, mols))
}

macro_rules! str_append {
    ($first:expr, $($x:expr),*) => {
        {
            $($first.push_str($x);)*
        }
    };
}

pub fn open_output_file(
    output_dir: &path::Path,
    activity: &str,
    summary_name: &str,
    ) -> Result<fs::File, io::Error>{
    // let mut output_file = None;
    let writer_path = output_dir.
                                join(format!("{}_{}.sdf", summary_name, activity));
    fs::create_dir_all(writer_path.parent().unwrap()).unwrap();
    let output_file = fs::File::create(writer_path).unwrap();

    Ok(output_file)
}

pub fn start_decompressor_threads(
    best_mols: & Vec<MolEntry>,
    target: &str, 
    input_dir: &path::PathBuf, 
    header: &Vec<String>, 
    pool: &Arc<Mutex<ThreadPool>>,
    // output_file: &mut fs::File,
    ) -> mpsc::Receiver<(String, std::vec::IntoIter<String>)>{
    let mut from_each_file = HashMap::<String, Vec<MolEntry>>::new();

    for mol in best_mols.iter(){
        from_each_file.entry(mol.file.clone()).or_insert(Vec::new()).push(mol.clone());
    }

    // let mut iterators = HashMap::<String, _>::new();
    let (tx, rx) = mpsc::channel();

    for (file, molentries) in from_each_file{
        let _target_name = target.to_string(); 
        let input_dir = input_dir.clone();
        let header = header.clone();
        let tx = tx.clone();
        let pool = pool.lock().unwrap();
        (*pool).execute(move || {
            // println!("Loading {} in target {}", file, target_name);

            let (_sdf_contents, mols) = load_mols(&input_dir, &file).unwrap();
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
            // println!("Done {} in target {}", file, target_name);
            tx.send((file, this_files_structures.into_iter())).unwrap();
        });
        // drop(pool);
    }
    drop(tx);
    rx
}

pub fn gather_work_and_write(rx: mpsc::Receiver<(String, std::vec::IntoIter<String>)>, best_mols: &Vec<MolEntry>, output_file: &mut fs::File, target: &str) -> Result<(), io::Error>{
    // println!("Writing block {} in target {}", best_mols[0].index_from_file.unwrap_or(-1), target);
    let mut best_structures = HashMap::<String, _>::new();
    // let mut waiting = false;
    // match rx.try_recv(){
    //     Ok((file, structures)) => {
    //         best_structures.insert(file, structures);
    //     }
    //     Err(mpsc::TryRecvError::Empty) => {
    //         println!("{} started waiting for work", target);
    //         waiting = true;
    //     }
    //     Err(mpsc::TryRecvError::Disconnected) => {
    //         panic!("{} disconnected", target);
    //     }
    // };
    let now = time::Instant::now();
    for received in rx{
        // if waiting{
        //     println!("{} stopped waiting for work", target);
        //     waiting = false;
        // }
        let (file, structures) = received;
        best_structures.insert(file, structures);
    }
    println!("{} waited for {}s", target, (now.elapsed().as_millis() as f32)/1000.);

    let mut chunk_to_write = String::new();


    for mol in best_mols.iter(){
        chunk_to_write.push_str(&(best_structures.get_mut(&mol.file).unwrap().next().unwrap()));
    }

    output_file.write_all(chunk_to_write.as_bytes())?;
    Ok(())
}

pub fn read_block(rest_of_lines: &mut io::Lines<io::BufReader<GzDecoder<fs::File>>>,
    leftover: Option<MolEntry>,
    _target: &str,
    equibind_best: &Option<HashMap<String, String>>,
    compare_type: &str
    ) -> Result<(
        // HashMap<String, usize>,
        Vec<MolEntry>,
        Option<MolEntry>,
        Option<String>,
    ), io::Error> {
    // let block = match leftover.as_ref(){
    //     Some(mol) => mol.index_from_file.unwrap_or(-1),
    //     None => 0,
    // };
    // println!("Reading block {} in {}", block, target);
    
    // println!("Reading block {}", );
    let mut name_to_ind = HashMap::<String, usize>::new();
    let mut best_mols  = Vec::<MolEntry>::new();
    let mut to_read_first_line = leftover.is_none();
    
    let mut prev_mol = match leftover{
        Some(mol) => {
            name_to_ind.insert(mol.name.clone(), 0);
            best_mols.push(mol.clone());
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
            if to_read_first_line{
                to_read_first_line = false;
            }
            else{
                leftover = Some(new_mol);
                // prev_mol = new_mol;
                break;
            }
        }
        
        match name_to_ind.get(&new_mol.name){
            Some(index) => {(&mut best_mols[*index]).update_to_best(&new_mol, compare_type, equibind_best);}
            None => {
                name_to_ind.insert(new_mol.name.clone(), best_mols.len());
                best_mols.push(new_mol.clone());
            }
        };
        
        prev_mol = new_mol;
    }
    // println!("Finished block {} in {}", block, target);
    Ok((best_mols, leftover, need_new_file))
}

pub fn load_best_equibind(input_dir: &path::Path, load: bool) -> Option<HashMap<String, String>>{
    if !load{
        return None;
    }
    let mut best_equibind = HashMap::<String, String>::new();
    let mut best_equibind_file = fs::File::open(input_dir.join("best_from_equibind.txt")).unwrap();
    let mut contents = String::new();
    best_equibind_file.read_to_string(&mut contents).unwrap();
    for line in contents.lines(){
        let line: Vec<&str> = line.split_whitespace().collect();
        best_equibind.insert(line[0].to_string(), line[1].to_string());
    }
    Some(best_equibind)
}


pub fn parse_target(target: &str,
    input_dir: path::PathBuf,
    output_dir: path::PathBuf,
    pool: Arc<Mutex<ThreadPool>>,
    summary_name: &str,
    compare_type: &str,
    same_as_equibind: bool,
    ) -> Result<(), io::Error> {
    println!("Starting parse of {}", target);
    let input_dir = input_dir.join(target);
    let output_dir = output_dir.join(target);
    
    let best_equibind = load_best_equibind(&input_dir, same_as_equibind);

    let summary_file = fs::File::open(input_dir.join(summary_name))?;
    // let summary_file = fs::File::open(input_dir.join("IDH1_actives.summary.gz"))?;
    let summary_decoder = GzDecoder::new(summary_file);
    let mut summary_bufreader = io::BufReader::with_capacity(8*1000, summary_decoder);
    
    let mut header = String::new();
    summary_bufreader.read_line(&mut header)?;
    let header = header.split_whitespace().map(|s| s.to_string()).collect::<Vec<String>>();

    let mut all_lines = summary_bufreader.lines();

    let leftover = None;
    let (best_mols, leftover, need_new_file) = read_block(&mut all_lines, leftover, &target, &best_equibind, compare_type)?;


    let mut prev_leftover = leftover;
    let mut prev_best_mols = best_mols;
    let mut output_file = open_output_file(&output_dir, &need_new_file.as_ref().unwrap(), summary_name)?;
    let mut prev_need_new_file: Option<String> = None;
    loop{
        let rx = start_decompressor_threads(&prev_best_mols, target, &input_dir, &header, &pool);
        let (best_mols, leftover, need_new_file) = read_block(&mut all_lines, prev_leftover, &target, &best_equibind, compare_type)?;
        

        gather_work_and_write(rx, &prev_best_mols, &mut output_file, &target)?;
        if prev_need_new_file.is_some(){
            let activity = prev_need_new_file.as_ref().unwrap();
            output_file = open_output_file(&output_dir, &activity, summary_name)?;
        }
        prev_leftover = leftover;
        prev_best_mols = best_mols;
        prev_need_new_file = need_new_file;
        if prev_leftover.is_none(){
            break;
        }
    }
    if prev_need_new_file.is_some(){
        let activity = prev_need_new_file.as_ref().unwrap();
        output_file = open_output_file(&output_dir, &activity, summary_name)?;
    }
    let rx = start_decompressor_threads(&prev_best_mols, target, &input_dir, &header, &pool);
    gather_work_and_write(rx, &prev_best_mols, &mut output_file, &target)?;

    Ok(())
}

pub fn gather_files(target_dir: &path::PathBuf, re: &Regex) -> Result<(), io::Error>{
    // let re = Regex::new(r".*newdefault(?:_)?(\d+)?.summary.gz_(.*.sdf)$").unwrap();
    // if let Ok(results) = fs::read_dir(&output_dir){
    //     for target in results{
    //         if let Ok(target) = target{
    //             if let Ok(filetype) = target.file_type(){
    //                 if !filetype.is_dir(){
    //                     continue;
    //                 }
    //             }
    if let Ok(files) = fs::read_dir(target_dir){
        let file_names = files.
            filter_map(|file| file.ok()).
            filter_map(|file| file.path().to_str().map(|name| name.to_string())).collect::<Vec<_>>();
        let (mut inactives, mut actives) : (Vec<_>, Vec<_>) = file_names.
            iter().
            filter_map(|name| re.captures(name)).
            partition(|cap| cap.get(2).unwrap().as_str() == "inactive.sdf");
        inactives.sort_by_key(|cap| cap.get(1).map(|m| m.as_str()).unwrap_or("0").parse::<i32>().unwrap());
        actives.sort_by_key(|cap| cap.get(1).map(|m| m.as_str()).unwrap_or("0").parse::<i32>().unwrap());
        for types in vec![inactives, actives]{
            let mut t = types.into_iter();
            // target.path().join(t.next().unwrap().get(0).unwrap().as_str());
            let first_match;
            match t.next(){
                Some(m) => first_match = m,
                None => continue
            }
            let (current_name, new_name) = (&first_match[0], target_dir.join(&first_match[2]));
            // let write_to = target.path().join(&t.next().unwrap()[0]);
            fs::rename(current_name, &new_name)?;
            let mut output = fs::OpenOptions::new().append(true).open(&new_name)?;
            for cap in t{
                // println!("{}", &cap[0]);
                let mut input = fs::File::open(&cap[0])?;
                io::copy(&mut input, &mut output)?;
                fs::remove_file(&cap[0])?;
            }
        }
        // println!("{:?}", &actives);
    }
            // }
        // }
    // }
    return Ok(())
}

pub fn parse_all_targets_in_dir(input_dir: path::PathBuf,
    output_dir: path::PathBuf,
    inner_threads: usize,
    outer_threads: usize,
    targets: &Option<Vec<String>>,
    compare_type: &str,
    same_as_equibind: bool
    ) -> Result<(), io::Error> {
    let only_parse_targets = targets.as_ref();
    let pool_outer = ThreadPool::new(outer_threads);
    let pool_inner = Arc::new(Mutex::new(ThreadPool::new(inner_threads)));
    let targets = fs::read_dir(&input_dir)?;

    for target in targets{
        let target = target?;
        if target.file_type()?.is_dir(){
            let curdir = target.path();
            for summary_file in glob(curdir.join("newdefault*.summary.gz").to_str().unwrap()).expect("Glob failed"){
                let summary_file = summary_file.unwrap().file_name().unwrap().to_str().unwrap().to_string();
                let target_name = curdir.file_name().unwrap().to_str().unwrap().to_string();
                match only_parse_targets{
                    Some(list) => {
                        if !list.contains(&target_name){
                            continue;
                        }
                    }
                    None => {}
                }
                let input_dir = input_dir.clone();
                let output_dir = output_dir.clone();
                let compare_type = compare_type.to_string();
                // println!("{}", summary_file.to_str().unwrap())
                let pool_inner = pool_inner.clone();
                pool_outer.execute(move || {
                    parse_target(&target_name, input_dir, output_dir, pool_inner, &summary_file, &compare_type, same_as_equibind).unwrap();
                });

            }
            // break;
        }
    }
    pool_outer.join();

    let re = Regex::new(r".*newdefault(?:_)?(\d+)?.summary.gz_(.*.sdf)$").unwrap();
    if let Ok(results) = fs::read_dir(&output_dir){
        for target in results{
            if let Ok(target) = target{
                if let Ok(filetype) = target.file_type(){
                    if !filetype.is_dir(){
                        continue;
                    }
                }
                let target_dir = target.path();
                let re = re.clone();
                pool_outer.execute(move || {
                    gather_files(&target_dir, &re).unwrap();
                });
            }
        }
    }
    pool_outer.join();
    Ok(())
}
