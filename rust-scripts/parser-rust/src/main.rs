use std::{path, fs, io, slice, str};
use flate2::read::GzDecoder;
use std::io::BufRead;
use std::io::prelude::*;
use regex::Regex;
use std::collections::HashSet;
// use std::time;
// use std::thread;
use threadpool::ThreadPool;

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
    file_name: &str,
    re: &Regex,
    already_opened_files: &mut HashSet<String>,
    ) -> Result<Option<fs::File>, io::Error>{
    let properties = match re.captures(file_name){
        Some(res) => res,
        None => {
            panic!("{}", file_name);
        }
    };
    let activity = match properties.name("activity"){
        Some(res) => &res.as_str()[1..],
        None => "active",
    };
    let writer_path = output_dir.
                                    join(activity).
                                    join(format!("{}_protein.sdf", &properties["protein"]));
    
    let mut output_file = None;

    let is_new = already_opened_files.insert(writer_path.to_str().unwrap().to_string());
    
    if is_new {
        fs::create_dir_all(writer_path.parent().unwrap()).unwrap();
        output_file = Some(fs::File::create(writer_path).unwrap());
    }
    else {
        output_file = Some(fs::OpenOptions::new().append(true).open(writer_path).unwrap());
    }
    Ok(output_file)
}

fn parse_target(target: &str, input_dir: path::PathBuf, output_dir: path::PathBuf) -> Result<(), io::Error> {
    
    let re = Regex::new(r"(?P<prefix>AID\d+)(?P<activity>_(?:in)?active)?(?:_(?P<idx>\d+))?-lig_(?P<protein>\w{4})-rec_docked_rescore$").unwrap();

    let input_dir = input_dir.join(target);
    let output_dir = output_dir.join(target);
    
    let summary_file = fs::File::open(input_dir.join("newdefault.summary.gz"))?;
    let summary_decoder = GzDecoder::new(summary_file);
    let mut summary_bufreader = io::BufReader::with_capacity(8*1000, summary_decoder);
    
    let mut header = String::new();
    summary_bufreader.read_line(&mut header)?;
    let header = header.split_whitespace().collect::<Vec<&str>>();

    let mut sdf_name = String::new();

    let mut sdf_contents = String::new();
    
    let mut mols = Vec::new();

    let mut prev_index = 1;

    let mut prev_name = String::new();
    
    // let mut i = 0;

    let mut output_file = None;

    let mut already_opened_files = HashSet::<String>::new();

    let mut skip_due_to_duplication = false;

    for line in summary_bufreader.lines(){
        let line: Vec<&str> = line.as_ref().unwrap().split_whitespace().collect();

        let index = line[0].parse::<usize>().unwrap();

        // Header: {0: "Rank", 1: "Title", 2: "Vina", 3: "Target", 4: "File", 5: "newdefault_CNNaffinity", 6: "newdefault_CNNscore"}
        
        let at_new_file = index <= prev_index;
        let file_actually_different = sdf_name != line[4];

        if at_new_file {
            // println!("Loading new file {}", line[4]);
            // let now = time::Instant::now();
            skip_due_to_duplication = !file_actually_different;
            if skip_due_to_duplication {
                prev_index = index;
                continue;
            }
            sdf_name = line[4].to_string();
            load_mols(&input_dir, &sdf_name, &mut sdf_contents, &mut mols)?;
            output_file = open_output_file(&output_dir, &sdf_name, &re, &mut already_opened_files)?;

            // println!("Loaded {} in {:?}", sdf_name, now.elapsed());
        }

        if skip_due_to_duplication{
            prev_index = index;
            continue;
        }

        if prev_name != line[1]{
            prev_name = line[1].to_string();
            let summary_aff = line[2].parse::<f32>().unwrap();
            let mut out_mol = mols[index-1].to_string();
            let sdf_aff = out_mol.lines().nth_back(1).unwrap().parse::<f32>().unwrap();
            let diff = (summary_aff - sdf_aff).abs();
            assert!(diff < 0.001);
            str_append!(out_mol, "> <", header[5], ">\n", line[5], "\n\n> <", header[6], ">\n", line[6], "\n\n$$$$\n");
            if let Some(writer) = output_file.as_mut(){
                writer.write(out_mol.as_bytes()).unwrap();
            }
            // i += 1;
        }

        prev_index = index;
        
        
    }
    Ok(())
}


fn parse_all_targets_in_dir(input_dir: path::PathBuf, output_dir: path::PathBuf, n_threads: usize) -> Result<(), io::Error> {
    let pool = ThreadPool::new(n_threads);

    let targets = fs::read_dir(&input_dir)?;
    for target in targets{
        let target = target?;
        if target.file_type()?.is_dir(){
            let target = target.file_name().to_str().unwrap().to_string();
            let input_dir = input_dir.clone();
            let output_dir = output_dir.clone();
            // println!("Spawning thread for {}", target);
            pool.execute(move || {
                parse_target(&target, input_dir, output_dir).unwrap();
            });
        }
    }
    pool.join();
    Ok(())
}
fn main() -> Result<(), io::Error> {

    let input_dir = path::PathBuf::from("/home/qzj517/POR-DD/data/gnina_LIT-PCBA_VS_data/lit-pcba");
    let output_dir = path::PathBuf::from("/home/qzj517/POR-DD/data/gnina_LIT-PCBA_VS_data/testing");
    parse_all_targets_in_dir(input_dir, output_dir, 15)?;
    Ok(())
}
