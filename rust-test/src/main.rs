use flate2::read::GzDecoder;
use std::fs;
use std::io::prelude::*;
// use std::io::Seek;
use std::str;
use std;
use std::path;
// use std::env;
use regex::Regex;
use std::collections::HashSet;
use std::io::Write;
// #[macro_use]
// extern crate fstrings;

fn _print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}

fn read_mol<T>(reader: &mut std::io::BufReader<T>, mol: &mut Vec<u8>) -> (String, String) where T: std::io::Read{
    // mol.clear();
    reader.read_until(b'$', mol).unwrap();
    reader.read_until(b'$', mol).unwrap();
    reader.read_until(b'$', mol).unwrap();
    reader.read_until(b'$', mol).unwrap();
    reader.read_until(b'\n', mol).unwrap();
    (String::from_utf8(mol.split(|x| *x == b'\n').next().unwrap().to_owned()).unwrap(),
     String::from_utf8(mol.split(|x| *x == b'\n').nth_back(3).unwrap().to_owned()).unwrap())
}

fn _look_for_mol_name<T>(reader: &mut std::io::BufReader<T>, mol: &mut Vec<u8>, _name: &str, _aff: &str) where T: std::io::Read {
    read_mol(reader, mol);
    // let old_name = mol.split(|x| *x == b'\n').next().unwrap().clone();
    // {
    let copy = mol.clone();
    let mut split_iter = copy.split(|x| *x == b'\n');
    let old_name = split_iter.next().unwrap();
    // }
    let mut new_name;
    loop {
        read_mol(reader, mol);
        let mut split_iter = mol.split(|x| *x == b'\n');
        new_name = split_iter.next().unwrap();
        if new_name != old_name {
            break;
        }
    }
    let sdf_aff = str::from_utf8(mol).unwrap().lines().nth_back(2).unwrap();
    // assert_eq!(format!("{:.3}", _aff), format!("{:.3}", sdf_aff));

}


fn main() {
    // let args: Vec<String> = env::args().collect();
    // let head = args[1].parse::<i32>().unwrap();
    // let head = 100;

    let dirname = path::PathBuf::from("/home/qzj517/POR-DD/data/gnina_LIT-PCBA_VS_data/lit-pcba");
    let output_dir = path::PathBuf::from("/home/qzj517/POR-DD/data/gnina_LIT-PCBA_VS_data/best_scoring");

    let summary_file_name = "newdefault.summary.gz";
    
    // let mut mol = vec![0u8];
    // let mut line = String::new();
    
    
    let mut current_name = String::new();
    let mut _current_aff = String::new();

    let mut current_sdf_file = String::new();
    let mut maybe_sdf_bufreader = None;
    let re = Regex::new(r"(?P<prefix>AID\d+)(?P<activity>_(?:in)?active)?(?:_(?P<idx>\d+))?-lig_(?P<protein>\w{4})-rec_docked_rescore$").unwrap();
    let mut already_opened_files = HashSet::<String>::new();
    let mut maybe_current_writer_file = None;
    
    for dir in dirname.read_dir().unwrap() {
        let dir = dir.unwrap().path();
        let summary_file_path = dir.join(summary_file_name);
        let summary_file = fs::File::open(summary_file_path).unwrap();
        let summary_decoder = GzDecoder::new(summary_file);
        let mut summary_bufreader = std::io::BufReader::new(summary_decoder);
        let mut header = String::new();
        summary_bufreader.read_line(&mut header).unwrap();
        let header = header.split_whitespace().collect::<Vec<&str>>();
        // println!("{:?}", header);
    
        // let mut i = 0;
        for line in summary_bufreader.lines(){
            
            let line: Vec<&str> = line.as_ref().unwrap().split_whitespace().collect();
            
            let mut mol: Vec<u8> = Vec::new();
            if current_sdf_file != line[4]{
                
                current_name = line[1].to_string();
                _current_aff = line[2].to_string();
                current_sdf_file = line[4].to_string();
                let sdf_file_path = dir.join(current_sdf_file.replace("_rescore", ".sdf.gz"));
                let sdf_file = fs::File::open(sdf_file_path).unwrap();
                let sdf_decoder = GzDecoder::new(sdf_file);
                maybe_sdf_bufreader = Some(std::io::BufReader::new(sdf_decoder));
                
                let properties = match re.captures(&current_sdf_file){
                    Some(res) => res,
                    None => {
                        panic!("{}", current_sdf_file);
                    }
                };
                let activity = match properties.name("activity"){
                    Some(res) => &res.as_str()[1..],
                    None => "active",
                };
                let new_writer_path = output_dir.join(dir.file_name().unwrap()).join(activity).join(format!("{}_protein.sdf", &properties["protein"]));
                if already_opened_files.insert(new_writer_path.to_str().unwrap().to_string()){ //Insert returns true if the value was not already present
                    fs::create_dir_all(new_writer_path.parent().unwrap()).unwrap();
                    maybe_current_writer_file = Some(fs::File::create(new_writer_path).unwrap());
                }
                else {
                    maybe_current_writer_file = Some(fs::OpenOptions::new().append(true).open(new_writer_path).unwrap());
                }
                read_mol(maybe_sdf_bufreader.as_mut().unwrap(), &mut mol);
            }
            else if current_name == line[1]{
                if let Some(reader) = maybe_sdf_bufreader.as_mut(){
                    read_mol(reader, &mut mol);
                    continue;
                }
            }
            else{
                current_name = line[1].to_string();
                _current_aff = line[2].to_string();
                if let Some(reader) = maybe_sdf_bufreader.as_mut(){
                    let (sdf_name, sdf_aff) = read_mol(reader, &mut mol);
                    assert_eq!(current_name, sdf_name);
                    // assert_eq!(format!("{:.3}", _current_aff), format!("{:.3}", sdf_aff));
                    // , &current_name, &_current_aff
                }
            }
            let mut mol = String::from_utf8(mol).unwrap();
            
            mol.insert_str(mol.len() - 5, &format!("> <{}>\n{}\n\n> <{}>\n{}\n\n", header[5], line[5], header[6], line[6]));
            if let Some(writer) = maybe_current_writer_file.as_mut(){
                writer.write(mol.as_bytes()).unwrap();
            }
            // i += 1;
            // if i > head {
                // break;
            // }
        }
    }
}
