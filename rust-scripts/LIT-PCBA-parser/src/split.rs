use crate::parse_gnina_data::parse_filename;
use std::{path, fs, io::{self, Seek, Read, BufReader, BufRead, Write}};

pub fn split_middle_at_different_file_idx(input_dir: &path::PathBuf, target: &str, file: &str) {

    let mut file_name = input_dir.join(target).join(file);
    std::process::Command::new("gunzip").args(vec![&file_name]).status().unwrap();
    file_name.set_extension("");
    println!("{:?}", file_name);
    let mut file = fs::File::open(&file_name).expect("Can't find file");
    let byte_len = file.seek(std::io::SeekFrom::End(0)).unwrap();
    let halfway = (byte_len as f64 * 0.45 ) as u64;
    println!("{}", halfway);
    
    
    // return ();
    
    // let halfway = 4980722428 - 1000;
    file.seek(std::io::SeekFrom::Start(halfway)).expect("seeking failed");
    let mut bufreader = BufReader::new(file);
    let mut buf = String::new();
    bufreader.read_line(&mut buf).expect("reading line failed");
    let mut actual_start = halfway + buf.len() as u64;
    buf.clear();
    bufreader.read_line(&mut buf).expect("reading line failed");
    let idx: Option<i32>;
    {
        let fname = buf.split_whitespace().nth(4).expect(&format!("Did not find 5 fields in line;\n{}\nAt {}", buf, actual_start));
        (_, idx, _) = parse_filename(fname)
    }
    // println!("{}", buf);
    actual_start += buf.len() as u64;
    loop{
        buf.clear();
        bufreader.read_line(&mut buf).expect("reading line failed");
        let new_idx: Option<i32>;
        {
            let fname = buf.split_whitespace().nth(4).expect(&format!("Did not find 5 fields in line;\n{}\nAt {}", buf, actual_start));
            (_, new_idx, _) = parse_filename(fname)
        }
        if new_idx != idx{
            break;
        }
        actual_start += buf.len() as u64;
    }
    println!("{}", actual_start);
    fn make_output_file(mut file_name: path::PathBuf, idx: u8) -> (path::PathBuf, fs::File){
        let name = file_name.file_name().unwrap().to_str().unwrap().to_string();
        let name = name.replacen(".", &format!("_{}.", idx), 1);
        file_name.set_file_name(&name);
        let write_to_file = fs::OpenOptions::new().create(true).write(true).open(&file_name).unwrap();
        (file_name, write_to_file)
    }

    let mut file = fs::File::open(&file_name).expect("Can't find file");
    let mut header = io::Read::by_ref(&mut file).bytes().
            filter_map(|c| c.ok()).
            take_while(|c| *c != b'\n').
            collect::<Vec<_>>();
    header.push(b'\n');
    
    
    let mut reader = std::io::Read::by_ref(&mut file).take(actual_start - header.len() as u64);
    let (writer_name, mut writer) = make_output_file(file_name.clone(), 1);

    writer.write(&header).unwrap();
    io::copy(&mut reader, &mut writer).unwrap();
    std::process::Command::new("gzip").args(vec!["-f", writer_name.to_str().unwrap()]).status().unwrap();
    
    
    // let mut reader = fs::File::open(&file_name).expect("Can't find file");
    let mut reader = file;
    reader.seek(std::io::SeekFrom::Start(actual_start)).unwrap();
    let (writer_name, mut writer) = make_output_file(file_name.clone(), 2);
    
    writer.write(&header).unwrap();
    io::copy(&mut reader, &mut writer).unwrap();
    std::process::Command::new("gzip").args(vec!["-f", writer_name.to_str().unwrap()]).status().unwrap();

    std::process::Command::new("gzip").args(vec![&file_name]).status().unwrap();
    let mut name = file_name.file_name().unwrap().to_str().unwrap().to_string();
    name.push_str(".gz");
    file_name.set_file_name(&name);
    name.insert(0, '_');
    let mut new_filename = file_name.clone();
    new_filename.set_file_name(name);
    std::process::Command::new("mv").args(vec![&file_name, &new_filename]).status().unwrap();
}