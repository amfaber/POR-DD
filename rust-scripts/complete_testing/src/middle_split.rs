#![allow(unused_imports)]

use clap::Parser;
use glob;
// use std::io::BufReader;
use std::path;
use std::fs;
use std::io::{Seek, BufReader, BufRead};
use regex::Regex;

// BYTE LOCATION OF IDX CHANGE: 4980722428

pub fn split_middle_at_different_file_idx() {
    let re = Regex::new(r"(?:AID\d+)(?P<activity>_(?:in)?active)?(?:_(?P<idx>\d+))?-lig_(?P<protein>\w{4})-rec_docked_rescore$").unwrap();
    let dir = path::PathBuf::from("/home/qzj517/POR-DD/data/gnina_LIT-PCBA_VS_data/lit-pcba/IDH1");
    let file_name = dir.join("newdefault.summary");
    let mut file = fs::File::open(file_name).expect("Can't find file");
    let byte_len = file.seek(std::io::SeekFrom::End(0)).unwrap();
    // let halfway = byte_len / 2 - 17967235;

    let halfway = 4980722428 - 1000;
    println!("{}", halfway);
    file.seek(std::io::SeekFrom::Start(halfway)).expect("seeking failed");
    let mut bufreader = BufReader::new(file);
    let mut buf = String::new();
    bufreader.read_line(&mut buf).expect("reading line failed");
    let mut actual_start = halfway + buf.len() as u64;
    buf.clear();
    bufreader.read_line(&mut buf).expect("reading line failed");
    let idx: String;
    {
        let fname = buf.split_whitespace().nth(4).expect(&format!("Did not find 5 fields in line;\n{}\nAt {}", buf, actual_start));
        idx = re.captures(fname).expect("regex didnt match").name("idx").expect("didnt find index in regex").as_str().to_string();
    }
    // println!("{}", buf);
    actual_start += buf.len() as u64;
    loop{
        buf.clear();
        bufreader.read_line(&mut buf).expect("reading line failed");
        let new_idx: &str;
        {
            let fname = buf.split_whitespace().nth(4).expect(&format!("Did not find 5 fields in line;\n{}\nAt {}", buf, actual_start));
            new_idx = re.captures(fname).expect("regex didnt match").name("idx").expect("didnt find index in regex").as_str();
        }
        let newfname = buf.split_whitespace().nth(4).expect(&format!("Did not find 5 fields in line;\n{}\nAt {}", buf, actual_start));
        if new_idx != idx{
            break;
        }
        actual_start += buf.len() as u64;
    }


    println!("{}", actual_start)
}