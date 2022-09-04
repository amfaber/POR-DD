#![allow(warnings)]
use regex::Regex;

use std::path;
use std::fs;
use std::time::Instant;
use std::io::{Seek, BufReader, BufRead};

pub fn test(){
    let re = Regex::new(r"(?:AID\d+)(?P<activity>_(?:in)?active)?(?:_(?P<idx>\d+))?-lig_(?P<protein>\w{4})-rec_docked_rescore$").unwrap();
    // let re = Regex::new(r"(?:AID\d+)(_(?:in)?active)?(?:_(\d+))?-lig_(\w{4})-rec_docked_rescore$").unwrap();
    // let dir = path::PathBuf::from("/home/qzj517/POR-DD/data/gnina_LIT-PCBA_VS_data/lit-pcba/TP53");
    let dir = path::PathBuf::from("/home/qzj517/POR-DD");
    let file_name = dir.join("whynothere.txt");
    let mut file = fs::File::open(file_name).expect("Can't find file");
    let mut bufreader = BufReader::new(file);
    // let mut buf = String::new();
    let mut lines = bufreader.lines();
    let mut lines = lines.skip(499500);
    lines.next();
    // bufreader.read_line(&mut buf);
    // buf.clear();
    let mut i = 1;
    let now = Instant::now();
    loop{
        // bufreader.read_line(&mut buf);
        let orig_line = lines.next().unwrap().unwrap();
        let line: Vec<&str> = orig_line.split_whitespace().collect();
        let index = line[0].parse::<usize>().unwrap()-1;
        let name = line[1].to_string();
        let mut vina = line[2].parse::<f32>().unwrap();
        let file = line[4].to_string();

        let (activity, index_from_file, protein) = parse_filename2(&file);

        let cnnaffinity = line[5].parse::<f32>().unwrap();
        let cnnascore = line[6].parse::<f32>().unwrap();
        // println!("{}", orig_line);
        // println!("activity: {}, idx: {}, protein: {}", activity, index_from_file.unwrap_or(-1), &protein);
        // println!("");
        if i == 500{
            break;
        }
        i += 1;
    }
    let no_regex_time = now.elapsed().as_micros();
    
    let mut i = 1;
    let now = Instant::now();

    loop{
        let line = lines.next().unwrap().unwrap();
        let line: Vec<&str> = line.split_whitespace().collect();
        let index = line[0].parse::<usize>().unwrap()-1;
        let name = line[1].to_string();
        let mut vina = line[2].parse::<f32>().unwrap();
        let file = line[4].to_string();

        let properties = match re.captures(&file){
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
            Some(res) => Some(res.as_str().parse::<i32>().unwrap()),
            None => None,
        };
        

        let cnnaffinity = line[5].parse::<f32>().unwrap();
        let cnnascore = line[6].parse::<f32>().unwrap();
        if i == 500{
            break;
        }
        i += 1;
    }
    let regex_time = now.elapsed().as_micros();
    println!("no regex: {}\nregex: {}", no_regex_time, regex_time);
    
}

pub fn parse_filename(filename: &str) -> (String, Option<i32>, String){

    let mut iter = filename.split("-");
    let mut iter2 = iter.next().unwrap().split("_");
    iter2.next();
    let mut activity;
    let mut index_from_file;
    match iter2.next(){
        Some(res) => {
            match res.parse::<i32>() {
                Ok(index) => {index_from_file = Some(index);
                    activity = "active".to_string();
                },
                Err(_) => {
                    activity = res.to_string();
                    index_from_file = Some(iter2.next().unwrap().parse::<i32>().unwrap());
                }
            }
        },
        None => {
            activity = "active".to_string();
            index_from_file = None;
        }

    }
    drop(iter2);
    let protein = iter.next().unwrap().split("_").nth(1).unwrap().to_string();
    (activity, index_from_file, protein)
}

pub fn parse_filename2(filename: &str) -> (String, Option<i32>, String){

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

pub fn parse_filename3(filename: &str) -> (String, Option<i32>, String){

    let mut iter = filename.split("-");
    let mut iter2 = iter.next().unwrap().split("_");
    iter2.next();
    let mut activity;
    let mut index_from_file;
    match iter2.next(){
        Some(res) => {
            activity = res.to_string();
            index_from_file = iter2.next().map(|ele| ele.parse::<i32>().unwrap())
        },
        None => {
            activity = "active".to_string();
            index_from_file = None;
        }

    }
    drop(iter2);
    let protein = iter.next().unwrap().split("_").nth(1).unwrap().to_string();
    (activity, index_from_file, protein)
}