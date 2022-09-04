#![allow(warnings)]
use std::io::{BufRead, BufReader, Read};
use std::{io, fs, time};

pub fn test1(){
    let filepath = "/home/qzj517/POR-DD/data/gnina_processed/LIT-PCBA/ADRB2/AID492947_inactive_T/3p0g_protein/gnina_out.sdf";
    let mut f = fs::File::open(filepath).unwrap();
    let mut bufreader = BufReader::new(f);
    let mut buffer = String::new();
    let now = time::Instant::now();
    let mut i = 1;
    let n_mols = 1000;
    let lines_per_mol = 100;
    let iter = bufreader.lines();
    for line in iter{

        if i == n_mols*lines_per_mol{
            break
        }
        i += 1;
    }
    dbg!(now.elapsed().as_micros());
    println!("");
}

pub fn test2(){
    let filepath = "/home/qzj517/POR-DD/data/gnina_processed/LIT-PCBA/ADRB2/AID492947_inactive_T/3p0g_protein/gnina_out.sdf";
    let mut f = fs::File::open(filepath).unwrap();
    let mut bufreader = BufReader::new(f);
    let mut buffer = String::new();
    let now = time::Instant::now();
    let mut i = 1;
    let n_mols = 1000;
    let lines_per_mol = 100;
    // let iter = bufreader.lines();
    loop{
        bufreader.read_line(&mut buffer);
        if i == n_mols*lines_per_mol{
            break
        }
        i += 1;
        // buffer.clear();
    }
    dbg!(now.elapsed().as_micros());
    let n_mols_found = buffer.split("$$$$\n").map(|s| s.to_string()).collect::<Vec<String>>().len();
    dbg!(n_mols_found);
    println!("");
}

pub fn test3() -> u128{
    let filepath = "/home/qzj517/POR-DD/data/gnina_processed/LIT-PCBA/ADRB2/AID492947_inactive_T/3p0g_protein/gnina_out.sdf";
    let mut f = fs::File::open(filepath).unwrap();
    let mut bufreader = BufReader::new(f);
    let mut buffer = Vec::<u8>::new();
    let now = time::Instant::now();
    let mut rest = [0, 0, 0, 0];
    let mut i = 1;
    let n_mols = 1000;
    let lines_per_mol = 100;
    let mut mols: Vec<String> = Vec::new();
    // let iter = bufreader.lines();
    // let mut slice = 0;
    loop{
        bufreader.read_until(b'$', &mut buffer);
        bufreader.read_until(b'$', &mut buffer);
        bufreader.read_until(b'$', &mut buffer);
        bufreader.read_until(b'$', &mut buffer);
        bufreader.read_until(b'\n', &mut buffer);

        // let mut handle = bufreader.take(4);
        // handle.read(&mut buffer);
        // bufreader = handle.into_inner();

        // bufreader.read_exact(&mut rest);
        
        if i == n_mols{
            break
        }
        i += 1;
        // buffer.clear();
        // let buffer_len = buffer.len();
        unsafe {
        let mol = std::str::from_utf8_unchecked(&buffer).to_string();
        mols.push(mol);
        }
        buffer.clear();
        // slice = buffer_len;
    }
    now.elapsed().as_micros()
}
pub fn test4() -> u128{
    let filepath = "/home/qzj517/POR-DD/data/gnina_processed/LIT-PCBA/ADRB2/AID492947_inactive_T/3p0g_protein/gnina_out.sdf";
    let mut f = fs::File::open(filepath).unwrap();
    let mut bufreader = BufReader::new(f);
    let mut buffer = Vec::<u8>::new();
    let now = time::Instant::now();
    let mut rest = [0, 0, 0, 0];
    let mut i = 1;
    let n_mols = 1000;
    let lines_per_mol = 100;
    let mut mols: Vec<String> = Vec::new();
    // let iter = bufreader.lines();
    // let mut slice = 0;
    loop{
        bufreader.read_until(b'$', &mut buffer);

        // bufreader.read_until(b'$', &mut buffer);
        // bufreader.read_until(b'$', &mut buffer);
        // bufreader.read_until(b'$', &mut buffer);
        // bufreader.read_until(b'\n', &mut buffer);

        let mut handle = bufreader.take(4);
        handle.read(&mut buffer);
        bufreader = handle.into_inner();

        // bufreader.read_exact(&mut rest);
        
        if i == n_mols{
            break
        }
        i += 1;
        // buffer.clear();
        // let buffer_len = buffer.len();
        // if i % 20 == 0{
            unsafe {
            let mol = std::str::from_utf8_unchecked(&buffer).to_string();
            mols.push(mol);
            }
            // let mol = std::str::from_utf8(&buffer).unwrap().to_string();
            // mols.push(mol);
        // }
        buffer.clear();
        // slice = buffer_len;
    }
    let out = now.elapsed().as_micros();
    // dbg!(&(mols[0]));
    out
    // dbg!(mols.len());
    // println!("");
}
pub fn test5() -> u128{
    let filepath = "/home/qzj517/POR-DD/data/gnina_processed/LIT-PCBA/ADRB2/AID492947_inactive_T/3p0g_protein/gnina_out.sdf";
    let mut f = fs::File::open(filepath).unwrap();
    let mut bufreader = BufReader::new(f);
    let mut buffer = Vec::<u8>::new();
    let now = time::Instant::now();
    let mut rest = [0, 0, 0, 0];
    let mut i = 1;
    let n_mols = 1000;
    let lines_per_mol = 100;
    let mut mols: Vec<String> = Vec::new();
    // let iter = bufreader.lines();
    // let mut slice = 0;
    loop{
        bufreader.read_until(b'$', &mut buffer);

        // bufreader.read_until(b'$', &mut buffer);
        // bufreader.read_until(b'$', &mut buffer);
        // bufreader.read_until(b'$', &mut buffer);
        // bufreader.read_until(b'\n', &mut buffer);

        // let mut handle = bufreader.take(4);
        // handle.read(&mut buffer);
        // bufreader = handle.into_inner();

        bufreader.read_exact(&mut rest);
        buffer.extend_from_slice(&rest);
        
        if i == n_mols{
            break
        }
        i += 1;
        // buffer.clear();
        // let buffer_len = buffer.len();
        unsafe {
        let mol = std::str::from_utf8_unchecked(&buffer).to_string();
        mols.push(mol);
        }
        buffer.clear();
        // slice = buffer_len;
    }
    now.elapsed().as_micros()
}