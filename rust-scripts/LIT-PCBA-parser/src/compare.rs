use std::{path::PathBuf, fs};

pub fn compare<P: Into<PathBuf>>(input_equibind: P, input_gnina: P){
    let input_equibind: PathBuf = input_equibind.into();
    let input_gnina: PathBuf = input_gnina.into();
    // let input_gnina = PathBuf::from(input_gnina);
    let targets = fs::read_dir(&input_equibind).unwrap().
        filter_map(|entry_res| entry_res.map(|entry| {
            match entry.file_type(){
                Ok(ft) => {
                    if ft.is_dir() {
                        Some(entry.file_name())
                    }
                    else {
                        None
                    }
                },
                Err(_) => None,
            }
        }).ok().flatten());
    for target in targets{
        println!("{:?}", input_gnina.join(&target));
        println!("{:?}", input_equibind.join(&target));
    }
    
}