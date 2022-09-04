use std::{path::PathBuf, fs, collections::{HashMap}, time::Instant};

use polars::{prelude::{DataFrame, NamedFromOwned, NamedFrom, ParquetWriter}, series::Series};
use itertools::multiunzip;
#[allow(non_snake_case)]

// #[derive(Debug)]
// struct Mol{
//     coords: Vec<f32>,
//     properties: Vec<(String, String)>
// }

fn extract_mols<'a>(input: &'a str, activity: &'a str) -> impl Iterator<Item=(f32, f32, f32, f32, String, String, String)> + 'a{
    let mols = input.split("$$$$\n");
    let mols = mols.filter_map(|molblock| {
        if molblock == ""{
            return None
        }
        let mut lines = molblock.lines();
        let name = lines.next().unwrap().to_string();
        // let name = name.unwrap();
        let iter = lines.by_ref().skip(3).map_while(|line| {
            let coords: Vec<_> = line.split_whitespace().take(4).filter_map(|element| element.parse::<f32>().ok()).collect();
            if coords.len() == 3{
                Some(coords)
            }
            else {
                None
            }
        });
        let coord_block: Vec<_> = iter.flatten().collect();
        let mut rest = lines.skip_while(|line| *line != "M  END").skip(1);
        let mut properties = Vec::new();
        loop{
            if let Some(propname) = rest.next(){
                let propvalue = rest.next().unwrap();
                properties.push((propname[3..propname.len()-1].to_string(), propvalue.to_string()));
                rest.next().unwrap();
            }
            else {
                break;
            }
        }
        // let out = Mol{
        //     coords: coord_block,
        //     properties,
        // };
        properties.push(("name".to_string(), name));
        return Some((
            properties[0].1.parse::<f32>().unwrap(),
            properties[1].1.parse::<f32>().unwrap(),
            properties[2].1.parse::<f32>().unwrap(),
            properties[3].1.parse::<f32>().unwrap(),
            properties[4].1.to_string(),
            properties[5].1.to_string(),
            activity.to_string(),
        ))
    });
    mols
}

pub fn compare<P: Into<PathBuf>>(input_gnina: P){
    let input_gnina: PathBuf = input_gnina.into();
    // let input_gnina = PathBuf::from(input_gnina);
    let targets = fs::read_dir(&input_gnina).unwrap().
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
    let mut big_df = DataFrame::empty();
    for target in targets{
        // let mut f = fs::File::open(input_gnina.join(&target).join("active.sdf")).unwrap();
        let now = Instant::now();
        let inactives_g = fs::read_to_string(&mut input_gnina.join(&target).join("inactive.sdf")).unwrap();
        let actives_g = fs::read_to_string(&input_gnina.join(&target).join("active.sdf")).unwrap();
        let reading = now.elapsed().as_millis();
        dbg!(reading);
        
        
        let gnina_mols = extract_mols(&actives_g, "active");
        let gnina_mols = gnina_mols.chain(extract_mols(&inactives_g, "inactive"));
        let now = Instant::now();
        let (
            minimizedAffinity,
            newdefault_CNNaffinity,
            newdefault_CNNscore,
            CNNaffinity_times_score,
            protein,
            name,
            activity,
            ): (
            Vec<f32>,
            Vec<f32>,
            Vec<f32>,
            Vec<f32>,
            Vec<String>,
            Vec<String>,
            Vec<String>,
            ) = multiunzip(gnina_mols);
        
        let len = minimizedAffinity.len();

        let minimizedAffinity = Series::from_vec("minimizedAffinity", minimizedAffinity);
        let newdefault_CNNaffinity = Series::from_vec("newdefault_CNNaffinity", newdefault_CNNaffinity);
        let newdefault_CNNscore = Series::from_vec("newdefault_CNNscore", newdefault_CNNscore);
        let CNNaffinity_times_score = Series::from_vec("CNNaffinity_times_score", CNNaffinity_times_score);
        let repeated_target: Vec<&str> = std::iter::repeat(&target.to_str().unwrap()).copied().take(len).collect();
        let target_series = Series::new("target", repeated_target);
        let protein = Series::new("protein", protein);
        let name = Series::new("name", name);
        let activity = Series::new("activity", activity);

        let df = DataFrame::new(vec![
            target_series,
            name,
            activity,
            protein,
            minimizedAffinity,
            newdefault_CNNaffinity,
            newdefault_CNNscore,
            CNNaffinity_times_score,
            ]).unwrap();
        big_df = big_df.vstack(&df).unwrap();
        let match_and_math = now.elapsed().as_millis();
        dbg!(match_and_math);

        // dbg!(df);
        
        // println!("Completed {}", &target.to_str().unwrap());
        // // break;
    }
    let file = fs::OpenOptions::new().create(true).write(true).open("df.parquet").unwrap();
    let w = ParquetWriter::new(file);
    w.finish(&mut big_df).unwrap();
}

fn main(){
    compare("/home/qzj517/POR-DD/data/gnina_LIT-PCBA_VS_data/parsed/vsscore_best")
}