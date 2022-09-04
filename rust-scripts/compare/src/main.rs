use std::{path::PathBuf, fs, collections::{HashMap}, time::Instant};

use polars::{prelude::{DataFrame, NamedFromOwned, NamedFrom, ParquetWriter}, series::Series};
use itertools::multiunzip;


#[derive(Debug)]
struct Mol{
    coords: Vec<f32>,
    properties: Vec<(String, String)>
}

fn extract_mols<'a>(actives: &'a str, inactives: &'a str) -> impl Iterator<Item=(String, Mol)> + 'a{
    let mols = inactives.split("$$$$\n").chain(actives.split("$$$$\n"));
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
        let out = Mol{
            coords: coord_block,
            properties,
        };
        return Some((name, out))
    });
    mols
}

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
    let mut big_df = DataFrame::empty();
    for target in targets{
        // let mut f = fs::File::open(input_gnina.join(&target).join("active.sdf")).unwrap();
        let now = Instant::now();
        let inactives_e = fs::read_to_string(&mut input_equibind.join(&target).join("inactive.sdf")).unwrap();
        let actives_e = fs::read_to_string(&input_equibind.join(&target).join("active.sdf")).unwrap();
        let inactives_g = fs::read_to_string(&mut input_gnina.join(&target).join("inactive.sdf")).unwrap();
        let actives_g = fs::read_to_string(&input_gnina.join(&target).join("active.sdf")).unwrap();
        let reading = now.elapsed().as_millis();
        dbg!(reading);
        
        
        
        let now = Instant::now();
        let equibind_mols = extract_mols(&actives_e, &inactives_e);
        let equibind_mols: HashMap<_,_> = equibind_mols.collect();
        let collect_equibind = now.elapsed().as_millis();
        dbg!(collect_equibind);
        
        let gnina_mols = extract_mols(&actives_g, &inactives_g);
        let now = Instant::now();

        let mut diffs = gnina_mols.filter_map(|(name, gnina_mol)| {
            equibind_mols.get(&name).map(
                |equibind_mol| {
                    let mut diff = 0.0;
                    let len = equibind_mol.coords.len() as u32;
                    assert_eq!(len, gnina_mol.coords.len() as u32);
                    for (coord_e, coord_g) in equibind_mol.coords.iter().zip(gnina_mol.coords.iter()){
                        diff += (coord_e - coord_g).powi(2);
                    }
                    
                    // diff /= 3 as f32;
                    // println!("{} {}", name, diff);
                    (
                        name,
                        diff,
                        len,
                        equibind_mol.properties[3].1.parse::<f32>().unwrap(),
                        equibind_mol.properties[4].1.parse::<f32>().unwrap(),
                        gnina_mol.properties[1].1.parse::<f32>().unwrap(),
                        gnina_mol.properties[3].1.parse::<f32>().unwrap()
                    )
                }
            )
        });
        let (names,
            diffs,
            lens,
            equibind_aff,
            equibind_VS,
            gnina_aff,
            gnina_VS): (Vec<String>,
                Vec<f32>,
                Vec<u32>,
                Vec<f32>,
                Vec<f32>,
                Vec<f32>,
                Vec<f32>) = multiunzip(diffs);
        
        // let names_series = ChunkedArray::<>::from_vec("names", names);
        // let diffs_series = ChunkedArray::<Float32Type>::from_vec("diffs", diffs);
        // let diffs_series: Series = diffs_series.into();
        let len = names.len();
        let names_series = Series::new("names", names);
        let diffs_series = Series::from_vec("diffs", diffs);
        let lens_series = Series::from_vec("lens", lens);
        let equibind_aff_series = Series::from_vec("equibind_aff", equibind_aff);
        let gnina_aff_series = Series::from_vec("gnina_aff", gnina_aff);
        let repeated_target: Vec<&str> = std::iter::repeat(&target.to_str().unwrap()).copied().take(len).collect();
        let target_series = Series::new("target", repeated_target);
        let equibind_VS_series = Series::from_vec("equibind_VS", equibind_VS);
        let gnina_VS_series = Series::from_vec("gnina_VS", gnina_VS);

        let df = DataFrame::new(vec![
            names_series,
            diffs_series,
            lens_series,
            equibind_aff_series,
            equibind_VS_series,
            gnina_aff_series,
            gnina_VS_series,
            target_series,
            ]).unwrap();
        big_df = big_df.vstack(&df).unwrap();
        let match_and_math = now.elapsed().as_millis();
        dbg!(match_and_math);

        dbg!(df);
        
        println!("Completed {}", &target.to_str().unwrap());
        // break;
    }
    let file = fs::OpenOptions::new().create(true).write(true).open("df.parquet").unwrap();
    let w = ParquetWriter::new(file);
    w.finish(&mut big_df).unwrap();
    // big_df.as_csv
    // ParquetWriter::<DataFrame>::to_parquet(&big_df, "./output.parquet").unwrap();
    // big_df
}

fn main(){
    compare("/home/qzj517/POR-DD/data/gnina_processed/LIT-PCBA_protein_best",
    "/home/qzj517/POR-DD/data/gnina_LIT-PCBA_VS_data/parsed/vsscore_as_equibind")
}