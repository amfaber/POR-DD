use std::{path::PathBuf, fs, collections::{HashMap}, time::Instant};

use polars::{prelude::{DataFrame, NamedFromOwned, NamedFrom, ParquetWriter}, series::Series};
#[allow(non_snake_case)]

// #[derive(Debug)]
// struct Mol{
//     coords: Vec<f32>,
//     properties: Vec<(String, String)>
// }

fn extract_mols<'a>(input: &'a str, activity: &'a str) -> impl Iterator<Item=Vec<(String, String)>> + 'a{
    let mols = input.split("$$$$\n");
    let mols = mols.filter_map(|molblock| {
        if molblock == ""{
            return None
        }
        let mut lines = molblock.lines();
        let name = lines.next().unwrap().to_string();
        // let name = name.unwrap();
        // let iter = lines.by_ref().skip(3).map_while(|line| {
        //     let coords: Vec<_> = line.split_whitespace().take(4).filter_map(|element| element.parse::<f32>().ok()).collect();
        //     if coords.len() == 3{
        //         Some(coords)
        //     }
        //     else {
        //         None
        //     }
        // });
        // let coord_block: Vec<_> = iter.flatten().collect();
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
        properties.push(("activity".to_string(), activity.to_string()));
        return Some(properties)
    });
    mols
}



pub fn extract_single_file<P: Into<PathBuf>>(sdfpath: P, outpath: P, activity: &str){
    let sdfpath: PathBuf = sdfpath.into();
    // let sdfpath = PathBuf::from(sdfpath);
    // let mut f = fs::File::open(sdfpath.join(&target).join("active.sdf")).unwrap();
    let now = Instant::now();
    let file_contents = fs::read_to_string(&sdfpath).unwrap();
    let reading = now.elapsed().as_millis();
    dbg!(reading);
    
    
    let mols = extract_mols(&file_contents, "N/A");
    let now = Instant::now();
    // mols.collect::<Vec<_>>();
    let mut properties: HashMap<_, _> = HashMap::new();
    mols.flatten().map(|ent| {
        let name = ent.0;
        let value = ent.1;
        properties.entry(name).or_insert(Vec::new()).push(value);
    }).fold(0, |_dummy_acc, _dummy_i| 0);
    // dbg!(&(test));
    // dbg!(&(properties));
    let iterating = now.elapsed().as_millis();
    dbg!(&iterating);
    // return ();
    // let (names, values): (Vec<_>, Vec<_>) = mols.flat_map(|mol| mol.0).unzip();
    // return ();
    // let now = Instant::now();
    // let mut properties: HashMap<_, _> = HashMap::new();
    // for (name, value) in names.into_iter().zip(values.into_iter()){
    //     properties.entry(name).or_insert(Vec::new()).push(value);
    // }
    // let hashmap_fill = now.elapsed().as_millis();
    // dbg!(&hashmap_fill);

    let now = Instant::now();
    let mut all_series = Vec::new();
    for (name, values) in properties.iter_mut(){
        all_series.push(Series::new(name, values));
    }
    let series_creation = now.elapsed().as_millis();
    dbg!(&series_creation);

    let mut df = DataFrame::new(all_series).unwrap();
    dbg!(&df);

    // let file = fs::OpenOptions::new().create(true).write(true).open(&outpath).unwrap();
    // let w = ParquetWriter::new(file);
    // w.finish(&mut df).unwrap();
}

