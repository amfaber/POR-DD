use pyo3::prelude::*;

use std::{path::PathBuf, fs, collections::{HashMap}, time::Instant};

use polars::{prelude::{DataFrame, NamedFrom, ParquetWriter}, series::Series};
#[allow(non_snake_case)]

fn extract_mols<'a>(input: &'a str, activity: Option<&'a str>) -> impl Iterator<Item=Vec<(String, String)>> + 'a{
    let activity = activity.map(|s| s.to_string());
    let mols = input.split("$$$$\n");
    let mols = mols.filter_map(move |molblock| {
        if molblock == ""{
            return None
        }
        let mut lines = molblock.lines();
        let name = lines.next().unwrap().to_string();
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
        properties.push(("name".to_string(), name));
        if let Some(activity) = &activity{
            properties.push(("activity".to_string(), activity.to_string()));
        }
        return Some(properties)
    });
    mols
}


#[pyfunction]
fn extract_single_file(sdfpath: &str, outpath: &str, activity: Option<&str>) -> PyResult<()>{
    let sdfpath: PathBuf = sdfpath.into();
    let outpath: PathBuf = outpath.into();
    let now = Instant::now();
    let file_contents = fs::read_to_string(&sdfpath).unwrap();
    let reading = now.elapsed().as_millis();
    // dbg!(reading);
    
    
    let mols = extract_mols(&file_contents, activity);
    let now = Instant::now();
    let mut properties: HashMap<_, _> = HashMap::new();
    mols.flatten().map(|ent| {
        let name = ent.0;
        let value = ent.1;
        properties.entry(name).or_insert(Vec::new()).push(value);
    }).fold(0, |_dummy_acc, _dummy_i| 0);

    let iterating = now.elapsed().as_millis();
    // dbg!(&iterating);

    let now = Instant::now();
    let mut all_series = Vec::new();
    for (name, values) in properties.iter_mut(){
        all_series.push(Series::new(name, values));
    }
    let series_creation = now.elapsed().as_millis();
    // dbg!(&series_creation);

    let mut df = DataFrame::new(all_series).unwrap();
    // // dbg!(&df);

    let file = fs::OpenOptions::new().create(true).write(true).open(&outpath).unwrap();
    let w = ParquetWriter::new(file);
    w.finish(&mut df).unwrap();
    Ok(())
}



/// A Python module implemented in Rust.
#[pymodule]
fn sdf_to_parquet(_py: Python, module: &PyModule) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(extract_single_file, module)?)?;
    Ok(())
}