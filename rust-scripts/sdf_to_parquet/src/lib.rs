use pyo3::prelude::*;
use std::any::Any;

use std::{path::PathBuf, fs, collections::{HashMap, HashSet}};
use polars::{prelude::{DataFrame, NamedFrom, ParquetWriter}, series::Series};

#[allow(non_snake_case)]

pub fn extract_mols<'a>(input: &'a str, activity: Option<&'a str>) -> impl Iterator<Item=Vec<(String, String)>> + 'a{
    let activity = activity.map(|s| s.to_string());
    let mols = input.split("$$$$");
    let mols = mols.filter_map(move |molblock| {
        if molblock == "" || molblock == "\n" || molblock == "\r\n"{
            return None
        }
        let mut lines = molblock.lines().enumerate();
        lines.next();
        let name = lines.next().unwrap().1.to_string();
        // dbg!(&name);
        let mut rest = lines.skip_while(|(_, line)| *line != "M  END").skip(1);
        let mut properties = Vec::new();
        loop{
            if let Some((i, mut propname)) = rest.next(){
                let mut this_it = propname.split('>');
                this_it.next();
                propname = this_it.next().unwrap();
                match rest.next(){
                    Some((j, propvalue)) => {
                        match propname.get(3..){
                            Some(name) => {
                                // dbg!(name);
                                // dbg!(propvalue);
                                properties.push((name.to_string(), propvalue.to_string()))
                            },
                            None => {
                                dbg!(propname);
                                dbg!(propvalue);
                                dbg!(i);
                                dbg!(j);
                                panic!("we shouldn't get here")
                            }
                        }
                    },
                    None => {
                        println!("Error encountered in sdf, saving properties only for the molecules before the error");
                        return None
                    }
                }
                rest.next();
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

trait MyDefault{
    fn mydefault() -> Self;
}

trait FillDefault: Any + 'static{
    // type Item: MyDefault;
    fn default_with_len(length: usize) -> Self where Self: Sized;
    fn push_default(&mut self);
    fn as_any(&self) -> &dyn Any;
    fn as_any_mut(&mut self) -> &mut dyn Any;
    fn into_any(self: Box<Self>) -> Box<dyn Any>;
}

impl MyDefault for String{
    fn mydefault() -> Self{
        "".to_string()
    }
}

impl MyDefault for f64{
    fn mydefault() -> Self{
        f64::NAN
    }
}

impl<T: MyDefault + Clone + 'static> FillDefault for Vec<T>{

    fn default_with_len(length: usize) -> Self {
        vec![T::mydefault(); length]
    }

    fn push_default(&mut self) {
        self.push(T::mydefault())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
    
    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }

    fn into_any(self: Box<Self>) -> Box<dyn Any>{
        self as Box<dyn Any>
    }
}


macro_rules! convert_to_series{
    ($vec:expr, $name:expr, $this_type: ty) => {
        if $vec.is::<Vec<$this_type>>(){
            Series::new($name, *$vec.downcast::<Vec<$this_type>>().unwrap())
        } else {
            panic!("we shouldn't be here")
        }
    };
    ($vec:expr, $name:expr, $this_type:ty, $($rest:ty),+) => {
        if $vec.is::<Vec<$this_type>>(){
            Series::new($name, *$vec.downcast::<Vec<$this_type>>().unwrap())
        } else {
            convert_to_series!($vec, $name, $($rest),+)
        }
    }
}

pub fn extract_single_file(sdfpath: &str, outpath: Option<&str>, float_cols: Option<Vec<String>>, activity: Option<&str>){
    let sdfpath: PathBuf = sdfpath.into();
    let default_out = sdfpath.parent().unwrap().join("props.parquet");
    let outpath = outpath.unwrap_or(default_out.to_str().unwrap());
    let outpath: PathBuf = outpath.into();
    let file_contents = fs::read_to_string(&sdfpath).unwrap();
    let mut float_set = HashSet::new();
    float_set.extend(float_cols.unwrap_or(
        ["minimizedAffinity", "minimizedRMSD", "CNNscore", "CNNaffinity", "CNN_VS", "CNNaffinity_variance"].into_iter()
        .map(|s| s.to_string()).collect::<Vec<_>>()
    ).into_iter());
    
    let mols = extract_mols(&file_contents, activity);
    let mut properties: HashMap<_, Box<dyn FillDefault + 'static>> = HashMap::new();
    let mut not_visited: HashSet<_> = HashSet::new();
    for (i, mol) in mols.enumerate(){
        // dbg!(i);
        for ent in mol{
            let name = ent.0;
            let value = ent.1;
            not_visited.remove(&name);
            if float_set.contains(&name){
                let value: f64 = value.parse().unwrap();
                properties.entry(name).or_insert(Box::new(Vec::<f64>::default_with_len(i))).as_any_mut()
                    .downcast_mut::<Vec<f64>>().unwrap().push(value);
            } else {
                properties.entry(name).or_insert(Box::new(Vec::<String>::default_with_len(i))).as_any_mut()
                    .downcast_mut::<Vec<String>>().unwrap().push(value);
            }
        }

        for missing_entry in not_visited.iter(){
            properties.get_mut(missing_entry).unwrap().push_default();
        }
        not_visited.extend(properties.keys().cloned());
    }

    let mut all_series = Vec::new();
    for (name, values) in properties.into_iter(){
        let values = values.into_any();
        all_series.push(convert_to_series!(values, &name, f64, String));
    }
    let mut df = DataFrame::new(all_series).unwrap();

    let file = fs::OpenOptions::new().create(true).write(true).truncate(true).open(&outpath).unwrap();
    let w = ParquetWriter::new(file);
    w.finish(&mut df).unwrap();
}

#[pyfunction]
#[pyo3(name = "extract_single_file")]
pub fn py_extract_single_file(sdfpath: &str, outpath: Option<&str>, float_cols: Option<Vec<String>>, activity: Option<&str>) -> PyResult<()>{
    Ok(extract_single_file(sdfpath, outpath, float_cols, activity))
}


/// A Python module implemented in Rust.
#[pymodule]
fn sdf_to_parquet(_py: Python, module: &PyModule) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(py_extract_single_file, module)?)?;
    Ok(())
}