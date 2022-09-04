use std::{path::PathBuf, fs, io::{BufReader, Seek}, io::Read, io::{BufRead, self, Write}, sync::mpsc};
use threadpool::ThreadPool;
use polars::prelude::*;
use std::collections::HashMap;


pub fn test(){
    let pool = ThreadPool::new(8);
    let input_dir = "/home/qzj517/POR-DD/data/gnina_processed/LIT-PCBA";
    let output_dir = PathBuf::from("/home/qzj517/POR-DD/data/gnina_processed/LIT-PCBA_protein_best");
    let lf = LazyFrame::scan_parquet(
        "/home/qzj517/POR-DD/complete_testing/full_filtered.parquet".into(),
        ScanArgsParquet::default()).unwrap();
    
    let df = lf.
    select(
        vec![cols(["directory", "ligands", "index", "protein"])]
    ).
    sort_by_exprs(vec![col("directory"), col("ligands"), col("index")], vec![false, false, false], false).
    collect().unwrap();
    
    let directory_ligands_dfs = df.partition_by(["directory", "ligands"]).unwrap();
    for directory_ligands_df in directory_ligands_dfs{
        let protein_dfs = directory_ligands_df.partition_by(["protein"]).unwrap();
        let first_row = directory_ligands_df.get(0).unwrap();
        let (target, ligands) = match (&first_row[0], &first_row[1]){
            (AnyValue::Utf8(target), AnyValue::Utf8(ligands)) => (target.to_string(), ligands.to_string()),
            _ => panic!("bah")
        };
        dbg!(&target);
        dbg!(&ligands);
        // dbg!(protein);
        let mut map = HashMap::new();
        // let mut file_positions = vec![0];
        let (tx, rx) = mpsc::channel();
        
        for protein_df in protein_dfs{
            // dbg!(protein);
            let protein = match protein_df.get(0).unwrap()[3]{
                AnyValue::Utf8(protein) => protein.to_string(),
                _ => panic!("bah")
            };
            let idx = protein_df["index"].clone();
            let tx = tx.clone();
            let input_dir = input_dir.clone();
            let target = target.clone();
            let ligands = ligands.clone();
            let protein = protein.clone();
            pool.execute(move || {
                // protein_df;
                let mols = load_mols(input_dir, &target, &ligands, &protein, &idx);
                tx.send((protein, mols));
            });
            // let test: Vec<Option<i64>>= protein_df["index"].i64().unwrap().into_iter().collect();
        }
        drop(tx);
        for received in rx{
            let (protein, mols) = received;
            dbg!(&protein);
            map.insert(protein, mols);
        }
        let protein_names = directory_ligands_df["protein"].utf8().unwrap().into_iter();
        let mut full_mols = String::new();
        for protein_name in protein_names{
            let protein_name = protein_name.unwrap();
            full_mols.push_str(&map.get_mut(protein_name).unwrap().next().unwrap());
        }
        let name = ligands.splitn(2, "_").nth(1).unwrap();
        // name.push_str(".sdf");
        let mut write_path = output_dir.join(target).join(name);
        write_path.set_extension("sdf");
        dbg!(&write_path);
        fs::create_dir(write_path.parent().unwrap());
        let mut writer = fs::OpenOptions::new().write(true).create(true).open(write_path).unwrap();
        writer.write_all(full_mols.as_bytes());
        // open(write_path);
        // println!("{}", full_mols);
        // break;
        // println!("Target: {}, ligands: {} done", target, ligands); 
    }
    
    // let df = df!("a" => &[1, 2, 3, 4], "b" => &[0, -1, 3, 5]).unwrap();
    // df.map

    // dbg!(directory_ligands_dfs);
}

fn next_mol(bufreader: &mut BufReader::<fs::File>, buffer: &mut Vec<u8>, protein: &str) -> u64{
    let bytes_read = bufreader.read_until(b'$', buffer).unwrap();
    if bytes_read == 0{
        panic!("EOF reached before finding all molecules")
    }
    buffer.pop();
    buffer.extend_from_slice(b"> <protein>\n");
    buffer.extend_from_slice(protein[0..4].as_bytes());
    buffer.extend_from_slice(b"\n\n$0000");
    let buffer_len = buffer.len();
    // buffer.resize(buffer_len + 4, 0);
    bufreader.read_exact(&mut buffer[buffer_len-4..buffer_len]).unwrap();
    (bytes_read + 4) as u64
}

fn load_mols<P: Into<PathBuf>>(input_dir: P,
    target: &str,
    ligands: &str,
    protein: &str,
    idx: &Series,
    // file_positions: &mut Vec<u64>,
    ) -> std::vec::IntoIter<String>{
    let input_dir: PathBuf = input_dir.into();
    let file_path = input_dir.join(target).join(ligands).join(protein).join("gnina_out.sdf");
    let idx_iter = idx.i64().unwrap().into_iter();
    let file = fs::File::open(file_path).unwrap();
    // let mut bufreader = BufReader::new(file);
    let mut buffer = Vec::new();
    let mut mols : Vec<String> = Vec::new();
    let mut i = 0;
    let mut bufreader = BufReader::new(file);
    for relevant_idx in idx_iter{
        let relevant_idx = relevant_idx.unwrap() as usize;
        // let mut position = *file_positions.last().unwrap();
        // file.seek(io::SeekFrom::Start(position));
        loop{
            next_mol(&mut bufreader, &mut buffer, protein);
            // file_positions.push(position);
            if i == relevant_idx{
                unsafe {
                    let mol = std::str::from_utf8_unchecked(&buffer).to_string();
                    mols.push(mol);
                }
                buffer.clear();
                i += 1;
                break;
            }
            buffer.clear();
            i += 1;
        }
        // file = bufreader.into_inner();

    }
    mols.into_iter()
}