// use std::sync::mpsc::Sender;
// use std::alloc::Global;

use std::{path, io};
use clap::Parser;
use any_metric_parser::utils;

#[derive(Parser, Debug)]
struct Args {
    input_dir: path::PathBuf,
    output_dir: path::PathBuf,

    #[clap(short, long, default_value = "15")]
    inner_threads: usize,
    
    #[clap(short, long, default_value = "8")]
    outer_threads: usize,
    
    #[clap(short, long, multiple = true)]
    targets: Option<Vec<String>>,
    
    #[clap(short, long, default_value = "vsscore", possible_values = vec!["vsscore", "cnnaffinity", "cnnscore", "vina"])]
    compare_type: String
}

fn main() -> Result<(), io::Error> {
    let args = Args::parse();
    println!("{:?}", args);
    let input_dir = args.input_dir;
    let output_dir = args.output_dir;
    let inner_threads = args.inner_threads;
    let outer_threads = args.outer_threads;
    let targets = args.targets;
    let compare_type = args.compare_type;
    utils::parse_all_targets_in_dir(input_dir, output_dir, inner_threads, outer_threads, &targets, &compare_type)?;
    Ok(())
}
