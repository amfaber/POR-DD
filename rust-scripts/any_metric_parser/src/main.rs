use std::{path, io};
use clap::{Parser, Subcommand};
use clap;
use any_metric_parser;

#[derive(Subcommand, Debug)]
#[clap(subcommand_negates_reqs = true)]
enum Commands{
    Split{
        input_dir: path::PathBuf,
        #[clap(short, long, multiple = true, required = true)]
        targets: Vec<String>,
    },
    Equibind{
        input_dir: path::PathBuf,
        output_dir: path::PathBuf,
        df_path: String
    },
    Compare{
        input_equibind: path::PathBuf,
        input_gnina: path::PathBuf,
    }
}

#[derive(Parser, Debug)]
#[clap(subcommand_negates_reqs = true)]
struct Args {

    #[clap(subcommand)]
    command: Option<Commands>,

    #[clap(required = true)]
    input_dir: Option<path::PathBuf>,

    #[clap(required = true)]
    output_dir: Option<path::PathBuf>,

    #[clap(short, long, default_value = "15")]
    inner_threads: usize,
    
    #[clap(short, long, default_value = "8")]
    outer_threads: usize,
    
    #[clap(short, long, multiple = true)]
    targets: Option<Vec<String>>,
    
    #[clap(short, long, default_value = "vsscore", possible_values = vec!["vsscore", "cnnaffinity", "cnnscore", "vina"])]
    compare_type: String,

    #[clap(short, long, action = clap::ArgAction::Set, default_value = "true")]
    same_as_equibind: bool,
}



fn main() -> Result<(), io::Error> {
    let args = Args::parse();
    println!("{:?}", args);
    // return Ok(());
    match args.command{
        Some(Commands::Split { input_dir, targets }) => {
            // let input_dir = args.input_dir;
            any_metric_parser::split::split_middle_at_different_file_idx(&input_dir, &targets[0], "newdefault.summary.gz");
        },
        Some(Commands::Equibind { input_dir, output_dir, df_path }) => {
            // let input_dir = args.input_dir;
            any_metric_parser::parse_equibind_data::parse_df(input_dir, output_dir, &df_path)
        },
        Some(Commands::Compare { input_equibind, input_gnina } ) => {
            any_metric_parser::compare::compare(input_equibind, input_gnina)
        },
        None => {
            let input_dir = args.input_dir.unwrap();
            let output_dir = args.output_dir.unwrap();
            let inner_threads = args.inner_threads;
            let outer_threads = args.outer_threads;
            let compare_type = args.compare_type;
            let same_as_equibind = args.same_as_equibind;
            let targets = args.targets;
            any_metric_parser::parse_gnina_data::parse_all_targets_in_dir(input_dir, output_dir, inner_threads, outer_threads, &targets, &compare_type, same_as_equibind)?;
        }
    };

    Ok(())
}