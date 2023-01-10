use sdf_to_parquet;
fn main(){
    // let inp = "/home/qzj517/POR-DD/data/raw_data/FDA_drugs/test_3D_opt_1216.sdf";
    let inp = "/home/qzj517/POR-DD/Internal_assays/Enamine_3mil_screen/best.sdf";
    let out = Some("/home/qzj517/POR-DD/rust-scripts/sdf_to_parquet/idk.parquet");
    sdf_to_parquet::extract_single_file(inp, out, None, None);
}