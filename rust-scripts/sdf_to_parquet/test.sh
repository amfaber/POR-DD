maturin develop --release
python -c "import sdf_to_parquet;\
import pandas as pd;\
sdf_to_parquet.extract_single_file('/home/qzj517/POR-DD/data/raw_data/FDA_drugs/gnina.sdf',\
 '/home/qzj517/POR-DD/data/raw_data/FDA_drugs/props.parquet', float_cols = ['CNN_VS', 'CNNaffinity']);\
 df = pd.read_parquet('/home/qzj517/POR-DD/data/raw_data/FDA_drugs/props.parquet');\
 print(df.dtypes);"