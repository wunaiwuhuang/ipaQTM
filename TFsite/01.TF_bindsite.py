import pandas as pd
import os

work_dir = "/data1/wuguojia/data/IPA_QTM_tcga/TFsite/"
os.chdir(work_dir)
os.makedirs("./all_tf", exist_ok=True)

# ReMap
remap_path = "./remap2022/remap2022_nr_macs2_hg38_v1_0.bed"
remap = pd.read_csv(remap_path, sep='\t', header=None, usecols=[0,1,2,3], names=["chr","start","end","TF"])
remap = remap[remap["chr"].isin([f"chr{i}" for i in range(1, 23)])]
remap["TF"] = remap["TF"].str.replace(r':.*', '', regex=True)

# ENCODE
encode_path = "./encode_v2av3/00.TF.anno.bed"
encode = pd.read_csv(encode_path, sep='\t', header=None, usecols=[0,1,2,3], names=["chr","start","end","TF"])
encode["chr"] = "chr" + encode["chr"].astype(str)
encode = encode[encode["chr"].isin([f"chr{i}" for i in range(1, 23)])]

# combine and save
tfsite = pd.concat([remap, encode]).drop_duplicates()
for tf, df in tfsite.groupby("TF"):
    outfile = os.path.join("all_tf", f"{tf}.txt")
    df.to_csv(outfile, sep="\t", header=False, index=False)
# or: tfsite.to_parquet("TF_binding_sites.parquet")
