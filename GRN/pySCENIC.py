import loompy
import pyscenic
import subprocess
import pandas as pd
import numpy as np

# File loom
f_loom_path_scenic = "TimeCourse.loom"
ds = loompy.connect(f_loom_path_scenic, mode='r')

#STEP 1: Gene regulatory network inference, and generation of co-expression modules
# transcription factors list
f_tfs = "hs_hgnc_tfs.txt" # human
cmd=f"pyscenic grn {f_loom_path_scenic} {f_tfs} -o adj.csv --num_workers 20"
subprocess.run(cmd, shell=True, check=True)
adjacencies = pd.read_csv("adj.tsv", index_col=False, sep='\t')

adjacencies.head()






# Close the file
ds.close()
