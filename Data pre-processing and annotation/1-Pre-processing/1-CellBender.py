import os
import sys
import subprocess

# Day for argument
if len(sys.argv) != 2:
    print("Usage: python cellbender.py <day>")
    sys.exit(1)

day = sys.argv[1]  # ex : "1" or "4" etc.

input_dir = "/home/ablanc/work/NetBIO2/P1/"
output_dir = "/home/ablanc/work/NetBIO2/P1/"
run_dir = f"run_count_J{day}"

input_file = os.path.join(input_dir, run_dir, "outs", "raw_feature_bc_matrix.h5")
output_file = os.path.join(output_dir, run_dir, "outs", "cleanedCellBender.h5")

if os.path.exists(input_file):
    cmd = [
        "cellbender", "remove-background",
        "--input", input_file,
        "--output", output_file
    ]
    print(f"Execution : {' '.join(cmd)}")
    subprocess.run(cmd)
else:
    print(f"[!] No file : {input_file}")

