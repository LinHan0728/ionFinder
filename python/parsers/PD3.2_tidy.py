#!/usr/bin/env python3
#to use: python PD3.2_tidy.py -p alscsftmt_b5_01_PSMs.txt -i alscsftmt_b5_01_InputFiles.txt -o output.tsv
#!/usr/bin/env python3
import argparse
import os
import pandas as pd

def main():
    parser = argparse.ArgumentParser(
        description="Filter PSM to single‐protein and merge with input files."
    )
    parser.add_argument(
        "--psm",
        "-p",
        required=True,
        help="Path to the PSMs .txt file (tab‐delimited)",
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Path to the InputFiles .txt file (tab‐delimited)",
    )
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="Full path (including filename) for the merged TSV output",
    )

    args = parser.parse_args()

    # Ensure parent directory of output exists
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    # Read the PSM and input files
    psm_1 = pd.read_csv(args.psm, sep="\t", header=0)
    inputfile = pd.read_csv(args.input, sep="\t", header=0)

    # Filter to only those rows with exactly one protein
    psm_1_flt = psm_1[psm_1["X..Proteins"] == 1]

    # Left join on the File.ID column
    merged = psm_1_flt.merge(inputfile, on="File.ID", how="left")

    # Write out the result
    merged.to_csv(
        args.output,
        sep="\t",
        index=False,
        quoting=pd.io.common.csv.QUOTE_NONE
    )

    print(f"Written merged file to {args.output}")

if __name__ == "__main__":
    main()
