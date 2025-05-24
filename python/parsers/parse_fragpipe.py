#!/usr/bin/env python3
#to use: python parse_fragpipe.py -i path/psm.tsv -o /Users/lin/MS_data/psm_parsed.tsv --fileExt mzML
#!/usr/bin/env python3
import os
import re
import argparse
import pandas as pd

from modules.molecular_formula import MolecularFormula
from modules import atom_table, utils

# ------------------------------------------------------------------------------
# 1) Map observed mass → canonical mod name
# ------------------------------------------------------------------------------
MASS_TO_NAME = {
    57.0215: 'carbamidomethyl',
    15.9949: 'oxidation',
    0.9840:  'deamidation',
}

def parse_assigned_mods_with_pos(s):
    """
    From "2Q(0.9840),5R(0.9840)" → [(1,'Q','0.9840'), (4,'R','0.9840')]
    (zero-based positions)
    """
    out = []
    if not isinstance(s, str) or not s:
        return out
    for m in re.finditer(r'\b(\d+)\s*([A-Za-z])\s*\(\s*([0-9]+(?:\.[0-9]+)?)\s*\)', s):
        pos      = int(m.group(1)) - 1
        residue  = m.group(2).upper()
        mass_str = m.group(3)
        out.append((pos, residue, mass_str))
    return out

# ------------------------------------------------------------------------------
# 2) Compute molecular formula via existing routine
# ------------------------------------------------------------------------------
MODIFICATION_REGEX = re.compile(r'([A-Z_])\((.*?)(\(.*?\))?\)')

def extractModifications(seq, fixed_mods):
    # strip inline parentheses tokens
    clean = ''
    matches = list(MODIFICATION_REGEX.finditer(seq))
    keep = [True] * len(seq)
    for m in matches:
        for i in range(m.start()+1, m.end()):
            keep[i] = False
    for ch, k in zip(seq, keep):
        if k and ch != '_':
            clean += ch

    aas = utils.strToAminoAcids(clean)
    formula = MolecularFormula(clean)

    for pos, res, mass_str in fixed_mods:
        mass = float(mass_str)
        closest = min(MASS_TO_NAME, key=lambda k: abs(k - mass))
        name = MASS_TO_NAME[closest]
        # treat R-citrullination chemically as N-deamidation
        table_res = 'N' if (name == 'deamidation' and res == 'R') else res
        if name in atom_table.MODIFICATIONS and table_res in atom_table.MODIFICATIONS[name]:
            moddef = atom_table.get_mod(name, table_res)
            mval = atom_table.calc_mass(moddef)
            formula.add_mod(name, table_res)
            aas[pos].mod += mval

    return aas, formula

def extractAllModifications(seqs, mods_list):
    aa_all, form_all = [], []
    for seq, mods in zip(seqs, mods_list):
        aa, fm = extractModifications(seq, mods)
        aa_all.append(aa)
        form_all.append(fm)
    return aa_all, form_all

# ------------------------------------------------------------------------------
# 3) Main parser
# ------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        prog='parse_fragpipe',
        description='Convert FragPipe/MSFragger psm.tsv → IonFinder TSV'
    )
    parser.add_argument('-i','--input',   required=True, help='tab-delimited psm.tsv')
    parser.add_argument('-o','--output',  default='parsed.tsv', help='output TSV path')
    parser.add_argument('--fileExt',       default='mzML',    help='extension for precursorFile')
    args = parser.parse_args()

    # Load and normalize
    df = pd.read_csv(args.input, sep='\t', low_memory=False)
    df.columns = df.columns.str.strip()

    # 3.1 Strip [123] tags, uppercase sequences
    raw = df['Modified Peptide'].fillna(df['Peptide'])
    clean_seqs = raw.str.replace(r'\[\d+\]', '', regex=True).str.upper()

    # 3.2 Parse Assigned Modifications with positions
    mods_list = df['Assigned Modifications'].apply(parse_assigned_mods_with_pos).tolist()

    # 3.3 Compute formulas (for output) via legacy routine
    _, form_all = extractAllModifications(clean_seqs.tolist(), mods_list)

    # 3.4 Derive sampleName, scanNum, and precursorFile
    sample_name    = df['Spectrum'].str.split(r'\.', n=1).str[0]
    scanNum        = df['Spectrum'].str.extract(r'\.(\d+)\.')[0].astype(int)
    precursor_file = sample_name + f".{args.fileExt}"

    # 3.5 Build final sequence strings by injecting markers AFTER each letter
    sequences = []
    for idx, base_seq in enumerate(clean_seqs):
        # map position → list of (residue,mass)
        pos_map = {}
        for pos, res, mass in mods_list[idx]:
            pos_map.setdefault(pos, []).append((res, mass))

        s = ''
        for pos, letter in enumerate(base_seq):
            s += letter
            for res, mass in pos_map.get(pos, []):
                if res == 'R':
                    s += '*'
                else:
                    s += f'({mass})'
        sequences.append(s)

    # 3.6 Assemble and write
    out = pd.DataFrame({
        'sampleName':    sample_name,
        'precursorFile': precursor_file,
        'sequence':      sequences,
        'formula':       [str(f) for f in form_all],
        'scanNum':       scanNum
    })
    out.to_csv(args.output, sep='\t', index=False)
    print(f"Wrote {len(out)} rows to {args.output}")

if __name__ == '__main__':
    main()
