#!/usr/bin/env python3
#to use: python parse_fragpipe.py -i /Users/lin/MS_data/tran_set5_psm.tsv -o /Users/lin/MS_data/tran_set5_psm_parsed.tsv --fileExt mzML
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
    0.9840:  'deamidation',   # Q/N deamidation or R citrullination
}

def parse_assigned_mods_with_pos(s):
    """
    Parse strings like "2Q(0.9840),5R(0.9840)" into
    a list of (pos,residue,mass_str) tuples (zero-based).
    """
    out = []
    if not isinstance(s, str) or not s:
        return out
    for m in re.finditer(
        r'\b(\d+)\s*([A-Za-z])\s*\(\s*([0-9]+(?:\.[0-9]+)?)\s*\)', s
    ):
        pos      = int(m.group(1)) - 1
        residue  = m.group(2).upper()
        mass_str = m.group(3)
        out.append((pos, residue, mass_str))
    return out

# ------------------------------------------------------------------------------
# 2) Build molecular formula from sequence + mods
# ------------------------------------------------------------------------------
MODIFICATION_REGEX = re.compile(r'([A-Z_])\((.*?)(\(.*?\))?\)')

def extractModifications(seq, fixed_mods):
    """
    Given seq (pure letters) and a list of (pos,residue,mass_str),
    return (aa_list, formula_object).
    """
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
            mval   = atom_table.calc_mass(moddef)
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
    parser.add_argument('--fileExt',       default='mzML',     help='extension for precursorFile')
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep='\t', low_memory=False)
    df.columns = df.columns.str.strip()

    # 3.1 Use the unmodified 'Peptide' column for base sequence
    base_peptides = df['Peptide'].str.upper().tolist()

    # 3.2 Parse Assigned Modifications into (pos,residue,mass_str)
    mods_list = df['Assigned Modifications'].apply(parse_assigned_mods_with_pos).tolist()

    # 3.3 Compute formulas
    _, form_all = extractAllModifications(base_peptides, mods_list)

    # 3.4 Derive sampleName, scanNum, and precursorFile from 'Spectrum'
    sample_name    = df['Spectrum'].str.split(r'\.', n=1).str[0]
    scanNum        = df['Spectrum'].str.extract(r'\.(\d+)\.')[0].astype(int)
    precursor_file = sample_name + f".{args.fileExt}"

    # 3.5 Build sequences with modifications AFTER each letter:
    sequences = []
    for idx, pep in enumerate(base_peptides):
        pos_map = {}
        for pos, res, mass in mods_list[idx]:
            pos_map.setdefault(pos, []).append((res, mass))

        s = ''
        for pos, letter in enumerate(pep):
            s += letter
            for res, mass in pos_map.get(pos, []):
                if res == 'R':
                    s += '*'
                else:
                    s += f'({mass})'
        sequences.append(s)

    # 3.6 Gather parent columns and Observed M/Z from input
    parent_id_list          = df['Protein ID']
    parent_protein_list     = df['Entry Name']
    parent_description_list = df['Protein Description']
    parent_mz_list          = df['Observed M/Z']

    # 3.7 Assemble output DataFrame
    out = pd.DataFrame({
        'sampleName':         sample_name,
        'precursorFile':      precursor_file,
        'parentId':           parent_id_list,
        'parentProtein':      parent_protein_list,
        'parentDescription':  parent_description_list,
        'parent_mz':          parent_mz_list,
        'sequence':           sequences,
        'formula':            [str(f) for f in form_all],
        'scanNum':            scanNum
    })

    out.to_csv(args.output, sep='\t', index=False)
    print(f"Wrote {len(out)} rows to {args.output}")

if __name__ == '__main__':
    main()
