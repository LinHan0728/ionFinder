#!/usr/bin/env python3
import os
import re
import argparse
import pandas as pd

from modules.molecular_formula import MolecularFormula
from modules import atom_table, utils

# ------------------------------------------------------------------------------
# 0) Normalize variant mod names to our canonical keys
# ------------------------------------------------------------------------------
NORMALIZE_MOD = {
    'deamidated':       'deamidation',
    'citrullination':   'deamidation',
    'tmt6plex':         'tmt6plex',
    'tmt10plex':        'tmt10plex',
    'tmtpro16plex':     'tmt16plex',
    'tmt16plex':        'tmt16plex',
    'oxidation':        'oxidation',
    'carbamidomethyl':  'carbamidomethyl',
}

# ------------------------------------------------------------------------------
# 1) Parse the “Modifications” column into (pos,residue,mod_name)
# ------------------------------------------------------------------------------
def parse_pd_mods(s):
    out = []
    if not isinstance(s, str) or not s:
        return out

    for part in s.split(';'):
        part = part.strip()
        # N-term
        m = re.match(r'^(?:N-?Term)\(([^\)]+)\)$', part, flags=re.IGNORECASE)
        if m:
            out.append((0, 'N-TERM', m.group(1)))
            continue
        # residue-pos, e.g. "M3(Oxidation)"
        m = re.match(r'^([A-Za-z])(\d+)\(([^\)]+)\)$', part)
        if m:
            res, idx, mod = m.groups()
            out.append((int(idx) - 1, res.upper(), mod))
    return out

# ------------------------------------------------------------------------------
# 2) Build molecular formula
# ------------------------------------------------------------------------------
def extract_form(seq, mods):
    aas     = utils.strToAminoAcids(seq)
    formula = MolecularFormula(seq)

    for pos, res, mod in mods:
        raw = mod.lower()
        key = NORMALIZE_MOD.get(raw, raw)
        # remap R‐citrullination to N‐deamidation chemistry
        table_res = 'N' if (key == 'deamidation' and res == 'R') else res

        if key in atom_table.MODIFICATIONS and table_res in atom_table.MODIFICATIONS[key]:
            moddef = atom_table.get_mod(key, table_res)
            mval   = atom_table.calc_mass(moddef)
            formula.add_mod(key, table_res)
            if res == 'N-TERM':
                aas[0].mod += mval
            elif 0 <= pos < len(aas):
                aas[pos].mod += mval

    return formula

# ------------------------------------------------------------------------------
# 3) Main parser
# ------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        prog='parse_pd32',
        description='Parse PD3.2 PSM → IonFinder TSV'
    )
    parser.add_argument('-i','--input',   required=True, help='PD3.2 PSM TSV')
    parser.add_argument('-o','--output',  default='parsed.tsv', help='output TSV')
    parser.add_argument('--fileExt',      default='mzML',      help='precursorFile extension')
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep='\t', low_memory=False)
    df.columns = df.columns.str.strip()

    # 3.1 sampleName, scanNum, and parentId
    df['FileBase'] = (
        df['File.Name']
        .apply(lambda x: re.split(r'[\\/]', str(x))[-1])
        .str.split('.', n=1).str[0]
    )
    sampleName    = df['FileBase']
    scanNum       = df['First.Scan']
    parentId      = df['Master.Protein.Accessions']  # <--- grab this column

    # 3.2 precursorFile
    precursorFile = sampleName + f".{args.fileExt}"

    # 3.3 clean base sequence from Annotated.Sequence
    def clean_seq(a):
        m = re.match(r'^\[[^\]]+\]\.([A-Za-z]+)\.\[[^\]]+\]$', a)
        if m:
            return m.group(1).upper()
        return re.sub(r'[^A-Za-z]', '', a).upper()

    base_seqs = df['Annotated.Sequence'].apply(clean_seq)

    # 3.4 parse modifications
    mods_list = df['Modifications'].apply(parse_pd_mods).tolist()

    # 3.5 build sequence + formula with numeric mass shifts,
    #     collapsing all N-term mods into one block
    sequences = []
    formulas  = []

    for pep, mods in zip(base_seqs, mods_list):
        # build formula
        formula = extract_form(pep, mods)
        formulas.append(str(formula))

        # group mods by position
        pos_map = {}
        for pos, res, mod in mods:
            pos_map.setdefault(pos, []).append((res, mod))

        # sum all N-terminal masses
        nterm_mass = 0.0
        for res, mod in pos_map.get(0, []):
            if res == 'N-TERM':
                raw       = mod.lower()
                key       = NORMALIZE_MOD.get(raw, raw)
                table_res = 'N' if key=='deamidation' else res
                moddef    = atom_table.get_mod(key, table_res)
                nterm_mass += atom_table.calc_mass(moddef)

        # start building the annotated sequence
        s = ''
        if nterm_mass:
            s += f'({nterm_mass:.4f})'

        for i, aa in enumerate(pep):
            s += aa
            # append all mods at this residue
            for res, mod in pos_map.get(i, []):
                if res == 'N-TERM':
                    continue
                raw       = mod.lower()
                key       = NORMALIZE_MOD.get(raw, raw)
                table_res = 'N' if (key=='deamidation' and res=='R') else res
                moddef    = atom_table.get_mod(key, table_res)
                mval      = atom_table.calc_mass(moddef)

                if key=='deamidation' and res=='R':
                    s += '*'
                else:
                    s += f'({mval:.4f})'

        sequences.append(s)

    # 3.6 assemble and write
    out = pd.DataFrame({
        'sampleName':    sampleName,
        'precursorFile': precursorFile,
        'parentId':      parentId,      # <--- include here
        'sequence':      sequences,
        'formula':       formulas,
        'scanNum':       scanNum
    })
    out.to_csv(args.output, sep='\t', index=False)
    print(f"Wrote {len(out)} rows to {args.output}")

if __name__ == '__main__':
    main()
