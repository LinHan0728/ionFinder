#!/usr/bin/env python3
#to use: python parse_fragpipe.py -i /Users/lin/MS_data/Tran_PXD026581/tran_set5_psm.tsv -o /Users/lin/MS_data/tran_set5_psm_parsed_tmt.tsv --fileExt mzML
#!/usr/bin/env python3
#!/usr/bin/env python3
import os
import re
import argparse
import pandas as pd

from modules.molecular_formula import MolecularFormula
from modules import atom_table, utils

# ------------------------------------------------------------------------------
# 1) Map observed mass → canonical mod name (for formula-building)
# ------------------------------------------------------------------------------
MASS_TO_NAME = {
    57.0215:    'carbamidomethyl',
    15.9949:    'oxidation',
    0.9840:     'deamidation',    # Q/N deamidation or R citrullination
    229.16293:  'tmt10plex',      # TMT-10
    304.20715:  'tmt16plex',      # TMTpro-16
}

def parse_assigned_mods_with_pos(s):
    """
    Parse "N-term(229.1629),1C(57.0215),2C(57.0215),4K(229.1629)"
    into [(0,'N-TERM','229.1629'), (0,'C','57.0215'), ...]
    """
    out = []
    if not isinstance(s, str) or not s:
        return out

    # N-term first
    for m in re.finditer(r'\bN-?term\(\s*([0-9]+(?:\.[0-9]+)?)\s*\)', s, flags=re.IGNORECASE):
        out.append((0, 'N-TERM', m.group(1)))

    # then residue-position
    for m in re.finditer(r'\b(\d+)\s*([A-Za-z])\s*\(\s*([0-9]+(?:\.[0-9]+)?)\s*\)', s):
        pos      = int(m.group(1)) - 1
        residue  = m.group(2).upper()
        mass_str = m.group(3)
        out.append((pos, residue, mass_str))

    return out

# ------------------------------------------------------------------------------
# 2) Build formula & per‐residue mod masses
# ------------------------------------------------------------------------------
MODIFICATION_REGEX = re.compile(r'([A-Z_])\((.*?)(?:\(.*?\))?\)')

def extractModifications(seq, mods):
    # strip inline parentheses to get clean seq
    clean = ''
    keep  = [True]*len(seq)
    for m in MODIFICATION_REGEX.finditer(seq):
        for i in range(m.start()+1, m.end()):
            keep[i] = False
    for ch,k in zip(seq,keep):
        if k and ch!='_':
            clean += ch

    aas     = utils.strToAminoAcids(clean)
    formula = MolecularFormula(clean)

    for pos, res, mass_str in mods:
        m = float(mass_str)
        closest = min(MASS_TO_NAME, key=lambda k: abs(k-m))
        name    = MASS_TO_NAME[closest]

        # treat R→citrullination as deamidation chemically
        table_res = 'N' if (name=='deamidation' and res=='R') else res
        if name in atom_table.MODIFICATIONS and table_res in atom_table.MODIFICATIONS[name]:
            moddef = atom_table.get_mod(name, table_res)
            delta  = atom_table.calc_mass(moddef)
            formula.add_mod(name, table_res)
            if 0 <= pos < len(aas):
                aas[pos].mod += delta

    return aas, formula

def extractAllModifications(seqs, mods_list):
    aa_all, fm_all = [], []
    for seq, mods in zip(seqs, mods_list):
        aa, fm = extractModifications(seq, mods)
        aa_all.append(aa)
        fm_all.append(fm)
    return aa_all, fm_all

# ------------------------------------------------------------------------------
# 3) Main parser
# ------------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(
        prog='parse_fragpipe',
        description='Convert FragPipe/MSFragger psm.tsv → IonFinder TSV'
    )
    p.add_argument('-i','--input',   required=True, help='tab-delimited psm.tsv')
    p.add_argument('-o','--output',  default='parsed.tsv', help='output TSV path')
    p.add_argument('--fileExt',      default='mzML',     help='extension for precursorFile')
    args = p.parse_args()

    df = pd.read_csv(args.input, sep='\t', low_memory=False)
    df.columns = df.columns.str.strip()

    # 3.1 Unmodified base peptides & parse mods
    base_peptides = df['Peptide'].str.upper().tolist()
    mods_list     = df['Assigned Modifications']\
                        .apply(parse_assigned_mods_with_pos)\
                        .tolist()
    _, form_all   = extractAllModifications(base_peptides, mods_list)

    # 3.2 sampleName, scanNum, precursorFile
    sample_name    = df['Spectrum'].str.split(r'\.', n=1).str[0]
    scanNum        = df['Spectrum'].str.extract(r'\.(\d+)\.')[0].astype(int)
    precursor_file = sample_name + f".{args.fileExt}"

    # 3.3 Build modified sequences
    sequences = []
    for pep, mods in zip(base_peptides, mods_list):
        # group by position
        pos_map = {}
        for pos,res,mass in mods:
            pos_map.setdefault(pos, []).append((res,mass))

        # collapse all N-term masses into one
        nterm_mass = 0.0
        for res,mass in pos_map.get(0, []):
            if res=='N-TERM':
                nterm_mass += float(mass)

        s = ''
        if nterm_mass:
            s += f'({nterm_mass:.4f})'

        # now each residue
        for i, aa in enumerate(pep):
            s += aa
            for res,mass in pos_map.get(i, []):
                if res=='N-TERM':
                    continue
                m    = float(mass)
                closest = min(MASS_TO_NAME, key=lambda k: abs(k-m))
                name = MASS_TO_NAME[closest]

                # R‐mod → star
                if name=='deamidation' and res=='R':
                    s += '*'
                else:
                    # if TMT use the numeric mass, else the mass
                    s += f'({m:.4f})'

        sequences.append(s)

    # 3.4 Parent & M/Z columns
    out = pd.DataFrame({
        'sampleName':         sample_name,
        'precursorFile':      precursor_file,
        'parentId':           df['Protein ID'],
        'parentProtein':      df['Entry Name'],
        'parentDescription':  df['Protein Description'],
        'parent_mz':          df['Observed M/Z'],
        'sequence':           sequences,
        'formula':            [str(f) for f in form_all],
        'scanNum':            scanNum
    })

    out.to_csv(args.output, sep='\t', index=False)
    print(f"Wrote {len(out)} rows to {args.output}")

if __name__ == '__main__':
    main()
