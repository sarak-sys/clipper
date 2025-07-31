import pandas as pd
import re
import logging


MOD_MASS_MAP = {
    'n[28.0313]': 'Dimethyl [N-Term]',
    'n[34.0631]': 'Dimethyl [N-Term]',
    'K[28.0313]': 'Dimethyl [K]',
    'K[34.0631]': 'Dimethyl [K]',
    'n[42.0106]': 'Acetyl [N-Term]',
    'K[42.0106]': 'Acetyl [K]',
    'C[57.0215]': 'Carbamidomethyl [C]'
}

#check if this is a QUANT format
def is_quant_format(df):
    """
    the QUANT formait is detected by the presence of at least one column "Light MaxLFQ Intensity".
    - presence of light modified peptide and light modified peptide

    """

    has_light_col = "Light Modified Peptide" in df.columns
    has_heavy_col = "Heavy Modified Peptide" in df.columns

    return "Light MaxLFQ Intensity" in " ".join(df.columns) and has_heavy_col and has_light_col


def extract_mass_modifications(peptide):
    if not isinstance(peptide, str):
        return set()
    mods_found = set()
    for mass_tag, pd_format in MOD_MASS_MAP.items():
        if mass_tag in peptide:
            mods_found.add(pd_format)
    return mods_found

def extract_lysine_mod(seq, has_dimethyl_k):
    """
    Only adds NxDimethyl [K#; K#] if a mass-based Dimethyl [K] is detected in light/heavy sequences
    Otherwise, ignores lysines entirely
    Ensures the Modifications column reflects only what is chemically observed
    """
    if not isinstance(seq, str) or not has_dimethyl_k:
        return None
    # Extract core sequence from annotated format like [-].PEPTIDESEQ.[X]
    match = re.search(r'\.(.*?)\.', seq)
    if not match:
        return None
    core = match.group(1)
    lysine_positions = [f'K{i}' for i, aa in enumerate(core, start=1) if aa == 'K']
    if lysine_positions:
        count = len(lysine_positions)
        return f"{count}xDimethyl [{'; '.join(lysine_positions)}]"
    return None

def combined_modifications(row):
    light_pep = row['Light Modified Peptide']
    heavy_pep = row['Heavy Modified Peptide']
    seq = row['Annotated Sequence']

    light_mods = extract_mass_modifications(light_pep)
    heavy_mods = extract_mass_modifications(heavy_pep)

    base_mods = light_mods.union(heavy_mods)

    # Check if lysine site-mods should be added.
    add_k_mod = 'Dimethyl [K]' in base_mods
    k_mod_string = extract_lysine_mod(seq, add_k_mod)

    final_mods = set(base_mods)
    if k_mod_string:
        # Replace generic 'Dimethyl [K]' with detailed version
        final_mods.discard('Dimethyl [K]')
        final_mods.add(k_mod_string)

    return '; '.join(sorted(final_mods)) if final_mods else None

# find the light intensity and heavy intinsty columns 
# define the normalized abundances columns
def rename_maxlfq_columns(df):
    def rename(col):
        match = re.match(r"(.+?)\s+(Light|Heavy) MaxLFQ Intensity", col)
        if match:
            sample, label = match.groups()
            label = label.lower()
            new_col = f"Abundances (Normalized): {label}{sample}"
            #logging.info(f"Renamed column '{col}' to '{new_col}'")
            logging.info(f" '{new_col}' created successfully")
            return new_col
        return col  # keep unchanged if no match

    df.columns = [rename(col) for col in df.columns]
    return df


def convert_quant_to_pd(df):
    """
    Convert from QUANT format to PD format

    """
    # # Load the Excel file
    # file_path = 'combined_modified_peptide_label_quant.xlsx'
    # df = pd.read_excel(file_path)

    # define the normalized abundances columns
    try:
        df = rename_maxlfq_columns(df)
        logging.info("Successfully created Abundances columns.")
    except:
        logging.error("could Not create Abundances columns.")
        

    # define 'Master Protein Accessions'
    
    if 'Protein ID' in df.columns:
        df.rename(columns={'Protein ID': 'Master Protein Accessions'}, inplace=True)
        logging.info("Successfully created column: 'Master Protein Accessions'")
    else:
        if 'Protein' in df.columns:
            df['Master Protein Accessions'] = df['Protein'].str.extract(r'\|([^|]+)\|')
            logging.info("Successfully created column: 'Master Protein Accessions'")
        else:
            logging.error("Could not create 'Master Protein Accessions'. Ensure that 'Protein ID' or 'Protein' exists in your input file.")
            # return None


    
    #Create the Annotated Sequence column using the columns seq,end,start
    try:
        df['Annotated Sequence'] = '[-].' + df['Peptide Sequence'] + '.[-]'
        df = df.drop(columns=['Peptide Sequence'])
        logging.info("Successfully created column: 'Annotated Sequence'")
    except KeyError:
        logging.error("Could not create 'Annotated Sequence'. Please ensure that the column 'Peptide Sequence' is present in your input file.")
        # return None
    
    #extract modifications

    required_columns = ['Light Modified Peptide', 'Heavy Modified Peptide', 'Annotated Sequence']
    missing = [col for col in required_columns if col not in df.columns]
    if missing:
        logging.error(f"Missing required column(s): {', '.join(missing)}")
    
    try:
        df['Modifications'] = df.apply(combined_modifications, axis=1)
        logging.info("Successfully created column: 'Modifications'")
    except:
        logging.error("Could NOT create the column 'Modifications'")
        # return None



    # modified_file_path = 'Modified_Combined.xlsx'
    # df.to_excel(modified_file_path, index=False)
    return df