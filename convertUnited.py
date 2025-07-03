import pandas as pd
import logging

def is_united_format(df):
    """
    recognize if the datafile is a UNITED form,
    if both:
    -dataframe first column is "seq" 
    - also has columns that starts with "ratio" 

    """

    first_col_is_seq = (df.columns[0].lower() == 'seq')
    has_ratio_cols = any(
    col.lower().startswith('ratio') and not col.lower().startswith('logratio')
    for col in df.columns)
    return first_col_is_seq and has_ratio_cols




def generate_modifications_combined(seq, mod):
    """
    create the Modifications column by combining the given columns seq , Mod .

    """
    modifications = []

    # Map N-terminal Mod
    mod_map = {
        'acet': 'Acetyl [N-Term]',
        'dimet': 'Dimethyl [N-Term]',
        'methyl': 'Methyl [N-Term]'
    }
    
    nterm_mod = mod_map.get(mod.lower())
    if nterm_mod:
        modifications.append(f"1x{nterm_mod}")

    # Find lysine positions for Dimethyl [K]
    lysine_positions = [f'K{i}' for i, aa in enumerate(seq, start=1) if aa == 'K']
    if lysine_positions:
        count = len(lysine_positions)
        modifications.append(f"{count}xDimethyl [{'; '.join(lysine_positions)}]")
    
    #if no suitable mapping is detcted,  modification is defined to be empty.
    return '; '.join(modifications) if modifications else None



def convert_united_to_pd(df):
    # # Load the Excel file
    # file_path = 'UF.xlsx'
    # df = pd.read_excel(file_path)

    # Extract content between "|" in the 'protein_id' column and rename the column to 'Master Protein Accessions'
    try:
        df['Master Protein Accessions'] = df['protein'].str.extract(r'\|([^|]+)\|')
        logging.info("Successfully created column: 'Master Protein Accessions'")
    except KeyError:
        logging.error("Could not create 'Master Protein Accessions'. Please ensure that the column 'protein' is present in your input file.")

    #replace values by NULLs
    df.replace(['','0 dev 0','num dev 0'], pd.NA,inplace=True)

    #Create the Annotated Sequence column using the columns seq,end,start
    try:
        df['Annotated Sequence'] = '[' + df['start'].str.upper() + '].' + df['seq'] + '.[' + df['end'].str.upper() + ']'
        logging.info("Successfully created column: 'Annotated Sequence'")
    except KeyError:
        logging.error("Could not create 'Annotated Sequence'. Please ensure that the columns 'seq', 'start', and 'end' are present in your input file.")

    # Drop the peptide_start and peptide_end columns
    df = df.drop(columns=['start', 'end'])
  
    # Apply the function
    try:
        df['Modifications'] = df.apply(lambda row: generate_modifications_combined(row['seq'], row['Mod']), axis=1)
        logging.info("Successfully created column: 'Modifications'")
    except KeyError:
        logging.error("Could not create 'Annotated Sequence'. Please ensure that the column 'Mod' is present in your input file.")
   
    # Drop the original 'protein' column
    df = df.drop(columns=['protein'])


    # Drop the original 'Mod' column
    df = df.drop(columns=['Mod'])


    # Set  reference intensity
    reference_intensity = 1000

    # Identify ratio columns 
    ratio_cols = [col for col in df.columns if col.startswith("ratio") and not col.startswith("logratio")]
    
    if not(ratio_cols):
       logging.info("No 'ratio' columns found. Abundances will NOT be created.") 

    # Create new normalized abundance columns
    for i, col in enumerate(ratio_cols, start=1):
        new_col = f"Abundances (Normalized): Repeat{i}"
        df[new_col] = df[col] * reference_intensity
        logging.info(f'Column "{col}" has been converted to "{new_col}".')


    # Drop the ratio columns 
    df.drop(columns=ratio_cols, inplace=True)

    # Save the modified dataframe to a new file
    # modified_file_path = 'Modified_United.xlsx'
    # df.to_excel(modified_file_path, index=False)
    return df