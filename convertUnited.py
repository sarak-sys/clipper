import pandas as pd
import logging

def is_united_format(df):
    """
    recognize if the datafile is a UNITED form,
    if light_area or heavy_area column exists 

    """
    has_light_or_heavy = any(('light_area' in col or 'heavy_area' in col) for col in df.columns)
    return  has_light_or_heavy




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
# Load the Excel file
# file_path = 'try_united.xlsx'
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
    
    if  "start" in df.columns and "end" in df.columns:
        df['Annotated Sequence'] = '[' + df['start'].str.upper() + '].' + df['seq'] + '.[' + df['end'].str.upper() + ']'
        logging.info("Successfully created column: 'Annotated Sequence'")
    else:
        if "before" in df.columns and "after" in df.columns:
            df['Annotated Sequence'] = '[' + df['before'].str.upper() + '].' + df['seq'] + '.[' + df['after'].str.upper() + ']'
            logging.info("Successfully created column: 'Annotated Sequence'")
        else:
            logging.error("Could not create 'Annotated Sequence'. Please ensure that the columns 'seq' and ('start', and 'end') or ('before' anf 'after') are present in your input file.")

    # # Drop the peptide_start and peptide_end columns
    # df = df.drop(columns=['start', 'end'])
    

    #if the Mod column exists it maps the values as in mod_map above, if more mapping needed, it can be added to mod_map whenever.
    #if no Mod column exists it creates one with value "dimet" (only needed to run).
    if "Mod" in df.columns:
        df["Mod"] = df["Mod"].fillna("Unmodified")
        df['Modifications'] = df.apply(lambda row: generate_modifications_combined(row['seq'], row['Mod']), axis=1)
        # Drop the original 'Mod' column
        #df = df.drop(columns=['Mod'])
        logging.info("Successfully created column: 'Modifications', Mod column has been converted to Modifications")
    else:
        df["Mod"] = "dimet"
        df['Modifications'] = df.apply(lambda row: generate_modifications_combined(row['seq'], row['Mod']), axis=1)
        logging.info("Successfully created column: 'Modifications'. Values are ''")

    # Drop the original 'protein' column
    df = df.drop(columns=['protein'])


    # # Drop the ratio columns 
    # df.drop(columns=ratio_cols, inplace=True)

    #create the abundances columns.
    rename_map = {}

    for col in df.columns:
        if col.startswith("light_area "):
            suffix = col[len("light_area "):]  # remove the prefix
            new_col = f"Abundances (Normalized): light{suffix}"
            rename_map[col] = new_col
        elif col.startswith("heavy_area "):
            suffix = col[len("heavy_area "):]
            new_col = f"Abundances (Normalized): heavy{suffix}"
            rename_map[col] = new_col

    # Apply renaming
    df.rename(columns=rename_map, inplace=True)
    # Save the modified dataframe to a new file
    # modified_file_path = 'Modified_try_United.xlsx'
    # df.to_excel(modified_file_path, index=False)
    return df