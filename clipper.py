import os
import re
import time
import logging
import concurrent.futures
from ast import literal_eval
from itertools import combinations, permutations

from datetime import datetime
from pathlib import Path
from tqdm import tqdm

import numpy as np
import pandas as pd
from pandas.errors import ParserError
from scipy.stats import f_oneway, rv_histogram
from statsmodels.stats.weightstats import ttest_ind
from statsmodels.stats.multitest import multipletests

from . import annutils
from .entry import Entry
from .logo import create_logo_helper
from .visualize import Visualizer
from . import protease_prediction as pp
from .globals import *

from convertQuant import convert_quant_to_pd , is_quant_format
from convertUnited import convert_united_to_pd , is_united_format

class Clipper:
    """Annotator class for processing and analyzing peptide proteomics data. All
    arguments are passed as a dictionary to the class constructor.

    Attributes:
        conditions (dict): A dictionary of conditions and their respective files.
        annot (pd.DataFrame): The input data as a Pandas DataFrame.
        outfolder (str): The path to the output folder.
        output_filetype (str): The type of output file to generate.
        logfolder (str): The path to the log folder.
        logfile (str): The path to the log file.
        figures (dict): A dictionary of generated figures.
        separate (bool): A flag to indicate if annotation data should be separate.
        pseudocounts (float): The pseudocount value for sequence logo generation.
        logo (str): The type of logo to generate.
        stat (bool): A flag to indicate if statistical calculations should be performed.
        df (pd.DataFrame): The dataframe to store processed data.
    """

    def __init__(self, args):

        """
        Initialize the Annotator class and set the attributes from the provided arguments.

        Parameters
        ----------
        args : dict
            A dictionary of arguments for setting various internal attributes related to input, output, data handling, and visualization parameters.
        """ 

        # global variables
        self.result_folder_name = result_folder_name
        self.data_folder_name = data_folder_name
        self.annotation_prefix = annotation_prefix

        self.alphafold_models_filename = alphafold_models_filename
        self.merops_filename = merops_filename
        self.merops_name_filename = merops_name_filename
        self.merops_sub_filename = merops_sub_filename
        self.protein_atlas_filename = protein_atlas_filename
        
        # output folders
        self.plot_protein_folder = plot_protein_folder
        self.plot_general_folder = plot_general_folder
        self.plot_fold_change_folder = plot_fold_change_folder
        self.plot_volcano_folder = plot_volcano_folder
        self.plot_piechart_folder = plot_piechart_folder
        self.plot_logo_folder = plot_logo_folder
        self.plot_enrichement_folder = plot_enrichement_folder
        self.plot_pathway_folder = plot_pathway_folder

        # input attributes
        self.infile_type = args["infile_type"]
        self.infile = args["infile"]
        self.alpha = args["alpha"]
        self.software = args["software"]

        # filtering and sanitizing
        self.level = args["level"]
        self.dropna = args["dropna"]
        self.fillna = args["fillna"]

        # annotation attributes
        self.separate = args["separate"]
        self.sleeptime = args["sleeptime"]
        self.noexo = args["noexo"]
        self.merops = None
        self.nomerops = args["nomerops"]
        self.calcstructure = args["calcstructure"]

        self.conditionfile = args["conditionfile"]
        self.proteasefile = args["proteasefile"]
        self.stat = args["stat"]
        self.pairwise = args["stat_pairwise"]
        self.significance = args["significance"]
        self.multipletesting = args["multipletesting"]
        self.multipletestingmethod = args['multipletestingmethod']
        self.available_models = None

        self.plot = args["visualize"]
        self.pymol_verbose = args["pymol_verbose"]
        self.cleavagevis = args["cleavagevis"]
        self.logo = args["logo"]
        self.logo_fc = args["logo_fc"]
        self.cleavagesitesize = args["cleavagesitesize"]
        self.volcano_foldchange = args["volcano_foldchange"]
        self.pseudocounts = args["pseudocounts"]
        self.enrichment = args["enrichment"]
        self.pathway = args["pathway"]

        # logging and file handling
        self.timestamp = args["timestamp"]
        self.logfile = args["logfile"]
        self.outfile_type = args["output_filetype"]
        self.outname = args["output_name"]

        self.conditions = None          # when condition file is read, this will contain a dictionary with condition and identifyers for each replicate in a dict
        self.conditionpermutations = [] # all permutations of conditions
        self.conditioncombinations = [] # all combinations of conditions
        self.figures = {}
        self.entwarnings = {'retrieval': [],    #used to store warnings during annotation, and will present them later in terminal for cleaner output
                            'peptides not found': [],
                            'model not available': []}

        self.basefolder = Path.cwd().absolute()
        self.resultfolder = self.basefolder / self.result_folder_name
        self.datafolder = self.basefolder / self.data_folder_name

        logging.info("Startup complete!\n")

    def validate_input_output_formats(self):

        """
        Validate the input and output file formats based on their file extensions.
        """

        if self.infile_type == "infer":
            if self.infile.endswith(".csv") or self.infile.endswith(".tsv") or self.infile.endswith(".txt"):
                self.infile_type = "csv"
            elif self.infile.endswith(".xlsx") or self.infile.endswith(".xls"):
                self.infile_type = "excel"
            else:
                self.raise_invalid_file_format_error("input")

        valid_output_formats = ["xlsx", "csv", "tsv", "pkl", "json"]
        if self.outfile_type not in valid_output_formats:
            self.raise_invalid_file_format_error("output")

        if self.conditionfile is not None and not os.path.exists(self.conditionfile):
            error_message = "Invalid condition file. Check if the path exists, and try again."
            logging.critical("Exit with code 6.")
            raise TypeError(error_message)
        
        if self.proteasefile is not None and not os.path.exists(self.proteasefile):
            error_message = "Invalid protease file. Check if the path exists, and try again."
            logging.critical("Exit with code 7.")
            raise TypeError(error_message)
            
        if self.conditionfile and not self.conditionfile.endswith(".txt"):
            self.raise_invalid_file_format_error("condition file")
        if self.proteasefile and not self.proteasefile.endswith(".txt"):
            self.raise_invalid_file_format_error("protease file")

    def raise_invalid_file_format_error(self, file_type):

        """
        Raise an error if the input or output file format is invalid.

        Parameters
        ----------
        file_type : str
            The type of the file that caused the error, e.g., "input", "output".
        """

        error_message = f"Invalid {file_type} format. Please select a valid {file_type} format and try again."
        logging.critical(f"{error_message}. Exiting with code 1.")
        raise TypeError(error_message)
    
    def set_input_output_paths(self):
        """
        Set the paths for the input and output files, and for the folders where different types of plots will be saved.
        """

        if self.outname:
            self.outfolder = self.resultfolder / self.outname
            self.outname = self.outname + self.annotation_prefix + self.outfile_type
        else:
            self.outfolder = self.resultfolder / self.timestamp
            self.outname = Path(self.infile).name.rsplit(".", 1)[0] + self.annotation_prefix + self.outfile_type

        self.temp_folder = self.outfolder / "tmp"
        self.protein_folder = self.outfolder / self.plot_protein_folder
        self.general_folder = self.outfolder / self.plot_general_folder
        self.fold_change_folder = self.outfolder / self.plot_fold_change_folder
        self.volcano_folder = self.outfolder / self.plot_volcano_folder
        self.piechart_folder = self.outfolder / self.plot_piechart_folder
        self.logo_folder = self.outfolder / self.plot_logo_folder
        self.enrichment_folder = self.outfolder / self.plot_enrichement_folder
        self.pathway_folder = self.outfolder / self.plot_pathway_folder

        self.folders = {"out": self.resultfolder, "data": self.datafolder, "protein": self.protein_folder, "general": self.general_folder,
                        "volcano": self.volcano_folder, "fold": self.fold_change_folder, "piechart": self.piechart_folder, "logo": self.logo_folder,
                        "enrichment": self.enrichment_folder, "pathway": self.pathway_folder}

    def load_data(self):
        """
        Load the input data into a Pandas DataFrame.
        """

        logging.info("Reading file...")
        self.read_file()
        logging.info("Read dataframe.")
        logging.info(f"Read input with {len(self.df)} peptides.")

    def set_software(self):
        """
        Identify the software used to generate the input file, based on the structure of the data.
        Generate the indexing patterns according to the identified software.
        """

        if self.software == "infer":
            try:
                self.df.loc[0, "Master Protein Accessions"]
                self.software = "pd"
            except KeyError:
                try:
                    self.df.loc[0, "PG.ProteinAccessions"]
                    self.software = "sm"
                except KeyError:
                    try:
                        self.df.loc[0, "Leading razor protein"]
                        self.software = "mq"   
                    except KeyError:
                        error_message = (
                            f"Invalid input. Please make sure input format is correct and contains accession and sequence columns with default names, and try again."
                        )
                        logging.critical("Invalid input. Exiting with code 4.")
                        raise TypeError(error_message)
            except:
                logging.critical("Invalid input. Exiting with code 5.")
                raise TypeError("Invalid input")

        logging.info(f"Input software is {self.software}")
        self.patterns = self.get_patterns()
        logging.debug(f"Patterns are: {self.patterns}")
        logging.info(f"Successfully generated indexing patterns for {self.software} input. See logfile for the found patterns.")
        logging.info("Format check complete.\n")

    def remove_empty_accessions(self):

        """
        Remove rows from the dataframe that contain empty accession numbers and writes affected lines to log.
        """

        col_acc = self.patterns['acc']
        invalid_acc = self.df[col_acc].isna()

        if invalid_acc.any():
            self.df = self.df[~invalid_acc].reset_index(drop=True)
            logging.warning(f'There were rows with no accession numbers in the loaded file, please check the log file for further information')

        for i, k in enumerate(invalid_acc):
            if k:
                logging.debug(f"Empty accession rows: {str(i + 1)}")
        

    def remove_empty_sequences(self):

        """
        Remove rows from the dataframe that contain empty sequences.
        """

        col_seq = self.patterns['seq']
        invalid_seq = self.df[col_seq].isna()
        if invalid_seq.any():
            logging.info(f"Empty sequence rows: {', '.join(map(str, invalid_seq.index + 1))}")
            self.df = self.df[~invalid_seq].reset_index(drop=True)

    def remove_invalid_alphabets(self):
            
        """
        Remove rows from the dataframe that contain invalid amino acid characters.
        """

        col_seq = self.patterns['seq']
        pattern = self.patterns['amino']
        invalid_alphabet = self.df[col_seq].str.contains(pattern)
        if invalid_alphabet.any():
            logging.info(f"Invalid sequence character rows: {', '.join(map(str, invalid_alphabet.index + 1))}")
            self.df = self.df[~invalid_alphabet].reset_index(drop=True)

    def sanitize(self):

        """
        Apply all sanitization steps to the dataframe:
        - Remove rows with empty accessions
        - Remove rows with empty sequences
        - Remove rows with invalid amino acid characters
        """
        
        self.remove_empty_accessions()
        self.remove_empty_sequences()
        self.remove_invalid_alphabets()

    def filter_df(self):

        """
        Filters the dataframe based on the level (nterm/quant) provided by the user.
        """

        try:
            if self.level == "nterm":
                logging.info("Level is set to include all N-terms")
                pattern = self.patterns['nterm']
            elif self.level == "quant":
                logging.info("Level is set to include quantified N-terms")
                pattern = self.patterns['nterm_label']
            elif self.level == "all":
                logging.info("Level is set to include all peptides")
            else:
                logging.warning("Unrecognized level argument. Falling back to all")
                return

            self.df = self.df[
                self.df[self.patterns['mod']].str.contains(pat=pattern, na=False)
            ].reset_index(drop=True)

        except:
            logging.critical("Exit with code 3")
            raise TypeError(
                f"Could not filter dataframe. Make sure software \
                {self.software} is correct, and try again."
            )

        if self.dropna:
            columns = self.df.columns[
                self.df.columns.str.contains(pat=self.patterns['quant'])
            ]
            self.df = self.df.dropna(subset=columns, how="all")
            
    def prepare(self):

        """
        Perform a series of preparatory steps based on user input, including validation, data loading, software setting, and folder creation.
        """

        self.validate_input_output_formats()
        self.set_input_output_paths()
        self.load_data()
        self.set_software()
        self.make_folders()
        logging.info("Initialization complete!\n")

        logging.info("Filtering and sanitizing input...")

        if self.level != "all":
            self.filter_df()

        self.sanitize()

    def make_folders(self):

        """
        Creates necessary output folders in the specified directory.
        """
        if not os.path.exists(self.outfolder):
            os.mkdir(self.outfolder)
        if not os.path.exists(self.temp_folder):
            os.mkdir(self.temp_folder)

        if self.plot:
            if not os.path.exists(self.general_folder):
                os.mkdir(self.general_folder)
            if not os.path.exists(self.fold_change_folder):
                os.mkdir(self.fold_change_folder)
            if not os.path.exists(self.volcano_folder):
                os.mkdir(self.volcano_folder)
            if not os.path.exists(self.piechart_folder):
                os.mkdir(self.piechart_folder)
            if self.stat:
                if self.cleavagevis in ['seq', 'both']:
                    if not os.path.exists(self.protein_folder):
                        os.mkdir(self.protein_folder)
                elif self.cleavagevis is not None:
                    logging.warning(f"-clvis argument '{self.cleavagevis}' is unknown, falling back to 'None'")
                    self.cleavagevis = None
                if self.logo:
                    if not os.path.exists(self.logo_folder):
                        os.mkdir(self.logo_folder)
                if self.enrichment:
                    if not os.path.exists(self.enrichment_folder):
                        os.mkdir(self.enrichment_folder)
                if self.pathway:
                    if not os.path.exists(self.pathway_folder):
                        os.mkdir(self.pathway_folder)

    def read_condition_file(self):
            
        """
        Reads and parses the condition file, storing information about conditions and corresponding channels.
        """

        with open(self.conditionfile, "r") as fh:
            self.conditions = {
                line.split()[0]: line.split()[1:] for line in fh.readlines()
            }

        for pair in combinations(self.conditions.keys(), 2):
            self.conditioncombinations.append(f'{pair[0]} vs. {pair[1]}')
        for pair in permutations(self.conditions.keys(), 2):
            self.conditionpermutations.append(f'{pair[0]} vs. {pair[1]}')

    def read_MEROPS(self):

        """
        Read the MEROPS data from csv files and store them in dataframes.
        """

        logging.info("Reading MEROPS data..")
        self.merops = pd.read_csv(self.datafolder / self.merops_filename)

        self.merops_name = pd.read_csv(self.datafolder / self.merops_name_filename)
        self.merops_name = self.merops_name[self.merops_name.type == "real"]

        self.merops_sub = pd.read_csv(self.datafolder / self.merops_sub_filename)
        logging.info("Read MEROPS data")

    def read_protein_atlas(self):

        """
        Read the Protein Atlas data from a TSV file and store it in a dataframe.
        """

        logging.info("Reading Protein Atlas data..")
        self.protein_atlas = pd.read_csv(self.datafolder / self.protein_atlas_filename, sep='\t')
        logging.info("Read Protein Atlas data.")

    def initialize_annotation(self):

        """
        Initialize a dataframe with empty annotation columns, with the same number of rows as the input dataframe.
        
        Parameters
        ----------
        length : int
            The number of rows in the input dataframe.
        """

        length = len(self.df)
        annot_columns = [
                        "query_sequence",
                        "query_accession",
                        "name",
                        "full_sequence",
                        "description",
                        "keywords",
                        "go_codes",
                        "go_names",
                        "proteoform_certainty%",
                        "acc_length",
                        "start_pep",
                        "end_pep",
                        "p1_position",
                        "cleavage_site",
                        f"p{self.cleavagesitesize}_p{self.cleavagesitesize}prime",
                        "nterm_annot",
                        "protease_uniprot",
                        "protease_merops_code",
                        "protease_merops_name",
                        ]

        self.annot = pd.DataFrame(
            columns=annot_columns,
            index=range(length),
        )
        
        # remove Merops columns if not wanted
        if self.nomerops:
            self.annot.drop(["protease_merops_code", "protease_merops_name"], axis=1, inplace=True)

        logging.info("Initialized annotation dataframe.\n")



    def read_file(self):
            
        """
        Reads the input file and returns a dataframe object. Supports CSV and Excel file formats.
        """

        try:
            if self.infile_type == "csv":
                logging.info("Input is csv")
                try:
                    df = pd.read_csv(self.infile, sep=",")
                except pd.errors.ParserError:
                    logging.info("Failed to parse with ',', trying with ';'")
                    try:
                        df = pd.read_csv(self.infile, sep=";")
                    except pd.errors.ParserError:
                        logging.info("Failed to parse with ';', trying with '\t'")
                        try:
                            df = pd.read_csv(self.infile, sep="\t")
                        except pd.errors.ParserError:
                            raise ValueError(f"Unsupported csv type: {self.infile_type}")
                        
            elif self.infile_type == "excel":
                logging.info("Input is excel")
                try:
                    df = pd.read_excel(self.infile, engine="openpyxl")
                except pd.errors.ParserError:
                    logging.info("Failed to parse excel file")
                    raise ValueError(f"Unsupported excel type: {self.infile_type}")
                
            else:
                raise ValueError(f"Unsupported file type: {self.infile_type}")
            
            #detect if the input file is QUANT or UNITED format
            if is_quant_format(df):
                logging.info("QUANT format is detected  — applying custom preprocessing")
                df = convert_quant_to_pd(df)
            elif is_united_format(df):
                logging.info("UNITED format is detected — applying custom preprocessing")
                df = convert_united_to_pd(df)
            else:
                logging.info("File format is not recognized as UNITED or QUANT — proceeding as normal")
        except (OSError, Exception) as err:
            self.handle_file_error(err)

        self.df = df

    def read_protease_file(self):

        """
        Reads information about proteases from a specified file.
        """

        with open(self.proteasefile, "r") as fh:
            self.proteases = {
                line.strip() for line in fh.readlines()
            }

    def handle_file_error(self, err):

        """
        Handles file reading errors by logging the error and raising an appropriate exception.
        Args:
            err (Exception): The error raised during file reading.
        """

        logging.critical(f"Could not read file due to error: {err}")
        logging.critical("Exit with code 2")
        raise TypeError(
            f"Could not read file. Make sure the path {self.infile} is correct, "
            f"file type {self.infile_type} is supported, and try again."
        )

    def infer_infile_annotation_status(self, arguments):
        """
        will check if the input file has already been annotated, and then decide whether to annotate or not.
        """
        infileannotation = {'uniprot': False,
                            'merops': False,
                            'exopeptidase': False
                            }
        cols = self.df.columns
        
        if 'acc_length' in cols:
            infileannotation['uniprot'] = True
        if 'protease_merops_code' in cols:
            infileannotation['protease_merops_code'] = True
        if 'exopeptidase' in cols:
            infileannotation['exopeptidase'] = True
        
        arguments['infileannotation'] = infileannotation
        
        return arguments
    
    def process_columns(self, pat):

        """
        Processes a specified column in the dataframe based on a pattern. Converts strings to numeric values and fills NaNs.
        Args:
            pat (str): A string pattern to identify the column(s) of interest.
        Returns:
            str: The pattern if column(s) matching the pattern exist, else None.
        """
            
        quant = self.df.columns[self.df.columns.str.contains(pat=pat)]

        if len(quant) > 0:
            try:
                columns = self.df.columns[self.df.columns.str.contains(pat)]
                for col in columns:
                    self.df[col] = self.df[col].astype(str).str.replace(",", ".").str.strip()

                    self.df[col] = pd.to_numeric(self.df[col], errors='coerce')
                    logging.debug(f"Converted values in {col} to numeric type")

                    mask = (self.df[col].notna()) & pd.to_numeric(self.df[col], errors='coerce').isna()
                    if mask.any():
                        logging.warning(f"Could not convert the following values in column {col}")

                    if self.fillna is not None:
                        # Fill NaN values with the user-specified value
                        self.df[col].fillna(float(self.fillna), inplace=True)
                        logging.info(f"Filled NaN values in column '{col}' with {float(self.fillna)}")
                logging.info('Finished converting quantification values to floats.')
                    
            except:
                logging.critical("Invalid input. Exiting with code 4.")
                raise TypeError(
                    f"Invalid input. Could not convert values to float. Make sure input format {self.software} is correct and there are no string literals in quant columns, and try again."
                    )
            return pat
        else:
            return None

    def get_patterns_sm(self):
            
        """
        Returns column patterns specific to 'sm' software.
        Raises:
            TypeError: If the input does not meet certain criteria.
        Returns:
            dict: A dictionary of column name patterns.
        """
        
        patterns = {}

        patterns['acc'] = "PG.ProteinAccessions"
        patterns['seq'] = "PEP.StrippedSequence"
        patterns['amino'] = "B|J|O|U|X|Z"

        try:
            test_for_first_row_nan = annutils.parse_acc(self.df.loc[0, patterns['acc']])
            if test_for_first_row_nan is None:
                logging.critical("The first row in the input file does not have an accession number. Please delete this row (or rows) manually and try again.")
            self.df.loc[0, patterns['seq']]
        except:
            logging.critical("Invalid input. Exiting with code 4.")
            raise TypeError(
                f"Invalid input. Please make sure input format {self.software} \
                is correct and contains sequence columns with default names, \
                and try again."
            )

        # Identify modification column
        if "P.MoleculeID" in self.df.columns:
            patterns['mod'] = "P.MoleculeID"
        elif "EG.PrecursorId" in self.df.columns:
            patterns['mod'] = "EG.PrecursorId"
        else:
            logging.critical("Invalid input. Exiting with code 4.")
            raise TypeError(
                f"Invalid input. Please make sure input format {self.software} \
                contains valid modification column (P.MoleculeID or EG.PrecursorId), \
                and try again."
            )

        # Check for the presence of columns and process if present
        pat_tmt = r"PEP\.TMT"
        pat_tot = r"EG\.TotalQuantity"
        patterns['quant'] = self.process_columns(pat_tmt) or self.process_columns(pat_tot)

        # Check modification type
        pat_tmt_mod = r"\[TMT"
        pat_dim_mod = r"Dimeth"

        mod_tmt = len(self.df[self.df[patterns['mod']].str.contains(pat_tmt_mod, na=False)])
        mod_dim = len(self.df[self.df[patterns['mod']].str.contains(pat_dim_mod, na=False)])

        if mod_tmt > 0:
            patterns['label'] = r"\[TMT"
            patterns['nterm'] = r"N-?ter"
            patterns['nterm_label'] = r"TMT.*_Nter"
            patterns['lysine_label'] = r"K\[TMT.{0,3}_Lys\]"
        elif mod_dim > 0:
            patterns['label'] = pat_dim_mod
            patterns['nterm'] = r"DimethNter0"
            patterns['nterm_label'] = r"\[DimethNter0\]"
            patterns['lysine_label'] = r"K\[DimethLys0\]"

        return patterns

    def get_patterns_pd(self):

        """
        Returns column patterns specific to 'pd' software.
        Raises:
            TypeError: If the input does not meet certain criteria.
        Returns:
            dict: A dictionary of column name patterns.
        """

        patterns = {}

        patterns['acc'] = "Master Protein Accessions"
        patterns['mod'] = "Modifications"
        patterns['nterm'] = r"\[N-Term\]"
        patterns['seq'] = "Sequence"
        patterns['amino'] = "B|J|O|U|X|Z"
        

        # go here
        first_row_is_nan = True
        while first_row_is_nan:
            if self.df.loc[0, patterns['acc']] is np.nan:
                self.df = self.df.drop([0]).reset_index(drop=True)
                logging.warning("Removed the current first row of the input file as the accession was not available.")
            else:
                first_row_is_nan = False

        try:
            annutils.parse_acc(self.df.loc[0, patterns['acc']])
            self.df.loc[0, patterns['seq']]
        except KeyError:
            try:
                # if sequence column is not present, check if modifications column is present
                annutils.parse_acc(self.df.loc[0, patterns['acc']])
                patterns['seq'] = "Annotated Sequence"
                annutils.parse_sequence(self.df.loc[0, patterns['mod']])
                patterns['amino'] = "\.[A-Z]*(B|J|O|U|X|Z)[A-Z]*\."
            except KeyError:
                logging.critical("Invalid input. Exiting with code 4.")
                raise TypeError(
                    f"Invalid input. Please make sure input \
                    format {self.software} is correct and contains \
                    sequence columns with default names, and try again."
                )

        # find the quant columns
        pat_scale = r'Abundances \(?Scaled\):.*'
        pat_norm = r'Abundances \(?Normalized\):.*'
        pat_grouped = r'Abundances \(Grouped\):.*'
        pat_raw = r'Abundance: .*'
        pat_other = r'Abundance.*'
        
        # Check for the presence of columns and process if present
        patterns['quant'] = self.process_columns(pat_scale) or \
                            self.process_columns(pat_norm) or \
                            self.process_columns(pat_grouped) or \
                            self.process_columns(pat_raw) or \
                            self.process_columns(pat_other)
        
        pat_label = r"TMT"
        pat_alt_label = r"Dimethyl"

        quant_tmt = len(self.df[self.df['Modifications'].str.contains(pat=pat_label, na=False)])
        quant_dimethyl = len(self.df[self.df['Modifications'].str.contains(pat=pat_alt_label, na=False)])

        if quant_tmt == 0:
            if quant_dimethyl == 0:
                patterns['label'] = None
                patterns['nterm_label'] = None
            else:
                patterns['label'] = pat_alt_label
                patterns['nterm_label'] = r"Dimethyl \[N-Term\]"
                patterns['lysine_label'] = r"Dimethyl \[K"

        else:
            patterns['label'] = pat_label
            patterns['nterm_label'] = r"TMT.* \[N-Term\]"
            patterns['lysine_label'] = r"TMT.{0,5} \[K"

        return patterns

    def get_patterns(self):
        """
        Returns column patterns based on the software used. This is determined by the 'software' attribute of the class.
        Raises:
            TypeError: If the 'software' attribute is not recognized.
        Returns:
            dict: A dictionary of column name patterns.
        """

        patterns = {}

        if self.software == 'sm':
            patterns = self.get_patterns_sm()
        elif self.software == 'pd':
            patterns = self.get_patterns_pd()
        else:
            logging.critical("Invalid software input. Exiting with code 4.")
            raise TypeError(f"Invalid software input. Please provide a valid software format and try again.")
        
        print(f'\n{datetime.now().strftime("%Y-%m-%d %H:%M:%S,%f")[0:-3]} [INFO] The column name patterns which were detected are:')
        for pattern, value in patterns.items():
            print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S,%f")[0:-3]} [INFO] - {pattern}: {value}')
        print("")

        return patterns

    def proteoform_check(self):

        """
        Computes the probability of proteoforms based on the master accession column in the dataframe.
        The probabilities are stored in the dataframe under a new column "proteoform_certainty%".
        """
            
        self.annot["proteoform_certainty%"] = self.df[self.patterns['acc']].apply(
            lambda x: 100 / len([i.strip() for i in x.split(';')])
        )
    
    def general_conditions(self):
        """
        Generates general statistics for the conditions supplied. If no conditions are given, default
        conditions are used. The method computes mean, standard deviation and coefficient of variation (CV)
        for each condition and stores them in the dataframe. It also computes fold changes and log2 fold 
        changes between each pair of conditions.
        """
        if self.conditions:
            def calc_stats(df, cols):
                mean = df[cols].mean(axis=1)
                std = df[cols].std(axis=1)
                return mean, std, std / mean

            def generate_condition_columns(pair, mean0, mean1):
                column_name = f"Fold change: {pair[0]} vs. {pair[1]}"
                column_log = f"Log2 fold change: {pair[0]} vs. {pair[1]}"
                fold_change = mean0 / mean1
                log_fold_change = np.log2(fold_change)
                return column_name, column_log, fold_change, log_fold_change

            for condition, channels in self.conditions.items():
                cols = self.df.columns[self.df.columns.str.contains("|".join(channels)) & self.df.columns.str.contains(self.patterns['quant'])]
                mean, std, cv = calc_stats(self.df, cols)
                self.annot[f"{condition}_mean"] = mean
                self.annot[f"{condition}_deviation"] = std
                self.annot[f"{condition}_CV"] = cv

            if len(self.conditions) > 1:
                for pair in permutations(self.conditions.keys(), 2):
                    cols0 = self.df.columns[self.df.columns.str.contains("|".join(self.conditions[pair[0]])) & self.df.columns.str.contains(self.patterns['quant'])]
                    cols1 = self.df.columns[self.df.columns.str.contains("|".join(self.conditions[pair[1]])) & self.df.columns.str.contains(self.patterns['quant'])]
                    mean0, _, _ = calc_stats(self.df, cols0)
                    mean1, _, _ = calc_stats(self.df, cols1)
                    column_name, column_log, fold_change, log_fold_change = generate_condition_columns(pair, mean0, mean1)
                    self.annot[column_name] = fold_change
                    self.annot[column_log] = log_fold_change


    def percentile_fold(self, percentile):

        """
        Checks fold change distribution and marks rows above a certain percentile.
        Args:
            percentile (float): The percentile above which rows will be marked.
        """

        left_cutoff = percentile
        right_cutoff = 1 - percentile

        if len(self.conditions) > 1:
            conditions_iter = permutations(self.conditions.keys(), 2)

            for pair in conditions_iter:
                column = f"Fold change: {pair[0]} vs. {pair[1]}"
                self.annot[f"Fold {pair[0]} vs. {pair[1]} significance"] = np.nan

                if self.significance == 'all':
                    subframe = self.annot
                elif self.significance == 'nterm':
                    is_internal = self.annot["nterm_annot"] == "Internal"
                    subframe = self.annot[~is_internal] if column.endswith("low") else self.annot[is_internal]
                else:
                    logging.info(f"Significance invalid argument {self.significance}, skipping")
                    continue

                finite_vals = subframe[column].dropna()
                finite_vals = finite_vals[np.isfinite(finite_vals)]
                hist = np.histogram(finite_vals)
                hist_dist = rv_histogram(hist)

                def classify_fold_change(value):
                    cd = hist_dist.cdf(value)
                    if cd > right_cutoff:
                        return "significant high"
                    elif cd < left_cutoff:
                        return "significant low"
                    return np.nan

                self.annot[f"Fold {pair[0]} vs. {pair[1]} significance"] = subframe[column].apply(classify_fold_change)

    def condition_statistics(self):
        """
        Performs t-test or ANOVA statistical significance tests based on the number of conditions.
        The test results are stored in the dataframe.
        """

        def perform_test(test_func, column_name, column_log, cols_per_condition, vals_per_condition):
            result = test_func(*vals_per_condition)
            statistic, p_value = result[0], result[1]
            self.annot[column_name] = p_value
            self.annot[column_log] = -np.log10(p_value)
        
        conditions = list(self.conditions.keys())
        if len(conditions) >= 2:
            test_func = ttest_ind if len(conditions) == 2 else f_oneway
            stat_name = f"{'Independent T-test p-value:' if len(conditions) == 2 else 'ANOVA p-value:'}"
            column_name = f"{stat_name} {' vs. '.join(conditions)}"
            column_log = f"-Log10 {column_name}"
            cols_per_condition = []
            vals_per_condition = []
            for cond in conditions:
                replicates = []
                for replicate in self.conditions[cond]:
                    replicates.extend([column for column in self.df.columns if re.search(replicate, column) and re.search(self.patterns['quant'], column)])
                cols_per_condition.append(replicates)
                vals_per_condition.append(np.log2(self.df[replicates]).T)
            perform_test(test_func, column_name, column_log, cols_per_condition, vals_per_condition)
        
        else:
            logging.warning("-stat (condition statistics) was requested, but <2 conditions were supplied through a condition file. Please recheck conditions.")

        if self.pairwise and len(conditions) > 2:
            for pair in permutations(self.conditions.keys(), 2):
                stat_name = "Independent T-test p-value:"
                column_name = f"{stat_name} {pair[0]} vs. {pair[1]}"
                column_log = f"-Log10 {column_name}"
                cols_per_condition = []
                vals_per_condition = []
                for cond in pair:
                    replicates = []
                    for replicate in self.conditions[cond]:
                        replicates.extend([column for column in self.df.columns if re.search(replicate, column) and re.search(self.patterns['quant'], column)])
                    cols_per_condition.append(replicates)
                    vals_per_condition.append(np.log2(self.df[replicates]).T)
                perform_test(ttest_ind, column_name, column_log, cols_per_condition, vals_per_condition)

    def correct_multiple_testing(self):

        """
        Corrects p-values for multiple testing using the chosen method.
        The corrected p-values are stored in the dataframe.
        """
            
        if not self.multipletesting or len(self.conditions) < 2:
            logging.info("Multiple testing corrected values are NOT used for visualizations")
            return   

        def _perform_correction(pvals, column_name, column_log, original_column_name):
            # Identify NaN values
            not_nan = ~np.isnan(pvals)

            # Remove NaN values before performing multiple testing correction
            pvals_no_nan = pvals[not_nan]

            # Check if there are p-values to correct
            if len(pvals_no_nan) > 0:
                corrected_pvals_no_nan = multipletests(pvals_no_nan, alpha=self.alpha, method=self.multipletestingmethod)[1]

                # Save the corrected p-values back to their original index positions
                corrected_pvals = np.empty_like(pvals)
                corrected_pvals[:] = np.nan
                corrected_pvals[not_nan] = corrected_pvals_no_nan

                # Insert the corrected p-value columns next to the original p-value columns
                original_column_idx = self.annot.columns.get_loc(original_column_name)
                self.annot.insert(original_column_idx + 2, column_name, corrected_pvals)
                self.annot.insert(original_column_idx + 3, column_log, -np.log10(corrected_pvals))

                logging.debug(f"Corrected p-values for {column_name} using {self.multipletestingmethod} method")
            else:
                logging.warning(f"No p-values to correct for {column_name}")

        if len(self.conditions.keys()) >= 2:
            stat_name = 'Independent T-test p-value' if len(self.conditions) == 2 else 'ANOVA p-value'
            column_name = f"Corrected {stat_name}: {' vs. '.join(self.conditions.keys())}"
            column_log = f"-Log10 {column_name}"

            original_column_name = f"{stat_name}: {' vs. '.join(self.conditions.keys())}"
            pvals = self.annot[original_column_name].values
            _perform_correction(pvals, column_name, column_log, original_column_name)

        if self.pairwise and len(self.conditions.keys()) > 2:
            stat_name = "Independent T-test p-value:"
            for pair in permutations(self.conditions.keys(), 2):
                original_column_name = f"{stat_name} {pair[0]} vs. {pair[1]}"
                column_name = f"Corrected {original_column_name}"
                column_log = f"-Log10 {column_name}"
                
                pvals = self.annot[original_column_name].values
                _perform_correction(pvals, column_name, column_log, original_column_name)

        logging.info(f"Finished multiple testing correction using the {self.multipletestingmethod} method.\n")

    def process_entry(self, loc):

        """
        Process a single entry and return the annotation.
        Args:
            loc (int): The location (index) of the entry in the dataframe.
        Returns:
            dict: A dictionary of annotation results.
        """
        acc = annutils.parse_acc(self.df.loc[loc, self.patterns['acc']])
        if self.patterns['seq'] == 'Annotated Sequence':
            seq = annutils.parse_sequence(self.df.loc[loc, self.patterns['seq']])
        else:
            seq = self.df.loc[loc, self.patterns['seq']]
        ent = Entry(acc, seq)
        ent.get_record(self.sleeptime)
        if ent.record is not None:
            ent.parse_general()
            ent.parse_cleavage(self.cleavagesitesize)
            if ent.cleavage_site is not None:
                ent.parse_protease()
                if self.nomerops is False:
                    ent.merops_protease(self.merops, self.merops_name)
            return annutils.map_dict(self.annot.loc[loc], ent.annot), ent.warnings
        else:
            return {"name": "HTTPError, not found"}, ent.warnings
    
    def entry_annotate(self, loc):

        """
        Single entry annotation function used with multiple threading.
        Args:
            loc (int): The location (index) of the entry in the dataframe.
        """

        self.annot.loc[loc], warnings = self.process_entry(loc)
        for key, value in warnings.items():
            if warnings[key] != []:
                self.entwarnings[key].append(value)

    def threaded_annotate(self, cores):

        """
        Annotation function using multiple threading. This function utilizes all available cores to
        process and annotate entries in the dataframe.
        """

        length = len(self.df)
        batch_length = min(cores, length)

        # if self.nomerops is False:
        #     self.read_MEROPS()
        #     logging.info("Read MEROPS data")

        logging.info("Fetching information from Uniprot...")
        with tqdm(range(0, length, batch_length), leave = 0) as t:
            for i in t:
                batch = list(range(i, i + batch_length))

                with concurrent.futures.ThreadPoolExecutor(
                    max_workers=batch_length
                ) as executor:
                    executor.map(self.entry_annotate, batch)

                elapsed = t.format_dict['elapsed']

        for key, values in self.entwarnings.items():
            if self.entwarnings[key] != []:
                logging.warning(f'There were "{key}" errors')
                for value in values:
                    logging.warning(f'  - {value[0]}')
                self.entwarnings[key] = []
        time.sleep(0.5)
        logging.info(f"Gathering info from Uniprot took {annutils.format_seconds_to_time(elapsed)}")
        

    def annotate_protein_atlas(self):

        """
        Annotates specific columns from Protein Atlas.
        
        This function reads from the protein atlas dataframe, selects specific columns, and merges 
        this data with the existing annotations based on 'Uniprot' identifier. 
        """

        self.read_protein_atlas()
        
        protein_atlas_columns = ['Uniprot', 'RNA tissue specific nTPM', 'RNA single cell type specific nTPM', 'Chromosome', 'Position', 'Protein class', 'Biological process', 'Molecular function', 'Disease involvement']
        protein_atlas_sub = self.protein_atlas[protein_atlas_columns].add_prefix('ProteinAtlas_')
        protein_atlas_sub.rename(columns={"ProteinAtlas_Uniprot": "Uniprot"}, inplace=True)

        # Check for duplicates in Uniprot column of Protein Atlas dataframe and drop them
        if protein_atlas_sub['Uniprot'].duplicated().any():
            logging.warning("Duplicates found in 'Uniprot' column of 'protein_atlas_sub'. Dropping duplicates.")
            protein_atlas_sub = protein_atlas_sub.drop_duplicates(subset='Uniprot')

        self.annot = pd.merge(self.annot, protein_atlas_sub, left_on='query_accession', right_on='Uniprot', how='left')
        self.annot.drop(columns='Uniprot', inplace=True)

    def exopeptidase(self):

        """
        Annotates dipeptidase and aminopeptidase activity by checking sequences for ragged patterns.
        
        It checks for sequences that end with a ragged pattern and have a similar sequence in the 
        dataframe with the last amino acid removed. Such sequences are annotated for aminopeptidase 
        and dipeptidase activity.
        """

        # keep track of sequences that have been checked
        cleared = set()
        # initialize column in the annotation dataframe
        self.annot["exopeptidase"] = np.nan
        # remove nan values and sort by length, so that longer sequences are checked first and cleared
        sequences = self.annot["query_sequence"].dropna()
        sequences.sort_values(key=lambda x: x.str.len(), kind="mergesort", ascending=False, inplace=True)

        exopeptidase_activity = False
        with tqdm(sequences, leave = 0) as t:
            for seq in t:
                if seq not in cleared:
                    # check if sequence ends with ragged pattern. If so, check if the same sequence
                    # with the last amino acid removed is also in the dataframe. If so, annotate.
                    rag_flag = False
                    cleared.add(seq)
                    char_match = seq[-5:]
                    matching_peptides = sequences[sequences.str.endswith(char_match)]

                    # if there are matches with 5 amino acids from the carboxy terminus, check if there is a ragged pattern
                    if len(matching_peptides) > 1:
                        compare = seq
                        lpep = len(seq)

                        # for each matching peptide, check if the same sequence with the last amino acid removed is also in the dataframe
                        for ind in matching_peptides.index:
                            pep = matching_peptides[ind]
                            same_seq_indices = self.annot[self.annot["query_sequence"] == pep].index
                            
                            # if there is a ragged pattern, annotate as such
                            if pep == compare[1:]:
                                cleared.add(pep)
                                compare = pep
                                logging.debug(f"1 {seq} {pep}")
                                exopeptidase_activity = True
                                # if the ragged pattern originated from dipetidase activity, annotate as such
                                if rag_flag and len(pep) == lpep - 1:
                                    for i in same_seq_indices:
                                        self.annot.loc[i, "exopeptidase"] = "Dipeptidase_seed_Aminopeptidase_activity"
                                # if the ragged pattern exists, annotate as aminopeptidase
                                else:
                                    for i in same_seq_indices:
                                        self.annot.loc[i, "exopeptidase"] = "Aminopeptidase_activity"
                                    rag_flag = False

                            elif pep == compare[2:]:
                                cleared.add(pep)
                                compare = pep
                                rag_flag = True
                                lpep = len(pep)
                                logging.debug(f"2 {seq} {pep}")
                                exopeptidase_activity = True
                                for i in same_seq_indices:
                                    self.annot.loc[i, "exopeptidase"] = "Dipeptidase_activity"
            elapsed = t.format_dict['elapsed']
            logging.info(f"Exopeptidase activity check took {annutils.format_seconds_to_time(elapsed)}")

        if exopeptidase_activity:
            logging.info("Exopeptidase activity was found, check logfile or output for more details on exact peptide sequences")

    def predict_protease_activity(self):
        """
        Predicts the protease activity for a given peptide sequence.

        This function reads the list of proteases from a file, constructs Position Specific Scoring
        Matrices (PSSMs) for each protease, and scores a given peptide sequence based on these PSSMs.
        The output is a string of protease code and scores.
        """

        self.annot["predicted_protease_activity"] = np.nan

        # Read the list of proteases
        protease_codes = pp.read_protease_file(self.proteasefile)
        logging.info(f"Protease codes: {protease_codes}")

        # Create PSSM for each protease and store in a dictionary
        pssms = pp.construct_pssms(protease_codes, self.merops, self.merops_sub, self.cleavagesitesize)
        
        with tqdm(range(len(self.annot)), leave = 0) as t:
            for i in t:
                # Get the peptide sequence
                cleavage = self.annot.loc[i, f"p{self.cleavagesitesize}_p{self.cleavagesitesize}prime"]
                # Score the peptide with each PSSM
                if isinstance(cleavage, str) and len(cleavage) == self.cleavagesitesize*2 and not '-' in cleavage:
                    scores = pp.score_proteases(pssms, cleavage)
                    self.annot.loc[i, "predicted_protease_activity"] = scores

            elapsed = t.format_dict['elapsed']
            logging.info(f"Protease activity prediction took {annutils.format_seconds_to_time(elapsed)}")

    def annotate_structure(self):

        """
        Annotates secondary structure and solvent accessibility for the cleavage site.
        
        This function checks if the secondary structure and solvent accessibility of the cleavage 
        site are available in the precalculated model data. If not, it calculates these properties 
        either for all sequences or for those sequences that are significantly different depending on 
        the 'calcstructure' attribute of the object. The significant difference is determined based 
        on a t-test or ANOVA.
        
        Args:
            alpha (float): The cutoff value for determining significant difference. Default is 0.05.
        """

        if self.available_models is None:
            available_models = annutils.read_alphafold_accessions(self.datafolder / self.alphafold_models_filename)
            self.available_models = available_models

        self.annot[f"secondary_structure p{self.cleavagesitesize}_p{self.cleavagesitesize}prime"] = np.nan
        self.annot[f"solvent_accessibility p{self.cleavagesitesize}_p{self.cleavagesitesize}prime"] = np.nan

        if self.calcstructure == "all":
            # Create a dictionary where each accession is a key and the value is a list of tuples
            # Each tuple contains the index of the cleavage site in the annot dataframe and the cleavage site position

            acc_cleavage_sites = {}
            for i in range(len(self.annot)):
                # get the accession and cleavage site
                acc = self.annot.loc[i, "query_accession"]
                cleavage_site = self.annot.loc[i, "p1_position"]

                # if the cleavage site is not nan, add it to the dictionary
                if float(cleavage_site).is_integer():
                    # Append the index and cleavage site to the list of the corresponding accession
                    acc_cleavage_sites.setdefault(acc, []).append((i, cleavage_site))

            # get the secondary structure and solvent accessibility of all cleavage sites

            structure_tmp_filepath = self.temp_folder / "structure_properties.txt"
            structure_properties = annutils.get_structure_properties(acc_cleavage_sites, structure_tmp_filepath, self.pymol_verbose, self.available_models)

            with open(structure_tmp_filepath, 'r') as f:
                lines = f.readlines()[0]
                structure_properties = literal_eval(lines)
            os.remove(structure_tmp_filepath)

            # assign to the annotation dataframe
            for (acc, cleavage_site), (index, ss, sa) in structure_properties.items():
                self.annot.loc[index, f"secondary_structure p{self.cleavagesitesize}_p{self.cleavagesitesize}prime"] = ss
                self.annot.loc[index, f"solvent_accessibility p{self.cleavagesitesize}_p{self.cleavagesitesize}prime"] = sa

        elif self.calcstructure == "sig":
            cols = []
            if len(self.conditions) > 1:
                if self.stat:
                    if self.pairwise or len(self.conditions) == 2:
                        conditions_iter = combinations(self.conditions.keys(), 2)
                        for pair in conditions_iter:
                            if self.multipletesting:
                                column_name = f"Corrected Independent T-test p-value: {pair[0]} vs. {pair[1]}"
                            else:
                                column_name = f"Independent T-test p-value: {pair[0]} vs. {pair[1]}"
                            cols.append(column_name)
                    else:
                        column_name = "ANOVA p-value: " + " vs. ".join(self.conditions.keys())
                        cols.append(column_name)

                    # iterate over all columns and entries in the dataframe
                    for column_name in cols:
                        subframe = self.annot[self.annot[column_name] <= self.alpha]
                        
                        acc_cleavage_sites = {}
                        for i, row in subframe.iterrows():  # use .iterrows() to keep track of the original index
                            # get the accession and cleavage site
                            acc = row["query_accession"]
                            cleavage_site = row["p1_position"]

                            # if the cleavage site is not nan, add it to the dictionary
                            if float(cleavage_site).is_integer():
                                # Append the original index and cleavage site to the list of the corresponding accession
                                acc_cleavage_sites.setdefault(acc, []).append((i, cleavage_site))

                        # get the secondary structure and solvent accessibility of all cleavage sites
                        structure_tmp_filepath = self.temp_folder / "structure_properties.txt"
                        structure_properties = annutils.get_structure_properties(acc_cleavage_sites, structure_tmp_filepath, self.pymol_verbose, self.available_models)

                        with open(structure_tmp_filepath, 'r') as f:
                            lines = f.readlines()[0]
                            structure_properties = literal_eval(lines)
                        os.remove(structure_tmp_filepath)

                        for (acc, cleavage_site), (index, ss, sa) in structure_properties.items():
                            # if the cleavage site is not nan (has not been annotated in previous iterations), annotate
                            if isinstance(cleavage_site, int) and pd.isnull(self.annot.loc[i, f"secondary_structure p{self.cleavagesitesize}_p{self.cleavagesitesize}prime"]) and pd.isnull(self.annot.loc[i, f"solvent_accessibility p{self.cleavagesitesize}_p{self.cleavagesitesize}prime"]):
                                self.annot.loc[index, f"secondary_structure p{self.cleavagesitesize}_p{self.cleavagesitesize}prime"] = ss
                                self.annot.loc[index, f"solvent_accessibility p{self.cleavagesitesize}_p{self.cleavagesitesize}prime"] = sa

                else:
                    logging.warning("cannot do structure calculations on significant peptides, as no statistics has been done. Please consider adding the -stat flag to the query!")
            else:
                logging.warning("There are less than 2 conditions provided. Peptide plotting on protein structure cannot be done for significant peptides as no significant peptides can be found!")

        else:
            raise ValueError(f"calcstructure argument must be either None, 'all', or 'sig'. It was '{self.calcstructure}' and of type {type(self.calcstructure)}")

    def visualize(self):

        """
        Calls the Visualizer class to create various plots, and stores these as figure objects.
        
        The Visualizer class is used to create a series of visualizations such as general statistics,
        CV plot, pie charts, heatmap, clustermap, PCA, and UMAP visualizations. These figures are 
        stored for later usage. Furthermore, if statistical and enrichment/pathway analyses are 
        carried out, additional figures are created. For multiple conditions, the method also generates
        volcano plots, fold change plots, galleries of significant peptides, and protein plots.
        """

        # initialize Visualizer class, and generate figures for general statistics, CV plot, pie charts, heatmap and clustermap
        vis = Visualizer(self.df, self.annot, self.conditions, self.software, self.patterns, self.temp_folder, self.pairwise, self.multipletesting)
        self.figures["General"] = vis.general()
        self.figures["CV"] = vis.cv_plot()
        self.figures["Piechart"] = vis.generate_pie_charts()
        self.figures["Heatmap"] = vis.heatmap()
        clustmap = vis.clustermap()
        if clustmap is not None:
            self.figures["Clustermap"] = clustmap
        self.figures["PCA"]  = vis.pca_visualization()
        self.figures["UMAP"]  = vis.umap_visualization()

        if self.stat and self.enrichment:
            self.figures["Enrichment"] = vis.plot_functional_enrichment(self.conditioncombinations, self.alpha)
        elif self.enrichment:
            logging.warning("Functional enrichment could not be performed as -stat was False (enrichment requires statistics)")

        if self.stat and self.pathway:
            if len(self.conditions) == 2 or self.pairwise:
                vis.plot_pathway_enrichment(self.conditioncombinations, self.alpha, folder=self.pathway_folder)
            else:
                logging.warning("Plotting of enriched pathways requires either 2 conditions or pairwise statistics to be performed. Please check that 2 conditions are provided, or if more than 2 conditions are provided that the flag -spw is given.")
        elif self.pathway:
            logging.warning("Enriched pathways could not be plotted as -stat was False (pathway enrichment requires statistics)")

        # if there are more than one condition, generate volcano, fold change and fold change at termini plots, and gallery of significant peptides
        if len(self.conditions) > 1:
            
            self.figures["Volcano"] = vis.volcano(self.volcano_foldchange, self.alpha)
            self.figures["Fold"] = vis.fold_plot()
            self.figures["Fold_nterm"] = vis.fold_termini()

            logging.info("Generating a PDF gallery of significant peptides...")
            vis.gallery(self.alpha, stat=self.stat, folder=self.general_folder)
            logging.info("Finished gallery generation.")
            
            if self.stat and self.cleavagevis and (self.pairwise or len(self.conditions) == 2):
                if self.available_models is None:
                    available_models = annutils.read_alphafold_accessions(self.datafolder / self.alphafold_models_filename)
                    self.available_models = available_models

                logging.info("Starting protein plotting...")
                if self.nomerops is False:
                    vis.plot_protein(self.alpha, pymol_verbose=self.pymol_verbose, folder=self.protein_folder, merops=self.merops, alphafold=self.available_models, level=self.cleavagevis)
                    logging.info("Finished protein plotting.")
                else:
                    vis.plot_protein(self.alpha, pymol_verbose=self.pymol_verbose, folder=self.protein_folder, alphafold=self.available_models, level=self.cleavagevis)
                    logging.info("Finished protein plotting.")


    def create_logos(self):

        """
        Creates sequence logos using the logomaker package.
        
        This method constructs sequence logos for each condition or pair of conditions based on
        statistical results and fold changes. These logos are created with the help of a helper 
        function that calculates Position Specific Scoring Matrices (PSSMs) using the logomaker 
        package, Bioconductor matrix calculation principles, and material from the DTU course 
        Algorithms in Bioinformatics. The logos are then stored and later written in figure files.
        """

        if len(self.conditions) > 1:
            try:
                if self.stat:
                    if len(self.conditions.keys()) == 2 or (self.pairwise and len(self.conditions.keys()) > 2):
                        # select columns to use based on whether to use multiple testing corrected values or not
                        if self.multipletesting:
                            column_log = "Corrected Independent T-test p-value:"
                        else:
                            column_log = "Independent T-test p-value:"

                        # create logo plots
                        for comparison in self.conditionpermutations:
                            column_name_test = f"{column_log} {comparison}"
                            column_name_fold = f"Log2 fold change: {comparison}"
                            column_test = self.annot[column_name_test]
                            column_fold = self.annot[column_name_fold]
                            data = self.annot[(column_test < self.alpha) & (column_fold > self.logo_fc)]
                            if len(data) == 0:
                                logging.info(f"No significant peptides found when doing logo plots for {comparison} comparison with a fold change > {self.logo_fc}, plot will not be made")
                            else:
                                self.figures[f"Logo {comparison} high"] = create_logo_helper(data, comparison, self.pseudocounts, self.logo, self.cleavagesitesize)

                            data = self.annot[(column_test < self.alpha) & (column_fold < -self.logo_fc)]
                            if len(data) == 0:
                                logging.info(f"No significant peptides found when doing logo plots for {comparison} comparison with a fold change < -{self.logo_fc}, plot will not be made")
                            else:
                                self.figures[f"Logo {comparison} low"] = create_logo_helper(data, comparison, self.pseudocounts, self.logo, self.cleavagesitesize)
                        logging.info("Created logo plots using peptides with statistically significant changed abundance between conditions")
                    
                    if len(self.conditions.keys()) > 2:
                        # select columns to use based on whether to use multiple testing corrected values or not
                        if self.multipletesting:
                            column_log = "Corrected ANOVA p-value:"
                        else:
                            column_log = f"ANOVA p-value:"
                        column_log = self.annot.columns[self.annot.columns.str.startswith(column_log)]

                        if len(column_log) == 1:
                            # create logo plots
                            column_name_test = column_log[0]
                            column_test = self.annot[column_name_test]
                            data = self.annot[(column_test < self.alpha)]
                            if len(data) == 0:
                                logging.info("No significant peptides found when doing logo plots for ANOVA comparison, plot will not be made")
                            else:
                                self.figures[f"Logo ANOVA only p-values"] = create_logo_helper(data, 'ANOVA', self.pseudocounts, self.logo, self.cleavagesitesize)
                                logging.info("Created logo plots for peptides with significant changes according to ANOVA alone")

                        else:
                            logging.info("Could not create logo plots with ANOVA p-values, as multiple similar columns were found")

                elif self.significance:
                    for comparison in self.conditionpermutations:
                        column = f"Fold {comparison} significance"
                        data = self.annot[self.annot[column] == "significant high"]
                        self.figures[f"Logo {comparison} high"] = create_logo_helper(data, comparison, self.pseudocounts, self.logo, self.cleavagesitesize)

                        data = self.annot[self.annot[column] == "significant low"]
                        self.figures[f"Logo {comparison} low"] = create_logo_helper(data, comparison, self.pseudocounts, self.logo, self.cleavagesitesize)

                    logging.info("Created logo plots using peptides with the highest and lowest percentile (deafault 5%) abundance fold change between conditions")

            except KeyError as err:
                logging.debug(f'ERROR in create_logos(): {err}')

        elif len(self.conditions) == 1:
            condition = list(self.conditions.keys())[0]
            self.figures[f"Logo {condition}"] = create_logo_helper(self.annot, condition, self.pseudocounts, self.logo, self.cleavagesitesize)
            logging.info("created logo plots with all provided peptides, as statistics is not possible with a single condition")

    def write_files(self):

        """
        Writes the output files to the specified location.
        
        This function exports the final data, which can be either the annotations alone or the merged 
        dataset of original data and annotations, to the desired format (csv, tsv, xlsx, json, or pkl). 
        It also saves any generated figures to the corresponding folders. At the end of the execution, 
        the logfile is copied to the output folder, and the output folder is compressed into a zip file.
        """
        import shutil

        outfile = self.outfolder / self.outname

        if self.separate:
            final_df = self.annot
        else:
            final_df = self.df.join(self.annot, rsuffix="_duplicatecol")
            final_df.drop([i for i in final_df.columns if '_duplicatecol' in i], axis=1, inplace=True)

        # Define a mapping of output types to saving methods
        saving_methods = {
            "csv": lambda df, path: df.to_csv(path, sep=",", index=False),
            "tsv": lambda df, path: df.to_csv(path, sep="\t", index=False),
            "xlsx": lambda df, path: df.to_excel(path, engine="openpyxl", index=False),
            "json": lambda df, path: df.to_json(path),
            "pkl": lambda df, path: df.to_pickle(path, compression="infer"),
        }

        # Save the dataframe using the appropriate method
        saving_methods[self.outfile_type](final_df, outfile)

        # Save figures
        if len(self.figures) > 0:
            annutils.save_figures(self.figures, self.folders)

        logging.info(f"Finished. Wrote results to {self.outfolder}.")

        log_basename = os.path.basename(self.logfile)
        shutil.copy(self.logfile, self.outfolder / log_basename)

        shutil.make_archive(self.outfolder, "zip", self.outfolder)

        # if temp_folder is there, delete even if it contains files
        if self.temp_folder.exists():
            shutil.rmtree(self.temp_folder)


        # if self.temp_folder is not None:
        # os.rmdir(self.temp_folder)

        return None