from multiprocessing import cpu_count

class Config:
    """Central configuration settings for the ddPrimer pipeline."""
    
    #############################################################################
    #                           Pipeline Mode Options
    #############################################################################
    ONLY_LASTZ_MULTIZ = False           # Run only LastZ/MultiZ alignments
    RUN_MATCH_CHECKER = False           # Run only the Match Checker
    RUN_MATCH_CHECKER_ON_OUTPUT = False # Run Match Checker on pipeline output
    DEBUG_MODE = False                   # Debug logging mode (enable with --debug flag)
    DISABLE_INTERNAL_OLIGO = False      # Disable internal oligo design
    
    #############################################################################
    #                           Performance Settings
    #############################################################################
    NUM_PROCESSES = max(1, int(cpu_count() * 0.75))  # Use 75% of cores
    BATCH_SIZE = 100
    MAF_CHUNK_SIZE = 10000
    SHOW_PROGRESS = True
    
    #############################################################################
    #                           Design Parameters
    #############################################################################
    # Basic primer design constraints
    PRIMER_MIN_SIZE = 18
    PRIMER_OPT_SIZE = 20
    PRIMER_MAX_SIZE = 23
    PRIMER_MIN_TM = 50.0
    PRIMER_OPT_TM = 57.5
    PRIMER_MAX_TM = 65.0
    PRIMER_MIN_GC = 50.0
    PRIMER_MAX_GC = 60.0
    
    # Product size constraints
    PRIMER_PRODUCT_SIZE_RANGE = [[90, 200]]
    
    # Pipeline parameters
    MIN_SEGMENT_LENGTH = 100
    RETAIN_TYPES = "gene"  # gff filtering: "gene", "mRNA", "CDS", "exon", etc.
    FILTER_MEANINGFUL_NAMES = True  # only use named genes from gff
    COUNT_AMBIGUOUS_AS_MISMATCH = False
    GENE_OVERLAP_MARGIN = 25
    RESTRICTION_SITE = "GGCC"
    PENALTY_MAX = 5.0
    MAX_PRIMER_PAIRS_PER_SEGMENT = 3
    PREFER_PROBE_MORE_C_THAN_G = True  # Set to False to disable
    
    #############################################################################
    #                           Validation Options
    #############################################################################
    VALIDATION_MODE = "TOLERANT"  # "STRICT" or "TOLERANT"
    ALLOW_AMP_MISMATCHES = 2  # Number of mismatches allowed in amplicon for TOLERANT mode
    ALLOW_AMP_MISMATCH_PERCENT = 5  # Percentage of mismatches allowed (0 = use absolute count)
    MAX_SEARCH_LENGTH = 1000000  # Limit search in large chromosomes to improve performance
    
    #############################################################################
    #                           BLAST Database Options
    #############################################################################
    DB_FASTA = None            # Path to a FASTA file to create a BLAST database from
    DB_OUTPUT_DIR = None       # Custom output directory for the BLAST database
    DB_NAME = None             # Custom name for the BLAST database
    USE_CUSTOM_DB = False      # Whether to use a custom database or the default
    DB_PATH = "/Library/Application Support/ddPrimer/Tair DB/TAIR10_db"
    
    # BLAST+ parameters
    BLAST_WORD_SIZE = 7
    BLAST_EVALUE = 10
    BLAST_MAX_TARGET_SEQS = 100
    BLAST_REWARD = 2
    BLAST_PENALTY = -3
    BLAST_GAPOPEN = 5
    BLAST_GAPEXTEND = 2
    BLAST_FILTER_FACTOR = 100  # E-value filtering factor
    
    # NUPACK parameters
    NUPACK_TEMPERATURE = 37  # Celsius
    NUPACK_SODIUM = 0.05     # Molar
    NUPACK_MAGNESIUM = 0.0   # Molar
    
    #############################################################################
    #                           Primer3 Settings
    #############################################################################
    # Complete Primer3 settings dictionary
    PRIMER3_SETTINGS = {
        # General settings
        "P3_FILE_TYPE": "settings", 
        "P3_FILE_ID": "User settings", 
        "P3P_DEBUG_MODE": 0, 
        "P3P_GB_ORIENTATION": "+",
        "P3P_PRIMER_NAME_ACRONYM_INTERNAL": "IN", 
        "P3P_PRIMER_NAME_ACRONYM_LEFT": "F", 
        "P3P_PRIMER_NAME_ACRONYM_RIGHT": "R",
        "P3P_PRIMER_NAME_ACRONYM_SPACER": "_", 
        
        # Primer conditions
        "PRIMER_ANNEALING_TEMP": 52.0, 
        "PRIMER_DMSO_CONC": 0.0,
        "PRIMER_DMSO_FACTOR": 0.6, 
        "PRIMER_DNA_CONC": 50.0, 
        "PRIMER_DNTP_CONC": 0.8, 
        "PRIMER_FIRST_BASE_INDEX": 1,
        "PRIMER_FORMAMIDE_CONC": 0.0, 
        "PRIMER_GC_CLAMP": 1, 
        "PRIMER_INSIDE_PENALTY": -1.0, 
        
        # Internal oligo parameters
        "PRIMER_INTERNAL_DMSO_CONC": 0.0,
        "PRIMER_INTERNAL_DMSO_FACTOR": 0.6, 
        "PRIMER_INTERNAL_DNA_CONC": 50.0, 
        "PRIMER_INTERNAL_DNTP_CONC": 0.0,
        "PRIMER_INTERNAL_FORMAMIDE_CONC": 0.0, 
        "PRIMER_INTERNAL_MAX_BOUND": 110.0, 
        "PRIMER_INTERNAL_MAX_GC": 80.0,
        "PRIMER_INTERNAL_MAX_HAIRPIN_TH": 47.0, 
        "PRIMER_INTERNAL_MAX_LIBRARY_MISHYB": 12.0, 
        "PRIMER_INTERNAL_MAX_NS_ACCEPTED": 0,
        "PRIMER_INTERNAL_MAX_POLY_X": 4, 
        "PRIMER_INTERNAL_MAX_SELF_ANY": 12.0, 
        "PRIMER_INTERNAL_MAX_SELF_ANY_TH": 47.0,
        "PRIMER_INTERNAL_MAX_SELF_END": 12.0, 
        "PRIMER_INTERNAL_MAX_SELF_END_TH": 47.0, 
        "PRIMER_INTERNAL_MAX_SIZE": 27,
        "PRIMER_INTERNAL_MAX_TM": 70, 
        "PRIMER_INTERNAL_MIN_3_PRIME_OVERLAP_OF_JUNCTION": 4,
        "PRIMER_INTERNAL_MIN_5_PRIME_OVERLAP_OF_JUNCTION": 7, 
        "PRIMER_INTERNAL_MIN_BOUND": -10.0, 
        "PRIMER_INTERNAL_MIN_GC": 30.0,
        "PRIMER_INTERNAL_MIN_QUALITY": 0, 
        "PRIMER_INTERNAL_MIN_SIZE": 15, 
        "PRIMER_INTERNAL_MIN_THREE_PRIME_DISTANCE": -1,
        "PRIMER_INTERNAL_MIN_TM": 64, 
        "PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME": "hnnnn",
        "PRIMER_INTERNAL_OPT_BOUND": 97.0,
        "PRIMER_INTERNAL_OPT_GC_PERCENT": 50.0, 
        "PRIMER_INTERNAL_OPT_SIZE": 20, 
        "PRIMER_INTERNAL_OPT_TM": 65,
        "PRIMER_INTERNAL_SALT_DIVALENT": 0.0, 
        "PRIMER_INTERNAL_SALT_MONOVALENT": 50.0, 
        
        # General primer constraints
        "PRIMER_LIBERAL_BASE": 1,
        "PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS": 0, 
        "PRIMER_LOWERCASE_MASKING": 0, 
        "PRIMER_MAX_BOUND": 110.0,
        "PRIMER_MAX_END_GC": 3, 
        "PRIMER_MAX_END_STABILITY": 9.0, 
        "PRIMER_MAX_GC": 60.0, 
        "PRIMER_MAX_HAIRPIN_TH": 47.0,
        "PRIMER_MAX_LIBRARY_MISPRIMING": 12.0, 
        "PRIMER_MAX_NS_ACCEPTED": 0, 
        "PRIMER_MAX_POLY_X": 4, 
        "PRIMER_MAX_SELF_ANY": 8.0,
        "PRIMER_MAX_SELF_ANY_TH": 47.0, 
        "PRIMER_MAX_SELF_END": 3.0, 
        "PRIMER_MAX_SELF_END_TH": 47.0, 
        "PRIMER_MAX_SIZE": 23,
        "PRIMER_MAX_TEMPLATE_MISPRIMING": 12.0, 
        "PRIMER_MAX_TEMPLATE_MISPRIMING_TH": 47.0, 
        "PRIMER_MAX_TM": 65.0,
        "PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION": 4, 
        "PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION": 7, 
        "PRIMER_MIN_BOUND": -10.0,
        "PRIMER_MIN_END_QUALITY": 0, 
        "PRIMER_MIN_GC": 50.0, 
        "PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE": 3,
        "PRIMER_MIN_QUALITY": 0, 
        "PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE": 3, 
        "PRIMER_MIN_SIZE": 18, 
        "PRIMER_MIN_TM": 50.0,
        "PRIMER_NUM_RETURN": 10, 
        "PRIMER_OPT_BOUND": 97.0, 
        "PRIMER_OPT_GC_PERCENT": 52.5, 
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_OPT_TM": 57.5, 
        "PRIMER_OUTSIDE_PENALTY": 0.0, 
        
        # Primer pair parameters
        "PRIMER_PAIR_MAX_COMPL_ANY": 8.0,
        "PRIMER_PAIR_MAX_COMPL_ANY_TH": 47.0, 
        "PRIMER_PAIR_MAX_COMPL_END": 3.0, 
        "PRIMER_PAIR_MAX_COMPL_END_TH": 47.0,
        "PRIMER_PAIR_MAX_DIFF_TM": 1, 
        "PRIMER_PAIR_MAX_LIBRARY_MISPRIMING": 24.0, 
        "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING": 24.0,
        "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH": 47.0, 
        "PRIMER_PAIR_WT_COMPL_ANY": 0.0, 
        "PRIMER_PAIR_WT_COMPL_ANY_TH": 0.0,
        "PRIMER_PAIR_WT_COMPL_END": 0.0, 
        "PRIMER_PAIR_WT_COMPL_END_TH": 0.0, 
        "PRIMER_PAIR_WT_DIFF_TM": 0.0,
        "PRIMER_PAIR_WT_IO_PENALTY": 0.0, 
        "PRIMER_PAIR_WT_LIBRARY_MISPRIMING": 0.0, 
        "PRIMER_PAIR_WT_PRODUCT_SIZE_GT": 0.0,
        "PRIMER_PAIR_WT_PRODUCT_SIZE_LT": 0.0, 
        "PRIMER_PAIR_WT_PRODUCT_TM_GT": 0.0, 
        "PRIMER_PAIR_WT_PRODUCT_TM_LT": 0.0,
        "PRIMER_PAIR_WT_PR_PENALTY": 1.0, 
        "PRIMER_PAIR_WT_TEMPLATE_MISPRIMING": 0.0, 
        "PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH": 0.0,
        
        # Primer selection and product parameters
        "PRIMER_PICK_ANYWAY": 0, 
        "PRIMER_PICK_INTERNAL_OLIGO": 1, 
        "PRIMER_PICK_LEFT_PRIMER": 1, 
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_PRODUCT_MAX_TM": 1000000.0, 
        "PRIMER_PRODUCT_MIN_TM": -1000000.0, 
        "PRIMER_PRODUCT_OPT_SIZE": 0,
        "PRIMER_PRODUCT_OPT_TM": 0.0, 
        "PRIMER_PRODUCT_SIZE_RANGE": "90-200", 
        "PRIMER_QUALITY_RANGE_MAX": 100,
        "PRIMER_QUALITY_RANGE_MIN": 0, 
        
        # Chemistry and thermodynamic parameters
        "PRIMER_SALT_CORRECTIONS": 1, 
        "PRIMER_SALT_DIVALENT": 3.8, 
        "PRIMER_SALT_MONOVALENT": 50.0,
        "PRIMER_SECONDARY_STRUCTURE_ALIGNMENT": 1, 
        "PRIMER_SEQUENCING_ACCURACY": 20, 
        "PRIMER_SEQUENCING_INTERVAL": 250,
        "PRIMER_SEQUENCING_LEAD": 50, 
        "PRIMER_SEQUENCING_SPACING": 500, 
        "PRIMER_TASK": "generic",
        "PRIMER_THERMODYNAMIC_PARAMETERS_PATH": "/opt/homebrew/Cellar/primer3/2.4.0/share/primer3/primer3_config/",
        "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1, 
        "PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT": 0, 
        "PRIMER_TM_FORMULA": 1,
        
        # Penalty weights
        "PRIMER_WT_BOUND_GT": 0.0, 
        "PRIMER_WT_BOUND_LT": 0.0, 
        "PRIMER_WT_END_QUAL": 0.0, 
        "PRIMER_WT_END_STABILITY": 0.0,
        "PRIMER_WT_GC_PERCENT_GT": 0.5, 
        "PRIMER_WT_GC_PERCENT_LT": 0.5, 
        "PRIMER_WT_HAIRPIN_TH": 0.0,
        "PRIMER_WT_LIBRARY_MISPRIMING": 0.0, 
        "PRIMER_WT_NUM_NS": 0.0, 
        "PRIMER_WT_POS_PENALTY": 0.0, 
        "PRIMER_WT_SELF_ANY": 0.0,
        "PRIMER_WT_SELF_ANY_TH": 0.0, 
        "PRIMER_WT_SELF_END": 0.0, 
        "PRIMER_WT_SELF_END_TH": 0.0, 
        "PRIMER_WT_SEQ_QUAL": 0.0,
        "PRIMER_WT_SIZE_GT": 1.0, 
        "PRIMER_WT_SIZE_LT": 1.0, 
        "PRIMER_WT_TEMPLATE_MISPRIMING": 0.0,
        "PRIMER_WT_TEMPLATE_MISPRIMING_TH": 0.0, 
        "PRIMER_WT_TM_GT": 1.0, 
        "PRIMER_WT_TM_LT": 1.0
    }
    
    @classmethod
    def get_primer3_global_args(cls):
        """
        Get global primer3 arguments as a dictionary.
        This combines the simplified settings with the complete settings.
        
        Returns:
            dict: Dictionary of primer3 settings
        """
        # Start with complete settings
        settings = cls.PRIMER3_SETTINGS.copy()
        
        # Override with simplified settings for commonly adjusted parameters
        settings.update({
            "PRIMER_MIN_SIZE": cls.PRIMER_MIN_SIZE,
            "PRIMER_OPT_SIZE": cls.PRIMER_OPT_SIZE,
            "PRIMER_MAX_SIZE": cls.PRIMER_MAX_SIZE,
            "PRIMER_MIN_TM": cls.PRIMER_MIN_TM,
            "PRIMER_OPT_TM": cls.PRIMER_OPT_TM,
            "PRIMER_MAX_TM": cls.PRIMER_MAX_TM,
            "PRIMER_MIN_GC": cls.PRIMER_MIN_GC,
            "PRIMER_MAX_GC": cls.PRIMER_MAX_GC,
            "PRIMER_NUM_RETURN": cls.MAX_PRIMER_PAIRS_PER_SEGMENT,
        })
        
        # Convert the product size range from list to string format if needed
        if isinstance(cls.PRIMER_PRODUCT_SIZE_RANGE, list):
            # Convert [[min1, max1], [min2, max2], ...] to "min1-max1 min2-max2 ..."
            size_range_str = " ".join([f"{r[0]}-{r[1]}" for r in cls.PRIMER_PRODUCT_SIZE_RANGE])
            settings["PRIMER_PRODUCT_SIZE_RANGE"] = size_range_str
        
        return settings
    
    @classmethod
    def format_settings_for_file(cls):
        """
        Format settings for writing to a primer3 settings file.
        
        Returns:
            str: Formatted settings string
        """
        settings = cls.get_primer3_global_args()
        return "\n".join(f"{key}={value}" for key, value in settings.items()) + "\n"
    
    @classmethod
    def load_from_file(cls, filepath):
        """
        Load settings from a primer3 settings file.
        
        Args:
            filepath (str): Path to the settings file
            
        Returns:
            bool: True if settings were loaded successfully
        """
        try:
            with open(filepath, 'r') as f:
                settings_text = f.read()
                
            # Parse settings
            settings = {}
            for line in settings_text.strip().split('\n'):
                if '=' in line and not line.startswith('#'):
                    key, value = line.split('=', 1)
                    
                    # Try to convert to appropriate type
                    if value.replace('.', '', 1).isdigit():
                        if '.' in value:
                            value = float(value)
                        else:
                            value = int(value)
                    elif value.lower() in ('true', 'false'):
                        value = value.lower() == 'true'
                        
                    settings[key.strip()] = value
            
            # Update settings
            cls.PRIMER3_SETTINGS.update(settings)
            
            # Update simplified settings if corresponding keys exist
            if "PRIMER_MIN_SIZE" in settings:
                cls.PRIMER_MIN_SIZE = settings["PRIMER_MIN_SIZE"]
            if "PRIMER_OPT_SIZE" in settings:
                cls.PRIMER_OPT_SIZE = settings["PRIMER_OPT_SIZE"]
            if "PRIMER_MAX_SIZE" in settings:
                cls.PRIMER_MAX_SIZE = settings["PRIMER_MAX_SIZE"]
            if "PRIMER_MIN_TM" in settings:
                cls.PRIMER_MIN_TM = settings["PRIMER_MIN_TM"]
            if "PRIMER_OPT_TM" in settings:
                cls.PRIMER_OPT_TM = settings["PRIMER_OPT_TM"]
            if "PRIMER_MAX_TM" in settings:
                cls.PRIMER_MAX_TM = settings["PRIMER_MAX_TM"]
            if "PRIMER_MIN_GC" in settings:
                cls.PRIMER_MIN_GC = settings["PRIMER_MIN_GC"]
            if "PRIMER_MAX_GC" in settings:
                cls.PRIMER_MAX_GC = settings["PRIMER_MAX_GC"]
            if "PRIMER_NUM_RETURN" in settings:
                cls.MAX_PRIMER_PAIRS_PER_SEGMENT = settings["PRIMER_NUM_RETURN"]
            
            # Handle product size range
            if "PRIMER_PRODUCT_SIZE_RANGE" in settings:
                value = settings["PRIMER_PRODUCT_SIZE_RANGE"]
                if isinstance(value, str):
                    # Parse "min1-max1 min2-max2 ..." format
                    ranges = []
                    for r in value.split():
                        if '-' in r:
                            min_val, max_val = r.split('-')
                            ranges.append([int(min_val), int(max_val)])
                    cls.PRIMER_PRODUCT_SIZE_RANGE = ranges
            
            # Handle BLAST database options
            if "DB_FASTA" in settings:
                cls.DB_FASTA = settings["DB_FASTA"]
                cls.USE_CUSTOM_DB = True
            
            if "DB_OUTPUT_DIR" in settings:
                cls.DB_OUTPUT_DIR = settings["DB_OUTPUT_DIR"]
                
            if "DB_NAME" in settings:
                cls.DB_NAME = settings["DB_NAME"]
                
            if "DB_PATH" in settings:
                cls.DB_PATH = settings["DB_PATH"]
            
            return True
        except Exception as e:
            print(f"Error loading settings from {filepath}: {e}")
            return False
    
    @classmethod
    def save_to_file(cls, filepath):
        """
        Save current settings to a file.
        
        Args:
            filepath (str): Path to save the settings
            
        Returns:
            bool: True if settings were saved successfully
        """
        try:
            with open(filepath, 'w') as f:
                f.write(cls.format_settings_for_file())
            return True
        except Exception as e:
            print(f"Error saving settings to {filepath}: {e}")
            return False
    
    @staticmethod
    def debug(message):
        """
        Print debug messages if debug mode is enabled.
        
        Args:
            message (str): The debug message to print
        """
        if Config.DEBUG_MODE:
            print(f"[DEBUG] {message}")
