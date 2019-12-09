# Input Files
IDENTIFICATION_MZML_FILE = "C:\\Users\\Tom\\PycharmProjects\\TheDeuteriumCalculator\\mzMLFiles\\IDRun.mzML"  # must be .mzML
IDENTIFICATION_CSV_FILE = "C:\\Users\\Tom\\PycharmProjects\\DeuteriumCalculator\\Input.csv"  # Must be .csv
PROTEIN_SEQUENCE_FILE = "C:\\Users\\Tom\\PycharmProjects\\TheDeuteriumCalculator\\Protein.txt"  # Must be .txt

# Output Files
FULL_HDX_OUTPUT = "HdxOutput"  # no file extension
SUMMARY_HDX_OUTPUT = "OutputTable.csv"  # Must be .csv
SUMMARY_TABLE_1 = "Summary_Table_1.csv"  # Must be .csv
SUMMARY_TABLE_2 = "Summary_Table_2.csv"  # Must be .csv

# Output Figures
WOODS_PLOT_NAME = "WoodsPlot"  # no file extension

# Parameters (Can be changed, defaults given in parenthesis)
NOISE_LIMIT = 10000  # All individual peaks with less than this intensity ignored (10000)
PPM_MATCH_TOLERANCE = 10  # Peaks with less difference than this value matched to sequence (10)
SLIDING_WINDOW_PPM_TOLERANCE = 1  # Peaks with less difference than this value combined within each sliding window (1)
SLIDING_WINDOW_SIZE = 30  # width of sliding window in seconds, should be integer divisible by SLIDE_FRACTION (60)
SLIDE_FRACTION = 3  # Fraction of the window that the window moves each each slide (3)
RETENTION_TOLERANCE = 30  # window of retention times to search for given peptide (+-) (30)
WOODS_PLOT_CONFIDENCE = 0.99  # Use to calculate confidence interval for differential woods plot (0-1)

# Constants (Do not change)
DEUTERIUM_MASS = 1.00627
MINUTES_TO_SECONDS = 60
SLIDE_AMOUNT = SLIDING_WINDOW_SIZE / SLIDE_FRACTION
PEPTIDE_MASS_DICTIONARY = {
    'A': 71.0779,
    'C': 103.1429,
    'D': 115.0874,
    'E': 129.114,
    'F': 147.1739,
    'G': 57.0513,
    'H': 137.1393,
    'I': 113.1576,
    'K': 128.1723,
    'L': 113.1576,
    'M': 131.1961,
    'N': 114.1026,
    'P': 97.1152,
    'Q': 128.1292,
    'R': 156.1857,
    'S': 87.0773,
    'T': 101.1039,
    'U': 150.0379,
    'W': 186.2099,
    'Y': 163.1733,
    'V': 99.1311
}
MASS_OF_WATER = 18.01528
MASS_OF_HYDROGEN = 1.007276
