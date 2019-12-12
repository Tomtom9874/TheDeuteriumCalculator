import numpy as np
import pandas as pd
import csv
from scipy.optimize import curve_fit
from scipy import stats
import PARAMETERS as CON
from pyteomics import mzml
from datetime import datetime
from os import path
from matplotlib import pyplot as plt


# fit_gaussian helper function
def gaussian_trend(x, a, sigma, mean):
    return a * np.exp(-(x - mean) ** 2 / (2 * sigma ** 2))


# fit_gaussian helper function
def calculate_r_squared(x_data, y_data, *p):
    y_predicted = gaussian_trend(x_data, *p)
    residuals = (y_predicted - y_data) ** 2
    ss_res = sum(residuals)
    ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    r_squared = round(r_squared, 3)
    return r_squared


# fits a Guassian curve to ordered list parameter, returns r^2
def fit_gaussian(y_data):
    x_data = list(range(len(y_data)))
    try:
        popt, pcov = curve_fit(gaussian_trend, x_data, y_data)
        r_squared = calculate_r_squared(x_data, y_data, *popt)
        return r_squared
    except RuntimeError:
        return "N/A"


# returns the maximum possible deuterium for a given sequence
def sequence_to_max_deuterium(sequence: str):
    max_deuterium = len(sequence) - 1
    for letter in sequence:
        if letter.lower() == 'p':
            max_deuterium -= 1
    return max_deuterium


# Checks that user PARAMETER configuration is (relatively) correct
def check_parameters():
    if not check_extension(CON.IDENTIFICATION_MZML_FILE, "mzml"):
        print("ERROR: IDENTIFICATION_MZML_FILE must be a .mzml")
    if not check_extension(CON.IDENTIFICATION_CSV_FILE, "csv"):
        print("ERROR: IDENTIFICATION_CSV_FILE must be a .csv")
    if not check_extension(CON.PROTEIN_SEQUENCE_FILE, "txt"):
        print("ERROR: PROTEIN_SEQUENCE_FILE must be a .txt")
    for letter in CON.FULL_HDX_OUTPUT:
        if letter == '.':
            print("ERROR: FULL_HDX_OUTPUT should not have file extension")
    if not check_extension(CON.RECOMMENDATION_TABLE_1, "csv"):
        print("ERROR: SUMMARY_TABLE_1 must be a .csv")
    if CON.PPM_MATCH_TOLERANCE > 100:
        print("WARNING: PPM_MATCH_TOLERANCE is set to {}, "
              "this is high".format(CON.PPM_MATCH_TOLERANCE))
    if CON.SLIDING_WINDOW_PPM_TOLERANCE > 10:
        print("WARNING: SLIDING_WINDOW_PPM_TOLERANCE is set to {}, "
              "recommendation is under 10.".format(CON.SLIDING_WINDOW_PPM_TOLERANCE))
    if CON.SLIDING_WINDOW_SIZE % CON.SLIDE_AMOUNT != 0:
        print("ERROR: SLIDING_WINDOW_SIZE must be divisible by SLIDE_AMOUNT")
    if CON.RETENTION_TOLERANCE <= 0:
        print("ERROR: RETENTION_TOLERANCE must be a positive number, "
              "greater than 30 recommended.")
    if CON.RETENTION_TOLERANCE < 30:
        print("WARNING: A retention tolerance of greater than 30 is recommended")
    if not 0 < CON.WOODS_PLOT_CONFIDENCE < 1:
        print("confidence should be between 0 and 1")
        raise ValueError


# helper function for check_parameters()
def check_extension(string, extension):
    if string[-len(extension):].lower() != extension:
        return False
    else:
        return True


# Gets user input of path and changes it to usable string
def get_path_input():
    input_path = r"{}".format(input())
    compatible_path = ""
    for letter in input_path:
        if letter == '\\':
            compatible_path += '\\'
        if letter != '"':
            compatible_path += letter
    return compatible_path


# Determines the location of the sequence within the full protein
def find_start_end(peptide: str, protein: str):
    sequential_matches = 0
    start, end = 0, 0
    index = 0
    while index < len(protein):
        if peptide[sequential_matches] == protein[index]:
            sequential_matches += 1
            if sequential_matches == len(peptide):
                end = index + 1
                return start, end
            elif sequential_matches == 1:
                start = index + 1
            index += 1
        else:
            if sequential_matches != 0:
                index = start
            else:
                index += 1
            sequential_matches = 0
    return "NULL", "NULL"


# Removes any non-alpha characters from the protein sequence
def parse_protein(file: str):
    sequence = ""
    with open(file, 'r') as f:
        file_reader = f.readlines()
        for line in file_reader:
            for character in line:
                if character.isalpha():
                    sequence += character
    return sequence


# Returns the ppm difference between two m/z values
def get_ppm(mz1, mz2):
    return abs(mz1 - mz2) / mz1 * 1000000


# Takes in a list of tuples and combines each where first elements is within a ppm tolerance
def tuple_combine(some_list):
    changed = True
    just_changed = False
    start_list = sorted(some_list, key=lambda x: x[0])
    new_list = []
    while changed and len(start_list) != 0:
        changed = False
        for i in range(len(start_list) - 1):
            if just_changed:
                just_changed = False
                continue
            ppm = abs(get_ppm(start_list[i][0], start_list[i+1][0]))
            if ppm < CON.SLIDING_WINDOW_PPM_TOLERANCE:
                new_list.append((((start_list[i][0] + start_list[i+1][0]) / 2),
                                 start_list[i][1] + start_list[i+1][1]))
                changed = True
                just_changed = True
            else:
                new_list.append(start_list[i])
        if not just_changed:
            new_list.append(start_list[len(start_list) - 1])
        start_list = new_list
        new_list = []
        just_changed = False
    return start_list


# Modified binary search which returns the highest-intensity peak within a ppm tolerance
def compare(target, charge, array, full_array):
    midpoint = int(len(array) / 2)
    if abs(get_ppm(target, array[midpoint][0] * charge)) <= CON.PPM_MATCH_TOLERANCE:
        return_list = [(array[midpoint][0], array[midpoint][1])]
        offset = 1
        while offset != 0 and (midpoint - offset) > 0:
            ppm = abs(get_ppm(target, array[midpoint + offset][0] * charge))
            if ppm <= CON.PPM_MATCH_TOLERANCE:
                return_list.append((array[midpoint - offset][0], array[midpoint - offset][1]))
                offset += 1
            else:
                offset = 0
        offset = 1
        while offset != 0 and (midpoint + offset < len(array)):
            ppm = abs(get_ppm(target, array[midpoint + offset][0] * charge))
            if ppm <= CON.PPM_MATCH_TOLERANCE:
                return_list.append((array[midpoint + offset][0], array[midpoint + offset][1]))
                offset += 1
            else:
                offset = 0
        high_intensity = (0, 0)
        for key, value in return_list:
            if value > high_intensity[1]:
                high_intensity = key, value

        ppm_error = get_ppm(target, high_intensity[0] * charge)
        return ppm_error, high_intensity[0],  high_intensity[1]

    elif len(array) == 1 or len(array) == 0:
        return 0, 0, 0
    elif array[midpoint][0] * charge <= target:
        return compare(target, charge, array[midpoint:], full_array)
    else:
        return compare(target, charge, array[0: midpoint], full_array)


# Converts scan number to retention time using the mzml file
def set_retention_times(file: str):
    retention_scan_dictionary = {}
    with mzml.read(file) as f:
        for scan in f:
            if scan["ms level"] == 2:
                scan_time = float(scan["scanList"]["scan"][0]["scan start time"])
                retention_scan_dictionary[scan["index"] + 1] = scan_time * CON.MINUTES_TO_SECONDS
    return retention_scan_dictionary


#################################################################################################
class FullExperiment:

    def __init__(self, time_points: list, differential: bool, replications: int):
        self.runs = {}
        self._is_differential = differential
        self._num_replicates = replications
        self._time_points = time_points
        self._peptides = []
        for time in self._time_points:
            self.runs[time] = {True: {}, False: {}}
        self._file_names = []
        self.deviations_by_time = {}
        self.fractional_deviations_by_time = {}
        for time in self._time_points:
            self.deviations_by_time[time] = []
            self.fractional_deviations_by_time[time] = []
        self.protein = parse_protein(CON.PROTEIN_SEQUENCE_FILE)
        self.add_file_names()

    def add_runs(self, time):
        count = 0
        for replication in range(self._num_replicates):
            self.add_run(time, False, replication, count)
            count += 1
        if self._is_differential:
            for replication in range(self._num_replicates):
                self.add_run(time, True, replication, count)
                count += 1

    def add_file_names(self):
        for time in self._time_points:
            complexity = False
            for replication in range(self._num_replicates):
                self.add_file(time, complexity, replication)
            if self._is_differential:
                complexity = True
                for replication in range(self._num_replicates):
                    self.add_file(time, complexity, replication)

    def get_file(self, index: int):
        return self._file_names[index]

    def add_deviation(self, time, deviation, pep):
        if pep.get_mass_shift() > 0:
            self.deviations_by_time[time].append(deviation)
            self.fractional_deviations_by_time[time].append(deviation / pep.get_max_deuterium())

    # adds the name of each file to the list
    def add_file(self, time, is_complex: bool, replication):
        complexity = "Free"
        if is_complex:
            complexity = "Complex"
        print("Enter path to mzML for time:", time, "replication:",
              replication + 1, "complexity:", complexity)
        file = ""
        while not path.exists(file):
            file = get_path_input()
            if not path.exists(file):
                print("Not a valid path name.")
        self._file_names.append(file)

    # Adds a new experimental run, adds the peptides on the first run
    def add_run(self, time, complexity: bool, replication, index: int):

        run = ExperimentalRun(run=replication, complexity=complexity, time=time)
        file = self._file_names[index]
        run.read_mzml(file)
        run.hydrogen_deuterium_exchange(CON.IDENTIFICATION_MZML_FILE)
        uptakes = []
        self._peptides = []
        for pep in run.get_peptides():
            self._peptides.append(pep)
            uptakes.append(pep.get_mass_shift())
        self.runs[time][complexity][replication] = uptakes

    def update_headers(self, top_header, bottom_header, header_name):
        top_header.append(header_name)
        for time in self._time_points:
            bottom_header.append(str(time) + " min")
            top_header.append("")
        top_header.pop()
    # Contains multiple subroutines that generate each type of output

    def average_replications(self, averages, deviations, pep, index, complexity):
        # TODO may need to be rearranged to support multiple time points
        for time in self._time_points:
            replications = []
            for replication in range(self._num_replicates):
                replications.append(self.runs[time][complexity][replication][index])
            average = sum(replications) / self._num_replicates
            averages.append(average)
            deviation = np.std(replications, ddof=1)
            self.add_deviation(time, deviation, pep)
            deviations.append(deviation)

    def generate_recommendation_table_1(self):
        top_header = ["Start", "End", "Sequence", "Peptide mass (Da)", "Retention time (min)"]
        bot_header = ["Start", "End", "Sequence", "Peptide mass (Da)", "Retention time (min)"]
        self.update_headers(top_header, bot_header, "Uptake FREE (D)")
        if self._is_differential:
            self.update_headers(top_header, bot_header, "Uptake COMPLEX (D)")
        self.update_headers(top_header, bot_header, "Uptake error (SD) - FREE (D)")
        if self._is_differential:
            self.update_headers(top_header, bot_header, "Uptake error (SD) - COMPLEX (D)")
        peptide_lines = []
        # Each loop generates one line
        for index, pep in enumerate(self._peptides):
            peptide = [pep.get_start(), pep.get_end(), pep.get_sequence(), pep.get_average_mass()]
            start, end = pep.get_rt_start_end()
            rt = (start + end) / 2 / CON.MINUTES_TO_SECONDS
            peptide.append(rt)
            averages = []
            deviations = []
            self.average_replications(averages, deviations, pep, index, False)
            if self._is_differential:
                self.average_replications(averages, deviations, pep, index, True)
            peptide.extend(averages)
            peptide.extend(deviations)
            peptide_lines.append(peptide)

        with open(CON.RECOMMENDATION_TABLE_1, "w+", newline='') as f:
            csv_writer = csv.writer(f)
            csv_writer.writerow(top_header)
            csv_writer.writerow(bot_header)
            for line in peptide_lines:
                csv_writer.writerow(line)
        print("Wrote to:", CON.RECOMMENDATION_TABLE_1)

    def generate_output(self):
        print("\nGenerating Output file(s)")
        self.generate_recommendation_table_1()
        # TODO Add Second Output
        if self._is_differential:
            self.generate_differential_woods_plot(CON.RECOMMENDATION_TABLE_1, CON.WOODS_PLOT_TITLE)

    # confidence is between 0 and 1, df is
    def calculate_confidence_limit(self, time, fractional):
        df = (self._num_replicates * 2) - 1
        num_replications = self._num_replicates * 2
        if fractional:
            deviations = self.fractional_deviations_by_time[time]
        else:
            deviations = self.deviations_by_time[time]
        if not deviations:
            raise ValueError("add deviations before calculating confidence limit.")
        stdev = sum(deviations) / len(deviations)
        alpha = 1 - CON.WOODS_PLOT_CONFIDENCE
        critical_value = stats.t.ppf(1 - (alpha / 2), df)
        print("Stdev:", stdev, "replications:", num_replications, "critical v.:", critical_value)
        return stdev / (num_replications ** 0.5) * critical_value

    # Saves a plot of name "file"_"time"s.png with a given title. can be fractional or absolute
    def generate_differential_woods_plot(self, file: str, title: str, fractional=True):
        plt.title(title)
        plt.xlabel("Sequence")
        plt.ylabel("Relative Uptake (Da)")
        for time in self._time_points:
            if fractional:
                plt.ylabel("Relative Fractional Uptake")
            confidence = self.calculate_confidence_limit(time, fractional)
            print("Confidence:", confidence)
            plt.plot((0, len(self.protein)), (0, 0), 'k:')
            plt.plot((0, len(self.protein)), (confidence, confidence), 'r--')
            plt.plot((0, len(self.protein)), (-confidence, -confidence), 'r--')
            df = pd.read_csv(file, header=[0, 1])
            time_col = str(time) + " min"
            for _, row in df.iterrows():
                x = (row["Start"]["Start"], row["End"]["End"])
                difference = row["Uptake COMPLEX (D)"][time_col] - row["Uptake FREE (D)"][time_col]
                if fractional:
                    difference /= sequence_to_max_deuterium(row['Sequence']['Sequence'])
                y = (difference, difference)
                plt.plot(x, y, 'k')
            plot_file_name = CON.WOODS_PLOT_NAME + "_" + str(time) + "s.png"
            plt.savefig(plot_file_name, dpi=600)

    def coverage_calculator(self):
        pass


##########################################################################
class ExperimentalRun:

    def __init__(self, run, complexity, time):
        self.all_peaks = []
        self.peptides = []
        self.windows = {}
        self._run = run
        self._complexity = complexity
        self._time = time

    # Getters
    def get_peptides(self):
        return self.peptides

    def get_run(self):
        return self._run

    def get_complexity(self):
        return self._complexity

    def get_time(self):
        return self._time

    def get_tuple_dictionary(self):
        return self.windows

    def get_retention_time(self, index: int):
        return self.all_peaks[index]["retention time"]

    # Setters
    def set_tuple_dictionary(self, tuple_dict: dict):
        self.windows = tuple_dict

    # converts peptide scan number to retention times
    def set_pep_retention_times(self, file: str):
        conversion_dictionary = set_retention_times(file)
        for pep in self.peptides:
            scan = pep.get_scan()
            rt = conversion_dictionary[scan]
            pep.set_retention_times(rt - CON.RETENTION_TOLERANCE, rt + CON.RETENTION_TOLERANCE)

    # adds peptides to peptide_list
    def read_input(self, file: str):
        with open(file, 'r') as f:
            csv_reader = csv.DictReader(f)
            for row in csv_reader:
                self.add_peptide(row['Peptide'], float(row['Precursor']),
                                 int(row['Charge']), int(row["ID"]), float(row["ScanNum"]))

    def process_scan(self, scan):
        retention_time = scan["scanList"]["scan"][0]["scan start time"] * CON.MINUTES_TO_SECONDS
        tuple_list = []
        for j in range(len(scan["m/z array"])):
            if scan["intensity array"][j] >= CON.NOISE_LIMIT:
                tuple_list.append((scan["m/z array"][j], scan["intensity array"][j]))
        self.all_peaks.append({"retention time": retention_time, "tuple list": tuple_list})

    # creates a dictionary for each scan with the RT and a
    # list with tuples containing the m/z and intensity
    def read_mzml(self, file: str):
        total = 0
        print("(initializing)")
        with mzml.read(file) as f:
            for scan in f:
                if scan["ms level"] == 1:
                    total += 1
        count = 0
        retention_time = None
        with mzml.read(file) as f:
            for scan in f:
                if scan["ms level"] == 1:
                    retention_time = scan["scanList"]["scan"][0]["scan start time"]
                    retention_time *= CON.MINUTES_TO_SECONDS
                    count += 1
                    if count % 200 == 0 or count == 1 or count == total:
                        print(count, "/", total)
                    self.process_scan(scan)
        self.set_tuple_dictionary(self.sliding_window(retention_time))

    # averages intensity of similar masses within RT ranges
    def sliding_window(self, retention_time: float):
        window_count = int(int((retention_time // CON.SLIDING_WINDOW_SIZE) + 1) *
                           (CON.SLIDING_WINDOW_SIZE / CON.SLIDE_AMOUNT)) - 1
        start_time = datetime.now()
        start = 0
        stop = CON.SLIDING_WINDOW_SIZE
        windows = []
        for i in range(window_count):
            windows.append((start, stop))
            start += CON.SLIDE_AMOUNT
            stop += CON.SLIDE_AMOUNT
        window_dictionary = {}
        for window in windows:
            key = window
            window_dictionary[key] = []
            for dictionary in self.all_peaks:
                if window[0] <= dictionary["retention time"] < window[1]:
                    window_dictionary[key].extend(dictionary["tuple list"])
        ret = ((retention_time // CON.SLIDING_WINDOW_SIZE) + 1) * CON.SLIDING_WINDOW_SIZE
        print("(Sliding Window)")
        print("RT {}s / {}s".format(0, int(ret)))
        for rt, tuple_list in window_dictionary.items():
            window_dictionary[rt] = tuple_combine(tuple_list)
            if int(float(rt[0])) != 0 and int(float(rt[0])) % 600 == 0:
                print("RT {}s / {}s".format(int(rt[0]), int(ret)))
        print("Sliding Window Complete")
        print("Time elapsed:", datetime.now() - start_time, "\n")
        return window_dictionary

    def add_peptide(self, sequence, mz, charge, pep_id, scan):
        self.peptides.append(Peptide(sequence, mz, charge, pep_id, scan))

    def iterlists(self, index):
        yield from self.all_peaks[index]["tuple list"]

    # outputs a .csv that has the data for a TIC scatter plot
    def generate_total_ion_chromatogram(self):
        mass_ratio = float(input("Enter desired m/z ratio: "))
        tolerance = float(input("Enter m/z tolerance: "))
        start_time = datetime.now()
        tic_list = []
        for i in self.all_peaks:
            retention_time = i["retention time"]
            for mz, intensity in i["tuple list"]:
                if abs(mass_ratio - mz) <= tolerance:
                    tic_list.append((retention_time, intensity))
        elapsed_time = datetime.now() - start_time
        print("\nTime to generate TIC: {}\n".format(elapsed_time))
        user_input = input("What would you like to name this file? ")
        with open(user_input + ".csv", 'w', newline='') as f:
            csv_writer = csv.writer(f)
            for i in tic_list:
                csv_writer.writerow(i)

    def match_peptides(self, pep):
        start, end = pep.get_rt_start_end()

        charge = pep.get_charge()
        pep_mass_over_charge = pep.get_mass_over_charge()
        pep_mass = pep_mass_over_charge * charge
        tuple_list = []
        # collects windows that match the rt of the sequence to user rt tolerance.
        for rt, tup_list in self.get_tuple_dictionary().items():
            if start <= rt[0] <= end or start <= rt[1] <= end:
                tuple_list.extend(tup_list)
        tuple_list.sort(key=lambda x: x[0])
        # searches for a match for each deuteration
        for det in range(pep.get_max_deuterium() + 1):
            ppm_error, mz, intensity = compare(pep_mass + det * CON.DEUTERIUM_MASS,
                                               charge, tuple_list, tuple_list)
            if ppm_error != 0:
                pep.set_deuterium(det, mz, intensity, ppm_error)
        pep.set_weighted_mass()

    def hydrogen_deuterium_exchange(self, identification_file: str):
        """This function reads the identification file and converts the
        scan #'s to retention times. Then it goes through each peptide and
        matches each deuterium if possible. This ends by writing to an output
        file which contains information on each potential deuteration. This
        should be run for each experimental run."""
        self.read_input(CON.IDENTIFICATION_CSV_FILE)
        self.set_pep_retention_times(identification_file)
        count = 0
        start_time = datetime.now()
        print("Begin matching")
        for pep in self.peptides:
            self.match_peptides(pep)
            count += 1
            if count % 20 == 0 or count == 1 or count == len(self.peptides):
                print(count, "/", len(self.peptides))
        print("Time to match:", datetime.now() - start_time)
        complexity = "Free"
        if self._complexity:
            complexity = "Complex"
        file = (CON.FULL_HDX_OUTPUT + "_" + str(self._time) + "s_"
                + str(complexity) + "_" + str(self._run + 1) + ".csv")
        self.write_hdx(file)

    # outputs tabular info of HDX results to .csv
    def write_hdx(self, file: str):
        output_exists = path.exists(CON.FULL_HDX_OUTPUT)

        with open(file, "w+", newline='') as f:
            csv_writer = csv.writer(f)
            if not output_exists:
                header = ["Start", "End", "Sequence", "Charge", "SequenceMz", "Complex",
                          "Deuterium", "RT", "Mz", "Intensity", "PpmError", "Average",
                          "Gaussian Fit"]
                csv_writer.writerow(header)
            lines = []
            for pep in self.peptides:
                line_list = pep.get_rows(self._complexity)
                for line in line_list:
                    lines.append(line)
            for line in lines:
                csv_writer.writerow(line)
            print("Wrote to:", file)


######################################################################
class Peptide:
    def __init__(self, sequence, mz, charge, pep_id, scan):
        self._windows = []
        self._sequence = sequence[2:len(sequence) - 2]
        self._charge = charge
        self._mass_over_charge = mz
        self._id = pep_id
        self._rt_start = 0
        self._rt_end = float("inf")
        self._scan = scan
        self._weighted_mass_to_charge = 0
        self._deuterium_dictionary = {}
        self._mass_shift = 0
        self.max_deuterium = sequence_to_max_deuterium(self._sequence)
        for det in range(self.max_deuterium + 1):
            self._deuterium_dictionary[det] = {"m/z": 0, "intensity": 0, "ppm": 0}
        self._average_mass = 0
        self.set_average_mass()
        self._protein = parse_protein(CON.PROTEIN_SEQUENCE_FILE)
        self._start, self._end = find_start_end(self._sequence, self._protein)
        if self._start == "NULL":
            print(self._sequence)
        self._fit = 0  # Gaussian fit

    # Getters
    def get_fit(self):
        return self._fit

    def get_mass_shift(self):
        return self._mass_shift

    def get_weighted_mass(self):
        return self._weighted_mass_to_charge

    def get_average_mass(self):
        return self._average_mass

    def get_start(self):
        return self._start

    def get_end(self):
        return self._end

    def get_sequence(self):
        return self._sequence

    def get_charge(self):
        return self._charge

    def get_mass_over_charge(self):
        return self._mass_over_charge

    def get_max_deuterium(self):
        return self.max_deuterium

    def get_rt_start_end(self):
        return self._rt_start, self._rt_end

    # returns peptide data for the detailed output
    def get_rows(self, complexity: str):
        line_list = []
        self.set_fit()
        for i in range(self.get_max_deuterium() + 1):

            mz, intensity, ppm = self.get_deuterium(i)
            line = ["", "", "", 0, 0, "", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            line[0] = self._start
            line[1] = self._end
            line[2] = self._sequence
            line[3] = self._charge
            line[4] = self._mass_over_charge
            line[5] = complexity
            line[6] = i
            line[7] = (self._rt_end + self._rt_start) / 2 / CON.MINUTES_TO_SECONDS
            line[8] = mz
            line[9] = intensity
            line[10] = ppm
            line[11] = self._weighted_mass_to_charge
            line[12] = self._fit
            line_list.append(line)
        return line_list

    # returns peptide data for the summary
    def get_table_row(self):
        line = ["", "", "", 0, 0, 0]
        line[0] = self._start
        line[1] = self._end
        line[2] = self._sequence
        line[3] = self._average_mass
        line[4] = (self._rt_end + self._rt_start) / 2 / CON.MINUTES_TO_SECONDS
        line[5] = self._mass_shift
        return line

    def get_scan(self):
        return self._scan

    # returns mz, intensity, ppm_error of of one match
    def get_deuterium(self, deuterium):
        mz = self._deuterium_dictionary[deuterium]["m/z"]
        intensity = self._deuterium_dictionary[deuterium]["intensity"]
        ppm = self._deuterium_dictionary[deuterium]["ppm"]
        return mz, intensity, ppm

    # Setters
    def set_fit(self):
        intensities = []
        for key in self._deuterium_dictionary.keys():
            intensities.append(self._deuterium_dictionary[key]["intensity"])
        self._fit = fit_gaussian(intensities)

    def set_windows(self, window):
        self._windows.append(window)

    def set_deuterium(self, det, mz, intensity, ppm_error):
        self._deuterium_dictionary[det]["m/z"] = mz
        self._deuterium_dictionary[det]["intensity"] = intensity
        self._deuterium_dictionary[det]["ppm"] = ppm_error

    def set_retention_times(self, start: float, end: float):
        self._rt_start = start
        self._rt_end = end

    def set_weighted_mass(self):
        total_intensity = 0
        mass = 0
        for det in self._deuterium_dictionary.keys():
            mz, intensity, ppm_error = self.get_deuterium(det)
            total_intensity += intensity
            mass += intensity * mz
        if total_intensity == 0:
            mass = 0
        else:
            mass /= total_intensity
        self._weighted_mass_to_charge = mass
        if self._weighted_mass_to_charge == 0:
            print("No matches found for:", self._sequence)
        weighted_mass = self._weighted_mass_to_charge - CON.MASS_OF_HYDROGEN
        self._mass_shift = weighted_mass * self._charge - self._average_mass

    # calculates the average mass from the sequence (Uses values in the PARAMETERS.py file)
    def set_average_mass(self):
        mass = 0
        for amino in self._sequence:
            mass += CON.PEPTIDE_MASS_DICTIONARY[amino]
        mass += CON.MASS_OF_WATER
        self._average_mass = mass


def get_time_points():
    num_time_points = 0
    while num_time_points < 1:
        num_time_points = int(input("How many time points? "))
        if num_time_points < 1:
            print("Please enter a positive integer.")
    print("Please enter each time point")
    time_points = []
    for i in range(num_time_points):
        usr_input = -1
        while usr_input < 0:
            usr_input = int(input("Time point #{} (s): ".format(i + 1)))
            if usr_input < 0:
                print("Please enter 0 or a positive integer.")
        time_points.append(usr_input)
    return time_points


def get_is_differential():
    is_differential = -1
    while is_differential == -1:
        differential_input = input("Was this a differential experiment (Y/N)? ").lower().strip()
        if differential_input[0] == 'y':
            is_differential = True
        elif differential_input[0] == 'n':
            is_differential = False
        else:
            print("Please enter 'Yes' or 'No'")
    return is_differential


def get_num_replications():
    num_replications = 0
    while num_replications < 1:
        num_replications = int(input("How many replications? "))
        if num_replications < 1:
            print("Please enter a positive integer.")
    return num_replications


def show_menu():
    print("Please enter the number of one of the following selections:")
    print("(1) Generate detailed output from mzML")
    print("(2) Generate Summary Output from detailed outputs")
    print("(3) Generate Figures from Summary Output")


def main():
    menu_input = None
    while menu_input != 'q':
        show_menu()
        menu_input = input().strip()[0]
        if menu_input == '1':
            start_time = datetime.now()
            check_parameters()
            # Get User Input
            time_points = get_time_points()
            is_differential = get_is_differential()
            num_replications = get_num_replications()
            # Perform peak matching and generate output
            experiment = FullExperiment(time_points, is_differential, num_replications)
            for time in time_points:
                experiment.add_runs(time)
            experiment.generate_output()
            print("Total Time Elapsed:", datetime.now() - start_time)


if __name__ == '__main__':
    main()
