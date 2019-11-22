import numpy as np
import pandas as pd
import csv
import PARAMETERS as CON
from pyteomics import mzml
from datetime import datetime
from os import path
from matplotlib import pyplot as plt


def sequence_to_max_deuterium(sequence: str):
    max_deuterium = len(sequence) - 1
    for letter in sequence:
        if letter.lower() == 'p':
            max_deuterium -= 1
    return max_deuterium


def generate_differential_woods_plot(file: str, time: int, title: str, fractional=True):
    # Format plot
    plt.title(title)
    plt.xlabel("Quasi-Sequence")
    plt.ylabel("Relative Uptake (Da)")
    if fractional:
        plt.ylabel("Relative Fractional Uptake")
    plt.plot((0, 95), (0, 0), 'k:')

    df = pd.read_csv(file, header=[0, 1])
    time_col = str(time) + " min"
    for _, row in df.iterrows():
        x = (row["Start"]["Start"], row["End"]["End"])
        difference = row["Uptake COMPLEX (D)"][time_col] - row["Uptake FREE (D)"][time_col]
        if fractional:
            difference /= sequence_to_max_deuterium(row['Sequence']['Sequence'])
        y = (difference, difference)
        if abs(difference) < 10:
            plt.plot(x, y, 'k')
    plot_file_name = CON.WOODS_PLOT_NAME + "_" + str(time) + "s.png"
    plt.savefig(plot_file_name, dpi=600)


def generate_woods_plot(file: str, time: int, is_complex: bool):
    df = pd.read_csv(file, header=[0, 1], skipinitialspace=True)
    time_col = str(time) + " min"
    if is_complex:
        complexity_col = "Uptake COMPLEX (D)"
    else:
        complexity_col = "Uptake FREE (D)"
    for index, row in df.iterrows():
        x = (row["Start"]["Start"], row["End"]["End"])
        y = (row[complexity_col][time_col], row[complexity_col][time_col])
        if y[0] > -100:
            plt.plot(x, y, 'k')
    plt.show()


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
    if not check_extension(CON.SUMMARY_TABLE_1, "csv"):
        print("ERROR: SUMMARY_TABLE_1 must be a .csv")
    if CON.PPM_MATCH_TOLERANCE > 100:
        print("WARNING: PPM_MATCH_TOLERANCE is set to {}, this is high".format(CON.PPM_MATCH_TOLERANCE))
    if CON.SLIDING_WINDOW_PPM_TOLERANCE > 10:
        print("WARNING: SLIDING_WINDOW_PPM_TOLERANCE is set to {}, "
              "recommendation is under 10.".format(CON.SLIDING_WINDOW_PPM_TOLERANCE))
    if CON.SLIDING_WINDOW_SIZE % CON.SLIDE_AMOUNT != 0:
        print("ERROR: SLIDING_WINDOW_SIZE must be divisible by SLIDE_AMOUNT")
    if CON.RETENTION_TOLERANCE <= 0:
        print("ERROR: RETENTION_TOLERANCE must be a positive number, greater than 30 recommended.")
    if CON.RETENTION_TOLERANCE < 30:
        print("WARNING: A retention tolerance of greater than 30 is recommended")


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


# Removes any non-alpha characters form the protein sequence
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


# Takes i0n a list of tuples and combines all where first elements is within a ppm tolerance
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
            if abs(get_ppm(start_list[i][0], start_list[i+1][0])) < CON.SLIDING_WINDOW_PPM_TOLERANCE:
                new_list.append((((start_list[i][0] + start_list[i+1][0]) / 2), start_list[i][1] + start_list[i+1][1]))
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
            if abs(get_ppm(target, array[midpoint - offset][0] * charge)) <= CON.PPM_MATCH_TOLERANCE:
                return_list.append((array[midpoint - offset][0], array[midpoint - offset][1]))
                offset += 1
            else:
                offset = 0

        offset = 1
        while offset != 0 and (midpoint + offset < len(array)):
            if abs(get_ppm(target, array[midpoint + offset][0] * charge)) <= CON.PPM_MATCH_TOLERANCE:
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
        for i in f:
            if i["ms level"] == 2:
                retention_scan_dictionary[i["index"] + 1] = (float(i["scanList"]["scan"][0]["scan start time"])
                                                             * CON.MINUTES_TO_SECONDS)

    return retention_scan_dictionary


###################################################################################################
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

    def get_file(self, index: int):
        return self._file_names[index]

    # adds the name of each file to the list
    def add_file(self, time, complexity: bool, replication):
        complexness = "Free"
        if complexity:
            complexness = "Complex"
        print("Enter path to mzML for time:", time, "replication:", replication + 1, "complexity:", complexness)
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

    def generate_output(self):
        print("")
        print("Generating Output file(s)")

        ############################
        # Generates Summary Output 1
        top_header = ["Start", "End", "Sequence", "Peptide mass (Da)", "Retention time (min)"]
        header = ["Start", "End", "Sequence", "Peptide mass (Da)", "Retention time (min)"]

        # Free Average
        top_header.append("Uptake FREE (D)")
        for time in self._time_points:
            header.append(str(time) + " min")
            top_header.append("")
        top_header.pop()

        # Complex Average
        if self._is_differential:
            top_header.append("Uptake COMPLEX (D)")
            for time in self._time_points:
                header.append(str(time) + " min")
                top_header.append("")
            top_header.pop()

        # Free Deviation
        top_header.append("Uptake error (SD) - Free (D)")
        for time in self._time_points:
            header.append(str(time) + " min")
            top_header.append("")
        top_header.pop()

        # Complex Deviation
        if self._is_differential:
            top_header.append("Uptake error (SD) - Complex (D)")
            for time in self._time_points:
                header.append(str(time) + " min")
                top_header.append("")
            top_header.pop()

        peptide_lines = []
        # Each loop generates one line
        for index, pep in enumerate(self._peptides):
            peptide = [pep.get_start(), pep.get_end(), pep.get_sequence(), pep.get_average_mass()]
            start, end = pep.get_rt_start_end()
            rt = (start + end) / 2 / CON.MINUTES_TO_SECONDS
            peptide.append(rt)
            averages = []
            deviations = []
            # Free portion
            for time in self._time_points:
                replications = []
                for replication in range(self._num_replicates):
                    replications.append(self.runs[time][False][replication][index])
                average = sum(replications) / self._num_replicates
                averages.append(average)
                deviation = np.std(replications)
                deviations.append(deviation)
            # Complex Portion
            if self._is_differential:
                for time in self._time_points:
                    replications = []
                    for replication in range(self._num_replicates):
                        replications.append(self.runs[time][True][replication][index])
                    average = sum(replications) / self._num_replicates
                    averages.append(average)
                    deviation = np.std(replications)
                    deviations.append(deviation)
            peptide.extend(averages)
            peptide.extend(deviations)
            peptide_lines.append(peptide)

        with open(CON.SUMMARY_TABLE_1, "w+", newline='') as f:
            csv_writer = csv.writer(f)
            csv_writer.writerow(top_header)
            csv_writer.writerow(header)
            for line in peptide_lines:
                csv_writer.writerow(line)
        print("Wrote to:", CON.SUMMARY_TABLE_1)

        ###################
        # Generate Output 2
        # TODO Add Second Output

        # Generates a Wood's plot for each time-point
        for time in self._time_points:
            generate_differential_woods_plot(CON.SUMMARY_TABLE_1, time, "Differential Woods' Plot")


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

    # creates a dictionary for each scan with the RT and a list with tuples containing the m/z and intensity
    def read_mzml(self, file: str):
        total = 0
        print("(initializing)")
        with mzml.read(file) as f:
            for i in f:
                if i["ms level"] == 1:
                    total += 1
        count = 0
        retention_time = None
        with mzml.read(file) as f:
            for i in f:
                if i["ms level"] == 1:
                    retention_time = i["scanList"]["scan"][0]["scan start time"] * CON.MINUTES_TO_SECONDS
                    count += 1
                    if count % 200 == 0 or count == 1 or count == total:
                        print(count, "/", total)
                    tuple_list = []
                    for j in range(len(i["m/z array"])):
                        if i["intensity array"][j] >= CON.NOISE_LIMIT:
                            tuple_list.append((i["m/z array"][j], i["intensity array"][j]))
                    self.all_peaks.append({"retention time": retention_time, "tuple list": tuple_list})
        print("")
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
            if int(float(rt[0])) != 0:
                if int(float(rt[0])) % 600 == 0:
                    print("RT {}s / {}s".format(int(rt[0]), int(ret)))
        print("Sliding Window Complete")
        print("Time elapsed:", datetime.now() - start_time)
        print("")
        return window_dictionary

    def add_peptide(self, sequence, mz, charge, pep_id, scan):
        self.peptides.append(Peptide(sequence, mz, charge, pep_id, scan))

    def iterlists(self, index):
        yield from self.all_peaks[index]["tuple list"]

    # outputs a .csv that has the data for a TIC scatter plot
    def generate_total_ion_chromatogram(self):
        # TODO Modify how windows are selected
        # TODO Use windows to generate retention time peptide was found.
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

    # main function of program
    def hydrogen_deuterium_exchange(self, identification_file: str):
        self.read_input(CON.IDENTIFICATION_CSV_FILE)
        self.set_pep_retention_times(identification_file)
        count = 0
        start_time = datetime.now()
        tuple_dictionary = self.get_tuple_dictionary()
        print("Begin matching")
        for pep in self.peptides:
            start, end = pep.get_rt_start_end()
            count += 1

            charge = pep.get_charge()
            pep_mass_over_charge = pep.get_mass_over_charge()
            pep_mass = pep_mass_over_charge * charge

            tuple_list = []
            for rt, tup_list in tuple_dictionary.items():
                if start <= rt[0] <= end or start <= rt[1] <= end:
                    tuple_list.extend(tup_list)
            tuple_list.sort(key=lambda x: x[0])

            for det in range(pep.get_max_deuterium() + 1):
                ppm_error, mz, intensity = compare(pep_mass + det * CON.DEUTERIUM_MASS,
                                                   charge, tuple_list, tuple_list)
                if ppm_error != 0:
                    pep.set_deuterium(det, mz, intensity, ppm_error)
            pep.set_weighted_mass()
            if count % 20 == 0 or count == 1 or count == len(self.peptides):
                print(count, "/", len(self.peptides))
        print("Time to match:", datetime.now() - start_time)
        complexity = "Free"
        if self._complexity:
            complexity = "Complex"
        file = CON.FULL_HDX_OUTPUT + "_" + str(self._time) + "s_" + str(complexity) + "_" + str(self._run) + ".csv"
        self.write_hdx(file)
        # self.write_table(CON.SUMMARY_HDX_OUTPUT)

    # outputs tabular info of HDX results to .csv
    def write_hdx(self, file: str):
        output_exists = path.exists(CON.FULL_HDX_OUTPUT)

        with open(file, "w+", newline='') as f:
            csv_writer = csv.writer(f)
            if not output_exists:
                header = ["Start", "End", "Sequence", "Charge", "SequenceMz", "Complex",
                          "Deuterium", "RT", "Mz", "Intensity", "PpmError", "Average"]
                csv_writer.writerow(header)
            lines = []
            for pep in self.peptides:
                line_list = pep.get_rows(self._complexity)
                for line in line_list:
                    lines.append(line)
            for line in lines:
                csv_writer.writerow(line)
            print("Wrote to:", file)

    # outputs peptide summary to .csv
    def write_table(self, file: str):
        with open(file, "w+", newline='') as f:
            csv_writer = csv.writer(f)

            header = ["Start", "End", "Peptide", "PeptideAverageMass",
                      "RT", "DeuteriumUptake"]
            csv_writer.writerow(header)
            lines = []
            for pep in self.peptides:
                lines.append(pep.get_table_row())
            for line in lines:
                csv_writer.writerow(line)
            print(file)


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
        self._weighted_mass = 0
        self._deuterium_dictionary = {}
        self._mass_shift = 0
        # Set max deuterium

        self.max_deuterium = sequence_to_max_deuterium(self._sequence)

        for det in range(self.max_deuterium + 1):
            self._deuterium_dictionary[det] = {"m/z": 0, "intensity": 0, "ppm": 0}
        self._average_mass = 0
        self.set_average_mass()
        self._protein = parse_protein(CON.PROTEIN_SEQUENCE_FILE)
        self._start, self._end = find_start_end(self._sequence, self._protein)
        if self._start == "NULL":
            print(self._sequence)

    # Getters
    def get_mass_shift(self):
        return self._mass_shift

    def get_weighted_mass(self):
        return self._weighted_mass

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

        for i in range(self.get_max_deuterium() + 1):
            mz, intensity, ppm = self.get_deuterium(i)
            line = ["", "", "", 0, 0, "", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
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
            line[11] = self._weighted_mass
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
        self._weighted_mass = mass
        if self._weighted_mass == 0:
            print("No matches found for:", self._sequence)
        self._mass_shift = (self._weighted_mass - CON.MASS_OF_HYDROGEN) * self._charge - self._average_mass

    # calculates the average mass from the sequence (Uses values in the PARAMETERS.py file)
    def set_average_mass(self):
        mass = 0
        for amino in self._sequence:
            mass += CON.PEPTIDE_MASS_DICTIONARY[amino]
        mass += CON.MASS_OF_WATER
        self._average_mass = mass


def main():
    start_time = datetime.now()
    check_parameters()

    # Gets user input on time points
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

    # Gets user input on whether the run is differential
    is_differential = -1
    while is_differential == -1:
        differential_input = input("Was this a differential experiment (Y/N)? ").lower().strip()
        if differential_input[0] == 'y':
            is_differential = True
        elif differential_input[0] == 'n':
            is_differential = False
        else:
            print("Please enter 'Yes' or 'No'")

    # Gets user input on the number of replications
    num_replications = 0
    while num_replications < 1:
        num_replications = int(input("How many replications? "))
        if num_replications < 1:
            print("Please enter a positive integer.")

    experiment = FullExperiment(time_points=time_points, differential=is_differential, replications=num_replications)

    # Gets all file names
    for time in time_points:
        complexity = False
        for replication in range(num_replications):
            experiment.add_file(time, complexity, replication)
        if is_differential:
            complexity = True
            for replication in range(num_replications):
                experiment.add_file(time, complexity, replication)

    # Adds each run
    count = 0
    for time in time_points:
        complexity = False
        for replication in range(num_replications):
            experiment.add_run(time, complexity, replication, count)
            count += 1
        if is_differential:
            complexity = True
            for replication in range(num_replications):
                experiment.add_run(time, complexity, replication, count)
                count += 1
    experiment.generate_output()

    print("Total Time Elapsed:", datetime.now() - start_time)


if __name__ == '__main__':
    main()