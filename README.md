# README
***
## Python Installation
Head over to <a href="https://www.python.org/downloads/">Python.org</a> to install Python. Currently the supported version is: 3.7.5.

The only change to the settings you need to make is checking the box which adds Python to the Path.

Packages are premade code usable by other software. This program depends on a number of them which can all be installed by using a single command within the command prompt (Open using windows+r, type cmd):

`python -m pip install lxml pyteomics numpy matplotlib`

_**note: there are no numbers, any "l" character is the letter. It is recommended to copy paste the above line directly.**_

## Preparation
***
Open the file named _**Parameters.py**_ or _**Parameters**_. This file is divided into three sections.

In the first you change the input files to the names of your input files **(including the file extension e.g. .txt, .mzML)**.

_**Identification mzML File-**_ The mzML file which you performed the database search on which found the sequences you intend to search.

_**Identification CSV File-**_ The CSV file created from the output of your data base search. It should have the following columns (with correct name): ID, ScanNum, Precursor, Charge, Peptide

**ID** is a unique key for each sequence and is simply a column of ascending integers. **ScanNum** is the scan where the peptide was identified. **Precursor** is the m/z of the identified sequence, 
**Charge** is the z integer value, and **Peptide** is the sequence of amino acids.


In the second are parameters you can change to tune your results. In detail these are:

* **PPM Match Tolerance:** Peaks with an m/z difference of less than this value are considered a match. Adjust based on the sensitivity of your Mass Spectrometer.
* **Noise Limit:** This filters all individual peaks with an intensity under the PPM Match Tolerance.
* **Sliding Window PPM Tolerance:** The sliding window will consider peaks with a PPM difference lower than this to be the same and sum their intensities.
* **Sliding Window Size:** This is the width of the sliding window used to merge peaks which of the same m/z and similar retention times. 
* **Slide Fraction:** Determines the overlap of the window as it slides. A value of 1 would give no overlap. Recommended to keep below 4. Make sure <strong>Sliding Window Size</strong> is divisible by this value.
* **Retention Tolerance:** extends the retention time of an identification by plus and minus this value. Increase this value if there is a significant difference between the elution time of peptides in the identification run compared to the experimental run. 
***
_**Modify the constants with care, this can lead to incorrect evaluation.**_
***
## Running Software

In the command prompt navigate to the folder containing the program and associated files. This can be done by repeatedly typing "cd" followed by a space and the folder you wish to enter. Then type run the program by typing "python ReadMzml.py"

**Example**

If the program was stored at "C:\Users\User\Programs\DeuterationCalculator\ReadMzml.py"


The command prompt begins in C:\Users\User

Type `"cd Programs"`

Type `"cd DueterationCalculator"`

In this directory


Type `"python ReadMzml.py"`


This will start the program.

You will be prompted for the details of your experiment, then you will be asked for the path of each each of you experimental files. On Windows you can copy the path by holding shift, right clicking the file, and clicking "copy as path."