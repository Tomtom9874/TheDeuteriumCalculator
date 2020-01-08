# README
## Python Installation
Head over to <a href="https://www.python.org/downloads/">Python.org</a> to install Python. Currently the supported version is: 3.8.1.

The only change to the settings you need to make is checking the box which adds Python to the Path.

Packages are premade code usable by other software. This program depends on a number of them which can all be installed by using a single command within the command prompt (Open using windows+r, type cmd):

`python -m pip install lxml pyteomics numpy scipy pandas matplotlib`

_**note: there are no numbers, any "l" character is the letter. It is recommended to copy paste the above line directly.**_

## Preparation
Open the file named _**Parameters.py**_ or _**Parameters**_. You can use a text editor or an IDE such as PyCharm. This file is divided into five sections.

In the first you change the input files to the names of your input files **(including the file extension e.g. .txt, .mzML)**.
Leave the 'r' before the path to your input files. On Windows you can copy the path to a file by holding shift while right clicking the file.

_**Identification mzML File-**_ The mzML file which you performed the database search on which found the sequences you intend to search.

_**Identification CSV File-**_ The CSV file created from the output of your data base search. It should have the following columns (with correct name): ID, ScanNum, Precursor, Charge, Peptide

**ID** is a unique key for each sequence and is simply a column of ascending integers. **ScanNum** is the scan where the peptide was identified. 
**Precursor** is the m/z of the identified sequence, 
**Charge** is the z integer value, and **Peptide** is the sequence of amino acids.

_**Protein Sequence File**_ This is a .txt file which simply contains the string of characters representing the amino acids in the protein being examined.

The second section is the output files. Here you can either give the path to the location of each output file or you can simply put the name of the file and it will be added to the working directory.
*Once again, do not remove the r before the path*

In the third are parameters you can change to tune your results. In detail these are:

* **Noise Limit:** This filters all individual peaks with an intensity under the PPM Match Tolerance.
* **PPM Match Tolerance:** Peaks with an m/z difference of less than this value are considered a match. Adjust based on the sensitivity of your Mass Spectrometer.
* **Sliding Window PPM Tolerance:** The sliding window will consider peaks with a PPM difference lower than this to be the same and sum their intensities.
* **Sliding Window Size:** This is the width of the sliding window used to merge peaks which of the same m/z and similar retention times. 
* **Slide Fraction:** Determines the overlap of the window as it slides. A value of 1 would give no overlap. Recommended to keep below 4. Make sure <strong>Sliding Window Size</strong> is divisible by this value.
* **Retention Tolerance:** extends the retention time of an identification by plus and minus this value. Increase this value if there is a significant difference between the elution time of peptides in the identification run compared to the experimental run. 
* **Woods' Plot Confidence:** Decides where to draw the confidence line on Woods' plot. A two tailed test is used in this calculation.

The forth section allows for choices on where to output the Woods' plot as well as details on the plot.
***
_**Modify the constants with care, this can lead to incorrect evaluation.**_
***
## Running Software

In the command prompt navigate to the folder containing the program and associated files. This can be done by repeatedly typing "cd" followed by a space and the folder you wish to enter. Then type run the program by typing "python TheDeuteriumCalculator.py"

**Example**

If the program was stored at "C:\Users\User\Programs\DeuterationCalculator\TheDeuteriumCalculator.py"


The command prompt begins in C:\Users\User

Type `"cd Programs"`

Type `"cd DueterationCalculator"`

In this directory


Type `"python TheDeuteriumCalculator"`

This will start the program.
Each time the program runs you will be prompted for a number of details on the program. They are:

* The number of time points
* The time at which data were collected for each time point
* Whether the experiment was differential or not
* The number of replications (DO NOT DOUBLE FOR DIFFERENTIAL EXPERIMENTS)

This will bring you to the main menu. If this is the first time running the program for the current experiment type "1." This will begin the portion of the program which reads the data from the experimental mzML files.
You will be prompted for the path to each mzML file. You can copy the path on Windows by holding shift and right clicking the file and selecting "Copy as Path." The program will then begin processing the data automatically until this step is finished.
At that point detailed outputs will be generated which contain information on every peak matched. 

Step two can either be selected immediately or one can exit the program and manually edit data. This can be acheived by changing the intensity of any points that are incorrect to zero. 
**If changes are made to the data, the program must be restarted before they can be recognized.** Step two will generate the final outputs, including summary data and the Woods' Plot. 

***
## Notes

* Step one only needs to be completed once per experiment, unless user error is noticed such as incorrect selection of mzML files.
* If you are working on multiple experiments simultaneously, one must be careful to change the names of any outputs or move them, as they will otherwise be overwritten.