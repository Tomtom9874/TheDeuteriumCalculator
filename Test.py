from matplotlib import pyplot as plt
import pandas as pd


def sequence_to_max_deuterium(sequence: str):
    max_deuterium = len(sequence) - 1
    for letter in sequence:
        if letter.lower() == 'p':
            max_deuterium -= 1
    return max_deuterium


def generate_woods_plot(file: str, time: int, is_complex: bool):
    df = pd.read_csv(file, header=[0, 1])
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
    plt.savefig("TestGraph.png", dpi=1200, antialiasing=False)


def main():
    generate_differential_woods_plot("Summary_Table_1_PARE.csv", 0, "ParE Woods Plot", False)


if __name__ == '__main__':
    main()
