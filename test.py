from os import path


def check_directory(file_path):
    end_path = 0
    for index, char in enumerate(file_path):
        if char == '\\' or char == '/':
            end_path = index
    if end_path > 0:
        return path.exists(file_path[0:end_path])


exists = check_directory(r"C:\Users\Tom\PycharmProjects\TheDeuteriumCalculator\Outputs\HdxOutput_0s_Complex_3.csv")
print(exists)
