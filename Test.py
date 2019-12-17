def read_sequence(string):
    has_period = ('.' in string)
    if not has_period:
        return string
    found_period = False
    return_string = ""
    for letter in string:
        if letter == '.' and found_period:
            return return_string
        elif letter == '.':
            found_period = True
        elif found_period:
            return_string += letter
    raise ValueError("Sequence formatted incorrectly")



def main():
    sequence = "m.fffffdd.dg"
    print(read_sequence(sequence))

if __name__ == '__main__':
    main()
