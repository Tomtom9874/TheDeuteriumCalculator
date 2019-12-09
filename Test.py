def doThing(item):
    return item * 2


def main():
    myList = [0, 2, 4, 6, 8]
    print([doThing(x) for x in myList])

if __name__ == '__main__':
    main()
