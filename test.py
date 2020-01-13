import PARAMETERS as CON


def get_ppm(mz1, mz2):
    return abs(mz1 - mz2) / mz1 * 1000000


def tuple_combine(some_list):
    changed = True
    just_changed = False
    start_list = sorted(some_list, key=lambda x: x[0])
    new_list = []
    for i in range(len(start_list)):
        tup = (start_list[i][0], start_list[i][1], 1)
        start_list[i] = tup

    while changed and len(start_list) != 0:
        changed = False
        for i in range(len(start_list) - 1):
            if just_changed:
                just_changed = False
                continue
            ppm = abs(get_ppm(start_list[i][0], start_list[i+1][0]))
            if ppm < CON.SLIDING_WINDOW_PPM_TOLERANCE:
                count1 = start_list[i][2]
                count2 = start_list[i + 1][2]
                mz1 = start_list[i][0] * count1
                mz2 = start_list[i + 1][0] * count2
                count = count1 + count2
                mz = (mz1 + mz2) / count
                intensity = start_list[i][1] + start_list[i + 1][1]
                new_list.append((mz, intensity, count))
                changed = True
                just_changed = True
            else:
                new_list.append(start_list[i])
        if not just_changed:
            new_list.append(start_list[len(start_list) - 1])
        start_list = new_list
        new_list = []
        just_changed = False
    for i in range(len(start_list)):
        tup = (start_list[i][0], start_list[i][1])
        start_list[i] = tup
    return start_list


my_list = [(400, 1000), (400.0001, 200), (400.0002, 30), (400.03, 4), (500, 9999)]
print(my_list)
print(tuple_combine(my_list))
