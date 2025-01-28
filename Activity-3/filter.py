with open("data.txt", "r") as f:
    header = f.readline()
    lines = f.readlines()[69716 - 1:69806 - 10]

with open("chr13.txt", "w") as f:
    f.write(header)
    f.writelines(lines)