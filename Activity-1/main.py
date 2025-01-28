col_hash = {}
row_hash = {}

with open("data.txt", "r") as file:
    headers = file.readline().split("\t")
    for line in file:
        # make the values of the row into a list using a tab as the delimiter
        values = line.split("\t")

        # row calculation
        chrom = values[0]
        start = values[1]
        stop = values[2]
        # create the string that will be the key for this GW
        window = str(chrom) + ": " + str(start) + " - " + str(stop)
        row_hash[window] = sum(map(int, values[3:]))

        # column calculation
        for i in range(3, len(headers)):
            key = headers[i]
            value = int(values[i])
            if key not in col_hash:
                col_hash[key] = value
            else:
                col_hash[key] += value

# Task  1
print("\nTask 1")
col_count = len(col_hash)
print("# of NPs:", col_count)

# Task 2
print("\nTask 2")
row_count = len(row_hash)
print("# of Genomic Windows:", row_count)

# Task 3
print("\nTask 3")
total_windows_present = sum(col_hash.values())
average = total_windows_present / col_count
print("Average # of Windows:", average)

# Task 4
print("\nTask 4")
min_key = min(col_hash, key=col_hash.get)
print("Min # of Windows:", col_hash[min_key], "(", min_key, ")")
max_key = max(col_hash, key=col_hash.get)
print("Max # of Windows:", col_hash[max_key], "(", max_key, ")")

# Task 5 
print("\nTask 5")
min_key = min(row_hash, key=row_hash.get)
print("Min # of Present Windows in GW:", row_hash[min_key], "(", min_key, ")")
max_key = max(row_hash, key=row_hash.get)
print("Max # of Present Windows in GW:", row_hash[max_key], "(", max_key, ")")
total = sum(row_hash.values())
print("Average # of Present Windows in GW:", total / row_count, "\n")

