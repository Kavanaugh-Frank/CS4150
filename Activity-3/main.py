import statistics
import numpy as np
import matplotlib.pyplot as plt
col_hash = {}
row_hash = {}

# calculate the z-score, radial pos, and compaction using the entire dataset, but only graph points that are present in chr13
# the basic statistics will be calculated using the chr13 dataset

with open("data.txt", "r") as file:
    headers = file.readline().split("\t")
    for line in file:
        # make the values of the row into a list using a tab as the delimiter
        values = line.split("\t")

        # row calculation
        chrom = values[0]
        start = values[1]
        stop = values[2]
        window = str(chrom) + " " + str(start) + " " + str(stop)

        # if chrom == "chr13" and start > "21690000" and stop < "24420000":
        # create the string that will be the key for this GW
        
        row_hash[window] = sum(map(int, values[3:]))

        # column calculation
        for i in range(3, len(headers)):
            key = headers[i]
            value = int(values[i])
            if key not in col_hash:
                col_hash[key] = value
            else:
                col_hash[key] += value

# getting the data values for the entire dataset so I can categorize chr13 data
row_values = list(row_hash.values())
row_average = statistics.mean(row_values)
row_std_dev = statistics.stdev(row_values)

col_values = list(col_hash.values())
col_average = statistics.mean(col_values)
col_std_dev = statistics.stdev(col_values)

print(f"Row Hash - Average: {row_average}, Standard Deviation: {row_std_dev}")
print(f"Col Hash - Average: {col_average}, Standard Deviation: {col_std_dev}")


# Graph Function
def graph(hash):
    categories = [v[0] for v in hash.values()]
    values = [v[1] for v in hash.values()]
    z_scores = [v[2] for v in hash.values()]

    category_color_map = {
        0: 'red',    
        1: 'blue',   
        2: 'green',   
        3: 'purple',  
        4: 'orange',  
        5: 'brown',   
        6: 'pink',    
        7: 'cyan',    
        8: 'magenta', 
        9: 'yellow',  
        10: 'black', 
        11: 'grey'   
    }

    plt.autoscale()
    for i in range(len(values)):
        jitter = np.random.normal(0, 0.02)
        plt.scatter(z_scores[i] + jitter, values[i] + jitter, color=category_color_map[categories[i]], alpha=0.5)

    plt.xlabel("Z-Score (X-axis)")
    plt.ylabel("Value (Y-axis)")
    plt.title("Scatter Plot of Values by Z-Score and Category")

    plt.grid(True, linestyle="--", alpha=.5)
    plt.show()

# Radial Position
def categorize_value_z_score_col(value, col_average, col_std_dev):
    z_score = (value - col_average) / col_std_dev
    if z_score < 0:
        return 0, z_score
    elif 0 <= z_score < 0.8:
        return 1, z_score
    elif 0.8 <= z_score < 1.6:
        return 2, z_score
    elif 1.6 <= z_score < 2.4:
        return 3, z_score
    elif 2.4 <= z_score < 3.2:
        return 4, z_score
    elif 3.2 <= z_score < 4.0:
        return 5, z_score
    else:
        return 6, z_score

# Compaction
def categorize_value_z_score_row(value, row_average, row_std_dev):
    z_score = (value - row_average) / row_std_dev
    if z_score < -2:
        return 0, z_score
    elif -2 <= z_score < -1.6:
        return 1, z_score
    elif -1.6 <= z_score < -1.2:
        return 2, z_score
    elif -1.2 <= z_score < -0.8:
        return 3, z_score
    elif -0.8 <= z_score < -0.4:
        return 4, z_score
    elif -0.4 <= z_score < 0.0:
        return 5, z_score
    elif 0.0 <= z_score < 0.4:
        return 6, z_score
    elif 0.4 <= z_score < 0.8:
        return 7, z_score
    elif 0.8 <= z_score < 1.2:
        return 8, z_score
    elif 1.2 <= z_score < 1.6:
        return 9, z_score
    elif 1.6 <= z_score < 2.0:
        return 10, z_score
    else:
        return 11, z_score

chr13_col_hash = {}
chr13_row_hash = {}
with open("chr13.txt", "r") as file:
    headers = file.readline().split("\t")
    for line in file:
        # make the values of the row into a list using a tab as the delimiter
        values = line.split("\t")

        # row calculation
        chrom = values[0]
        start = values[1]
        stop = values[2]
        window = str(chrom) + " " + str(start) + " " + str(stop)
        
        chr13_row_hash[window] = sum(map(int, values[3:]))

        # column calculation
        for i in range(3, len(headers)):
            key = headers[i]
            value = int(values[i])
            if key not in chr13_col_hash:
                chr13_col_hash[key] = value
            else:
                chr13_col_hash[key] += value

    keys_to_delete = [key for key, value in chr13_col_hash.items() if value == 0]
    for key in keys_to_delete:
        del chr13_col_hash[key]

    # Task  1
    print("\nTask 1")
    col_count = len(chr13_col_hash)
    print("# of NPs:", col_count)

    # Task 2
    print("\nTask 2")
    row_count = len(chr13_row_hash)
    print("# of Genomic Windows:", row_count)

    # Task 3
    print("\nTask 3")
    total_windows_present = sum(chr13_col_hash.values())
    average = round(total_windows_present / col_count, 2)
    print("Average # of Windows:", average)

    # Task 4
    print("\nTask 4")
    min_key = min(chr13_col_hash, key=chr13_col_hash.get)
    print("Min # of Windows:", chr13_col_hash[min_key], "(", min_key, ")")
    max_key = max(chr13_col_hash, key=chr13_col_hash.get)
    print("Max # of Windows:", chr13_col_hash[max_key], "(", max_key, ")")

    # Task 5 
    print("\nTask 5")
    min_key = min(chr13_row_hash, key=chr13_row_hash.get)
    print("Min # of Present Windows in GW:", chr13_row_hash[min_key], "(", min_key, ")")
    max_key = max(chr13_row_hash, key=chr13_row_hash.get)
    print("Max # of Present Windows in GW:", chr13_row_hash[max_key], "(", max_key, ")")
    total = sum(chr13_row_hash.values())
    print("Average # of Present Windows in GW:", round(total / row_count, 2), "\n")


# categorize the values in the chr13 dataset using the average and std deviation from the entire dataset
category_counts_row = {i: 0 for i in range(12)}

for key, value in chr13_row_hash.items():
    category, z_score = categorize_value_z_score_row(value, row_average, row_std_dev)
    chr13_row_hash[key] = [category, value, z_score]
    category_counts_row[category] += 1

graph(chr13_row_hash)
print("Category Row Counts for CHR13 against the entire dataset", category_counts_row)

# create a new hash that holds the union of every key that is present in both chr13_col_hash and col_hash
union_hash = {key: col_hash[key] for key in chr13_col_hash if key in col_hash}
category_counts_col = {i: 0 for i in range(7)}

for key, value in union_hash.items():
    category, z_score = categorize_value_z_score_col(value, col_average, col_std_dev)
    union_hash[key] = [category, value, z_score]
    category_counts_col[category] += 1

graph(union_hash)
print("Category Column Counts for CHR13 against the entire dataset", category_counts_col)