import statistics
import matplotlib.pyplot as plt
import numpy as np
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

def categorize_value_z_score_col(value, col_average, col_std_dev):
    z_score = (value - col_average) / col_std_dev
    if z_score < -2.5:
        return 0, z_score
    elif -2.5 <= z_score < -1.5:
        return 1, z_score
    elif -1.5 <= z_score < -0.5:
        return 2, z_score
    elif -0.5 <= z_score < 0.5:
        return 3, z_score
    elif 0.5 <= z_score < 1.5:
        return 4, z_score
    elif 1.5 <= z_score < 2.5:
        return 5, z_score
    else:
        return 6, z_score

def categorize_value_z_score_row(value, row_average, row_std_dev):
    z_score = (value - row_average) / row_std_dev
    if z_score < -2.5:
        return 0, z_score
    elif -2.5 <= z_score < -2.0:
        return 1, z_score
    elif -2.0 <= z_score < -1.5:
        return 2, z_score
    elif -1.5 <= z_score < -1.0:
        return 3, z_score
    elif -1.0 <= z_score < -0.5:
        return 4, z_score
    elif -0.5 <= z_score < 0.0:
        return 5, z_score
    elif 0.0 <= z_score < 0.5:
        return 6, z_score
    elif 0.5 <= z_score < 1.0:
        return 7, z_score
    elif 1.0 <= z_score < 1.5:
        return 8, z_score
    elif 1.5 <= z_score < 2.0:
        return 9, z_score
    elif 2.0 <= z_score < 2.5:
        return 10, z_score
    else:
        return 11, z_score

def graph(col_hash, down_sample=False):
    categories = [v[0] for v in col_hash.values()]
    values = [v[1] for v in col_hash.values()]
    z_scores = [v[2] for v in col_hash.values()]

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

    if down_sample:
        indices = np.random.choice(len(categories), size=10000, replace=False)
        categories = [categories[i] for i in indices]
        values = [values[i] for i in indices]
        z_scores = [z_scores[i] for i in indices]

    # Plot the data
    plt.figure(figsize=(10, 6))

    for i in range(len(categories)):
        plt.scatter(z_scores[i], values[i], color=category_color_map[categories[i]])

    plt.xlabel("Z-Score (X-axis)")
    plt.ylabel("Value (Y-axis)")
    plt.title("Scatter Plot of Values by Z-Score and Category")

    plt.grid(True, linestyle="--", alpha=0.3)
    plt.show()

# Radial Position
print("Radial Position")
# Average
col_count = len(col_hash)
values = list(col_hash.values())
total_windows_present = sum(values)
col_average = round(total_windows_present / col_count, 2)
print("Average # of Detected Windows:", col_average)

# Standard Deviation
col_std_dev = round(statistics.stdev(values), 2)
print("Standard Deviation of Detected Windows:", col_std_dev)

category_counts_col = {i: 0 for i in range(7)}

for key, value in col_hash.items():
    category, z_score = categorize_value_z_score_col(value, col_average, col_std_dev)
    col_hash[key] = [category, value, z_score]
    category_counts_col[category] += 1

# graph(col_hash)

print("# of items in each category (0 & 6 are candidate outliers)", category_counts_col)
# print("Candidate Outliers:", [[key, value[1]] for key, value in col_hash.items() if value[0] == 0 or value[0] == 6])

# Compaction
print("\nCompaction")
# Average
row_count = len(row_hash)
values_row = list(row_hash.values())
total_windows_present = sum(values_row)
row_average = round(total_windows_present / row_count, 2)
print("Average # of Detected Windows:", row_average)

# Standard Deviation
row_std_dev = round(statistics.stdev(values_row), 2)
print("Standard Deviation of Detected Windows:", row_std_dev)

category_counts_row = {i: 0 for i in range(12)}

for key, value in row_hash.items():
    category, z_score = categorize_value_z_score_row(value, row_average, row_std_dev)
    row_hash[key] = [category, value, z_score]
    category_counts_row[category] += 1

# graph(row_hash, True)

print("# of items in each category (0 and 11 are candidate outliers)",category_counts_row)
# print("Candidate Outliers:", [[key, value[1]] for key, value in row_hash.items() if value[0] == 0 or value[0] == 11])