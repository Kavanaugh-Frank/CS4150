import matplotlib.pyplot as plt

# filter out the rows of the full genome to only the HIST1 region
with open("data.txt", "r") as f:
    header = f.readline()
    all_lines = f.readlines()
    lines = all_lines[69714:69714 + 81]

# creating a new text file that is just the HIST1 Chr13 region
with open("chr13.txt", "w") as f:
    f.write(header)
    f.writelines(lines)

# now reading the information from the HIST1 Chr13 region
# and getting the header, and full information of columns that have
# at least one '1' in them
with open("chr13.txt", "r") as f:
    column_header = f.readline().split("\t")  
    all_lines = f.readlines()

    # get the headers and windows from the dataset
    extracted_headers = [row.split("\t")[0:2] for row in all_lines]
    extracted_windows = [row.split("\t")[2:] for row in all_lines]

    # find the columns that need to be deleted (those without any 1's)
    cols_to_delete = []
    for col in range(len(extracted_windows[0])):
        zeros = True
        for row in range(len(extracted_windows)):
            if int(extracted_windows[row][col]) == 1:
                zeros = False
        if zeros:
            cols_to_delete.append(col)

    # removing the columns that are in cols_to_delete
    for row in range(len(extracted_windows)):
        new_row = []
        for col in range(len(extracted_windows[row])):
            if col not in cols_to_delete:
                new_row.append(extracted_windows[row][col])
        extracted_windows[row] = new_row

num_windows = len(extracted_headers)

def calculate_normal_linkage(row_a, row_b):
    # co-segregation
    co_segregation = 0
    for i in range(len(row_a)):
        if int(row_a[i]) == 1 and int(row_b[i]) == 1:
            co_segregation += 1

    co_segregation = float(co_segregation) / len(row_a)

    # Detection Freq function
    def calculate_detection_freq(row):
        detection_freq = 0
        for i in range(len(row)):
            if int(row[i]) == 1:
                detection_freq += 1
        
        return float(detection_freq) / len(row)
    
    # Detection Freq of A and B
    freq_a = calculate_detection_freq(row_a)
    freq_b = calculate_detection_freq(row_b)


    # Linkage
    linkage = (co_segregation - (freq_a * freq_b))
    
    # Normalized Linkage
    max_linkage = None

    if linkage < 0:
        max_linkage = min(freq_a*freq_b, (1-freq_a)*(1-freq_b))
    elif linkage > 0:
        max_linkage = min(freq_a*(1-freq_b), (1-freq_a)*freq_b)
    else:
        # if linkage == 0
        return linkage

    return (linkage/max_linkage)

# making a blank table
normalized_linkage_table = [[1 for _ in range(num_windows)] for _ in range(num_windows)]

# for every combination of windows
for i in range(num_windows):
    for j in range(num_windows):
        normalized_linkage_table[i][j] = calculate_normal_linkage(extracted_windows[i], extracted_windows[j])


print("Number of Windows: ", len(extracted_windows))
print("Number of NPs: ", len(extracted_windows[0]))

def plot_heatmap(matrix, title):
    plt.autoscale()
    plt.imshow(matrix, cmap='winter')
    plt.colorbar()
    plt.xlabel("Genomic Windows")
    plt.ylabel("Genomic Windows")
    plt.title(title)
    plt.show()

plot_heatmap(normalized_linkage_table, "Normalized Linkage Table for Genomic Windows in the HIST1 region")