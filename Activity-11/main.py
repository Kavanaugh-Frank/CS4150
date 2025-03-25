import matplotlib.pyplot as plt

# filter out the rows of the full genome to only the HIST1 region
with open("data.txt", "r", encoding="utf-8") as f:
    header = f.readline()
    all_lines = f.readlines()
    lines = all_lines[69714 : 69714 + 81]

# creating a new text file that is just the HIST1 Chr13 region
with open("chr13.txt", "w", encoding="utf-8") as f:
    f.write(header)
    f.writelines(lines)

# now reading the information from the HIST1 Chr13 region
# and getting the header, and full information of columns that have
# at least one '1' in them
with open("chr13.txt", "r", encoding="utf-8") as f:
    column_header = f.readline().split("\t")
    all_lines = f.readlines()

    # get the headers and windows from the dataset
    extracted_headers = [row.split("\t")[0:2] for row in all_lines]
    extracted_windows = [row.split("\t")[2:] for row in all_lines]

    # find the columns that need to be deleted (those without any 1's)
    cols_to_delete = []
    for col, _ in enumerate(extracted_windows[0]):
        all_zero = True
        for row, _ in enumerate(extracted_windows):
            if int(extracted_windows[row][col]) == 1:
                all_zero = False
        if all_zero:
            cols_to_delete.append(col)

    # removing the columns that are in cols_to_delete
    for row_index, row in enumerate(extracted_windows):
        new_row = []
        for col_index, value in enumerate(row):
            if col_index not in cols_to_delete:
                new_row.append(value)
        extracted_windows[row_index] = new_row

num_windows = len(extracted_headers)


def calculate_normal_linkage(row_a, row_b):
    # co-segregation
    co_segregation = 0
    for idx, _ in enumerate(row_a):
        if int(row_a[idx]) == 1 and int(row_b[idx]) == 1:
            co_segregation += 1

    co_segregation = float(co_segregation) / len(row_a)

    # detection freq function
    def calculate_detection_freq(inner_row):
        detection_freq = 0
        for _, inner_value in enumerate(inner_row):
            if int(inner_value) == 1:
                detection_freq += 1
        return float(detection_freq) / len(inner_row)

    # detection freq of A and B
    freq_a = calculate_detection_freq(row_a)
    freq_b = calculate_detection_freq(row_b)

    # linkage
    linkage = co_segregation - (freq_a * freq_b)

    # normalized linkage
    max_linkage = None

    if linkage < 0:
        max_linkage = min(freq_a * freq_b, (1 - freq_a) * (1 - freq_b))
    elif linkage > 0:
        max_linkage = min(freq_a * (1 - freq_b), (1 - freq_a) * freq_b)
    else:
        # if linkage == 0
        return linkage

    return linkage / max_linkage


# making the normalized linkage table
normalized_linkage_table = [
    [
        calculate_normal_linkage(extracted_windows[i], extracted_windows[j]) for j in range(num_windows)
    ]
    for i in range(num_windows)
]


print("Number of Windows: ", len(extracted_windows))
print("Number of NPs: ", len(extracted_windows[0]))


def plot_heatmap(matrix, title):
    plt.autoscale()
    plt.imshow(matrix, cmap="winter")
    plt.colorbar()
    plt.xlabel("Genomic Windows")
    plt.ylabel("Genomic Windows")
    plt.title(title)
    plt.savefig(title)


plot_heatmap(
    normalized_linkage_table,
    "Normalized Linkage Table for Genomic Windows in the HIST1 region",
)
