import networkx as nx
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
    extracted_headers = [row.split("\t")[0:3] for row in all_lines]
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

# some constant-ish value for all the divisions
row_length = len(extracted_windows[0])

# precompute the detection frequency
detection_freq = [sum(map(int, row)) / row_length for row in extracted_windows]

def calculate_normal_linkage(row_a, row_b, idx_a, idx_b):
    # co-segregation
    co_segregation = 0
    for idx, _ in enumerate(row_a):
        if int(row_a[idx]) == 1 and int(row_b[idx]) == 1:
            co_segregation += 1

    co_segregation = float(co_segregation) / row_length

    # detection freq of A and B
    freq_a = detection_freq[idx_a]
    freq_b = detection_freq[idx_b]

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
        calculate_normal_linkage(extracted_windows[i], extracted_windows[j], i, j) for j in range(num_windows)
    ]
    for i in range(num_windows)
]

# getting the unique linkage values into a single dimension array
# so that I can get the 75th percentile value
oneD_Linkage = []

# looping through the 2D array of normalized linkage values
# and appending them to a 1D array for sorting
# unique values only, and nothing where A = B
for i in range(len(normalized_linkage_table)):
    for j in range(i + 1, len(normalized_linkage_table)):
        oneD_Linkage.append(normalized_linkage_table[i][j])

print("Number of Unique Combinations of Normalized Linkage Table: ", len(oneD_Linkage))

# sort the array so that I can easily find the 75th percentile value
oneD_Linkage.sort()
percentile_75 = oneD_Linkage[int(0.75 * len(oneD_Linkage))]
print("75th Percentile Value: ", percentile_75)

# create a graph object
network_graph = nx.Graph()
# list that will store tuples of edges
edges = []

# going and finding the edges in the normalized linkage matrix
# an edge is defined as any value over the 75th percentile
for i in range(len(normalized_linkage_table)):
    for j in range(i + 1, len(normalized_linkage_table)):
        if normalized_linkage_table[i][j] > percentile_75:
            edges.append((i, j))

# Degree of Centrality
# 81 zeros for each of the windows
degreeOfCentrality = [0 for _ in range(len(normalized_linkage_table))]

# since there is 2 nodes that connect an edge
# have to increment both nodes for a single edge
for edge in edges:
    degreeOfCentrality[edge[0]] += 1
    degreeOfCentrality[edge[1]] += 1

# dividing each total count of edges with the N - 1 number of nodes (80 windows)
degreeOfCentrality = [(extracted_headers[i], i, (degree / 80)) for i, degree in enumerate(degreeOfCentrality)]

# sorting using the degree as the key
degreeOfCentrality.sort(key=lambda x: x[2])

# Calculate min, max, and average degree of centrality
min_degree = min(degreeOfCentrality, key=lambda x: x[2])
max_degree = max(degreeOfCentrality, key=lambda x: x[2])
avg_degree = sum(degree[2] for degree in degreeOfCentrality) / len(degreeOfCentrality)

print("Minimum Degree of Centrality:", min_degree)
print("Maximum Degree of Centrality:", max_degree)
print("Average Degree of Centrality:", avg_degree)


# writing the degrees of centrality in ascending order in to a file for easier reading
with open("output.txt", "w", encoding="UTF-8") as f:
    for degree in degreeOfCentrality:
        f.write(f"{degree}\n")


# Graphing Section
network_graph.add_edges_from(edges)

# seed the graphs random positioning so that it is the same each time
# and can be reproduced by my teammates to check for accuracy
pos = nx.spring_layout(network_graph, seed=123)  

# Adjust the spring layout to force more spacing between nodes
pos = nx.spring_layout(network_graph, seed=123)  # Increase the `k` value for more spacing

# plotting the network figure, using NetworkX
plt.figure(figsize=(8, 8))
nx.draw(
    network_graph, 
    pos=pos, 
    with_labels=True, 
    node_color="blue", 
    edge_color="gray",
    node_size=250, 
    font_size=10, 
    font_color="white"
)
plt.title("Network Graph of the Windows in the Hist1 Region")
plt.savefig("Network")


