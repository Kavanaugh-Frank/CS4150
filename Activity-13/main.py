import networkx as nx
import matplotlib.pyplot as plt
import math
import csv

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
    for index, header in enumerate(extracted_headers):
        extracted_headers[index] = f"{header[0]}:{header[1]}-{header[2]}"
    extracted_windows = [row.split("\t")[3:] for row in all_lines]

    # find the columns that need to be deleted (those without any 1's)
    cols_to_delete = []
    for col, _ in enumerate(extracted_windows[0]):
        all_zero = True
        for row, _ in enumerate(extracted_windows):
            if int(extracted_windows[row][col].strip()) == 1:
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

# another constant-ish value for the number of rows in the linkage table
# for this data set it should be 81, but for my sake just call the len function this
normalized_table_row_length = len(normalized_linkage_table)

# looping through the upper triangle of the 2D array of normalized linkage values
# and appending them to a 1D array for sorting
# unique values only, and nothing where A = B
for i in range(normalized_table_row_length):
    for j in range(i + 1, normalized_table_row_length):
        oneD_Linkage.append(normalized_linkage_table[i][j])

print("Number of Unique Combinations of Normalized Linkage Table: ", len(oneD_Linkage))

# sort the array so that I can easily find the 75th percentile value
oneD_Linkage.sort()
percentile_75 = oneD_Linkage[math.ceil(0.75 * len(oneD_Linkage))]
print("75th Percentile Value: ", percentile_75)

# list that will store tuples of edges
edges = []

# going and finding the edges in the normalized linkage matrix
# an edge is defined as any value over the 75th percentile
for i in range(normalized_table_row_length):
    for j in range(i + 1, normalized_table_row_length):
        if normalized_linkage_table[i][j] > percentile_75:
            edges.append((i, j))

# Degree of Centrality
# row_length zeros for each of the windows
degreeOfCentrality = [0 for _ in range(normalized_table_row_length)]

# since there is 2 nodes that connect an edge
# have to increment both nodes for a single edge
for edge in edges:
    degreeOfCentrality[edge[0]] += 1
    degreeOfCentrality[edge[1]] += 1

# dividing each total count of edges with the row_length - 1 number of nodes
# creating a tuple for each index -> (name, index, degree of centrality)
degreeOfCentrality = [
    (extracted_headers[i], i, (degree / (normalized_table_row_length - 1))) for i, degree in enumerate(degreeOfCentrality)
]

# sorting using the degree as the key
degreeOfCentrality.sort(key=lambda x: x[2])

# calculating min, max, and average degree of centrality
min_degree = min(degreeOfCentrality, key=lambda x: x[2])
max_degree = max(degreeOfCentrality, key=lambda x: x[2])
avg_degree = sum(degree[2] for degree in degreeOfCentrality) / len(degreeOfCentrality)

print("Minimum Degree of Centrality:", min_degree)
print("Maximum Degree of Centrality:", max_degree)
print("Average Degree of Centrality:", avg_degree)

# grabbing the 5 nodes with the largest degree of centrality
node_hubs = degreeOfCentrality[-5:]

print("Hubs:", node_hubs)

# function to find and return the neighbors of each of the node hubs
def find_neighbors(hub):
    neighbors = [] # list of the neighbors
    for edge in edges:
        if hub[1] not in edge:
            continue
        # check if the first edge is within the pair and if so add it
        # if we made it this fair then the second node is in the pair so add that
        neighbors.append(edge[0] if edge[1] == hub[1] else edge[1])
    return neighbors

# creating a list that holds all the neighbors for each of the hubs
neighbor_list = []
for i, hub in enumerate(node_hubs):
    neighbor_list.append(find_neighbors(hub))
    print("Size of Community", i, ":", len(neighbor_list[i]))
    print(neighbor_list[i])
# function to gather the information about each of the windows
# and the hist1 features within them
# returns a dictionary with the results
def preprocess_hist1_features(file_path):
    hist1_dict = {}
    with open(file_path, "r", encoding="UTF-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            hist1_dict[row["name"]] = int(row["Hist1"])
    return hist1_dict

# function to gather the information about each of the windows
# and the LAD features within them
# returns a dictionary with the results
def preprocess_lad_features(file_path):
    hist1_dict = {}
    with open(file_path, "r", encoding="UTF-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            hist1_dict[row["name"]] = int(row["LAD"])
    return hist1_dict

# call the preprocessing functions to gather the feature information
hist1_features = preprocess_hist1_features("Hist1_region_features.csv")    
lad_features = preprocess_lad_features("Hist1_region_features.csv")

# look through each of the center hubs
for index, hub in enumerate(node_hubs):
    hist1_count = 0
    lad_count = 0
    
    # for each center hub go through their neighbors and count the number of HIST1 and LAD features
    for node in neighbor_list[index]:
        node_name = degreeOfCentrality[node][0]
        if hist1_features[node_name] == 1:
            hist1_count += 1
        if lad_features[node_name] == 1:
            lad_count += 1

    # the number of neighbors each node has
    # used to calculate the percentages
    community_size = len(neighbor_list[index])

    # calculating the percentage of neighbors that has either the HIST1 or LAD feature
    hist1_percentage = (hist1_count / community_size) * 100
    print(f"\nCommunity {index}: {hist1_percentage:.2f}% of nodes have the HIST1 feature.")

    lad_percentage = (lad_count / community_size) * 100
    print(f"Community {index}: {lad_percentage:.2f}% of nodes have the LAD feature.")

def graph_community_interactions(community, community_index):
    # helper function to get the interactions between nodes in the same community
    def get_community_interactions(community):
        community_edges = []
        for edge in edges:
            if edge[0] in community and edge[1] in community:
                community_edges.append(edge)
        return community_edges
    community_edges = get_community_interactions(community)

    # make the graphing object
    network_graph = nx.Graph()
    # add the edges to the graph, which in turn adds the nodes
    network_graph.add_edges_from(community_edges)

    # plotting the network figure, using NetworkX
    plt.figure(figsize=(8, 8))

    # calculate node sizes based on the number of connections (degree of centrality)
    # multiply by the constant of 2500 to for a larger size in the graph and better visuals
    node_sizes = {index: centrality * 2500 for _, index, centrality in degreeOfCentrality}

    # make the layout of the graph, seed it with 123 for constant positions
    # and set k to 1.8 for spacing between nodes and less overlap
    pos = nx.spring_layout(network_graph, seed=123, k=1.8)  

    # put the title on the graph
    plt.title(f"Network Graph of the Windows in Community {community_index}")

    # draw the nodes and edges
    nx.draw(
        network_graph,
        pos=pos,
        with_labels=True,
        node_color="yellow",
        edge_color="gray",
        node_size=[node_sizes[node] for node in network_graph.nodes()],
        font_size=10,
        font_color="black",
        linewidths=1,
        edgecolors="black"
    )

    # save the graph
    plt.savefig(f"Network {community_index}")

# actually making the network graphs for each of the communities
for index, community in enumerate(neighbor_list):
    graph_community_interactions(community, index)

def plot_heatmap(matrix, community, community_index):
    # create a new matrix with zeros for rows/cols not in the community
    filtered_matrix = [[0 for _ in range(len(matrix))] for _ in range(len(matrix))]
    for i in community:
        for j in community:
            filtered_matrix[i][j] = matrix[i][j]

    plt.figure(figsize=(10, 8))
    plt.imshow(filtered_matrix, cmap="winter", vmin=-1, vmax=1)
    plt.colorbar()
    plt.xlabel("Genomic Windows")
    plt.ylabel("Genomic Windows")
    plt.title(f"Heatmap for Community {community_index}")
    plt.savefig(f"Heatmap for Community {community_index}")
    plt.close()

# plot heatmap for each community
for index, community in enumerate(neighbor_list):
    plot_heatmap(normalized_linkage_table, community, index)
