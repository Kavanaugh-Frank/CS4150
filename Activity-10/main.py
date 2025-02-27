import random
import math
import sys
import matplotlib.pyplot as plt

random.seed()

# grabbing the CLI argument for the number of iterations the k means clustering will be going through
if len(sys.argv) != 2:
    print("Usage: python main.py <number_of_iterations>")
    sys.exit(1)

number_of_iterations = int(sys.argv[1])

# grab the indices of each of the features that need to be 
features = [
    "Hist1", "Vmn", "LAD", 
    "RNAPII-S2P", "RNAPII-S5P", "RNAPII-S7P", 
    "Enhancer", "H3K9me3", "H3K20me3", "h3k27me3", "H3K36me3", 
    "NANOG", "pou5f1", "sox2", "CTCF-7BWU"
]
# this is a map that has the feature name as well as the index it the data in its column will be in
feature_map = {}
with open("Hist1_region_features.csv", "r") as f:
    headers = f.readline().split(",")
    headers = [header.strip() for header in headers]
    print(headers)
    for feature in features:
        feature_map[feature] = headers.index(feature)

    lines = f.readlines()
    csv_data = [lines.split(",") for lines in lines]

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
    
    extracted_nps = []
    extracted_headers = []

    for col_index in range(3, len(column_header)):
        column = [line.split()[col_index] for line in all_lines]
        non_zero_found = False
        for value in column:
            if int(value) != 0:
                non_zero_found = True
                break
        if non_zero_found:
            extracted_nps.append(column)
            extracted_headers.append(column_header[col_index])

print("Number of NPs: ", len(extracted_nps))
print("Number of Windows: ", len(extracted_nps[0]))

# function to create a Jaccard Similarity 2D Matrix
# using the normalized denominator
def calculate_similarity_matrix(data):
        result_matrix = [[None for _ in range(len(data))] for _ in range(len(data))]
        
        for col1 in range(len(data)):
            for col2 in range(len(data)):
                W = 0 # A = 1, B = 1
                X = 0 # A = 1, B = 0
                Y = 0 # A = 0, B = 1

                for index in range(len(data[0])):
                    A = int(data[col1][index])
                    B = int(data[col2][index])
                    if A == 1 and B == 1:
                        W += 1
                    elif A == 1 and B == 0:
                        X += 1
                    elif A == 0 and B == 1:
                        Y += 1
            
                denominator = min(W + X, W + Y)
                if denominator == 0:
                    J = 0
                else:
                    J = W / denominator

                result_matrix[col1][col2] = J
        return result_matrix

# the k-medoids clustering functions
def k_means(extracted_nps, k_values, num_groups, max_iter=100):   
    # takes in the jaccard matrix, and the current k-values and assigns each
    # column index to a k-group 
    def assign_data_points_to_groups(jaccard_matrix, k_values, k_groups):
        for i in range(len(jaccard_matrix)):
            similarities = [jaccard_matrix[k][i] for k in k_values]
            max_value = max(similarities)
            if similarities.count(max_value) > 1:
                min_group_size = float('inf')
                min_group_size_index = -1
                for idx, value in enumerate(similarities):
                    if value == max_value:
                        group_size = len(k_groups[idx])
                        if group_size < min_group_size:
                            min_group_size = group_size
                            min_group_size_index = idx
                highest_similarity = min_group_size_index
            else:
                highest_similarity = similarities.index(max(similarities))
            k_groups[highest_similarity].append(i)
        
        return k_groups

    # takes in a group's jaccard matrix and finds the column that has the shortest total distance
    # it then returns the index of that column in the original group
    def calculate_medoid_on_columns(group_result_matrix, original_indices):
        matrix_size = len(group_result_matrix)
        lowest_col = float('inf')
        lowest_col_idx = -1 

        for col in range(matrix_size):
            column_total_distance = 0
            for other_col in range(matrix_size):
                if col != other_col:
                    column_total_distance += 1 - group_result_matrix[col][other_col]

            if column_total_distance < lowest_col:
                lowest_col = column_total_distance
                lowest_col_idx = col

        return original_indices[lowest_col_idx], lowest_col

    for _ in range(max_iter):
        jaccard_matrix = calculate_similarity_matrix(extracted_nps)

        k_groups = [[] for _ in range(num_groups)]

        k_groups = assign_data_points_to_groups(jaccard_matrix, k_values, k_groups)

        new_k_values = []
        for group in k_groups:
            original_indices = [i for i in group]
            k_group = [extracted_nps[i] for i in original_indices]
            group_jaccard_matrix = calculate_similarity_matrix(k_group)
            medoid, _ = calculate_medoid_on_columns(group_jaccard_matrix, original_indices)
            new_k_values.append(medoid)

        if sorted(new_k_values) == sorted(k_values):
            print("Converge")
            break
        k_values = new_k_values

    return new_k_values, k_groups

# calculate the total variance of the k-groups
# using a jaccard matrix for each of the groups
def calculate_variance(k_groups, data, headers):
    total_variation = 0
    total_average_variation = 0
    for _, k_group in enumerate(k_groups):
        # grabbing the rest of the columns that belong to the group, including the center
        result = []
        header = []
        header_idx_count = 0
        for index in k_group:
            header.append(headers[index])
            result.append(data[index])
            header_idx_count += 1

        group_result_matrix = calculate_similarity_matrix(result)
        
        size = len(group_result_matrix)
        total = 0
        for i in range(size):
            for j in range(size):
                total += group_result_matrix[i][j]
        average = total / (size * size)

        variation = 0
        for i in range(size):
            for j in range(size):
                variation += (group_result_matrix[i][j] - average) ** 2

        total_variation += variation
        average_variation = math.sqrt(variation / (size**2 - 1))
        total_average_variation += average_variation

    return total_average_variation


num_groups = 3
lowest_variance = float('inf')
best_initial_k_values = []
best_ending_k_values = []
best_k_groups = []
for _ in range(number_of_iterations):
    k_values = random.sample(range(len(extracted_nps)), num_groups)

    print("Initial Medoids: ", k_values)

    initial_medoid = k_values.copy()

    k_values, clusters = k_means(extracted_nps, k_values, num_groups)

    variance = calculate_variance(clusters, extracted_nps, extracted_headers)

    if variance < lowest_variance:
        lowest_variance = variance
        best_initial_k_values = initial_medoid
        best_ending_k_values = k_values
        best_k_groups = clusters

print("Final Medoids: ", k_values)
print("Lowest Variance: ", lowest_variance)
print("Cluster Sizes: ")
for i, cluster in enumerate(clusters):
    print(f"  Cluster {i+1}: {len(cluster)} elements")



for i, group in enumerate(best_k_groups):
    feature_counts = {feature: 0 for feature in features}
    group_data = [extracted_nps[i] for i in group]

    group_size = len(group)

    for feature in features:
        feature_index = feature_map[feature]
        feature_count = 0

        for row in group_data:
            for col_idx, col in enumerate(row):  
                if int(row[col_idx]) == 1 and int(csv_data[col_idx][feature_index]) == 1:
                    feature_count += 1

        feature_counts[feature] = feature_count / (len(group_data[0]) * len(group_data)) * 100

    print(feature_counts)
    # Create a spider chart
    labels = list(feature_counts.keys())
    stats = list(feature_counts.values())

    angles = [n / float(len(labels)) * 2 * math.pi for n in range(len(labels))]
    stats += stats[:1]
    angles += angles[:1]

    ax = plt.subplot(111, polar=True)
    plt.xticks(angles[:-1], labels)

    ax.plot(angles, stats)
    ax.fill(angles, stats, 'b', alpha=0.1)
    ax.set_xlabel(f"Cluster {i+1} with Medoid of {extracted_headers[best_ending_k_values[i]]}")

    # plt.title(f'Group {i+1} Feature Averages')
    plt.show()
