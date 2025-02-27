import random
import math
import sys
import matplotlib.pyplot as plt
import statistics
import numpy as np

random.seed()

# grabbing the CLI argument for the number of iterations the k means clustering will be going through
if len(sys.argv) != 2:
    print("Usage: python main.py <number_of_iterations>")
    sys.exit(1)

number_of_iterations = int(sys.argv[1])

col_hash = []

# grabbing data so that I can get the average and the standard deviation
with open("data.txt", "r") as file:
    full_headers = file.readline().split("\t")
    col_hash = [0 for _ in range(0, len(full_headers) - 3)]
    for line in file:
        values = line.split("\t")
        for col in range(3, len(values)):
            col_hash[col - 3] += int(values[col])

col_average = statistics.mean(col_hash)
col_std_dev = statistics.stdev(col_hash)
print("The Average Number of 1s for the Columns (NPs): ", col_average)
print("The Standard Deviation for the Columns (NPs): ", col_std_dev)

# Filter the Data.txt file to the required Chr13 data
with open("data.txt", "r") as f:
    header = f.readline()
    all_lines = f.readlines()
    lines = all_lines[69714:69714 + 81]

with open("chr13.txt", "w") as f:
    f.write(header)
    f.writelines(lines)

chr13_col_hash = {}

with open("chr13.txt", "r") as f:
    column_header = f.readline().split("\t")  
    all_lines = f.readlines()
    
    # extracting the NPs that detect a window
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

def k_means(extracted_nps, k_values, num_groups, max_iter=100):    
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

        return original_indices[lowest_col_idx], lowest_col  # Return actual column index

    for _ in range(max_iter):
        jaccard_matrix = calculate_similarity_matrix(extracted_nps)

        k_groups = [[] for _ in range(num_groups)]

        assign_data_points_to_groups(jaccard_matrix, k_values, k_groups)

        new_k_values = []
        for group in k_groups:
            original_indices = [i for i in group]
            # k_group = np.array(extracted_nps)[original_indices, :]
            k_group = [extracted_nps[i] for i in original_indices]
            group_jaccard_matrix = calculate_similarity_matrix(k_group)
            medoid, _ = calculate_medoid_on_columns(group_jaccard_matrix, original_indices)
            new_k_values.append(medoid)

        if sorted(new_k_values) == sorted(k_values):
            print("Converge")
            break
        k_values = new_k_values

    return new_k_values, k_groups

def calculate_variance(k_values, k_groups, data, headers):
    total_variation = 0
    total_average_variation = 0
    for iteration, k_group in enumerate(k_groups):
        # grabbing the rest of the columns that belong to the group, including the center
        result = []
        header = []
        header_idx_count = 0
        for index in k_group:
            # keep track of the position of the k-value header so that we know where it is in the similarity matrix
            # and can use it to find total variation later
            if index == k_values[iteration]:
                k_value_header = headers[index]
                k_value_header_position = header_idx_count
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
        if size == 1:
            # if the size of any of the groups is just 1, then I am assuming that
            # the variation will not be good enough so I sabotage with a high value
            # and size == 1 throws an error when calculating the average variation
            average_variation = 10000
        else:
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

    variance = calculate_variance(k_values, clusters, extracted_nps, extracted_headers)

    if variance < lowest_variance:
        lowest_variance = variance
        best_initial_k_values = initial_medoid
        best_ending_k_values = k_values
        best_k_groups = clusters

print(k_values, lowest_variance, len(clusters[0]), len(clusters[1]), len(clusters[2]))

def categorize_value_z_score_col(value, col_average, col_std_dev):
    z_score = (value - col_average) / col_std_dev
    if z_score < -2:
        return 0, z_score
    elif -2 <= z_score < -1:
        return 1, z_score
    elif -1 <= z_score < 0:
        return 2, z_score
    elif 0 <= z_score < 1:
        return 3, z_score
    else:
        return 4, z_score

category_counts_col = {i: [0] * len(best_k_groups) for i in range(5)}

for i, group in enumerate(best_k_groups):
    data_headers = [extracted_headers[index] for index in group]
    data_group = []
    for header in data_headers:
        data_group.append(col_hash[full_headers.index(header) - 3])
    print(len(data_group))
    for value in data_group:
        category, _ = categorize_value_z_score_col(value, col_average, col_std_dev)
        category_counts_col[category][i] += 1

# Plotting the bar graph for all groups
plt.figure(figsize=(10, 6))
bar_width = 0.2
for i in range(len(best_k_groups)):
    counts = [(category_counts_col[cat][i] / len(best_k_groups[i]) * 100) for cat in range(5)]
    plt.bar([j + bar_width * i for j in range(0, 5)], counts, width=bar_width, label=f'Group {i+1}')

plt.xlabel('Category')
plt.ylabel('Percentage')
plt.title('Category Counts for All Groups')
plt.xticks([i + bar_width for i in range(0, 5)], [str(i) for i in range(5)])
plt.legend()
# plt.show()
plt.savefig("category_counts.png")
