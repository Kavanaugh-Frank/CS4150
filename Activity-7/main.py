import random
import math
import sys

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

    data = extracted_nps
    headers = extracted_headers

def calculate_similarity_matrix(data):
    # making a None matrix # of Cols x # of Cols
    result_matrix = [[None for _ in range(len(data))] for _ in range(len(data))]
    
    for col1 in range(len(data)):
        for col2 in range(0, len(data)):
            W = 0 # A = 1, B = 1
            X = 0 # A = 1, B = 0
            Y = 0 # A = 0, B = 1

        #going down the rows to increment the 3 counters based on the values of A and B
            for index in range(len(data[0])):
                A = int(data[col1][index])
                B = int(data[col2][index])
                if A == 1 and B == 1:
                    W += 1
                elif A == 1 and B == 0:
                    X += 1
                elif A == 0 and B == 1:
                    Y += 1
        
            J = W / min(X + W, Y + W)

            # setting the Matrix with the similarity scores
            result_matrix[col1][col2] = J
    return result_matrix

def calculate_medoid_on_columns(group_result_matrix):
    matrix_size = len(group_result_matrix)
    lowest_col = float('inf')
    lowest_col_idx = -1 

    for col in range(matrix_size):
        column_total_distance = 0
        for other_col in range(col + 1, matrix_size):
            column_total_distance += (1 - group_result_matrix[col][other_col])

        if column_total_distance < lowest_col:
            lowest_col = column_total_distance
            lowest_col_idx = col

    return lowest_col_idx, lowest_col

def calculate_least_diff_col(group_result_matrix):
    lowest_col = float("inf")
    lowest_col_idx = -1

    matrix_size = len(group_result_matrix)

        # calculate the average of each row
    row_sums = []
    for i in range(matrix_size):
        row_sum = sum(group_result_matrix[i]) / matrix_size
        row_sums.append(row_sum)

        # calculate the total difference of each column compared to the rest of the columns
        # to find the column with the shortest distance from the rest of the columns
    for col in range(matrix_size):
        column_total_diff = 0
        for i in range(matrix_size):
                # calculating the total distance of the specific col/row against the row average
            row_total_diff = abs(float(group_result_matrix[i][col]) - row_sums[i])
            column_total_diff += row_total_diff
            
            # keeping track of the column with the lowest total difference against the averages
        if column_total_diff < lowest_col:
            lowest_col = column_total_diff
            lowest_col_idx = col
    return lowest_col_idx

def k_medoid_clustering(result_matrix, data, headers, number_of_groups, k_values):
    k_groups = [[] for _ in range(number_of_groups)]

    previous_k_groups = []
    past_lowest_col = float('inf')  # past lowest column cost

    for i in range(100):
        new_k_values = []
        current_lowest_col = float('inf')  # track the current lowest column cost

        # assigning columns to groups
        k_groups = [[] for _ in range(number_of_groups)]
        for i in range(len(result_matrix)):
            similarities = [result_matrix[k][i] for k in k_values]
            max_value = max(similarities)
            if similarities.count(max_value) > 1:
                # there is a tie, so just whichever group at random
                max_indices = [value for value in similarities if value == max_value]
                max_index = similarities.index(random.choice(max_indices))
            else:
                max_index = similarities.index(max(similarities))
            k_groups[max_index].append(i)

        # check for oscillation
        if k_groups in previous_k_groups:
            print("Oscillated")
            # get some new start values
            k_values = [random.randint(0, len(data) - 1) for _ in range(number_of_groups)]
            k_groups = [[] for _ in range(number_of_groups)]
            for i in range(len(result_matrix)):
                similarities = [result_matrix[k][i] for k in k_values]
                max_index = similarities.index(max(similarities))
                k_groups[max_index].append(i)

        # keep track of previous k-groups to detect oscillation
        previous_k_groups.append(k_groups.copy())
        if len(previous_k_groups) > int(number_of_groups * 1.5):
            previous_k_groups.pop(0)

        # handle empty groups
        for group_idx in range(number_of_groups):
            if k_groups[group_idx] == []:
                while True:
                    new_medoid = random.randint(0, len(data) - 1)
                    if new_medoid not in new_k_values:
                        break
                k_values[group_idx] = new_medoid
                k_groups[group_idx] = [new_medoid]

        for k_group in k_groups:
            result = [data[index] for index in k_group]
            group_result_matrix = calculate_similarity_matrix(result)

            lowest_col_idx, lowest_col = calculate_medoid_on_columns(group_result_matrix)

            # only accept the change if it brings down the cost of the medoid
            if lowest_col < current_lowest_col:
                current_lowest_col = lowest_col

            if lowest_col < past_lowest_col:
                new_k_values.append(k_group[lowest_col_idx])
            else:
                new_k_values.append(k_values[k_groups.index(k_group)])

        # update the past lowest column cost
        past_lowest_col = current_lowest_col

        # check for convergence
        if new_k_values == k_values:
            break

        k_values = new_k_values

    return k_values, k_groups

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
            average_variation = 1000
        else:
            average_variation = math.sqrt(variation / (size**2 - 1))
        total_average_variation += average_variation

    return total_average_variation


# main body of the file, runs the algorithm x number of times and keeps track of the best result

# check if the number of iterations is provided as a command line argument
if len(sys.argv) != 2:
    print("Usage: python main.py <number_of_iterations>")
    sys.exit(1)

number_of_groups = 3

# calculate the similarity matrix for the text file data
result_matrix = calculate_similarity_matrix(data)
past_lowest_variance = float('inf')
best_initial_k_values = []
best_ending_k_values = []

for _ in range(0, int(sys.argv[1])):
    # init the k_values with random NPs that exist inside the data
    k_values = [random.randint(0, len(data) - 1) for _ in range(number_of_groups)]

    initial_k_values = k_values.copy()
    
    k_values, k_groups = k_medoid_clustering(result_matrix, data, headers, number_of_groups, k_values)

    variance = calculate_variance(k_values, k_groups, data, headers)

    if variance < past_lowest_variance:
        past_lowest_variance = variance
        best_initial_k_values = initial_k_values
        best_ending_k_values = k_values
        best_k_groups = k_groups

print("Best Initial K Values: ", best_initial_k_values)
print("Best Ending K Values: ", best_ending_k_values)
print("Best Total Variance: ", past_lowest_variance)