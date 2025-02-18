import random
import math

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

    return lowest_col_idx

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

result_matrix = calculate_similarity_matrix(data)

# seed the pseudo random number generator
# random.seed(15)
random.seed()

# how many k-groups should there be
# this can be changed to any number
number_of_groups = 3

# init the k_values with 3 random NPs that exist inside the data
k_values = [random.randint(0, len(data) - 1) for _ in range(number_of_groups)]

# make empty k-groups
k_groups = [[] for _ in range(number_of_groups)]

print(f"Randomly selected NPs {k_values}: {', '.join([headers[k] for k in k_values])}")

previous_k_groups = []
# while True:
for i in range(100):
    new_k_values = []
    old_k_groups = k_groups

    k_groups = [[] for _ in range(number_of_groups)]
    for i in range(len(result_matrix)):
        similarities = [result_matrix[k][i] for k in k_values]
        max_index = similarities.index(max(similarities))
        k_groups[max_index].append(i)

    # if the current k-groups are the same as the last 6 k-groups
    # we determined the algorithm has oscillated between values and stop
    if k_groups in previous_k_groups:
        print("Oscillated\n")
        # get some more random k-values to replace the oscillated ones
        k_values = [random.randint(0, len(data) - 1) for _ in range(number_of_groups)]
        k_groups = [[] for _ in range(number_of_groups)]
        for i in range(len(result_matrix)):
            similarities = [result_matrix[k][i] for k in k_values]
            max_index = similarities.index(max(similarities))
            k_groups[max_index].append(i)
        # break

    # this is the list of previous k-groups that are compared to the current k-groups
    # keep the list to 6 total past k-groups
    previous_k_groups.append(k_groups.copy())
    if len(previous_k_groups) > int(number_of_groups * 1.5):
        previous_k_groups.pop(0)

    # just in case a group is empty make a new random medoid for that group, and make that medoid the only member of that group
    for group_idx in range(number_of_groups):
        if k_groups[group_idx] == []:
            print(f"Group {group_idx} is empty, grabbing new random medoid for this group")
            # so the new medoid is not the same as the other k-values
            while True:
                new_medoid = random.randint(0, len(data) - 1)
                if new_medoid not in new_k_values:
                    break
            k_values[group_idx] = new_medoid
            k_groups[group_idx] = [new_medoid]

    # grabbing the new medoid for each group
    for k_group in k_groups:
        # getting the data from each of the indices in the group
        result = []
        header = []
        for index in k_group:
            header.append(headers[index])
            result.append(data[index])

        group_result_matrix = calculate_similarity_matrix(result)

        # two separate ways to get the new medoid

        lowest_col_idx = calculate_medoid_on_columns(group_result_matrix)
        # lowest_col_idx = calculate_least_diff_col(group_result_matrix)

        new_k_values.append(k_group[lowest_col_idx])

    print("New K Values: ", new_k_values)
    
    # if the new k-values are the same as the old k-values
    # this means the algorithm has converged
    if new_k_values == k_values:
        break

    k_values = new_k_values


# Determining how well the algorithm did

total_variation = 0
total_average_variation = 0
print("Final K Values: ", k_values, "\n")
for iteration, k_group in enumerate(k_groups):
    # getting the data in the column of the center of the group
    k_value_data = data[k_values[iteration]]

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
    
    average_variation = math.sqrt(variation / (size**2 - 1))
    total_average_variation += average_variation

    print("K_Value_header: ", k_value_header)
    print(f"Variance (for group {iteration + 1}): ", round(average_variation, 5), "\n")
print("Total Variance (Lower the Better): ", round(total_average_variation, 5))