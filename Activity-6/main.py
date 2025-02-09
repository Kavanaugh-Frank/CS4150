import random

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

result_matrix = calculate_similarity_matrix(data)

# seed the pseudo random number generator
# random.seed(15)
random.seed()

# how many k-groups should there be
number_of_groups = 3
# Ensure k_values are within the range of the data
k_values = [random.randint(0, len(data) - 1) for _ in range(number_of_groups)]
# make "number_of_groups" number of empty lists to hold values that belong to each group
k_groups = [[] for _ in range(number_of_groups)]

print(f"Randomly selected NPs {k_values}: {', '.join([headers[k] for k in k_values])}")

previous_k_groups = []
while True:
    new_k_values = []
    old_k_groups = k_groups

    k_groups = [[] for _ in range(number_of_groups)]
    for i in range(len(result_matrix)):
        similarities = [result_matrix[k][i] for k in k_values]
        max_index = similarities.index(max(similarities))
        k_groups[max_index].append(i)

    if k_groups in previous_k_groups:
        print("Oscillated")
        break

    previous_k_groups.append(k_groups.copy())
    if len(previous_k_groups) > 4:
        previous_k_groups.pop(0)

    for group_idx in range(number_of_groups):
        if k_groups[group_idx] == []:
            print(f"Group {group_idx} is empty, grabbing new random medoid for this group")
            new_medoid = random.randint(0, len(data) - 1)
            while new_medoid in new_k_values:
                new_medoid = random.randint(0, len(data) - 1)
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

        # loop through each column and calculate the total difference of that columns rows compared to the rest of the columns
        lowest_col = float("inf")
        lowest_col_idx = 0
        for col in range(0, len(group_result_matrix)):
            # the total difference the column has with every row against every column
            column_total_diff = 0
            for i in range(0, len(group_result_matrix)):
                # the total difference the column has with this specific row against every column
                row_total_diff = 0
                for j in range(0, len(group_result_matrix)):
                    row_total_diff += abs(float(group_result_matrix[i][col]) - float(group_result_matrix[i][j]))
                column_total_diff += row_total_diff / len(group_result_matrix)
            if column_total_diff < lowest_col:
                lowest_col = column_total_diff
                lowest_col_idx = col

        new_k_values.append(k_group[lowest_col_idx])

    print("New K Values: ", new_k_values)
    if new_k_values == k_values:
        print("Converged 2")
        break

    k_values = new_k_values


total_variation = 0
print("Final K Values: ", k_values)
for iteration, k_group in enumerate(k_groups):
    # getting the data in the column of the center of the group
    k_value_data = data[k_values[iteration]]

    # grabbing the rest of the columns that belong to the group, including the center
    result = []
    header = []
    header_idx_count = 0
    for index in k_group:
        if index == k_values[iteration]:
            k_value_header = headers[index]
            k_value_header_position = header_idx_count
        header.append(headers[index])
        result.append(data[index])
        header_idx_count += 1

    print("k_value_header: ", k_value_header)
    # print("k_value_header_position: ", k_value_header_position)
    group_result_matrix = calculate_similarity_matrix(result)

    variation = 0
    for i in range(0, len(group_result_matrix)):
        for j in range(0, len(group_result_matrix)):
            variation += abs(group_result_matrix[i][k_value_header_position] - group_result_matrix[i][j])

    total_variation += variation
    print("Number of NPs in this group: ", len(group_result_matrix))
    print(f"Total Variation (for group {iteration + 1}): ", round(variation, 5))
    print(f"Average Difference from the K-Value (for group {iteration + 1}): ", round(variation / (len(group_result_matrix) * len(group_result_matrix)), 5), "\n")

print("Total Variation: ", round(total_variation, 5))