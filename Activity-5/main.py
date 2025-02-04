import random
import csv

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

    for col_index in range(3, len(column_header)):
        column = [line.split()[col_index] for line in all_lines]
        non_zero_found = False
        for value in column:
            if int(value) != 0:
                non_zero_found = True
                break
        if non_zero_found:
            extracted_nps.append(column)

    data = extracted_nps

# making a None matrix # of Cols x # of Cols
result_matrix = [[None for _ in range(len(data))] for _ in range(len(data))]

# doing the Jaccard calculation
for col1 in range(len(data)):
    for col2 in range(0, col1 + 1):
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


# seed the psuedo random number generator
random.seed()

# how many k-groups should there be
number_of_groups = 3
# generate radom k values for each of the k-groups
k_values = [random.uniform(0, 1) for _ in range(number_of_groups)]
# make "number_of_groups" number of empty lists to hold values that belong to each group
k_groups = [[] for _ in range(number_of_groups)]

# going through each and categoirizing them into the k-groups
for i in range(len(result_matrix)):
    for j in range(0, i+1):
        sim_score = [abs(k - result_matrix[i][j]) for k in k_values]
        min_index = sim_score.index(min(sim_score))
        k_groups[min_index].append((i, j, result_matrix[i][j], abs(k_values[min_index] - result_matrix[i][j])))


# sort the k_groups by the distance
k_groups = [sorted(group, key=lambda x: x[3]) for group in k_groups]

with open("k_groups.csv", "w", newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    headers = [f"Group {i+1}" for i in range(number_of_groups)]
    csvwriter.writerow(headers)
    
    max_length = max(len(k_groups[0]), len(k_groups[1]), len(k_groups[2]))
    for i in range(max_length):
        row = []
        for group in k_groups:
            if i < len(group):
                row.append([f"Value: {group[i][2]:.5f}", f"Distance: {group[i][3]:.5f}"])
            else:
                row.append("")
        csvwriter.writerow(row)

for i in range(number_of_groups):
    print(f"Group {i+1} has {len(k_groups[i])} pairs")