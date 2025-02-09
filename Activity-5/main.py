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

# making a None matrix # of Cols x # of Cols
result_matrix = [[None for _ in range(len(data))] for _ in range(len(data))]

# doing the Jaccard calculation
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


# seed the psuedo random number generator
random.seed()

# how many k-groups should there be
number_of_groups = 3
# random NPs
k_values = [random.randint(0,162) for _ in range(number_of_groups)]

# make "number_of_groups" number of empty lists to hold values that belong to each group
k_groups = [[] for _ in range(number_of_groups)]

print(f"Randomly selected NPs {k_values}: {', '.join([headers[k + 3] for k in k_values])}")
for i in range(len(result_matrix)):
    similarities = [result_matrix[k][i] for k in k_values]
    max_index = similarities.index(max(similarities))
    k_groups[max_index].append(i)


for i in range(number_of_groups):
    print(f"Group {i+1} has {len(k_groups[i])}")