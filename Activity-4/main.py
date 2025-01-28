import matplotlib.pyplot as plt

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

# making a blank matrix # of Cols x # of Cols
result_matrix = [[0 for _ in range(len(data))] for _ in range(len(data))]

# doing the Jaccard calculation
for col1 in range(len(data)):
    for col2 in range(len(data)):
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
        
        # just in case there is an attempt to divide by zero
        if W + X + Y == 0:
            J = 0
        else:
            J = W / (W + X + Y)

        # setting the Matrix with the similarity scores
        result_matrix[col1][col2] = J

# heatmap of the result matrix
plt.figure(figsize=(10, 8))
plt.imshow(result_matrix, cmap='plasma')
plt.colorbar()
plt.title('Heatmap from the Jaccard Calculations')
plt.show()


# dump the result matrix to a CSV file
with open("result_matrix.csv", "w") as f:
    for row in result_matrix:
        output = ",".join(f"{value:.5f}" for value in row) + "\n"
        f.write(output)
