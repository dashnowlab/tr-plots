import csv

def read_tsv(file_path):
    data = []
    with open(file_path, 'r') as file:
        tsv_reader = csv.reader(file, delimiter='\t')
        for row in tsv_reader:
            data.append(row)
    return data

# Example usage:
file_path = '/Users/annelisethorn/Documents/Anschutz/Datasets/sample_information.tsv'
tsv_data = read_tsv(file_path)

# # Print the data
# for row in tsv_data:
#     print(row)

# Print the first column
first_column = [row[0] for row in tsv_data[1:0]]  # Skip first row (header)
print("First column:", first_column)

second_column = [row[1] for row in tsv_data[1:0]]  # Skip first row (header)
print("Second column:", second_column)

third_column = [row[2] for row in tsv_data[1:0]]  # Skip first row (header)
print("Third column:", third_column)

fourth_column = [row[3] for row in tsv_data[1:]]  # Skip first row (header)
print("Fourth column:", fourth_column)