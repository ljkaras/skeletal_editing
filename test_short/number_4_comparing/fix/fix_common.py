def modify_results_file(input_file, output_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    modified_lines = []
    current_file1 = None

    for line in lines:
        if line.strip().startswith('Comparison between'):
            # Extract the first file name from the line
            parts = line.strip().split(' ')
            current_file1 = parts[2]  # Assuming "Comparison between file1 and file2:"

        if line.strip().startswith('Number of common molecules'):
            # Append the first file name before the colon to this line
            count = line.strip().split(':')[-1].strip()
            line = f'Number of common molecules in {current_file1}: {count}\n'

        modified_lines.append(line)

    with open(output_file, 'w') as file:
        file.writelines(modified_lines)

# Usage
input_file = 'results.txt'
output_file = 'modified_results.txt'
modify_results_file(input_file, output_file)
