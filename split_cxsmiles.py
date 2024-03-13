import os
import gzip

def split_and_compress_file(filename, max_size_mb=2000):
    """
    Splits a large file into smaller files and compresses them using gzip.

    Parameters:
    filename (str): The path to the input file.
    max_size_mb (int): The maximum size of each split file in megabytes.
    """
    max_size = max_size_mb * 1024 * 1024  # Convert MB to bytes
    part_number = 1
    current_size = 0
    current_file = None

    with open(filename, 'r') as file:
        for line in file:
            if current_file is None or current_size + len(line.encode('utf-8')) > max_size:
                if current_file is not None:
                    current_file.close()
                    compress_file(current_file_name)
                current_file_name = f"{filename}_part_{part_number}.cxsmiles"
                current_file = open(current_file_name, 'w')
                part_number += 1
                current_size = 0
                print(f"Creating and compressing new file: {current_file_name}")
            current_file.write(line)
            current_size += len(line.encode('utf-8'))

    if current_file is not None:
        current_file.close()
        compress_file(current_file_name)

def compress_file(file_name):
    """
    Compresses a file using gzip and deletes the original file.
    Parameters:
    file_name (str): The path to the file to compress.
    """
    with open(file_name, 'rb') as f_in:
        with gzip.open(f'{file_name}.gz', 'wb') as f_out:
            f_out.writelines(f_in)
    os.remove(file_name)  # Delete the original file after compression

# Example usage:
filename = "real350.cxsmiles"  # Replace with the path to your large file
split_and_compress_file(filename, max_size_mb=2000)  # Splits into ~2GB files before compression
