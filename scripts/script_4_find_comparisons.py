import sys
import os
import glob


def group_files_by_pattern(folder_path, output_file_path):
    files_by_pattern = {}

    # Find all CSV files in the folder
    csv_files = glob.glob(os.path.join(folder_path, "*2*.csv"))

    # Group files by the pattern in their filenames
    for file_path in csv_files:
        filename = os.path.basename(file_path)
        pattern_key = filename.split("2")[-1].split(".")[
            0
        ]  # Extract the pattern after "2" in the filename
        files_by_pattern.setdefault(pattern_key, []).append(file_path)

    # Write grouped file paths to a text file
    with open(output_file_path, "w") as outfile:
        for pattern_key, file_paths in files_by_pattern.items():
            for file_path in file_paths:
                comparison_name = pattern_key.split("_")[0]
                # print(comparison_name)
                comparison_file_path = os.path.join(
                    "frameworks/", "library_" + comparison_name + ".csv"
                )
                outfile.write(f"{file_path},{comparison_file_path}\n")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <folder_path>")
        sys.exit(1)

    folder_path = sys.argv[1]
    if not os.path.isdir(folder_path):
        print(f"Error: {folder_path} is not a valid directory.")
        sys.exit(1)

    # output_file_path = 'comparing.txt'
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_file_path = os.path.join(script_dir, "comparing.txt")
    group_files_by_pattern(folder_path, output_file_path)
