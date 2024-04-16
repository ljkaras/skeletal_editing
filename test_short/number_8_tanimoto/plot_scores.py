import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys

def plot_similarity_histogram(file_path):
    # Load the similarity matrix from the HDF5 file
    with h5py.File(file_path, 'r') as f:
        similarity_matrix = f['similarities'][:]

    # Flatten the matrix to get a list of all similarity scores
    similarity_scores = similarity_matrix.flatten()

    # Calculate and print the average similarity score
    average_similarity = np.mean(similarity_scores)
    print(f"Average similarity score: {average_similarity:.4f}")

    # Plot the histogram
    plt.hist(similarity_scores, bins=np.arange(0, 1.01, 0.01), edgecolor='black')
    plt.xlabel('Tanimoto Similarity')
    plt.ylabel('Frequency')

    # Set the title based on the filename without the extension
    title = file_path.rsplit('.', 1)[0]
    plt.title(f'Histogram of Tanimoto Similarities: {title}')

    # Save the plot as a PNG file
    plt.savefig(f'{title}_similarity_histogram.png')
    plt.show()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script.py <file_path>")
        sys.exit(1)

    file_path = sys.argv[1]
    plot_similarity_histogram(file_path)
