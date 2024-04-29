from sklearn.cluster import KMeans
import numpy as np

# Dummy data
X = np.array([[1, 2], [1, 4], [1, 0],
              [10, 2], [10, 4], [10, 0]])

kmeans = KMeans(n_clusters=2, random_state=0, n_init=10).fit(X)
print(kmeans.labels_)

