import argparse
import hashlib
import time
import csv

def read_fasta(filename):
    genomes = {}

    with open(filename, "r") as file:
        header = None
        sequence = []

        for line in file:
            line = line.strip()

            if line.startswith(">"):
                if header is not None:
                    genomes[header] = "".join(sequence)

                header = line[1:]
                sequence = []

            else:
                sequence.append(line)

        if header is not None:
            genomes[header] = "".join(sequence)

    return genomes


def get_kmers(sequence, k):
    kmers = set()

    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]

        if all(base in "ACGT" for base in kmer):
            kmers.add(kmer)

    return kmers


def hash_kmer(kmer):
    data = kmer.encode("ascii")

    digest = hashlib.blake2b(
        data,
        digest_size=8
    ).digest()

    return int.from_bytes(digest, byteorder="big")


def make_sketch(kmers, density):
    sketch = set()

    for kmer in kmers:
        h = hash_kmer(kmer)

        # Convert hash to a value between 0 and 1
        normalized_hash = h / (2**64)

        if normalized_hash < density:
            sketch.add(h)

    return sketch

def estimate_jaccard(sketch_a, sketch_b):
    intersection = sketch_a & sketch_b
    union = sketch_a | sketch_b

    if len(union) == 0:
        return 0.0

    return len(intersection) / len(union)


def pairwise_distances(sketches):
    names = list(sketches.keys())
    distances = {}

    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            name_a = names[i]
            name_b = names[j]

            jaccard = estimate_jaccard(
                sketches[name_a],
                sketches[name_b]
            )

            distance = 1 - jaccard

            distances[(name_a, name_b)] = distance

    return distances



# the real deal


parser = argparse.ArgumentParser()

parser.add_argument("filename")
parser.add_argument("k", type=int)
parser.add_argument("density", type=float)

args = parser.parse_args()

genomes = read_fasta(args.filename)

k = args.k
density = args.density

# Start measuring computational time
start_time = time.perf_counter()

sketches = {}

for name, sequence in genomes.items():
    kmers = get_kmers(sequence, k)
    sketches[name] = make_sketch(kmers, density)
    
after_sketch = time.perf_counter()
print(f"Sketching time: {after_sketch - start_time:.3f} s")

print(name, "number of k-mers:", len(kmers))
print(name, "sketch size:", len(sketches[name]))

distances = pairwise_distances(sketches)

after_distances = time.perf_counter()
print(f"Distance calculation time: {after_distances - after_sketch:.3f} s")

print("Pairwise distances:")

for pair, distance in distances.items():
    print(pair, "=", distance)

    #Task2
def make_distance_matrix(distances, names):
    matrix = []

    for name_a in names:
        row = []

        for name_b in names:
            if name_a == name_b:
                distance = 0.0

            elif (name_a, name_b) in distances:
                distance = distances[(name_a, name_b)]

            else:
                distance = distances[(name_b, name_a)]

            row.append(distance)

        matrix.append(row)

    return matrix

def save_distance_matrix(filename, names, matrix):
    with open(filename, "w", newline="") as file:
        writer = csv.writer(file)

        # Header row
        writer.writerow([""] + names)

        # Distance matrix
        for name, row in zip(names, matrix):
            writer.writerow([name] + row)


names = list(genomes.keys())

matrix = make_distance_matrix(distances, names)

# Save the distance matrix to a CSV file
save_distance_matrix("distances.csv", names, matrix)

print("Saved distance matrix: distances.csv")

print("Distance matrix:")

print("        ", names)

for name, row in zip(names, matrix):
    print(name, row)




"""def neighbor_joining(names, matrix):
    labels = list(names)

    # Store distances using pairs of labels
    distances = {}

    for i in range(len(labels)):
        for j in range(len(labels)):
            distances[(labels[i], labels[j])] = matrix[i][j]

    children = {}
    next_node = 1

    while len(labels) > 2:
        n = len(labels)

        # Calculate the total distance for each label
        total = {}

        for i in labels:
            total[i] = sum(
                distances[(i, j)]
                for j in labels
                if i != j
            )

        # Find the pair with the smallest Q value
        best_pair = None
        best_q = float("inf")

        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                a = labels[i]
                b = labels[j]

                q = (n - 2) * distances[(a, b)] - total[a] - total[b]

                if q < best_q:
                    best_q = q
                    best_pair = (a, b)

        a, b = best_pair

        # Create a new internal node
        new_node = "Node" + str(next_node)
        next_node += 1

        children[new_node] = (a, b)

        # Calculate distances from the new node
        for c in labels:
            if c != a and c != b:
                new_distance = (
                    distances[(a, c)]
                    + distances[(b, c)]
                    - distances[(a, b)]
                ) / 2

                distances[(new_node, c)] = new_distance
                distances[(c, new_node)] = new_distance

        # Remove a and b
        labels.remove(a)
        labels.remove(b)

        # Add the new node
        labels.append(new_node)

    # Connect the final two nodes
    root = "Root"
    children[root] = (labels[0], labels[1])

    return children, root
"""

def upgma(names, matrix):
    clusters = list(names)

    distances = {}

    for i in range(len(clusters)):
        for j in range(len(clusters)):
            distances[(clusters[i], clusters[j])] = matrix[i][j]

    # Number of original genomes in each cluster
    sizes = {}

    for name in names:
        sizes[name] = 1

    children = {}
    next_node = 1

    while len(clusters) > 1:

        # Find the two closest clusters
        best_pair = None
        best_distance = float("inf")

        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                a = clusters[i]
                b = clusters[j]

                if distances[(a, b)] < best_distance:
                    best_distance = distances[(a, b)]
                    best_pair = (a, b)

        a, b = best_pair

        # Create a new internal node
        new_node = "Node" + str(next_node)
        next_node += 1

        children[new_node] = (a, b)

        # Size of the new cluster
        new_size = sizes[a] + sizes[b]
        sizes[new_node] = new_size

        # Calculate distances from the new cluster
        for c in clusters:

            if c != a and c != b:

                new_distance = (
                    sizes[a] * distances[(a, c)]
                    + sizes[b] * distances[(b, c)]
                ) / new_size

                distances[(new_node, c)] = new_distance
                distances[(c, new_node)] = new_distance

        # Remove old clusters
        clusters.remove(a)
        clusters.remove(b)

        # Add the new cluster
        clusters.append(new_node)

    root = clusters[0]

    return children, root

def print_tree(node, children, level=0):
    if node not in children:
        print("  " * level + node)
        return

    print("  " * level + node)

    left, right = children[node]

    print_tree(left, children, level + 1)
    print_tree(right, children, level + 1)
    
children, root = upgma(names, matrix)

after_tree = time.perf_counter()
print(f"UPGMA time: {after_tree - after_distances:.3f} s")

print("\nEvolutionary tree:")

print_tree(root, children)

# Stop measuring computational time
elapsed_time = time.perf_counter() - start_time

#print(f"\nElapsed time: {elapsed_time:.3f} seconds")
print(f"Total computation time: {after_tree - start_time:.3f} s")