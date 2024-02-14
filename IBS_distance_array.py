import numpy as np
import argparse

# Step 1: Parse PLINK genome results file
def parse_plink_genome_file(genome_file):
    ibs_distances = []

    with open(genome_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            fields = line.strip('n')
            individual1 = fields[4:18]
            individual2 = fields[22:36]
            ibs_distance = float(fields[81:91])  # Assuming IBS distance is in column 10 (0-indexed)
            ibs_distances.append((individual1, individual2, ibs_distance))

    return ibs_distances

# Step 2: Compute IBS distance array
def compute_ibs_distance_array(ibs_distances):
    individuals = sorted(set([x[0] for x in ibs_distances]))
    num_individuals = len(individuals)
    ibs_distance_array = np.zeros((num_individuals, num_individuals))

    # Populate the distance array
    for i, ind1 in enumerate(individuals):
        for j, ind2 in enumerate(individuals):
            if i != j:
                # Find the IBS distance between the pair (ind1, ind2)
                pair_distance = next(d[2] for d in ibs_distances if (d[0] == ind1 and d[1] == ind2) or (d[0] == ind2 and d[1] == ind1))
                ibs_distance_array[i, j] = pair_distance
    return ibs_distance_array, individuals

# Access command-line arguments
parser = argparse.ArgumentParser(description='generate the IBS distance array')
# Add arguments
parser.add_argument('--geDist', type=str, help='genome dist file from PLINK')
parser.add_argument('--output', type=str, help='output directory and prefix')
# Parse arguments
args = parser.parse_args()
print(args.geDist)
ibs_distances=parse_plink_genome_file(args.geDist)
ibs_distance_array, individuals = compute_ibs_distance_array(ibs_distances)
np.save(args.output+".npy",np.array([ibs_distance_array, individuals]))
