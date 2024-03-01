import numpy as np
import argparse

# Step 1: Parse PLINK genome results file
def parse_grm_file(genome_file):
    grm_distances = []

    with open(genome_file, 'r') as f:
        for line in f:
            fields = line.strip("\n").split("\t")
            individual1 = fields[0]
            individual2 = fields[1]
            grm_distance = fields[3]  # Assuming IBS distance is in column 10 (0-indexed)
            grm_distances.append((individual1.strip(" "), individual2.strip(" "), grm_distance))
    return grm_distances

# Step 2: Compute GRM distance array
def compute_grm_distance_array(grm_distances):
    a=[]
    for i in grm_distances:
        a.append(i[0])
        a.append(i[1])
    individuals = sorted(set(a))
    num_individuals = len(individuals)
    grm_distance_array = np.zeros((num_individuals, num_individuals))

    # Populate the distance array
    for j in grm_distances:
        index1=int(j[0])-1
        index2=int(j[1])-1
        # Find the IBS distance between the pair (ind1, ind2)
        grm_distance_array[index1,index2] = float(j[2])
        grm_distance_array[index2,index1] = float(j[2])
    return grm_distance_array

# Access command-line arguments
parser = argparse.ArgumentParser(description='generate the GRM distance array')
# Add arguments
parser.add_argument('--geDist', type=str, help='genome dist file from PLINK')
parser.add_argument('--output', type=str, help='output directory and prefix')
# Parse arguments
args = parser.parse_args()
distances=parse_grm_file(args.geDist)
ibs_distance_array= compute_grm_distance_array(distances)
np.save(args.output+".npy",ibs_distance_array)