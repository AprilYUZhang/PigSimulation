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
            ibs_distances.append((individual1.strip(" "), individual2.strip(" "), ibs_distance))

    return ibs_distances

# Step 2: Compute IBS distance array
def compute_ibs_distance_array(ibs_distances):
    indi_dict={}
    individuals = sorted(set([x[0] for x in ibs_distances]))
    for i,j in enumerate(individuals):
        indi_dict[j]=i
    num_individuals = len(individuals)
    ibs_distance_array = np.zeros((num_individuals, num_individuals))

    # Populate the distance array
    for j in ibs_distances:
        index1=indi_dict[j[0]]
        index2=indi_dict[j[1]]
        # Find the IBS distance between the pair (ind1, ind2)
        ibs_distance_array[index1,index2] = j[2]
        ibs_distance_array[index2,index1] = j[2]
    return ibs_distance_array, individuals

# Access command-line arguments
parser = argparse.ArgumentParser(description='generate the IBS distance array')
# Add arguments
parser.add_argument('--geDist', type=str, help='genome dist file from PLINK')
parser.add_argument('--output', type=str, help='output directory and prefix')
# Parse arguments
args = parser.parse_args()
ibs_distances=parse_plink_genome_file(args.geDist)
ibs_distance_array, individuals = compute_ibs_distance_array(ibs_distances)
np.save(args.output+".npy",np.array([ibs_distance_array, individuals]))
