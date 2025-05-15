# Import needed modules
import os
import sys
import glob
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo


#Store the path to directory containing fasta file inputs
in_dir = "/shared/forsythe/BB485/Week06/Brass_CDS_seqs/"

#Store the path to a new folder where we want outputs written
out_dir = "/scratch/mendozde/BB485/Week06/"

all_files = glob.glob(in_dir + "*fasta")
#print(all_files)

# start a for-loop through each file name, inside the loop
# create a mafft command for the file of interest
# call that mafft command

for file in all_files:
    #print(file)
    new_file_path = file.replace(in_dir, out_dir)
    print(new_file_path)
    
    aln_cmd = f'mafft --auto --quiet {file} > {new_file_path}'
    print(aln_cmd)
    os.system(aln_cmd)
    
# Create a for loop to run iqtree

all_aln_files = glob.glob(out_dir + "*fasta")

# Loop through align files and loop through iqtree for each, see example code and alter for needs


# Read in the trees and test the topologies
all_aln_files = glob.glob(out_dir + "*fasta")

for aln in all_aln_files:
    tree_command = f"iqtree -s {aln} -m TEST -nt 2"
    print(tree_command)
    os.system(tree_command)  # Uncomment this when ready to run

# Step 3: Read and root trees, determine topology
topology_counts = {"12top": 0, "13top": 0, "23top": 0, "Unknown": 0}

all_tree_files = glob.glob(out_dir + "*.treefile")

for tree in all_tree_files:
    try:
        temp_tree = Phylo.read(tree, "newick")

        # Find Es_ tip (outgroup)
        for tip in temp_tree.get_terminals():
            if "Es_" in tip.name:
                es_tip = tip
                break

        temp_tree.root_with_outgroup(es_tip)
        all_terminal_branches = temp_tree.get_terminals()

        # Store tips
        Bs_temp = Cr_temp = At_temp = out_temp = None

        for t in all_terminal_branches:
            if "Bs_" in t.name:
                Bs_temp = t 
            elif "Cr_" in t.name:
                Cr_temp = t
            elif "At_" in t.name:
                At_temp = t
            else:
                out_temp = t

        # Define topology groups
        P1_and_P2 = [Bs_temp, Cr_temp]
        P1_and_P3 = [Bs_temp, At_temp]
        P2_and_P3 = [Cr_temp, At_temp]

        # Check monophyly
        if temp_tree.is_monophyletic(P1_and_P2):
            topo_str = "12top"
        elif temp_tree.is_monophyletic(P1_and_P3):
            topo_str = "13top"
        elif temp_tree.is_monophyletic(P2_and_P3):
            topo_str = "23top"
        else:
            topo_str = "Unknown"

        topology_counts[topo_str] += 1
        print(f"{tree} → {topo_str}")
    except Exception as e:
        print(f"Error processing {tree}: {e}")
        topology_counts["Unknown"] += 1

# Print final topology counts
print("\nTopology Counts:")
for topo, count in topology_counts.items():
    print(f"{topo}: {count}")