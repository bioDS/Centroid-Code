__author__ = 'Lars Berling'

import os
import re
import time

from subprocess import Popen
from rpy2 import robjects
from treeoclock.summary.annotate_centroid import annotate_centroid
from treeoclock.summary.centroid import Centroid
from treeoclock.trees.time_trees import TimeTreeSet

global re_tree
re_tree = re.compile('\t?tree .*$', re.I)
global work_folder
work_folder = f"{os.path.dirname(os.path.realpath(__file__))}/Simulations"


def sim_tree(n, working_dir):
    # Simulates tree via the sim.taxa function in R (only internally since we do not need this tree)
    # Convert that tree to a ranked tree (conv_true_*.tree)
    # This will be converted to a nexus tree (true_*.tree) which will be used for alignment simulation

    robjects.r(f'''
        library(maps)
        library(ape)
        library(phytools)
        library(e1071)
        # tree <- rcoal({n}, tip.label = c(1:{n}), br = rdiscrete({n - 1}, probs = range(1:{n - 1})))
        # tree <- rcoal({n}, tip.label = c(1:{n}), br = 'coalescent')
        tree <- rcoal({n}, tip.label = c(1:{n}))
        tree <- as.multiPhylo(tree)
        names(tree) <- c('TrueTree')
        write.nexus(tree, file ="{work_folder}/{working_dir}/true_{working_dir}.tree")
    ''')
    # Generate the Cluster format of the sim-taxa tree


def sim_alignment(n, l=800, r=0.005, nbr=''):
    # simulates a tree via the sim_tree function and then simulates an alignment along that tree via the simSeq funciton
    # from phangorn R package
    # Writes the output to a fasta file

    working_dir = f"{n}_{l}_{str(r - int(r)).split('.')[1]}"
    if nbr != '':
        working_dir += f'_{nbr}'

    try:
        os.mkdir(f'{work_folder}/{working_dir}')
    except FileExistsError:
        pass

    # Generating a tree and converting it so that the origin is a ranked tree ()
    sim_tree(n, working_dir)

    # Calling R phangorn to simulate a sequence along the true tree
    robjects.r(f'''
        library(phangorn)
        tree <- read.nexus("{work_folder}/{working_dir}/true_{working_dir}.tree")  
        s <- simSeq(tree, l = {l}, rate = {r})
        write.phyDat(s, "{work_folder}/{working_dir}/{working_dir}.fas", format = "fasta")
    ''')


def generate_beast_xml(n, l, r, nbr=''):
    # Function will generate a BEAST2 xml input file using the R beautier library
    # Additionally changing the substitution rate with a python regex afterwards

    working_dir = f"{n}_{l}_{str(r - int(r)).split('.')[1]}"
    if nbr != '':
        working_dir += f'_{nbr}'

    robjects.r(f'''
        library(beautier)

        create_beast2_input_file(input_filename = "{work_folder}/{working_dir}/{working_dir}.fas", 
                         output_filename = "{work_folder}/{working_dir}/{working_dir}.xml",
                         site_model = create_jc69_site_model(), 
                         clock_model = create_strict_clock_model(), 
                         tree_prior = create_ccp_tree_prior(), 
                         mcmc = create_mcmc(chain_length = 2000000,
                                            store_every = 2000,
                                            pre_burnin = 500000,
                                            treelog = create_treelog(log_every = 2000)),
                         beauti_options = create_beauti_options(beast2_version = "2.6"))
    ''')
    rate_regex = re.compile('(?P<Prefix>name="mutationRate">)(?P<rate>\d+\.\d+)(?P<Suffix></parameter>)')
    with open(f'{work_folder}/{working_dir}/{working_dir}.xml',
              'r+') as f:
        content = f.read()
        content = rate_regex.sub(f"\g<Prefix>{r}\g<Suffix>", content)
        f.seek(0)
        f.write(content)
        f.truncate()


def run_beast(n, l, r, nbr=''):
    # Will execute Beast on the xml file for the given simulation via the beastier R package

    working_dir = f"{n}_{l}_{str(r - int(r)).split('.')[1]}"
    if os.path.exists(f'/Applications/BEAST 2.6.3/bin/beast'):
        beast = f'/Applications/BEAST\ 2.6.3/bin/beast'  # Mac OS path
    else:
        beast = f'/home/lars/.local/share/beast/bin/beast'  # Dora path
    if nbr != '':
        working_dir += f'_{nbr}'
    
    cwd = os.getcwd()  # Getting current Working directory
    os.chdir(f"{work_folder}/{working_dir}/")
    os.system(f"{beast} -beagle_cpu -overwrite -threads -1 {working_dir}.xml > beast_output.txt 2>&1")
    os.chdir(cwd)  # Resetting the working directory


def run_treeannotator(n, l, r, nbr='', heights="ca"):
    # /Applications/BEAST\ 2.6.3/bin/treeannotator 10_800_004_0/10_800_004_0.trees MCC_1.tree
    # Runs the tree annotator on the given file producing the default MCC output

    if heights not in ["mean", "median", "ca"]:
        raise ValueError(f"Unrecognized heights option for MCC tree! {heights}")

    working_dir = f"{n}_{l}_{str(r - int(r)).split('.')[1]}"
    if nbr != '':
        working_dir += f'_{nbr}'

    if not os.path.exists(
            f'{work_folder}/{working_dir}/MCC_{heights}.tree'):

        if os.path.exists(f'/Applications/BEAST 2.6.3/bin/treeannotator'):
            path = f'/Applications/BEAST\ 2.6.3/bin/treeannotator'  # Mac OS path
        else:
            path = f'/home/lars/.local/share/beast/bin/treeannotator'  # Dora path

        Popen(f'{path}'
              f' -heights {heights}'
              f' {work_folder}/{working_dir}/{working_dir}.trees'
              f' {work_folder}/{working_dir}/MCC_{heights}.tree',
              shell=True)
        # Waiting for treeannotator to finish calculating the MCC tree, otherwise error since next part is faster
        while not os.path.exists(
                f'{work_folder}/{working_dir}/MCC_{heights}.tree'):
            time.sleep(1)


def run_centroid(n, l, r, nbr):
    # Converts the set of trees to cluster format and runs the centroid algorithm (increasing subsample variation)
    working_dir = f"{n}_{l}_{str(r - int(r)).split('.')[1]}"
    if nbr != '':
        working_dir += f'_{nbr}'

    try:
        myts = TimeTreeSet(f'{work_folder}/{working_dir}/{working_dir}.trees')
        mycen = Centroid(start="FM", variation="inc_sub")
        cen, sos = mycen.compute_centroid(myts)  # computing the centroid
        cen = annotate_centroid(cen, myts)  # annotation of the centroid with branch lengths
        # writing the sos of the centroid computed to a file
        file = open(f"{work_folder}/{working_dir}/centroid_sos.log", "w+")
        file.write(f"Sum of squared RNNI distances of centroid: {sos}\n")
        file.close()
        cen.write_nexus(myts.map, file_name=f"{work_folder}/{working_dir}/centroid.tree", name=f"Cen_{working_dir}")
    except FileExistsError:
        pass

def run_simulation(n, l, r, n_simulations=1):

    # creating the work_folder if not existing    
    try:
        os.mkdir(f"{work_folder}")
    except FileExistsError:
        pass

    # combines all functions above to create a folder containing the true tree, the beast output .trees and a MCC tree
    sim = 0
    for cur in range(n_simulations):
        working_dir = f"{n}_{l}_{str(r - int(r)).split('.')[1]}_{sim}"
        while os.path.exists(f'{work_folder}/{working_dir}/'):
            sim += 1
            working_dir = f"{n}_{l}_{str(r - int(r)).split('.')[1]}_{sim}"

        sim_alignment(n, l, r, sim)
        generate_beast_xml(n, l, r, sim)
        run_beast(n, l, r, sim)
        run_treeannotator(n, l, r, sim, heights="median")
        run_treeannotator(n, l, r, sim, heights="mean")
        run_treeannotator(n, l, r, sim, heights="ca")
        run_centroid(n, l, r, sim)

        print(f"Finished {cur + 1} / {n_simulations} simulations")
        sim += 1


if __name__ == '__main__':
    run_simulation(n=10, l=800, r=0.005, n_simulations=1)

