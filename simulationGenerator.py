__author__ = 'Lars Berling'

import os
import re
import time

from subprocess import Popen
from tetres.summary.annotate_centroid import annotate_centroid
from tetres.summary.centroid import Centroid
from tetres.trees.time_trees import TimeTreeSet
from rpy2.robjects.packages import importr
from error_measures import r_install_decorator
from apply_new_mapping import apply_new_mapping
from nwk_parser import nwk_parser


global re_tree
re_tree = re.compile('\t?tree .*$', re.I)
global work_folder
work_folder = f"{os.path.dirname(os.path.realpath(__file__))}/Simulations"

# The following path needs to be set to the correct locaiton of BEAST/treeannotator on your machine
global BEAST_PATH
BEAST_PATH = "/bin/beast2-mcmc"
global TREEANNOTATOR_PATH
TREEANNOTATOR_PATH = "/home/lbe74/beast/bin/treeannotator"


def beast_path_checker(func):
    if not os.path.exists(BEAST_PATH):
        raise ValueError("You need to set the BEAST_PATH varaiable to BEAST2.6 on your system! (l. 25 simulationGenerator.py)")
    if not os.path.exists(TREEANNOTATOR_PATH):
        raise ValueError("You need to set the TREEANNOTATOR_PATH varaiable to a treeannotator executable on your system! (l. 27 simulationGenerator.py)")
    return func


@r_install_decorator
def sim_tree(n, working_dir):
    # Simulates tree via the sim.taxa function in R (only internally since we do not need this tree)
    # Convert that tree to a ranked tree (conv_true_*.tree)
    # This will be converted to a nexus tree (true_*.tree) which will be used for alignment simulation

    ape = importr('ape')
    phytools = importr('phytools')

    tree = phytools.as_multiPhylo(ape.rcoal(n, tip_label = [i for i in range(1, n+1)]))
    tree.names = ("TrueTree")
    ape.write_nexus(tree, file=f"{work_folder}/{working_dir}/true_{working_dir}.tree")


@r_install_decorator
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
    ape = importr('ape')
    phangorn = importr('phangorn')

    tree = ape.read_nexus(f"{work_folder}/{working_dir}/true_{working_dir}.tree")
    seq = phangorn.simSeq(tree, l=l, rate=r)
    phangorn.write_phyDat(seq,
                          f"{work_folder}/{working_dir}/{working_dir}.fas",
                          format = "fasta")


@r_install_decorator
def generate_beast_xml(n, l, r, nbr=''):
    # Function will generate a BEAST2 xml input file using the R beautier library
    # Additionally changing the substitution rate with a python regex afterwards

    working_dir = f"{n}_{l}_{str(r - int(r)).split('.')[1]}"
    if nbr != '':
        working_dir += f'_{nbr}'

    beautier = importr("beautier")
    beautier.create_beast2_input_file(input_filename = f"{work_folder}/{working_dir}/{working_dir}.fas",
                                      output_filename = f"{work_folder}/{working_dir}/{working_dir}.xml",
                                      site_model=beautier.create_jc69_site_model(),
                                      clock_model = beautier.create_strict_clock_model(),
                                      tree_prior = beautier.create_ccp_tree_prior(),
                                      mcmc = beautier.create_mcmc(chain_length = 2000000,
                                                         store_every = 2000,
                                                         pre_burnin = 500000,
                                                         treelog = beautier.create_treelog(log_every = 2000)),
                                      beauti_options = beautier.create_beauti_options(beast2_version = "2.6"))

    rate_regex = re.compile('(?P<Prefix>name="mutationRate">)(?P<rate>\d+\.\d+)(?P<Suffix></parameter>)')
    with open(f'{work_folder}/{working_dir}/{working_dir}.xml',
              'r+') as f:
        content = f.read()
        content = rate_regex.sub(f"\g<Prefix>{r}\g<Suffix>", content)
        f.seek(0)
        f.write(content)
        f.truncate()


@beast_path_checker
def run_beast(n, l, r, nbr=''):
    # Will execute Beast on the xml file for the given simulation via the beastier R package

    working_dir = f"{n}_{l}_{str(r - int(r)).split('.')[1]}"

    if not os.path.exists(BEAST_PATH):
        raise ValueError("Need to set path to BEAST!")

    if nbr != '':
        working_dir += f'_{nbr}'
    
    cwd = os.getcwd()  # Getting current Working directory
    os.chdir(f"{work_folder}/{working_dir}/")
    os.system(f"{BEAST_PATH} -beagle_cpu -overwrite {working_dir}.xml > beast_output.txt 2>&1")
    os.chdir(cwd)  # Resetting the working directory


@beast_path_checker
def run_treeannotator(n, l, r, nbr='', heights="ca"):
    # Runs the tree annotator on the given file producing the default MCC output

    if heights not in ["mean", "median", "ca", "keep"]:
        raise ValueError(f"Unrecognized heights option for MCC tree! {heights}")

    working_dir = f"{n}_{l}_{str(r - int(r)).split('.')[1]}"
    if nbr != '':
        working_dir += f'_{nbr}'

    if not os.path.exists(
            f'{work_folder}/{working_dir}/MCC_{heights}.tree'):

        if not os.path.exists(TREEANNOTATOR_PATH):
            raise ValueError("Need to set path to BEAST!")

        Popen(f'{TREEANNOTATOR_PATH}'
              f' -heights {heights}'
              f' {work_folder}/{working_dir}/{working_dir}.trees'
              f' {work_folder}/{working_dir}/MCC_{heights}.tree',
              shell=True)
        # Waiting for treeannotator to finish calculating the MCC tree, otherwise error since next part is faster
        while not os.path.exists(
                f'{work_folder}/{working_dir}/MCC_{heights}.tree'):
            time.sleep(1)
        if not os.path.exists(
                f'{work_folder}/{working_dir}/conv_MCC_{heights}.tree'):
            apply_new_mapping(
                f'{work_folder}/{working_dir}/MCC_{heights}.tree',
                f'{work_folder}/{working_dir}/{working_dir}.trees',
                f'{work_folder}/{working_dir}/MCC_{heights}_remap.tree')
            nwk_parser(
                f'{work_folder}/{working_dir}/MCC_{heights}_remap.tree',
                f'{work_folder}/{working_dir}/conv_MCC_{heights}.tree',
                re_tree)


def run_centroid(n, l, r, nbr):
    # Converts the set of trees to cluster format and runs the centroid algorithm (increasing subsample variation)
    working_dir = f"{n}_{l}_{str(r - int(r)).split('.')[1]}"
    if nbr != '':
        working_dir += f'_{nbr}'

    if not os.path.exists(f'{work_folder}/{working_dir}/conv_true_{working_dir}.tree'):
        apply_new_mapping(
            f'{work_folder}/{working_dir}/true_{working_dir}.tree',
            f'{work_folder}/{working_dir}/{working_dir}.trees',
            f'{work_folder}/{working_dir}/true_remap_{working_dir}.tree')
        nwk_parser(
            f'{work_folder}/{working_dir}/true_remap_{working_dir}.tree',
            f'{work_folder}/{working_dir}/conv_true_{working_dir}.tree',
            re_tree)

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
        with open(f"{work_folder}/{working_dir}/tmp.tree", 'w+') as f:
            f.write(cen.get_newick())
        nwk_parser(f'{work_folder}/{working_dir}/tmp.tree',
                   f'{work_folder}/{working_dir}/conv_centroid.tree')
        os.remove(f'{work_folder}/{working_dir}/tmp.tree')  # Removing the tmp.tree
    except FileExistsError:
        pass


@beast_path_checker
@r_install_decorator
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
    run_simulation(n=5, l=800, r=0.005, n_simulations=1)

