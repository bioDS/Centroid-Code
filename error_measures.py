__author__ = "Lars Berling"

from Bio import Phylo
import re
import json
from rpy2.robjects.packages import importr, isinstalled
from rpy2.robjects.vectors import StrVector

"""
Defining Error Measures for the comparison of two trees
See the 'Looking for trees in the forest paper'
"""


def r_install_decorator(func):
    packnames = ("ape", "phangorn", "phytools", "beautier")
    names_to_install = [x for x in packnames if not isinstalled(name=x)]
    if any(names_to_install):
        utils = importr("utils")
        utils.install_packages(StrVector(names_to_install))
    return func


def rank_error(tree_reference, tree):
    # Calculates the sum of absolute differences between all clade ranks from tree to the tree_reference
    # less is better
    sum = 0
    for i, clade in enumerate(tree_reference):
        # i=height-1
        if not clade == tree[i]:
            # if clades at same age are different
            if clade in tree:
                # if the clade of tree_reference is also in tree
                sum += abs((i+1) - tree.index(clade))
            else:
                # clade is not in tree_reference, need to find MRCA
                for x in tree:
                    if clade.issubset(x):
                        mrca = x
                        break
                mrca_age = tree.index(mrca) + 1
                sum += abs((i+1) - mrca_age)
    return sum


def ca_error(T_reference, T):
    # CAE
    # Calculates the heights of all clades and then compares heights in tree1 to heights in tree2 as reference
    # Not really heights but rather the distance to the root!
    tree_ref = Phylo.read(T_reference, 'nexus')
    tree = Phylo.read(T, 'nexus')

    # Compare the clade heights
    rec_nameing(tree_ref.root, tree_ref)
    rec_nameing(tree.root, tree)

    t_ref_list = tree_ref.root.name.split(",{")
    t_ref_dict = {}
    for i in t_ref_list:
        k, v = i.split(":")
        # Delete all non digits characters from k and sort the integers in the cluster
        k = sorted(re.sub("[{}]", "", k).split(','))
        t_ref_dict[json.dumps(k)] = float(v)

    t_list = tree.root.name.split(",{")
    t_dict = {}
    for i in t_list:
        k, v = i.split(":")
        k = sorted(re.sub("[{}]", "", k).split(','))
        t_dict[json.dumps(k)] = float(v)
    error = 0
    for k, v in t_ref_dict.items():
        if k in t_dict.keys():
            error += abs(v - t_dict[k])
        else:
            t_clades = [set(json.loads(k)) for k in t_dict.keys()]
            cur_t_ref_clade = set(json.loads(k))
            mrca_in_t = t_dict[
                json.dumps(sorted(list(min([i for i in t_clades if cur_t_ref_clade.issubset(i)], key=len))))]
            error += abs(v - mrca_in_t)
    return error


def clades_missed_error(tree, tree_reference):
    # number of clades that are present in tree2 but not in tree1
    return len([x for x in tree_reference if x not in tree])


def clades_called_error(tree, tree_reference):
    # paper definition does not make sense
    # i will implement the text i.e.
    # -1 for clades that are called incorrectly, meaning present in tree1 that are not in tree2
    # +1 for clades that are called correctly, meaning present in both tree1 and tree2
    incorrect_clades = clades_missed_error(tree_reference, tree)
    correct_clades = len([x for x in tree if x in tree_reference])
    return correct_clades - incorrect_clades


def clade_credibility_score(tree, treeset):
    # product of clade frequencies, should be max for the MCC tree
    score = 1
    for clade in tree:
        clade_sum = 0
        for t in treeset:
            if clade in t:
                clade_sum += 1
        score *= (clade_sum / len(treeset))
    return score


@r_install_decorator
def log_likelihood(tree_file, alignment_file):
    """
    Computes the log likelihood of a given tree and alignment pair via the R package phangorn

    :param tree_file: Path to tree file in nexus format
    :param alignment_file: Path to alignment in fasta format
    :return: Log likelihood of the tree given the alignment, float
    """

    ape = importr('ape')
    phangorn = importr('phangorn')

    alignment = phangorn.read_phyDat(alignment_file, format="fasta")
    tree = ape.read_nexus(tree_file)
    likelihood = phangorn.pml(tree, alignment, rate=0.005)
    return likelihood.rx2("logLik")[0]

@r_install_decorator
def tree_distances(tree1, tree2):
    """
    Computes different distances between tree1 and tree2, returns a dict with keys name of distances and value distances

    :param tree1: file for tree1 NEXUS format
    :param tree2: file for tree2 NEXUS format
    :return: dictionary with distances
    """

    out = {}

    ape = importr('ape')
    phangorn = importr('phangorn')

    t1 = ape.read_nexus(tree1)
    t2 = ape.read_nexus(tree2)

    distances = phangorn.treedist(t1, t2)
    out['RF'] = distances[0]  # symmetric difference -- RF distance
    out['KF'] = distances[1]  # branch score difference -- Kuhner Felsenstein
    out['PD'] = distances[2]  # path difference
    out['wPD'] = distances[3]  # quadratic path difference

    wRF = phangorn.wRF_dist(t1, t2)
    out['wRF'] = wRF[0]

    nRF = phangorn.RF_dist(t1, t2, normalize=True)
    out['nRF'] = nRF[0]

    nwRF = phangorn.wRF_dist(t1, t2, normalize=True)
    out['nwRF'] = nwRF[0]

    return out
