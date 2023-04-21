__author__ = 'Lars Berling'

from Bio import Phylo
import io
import re as re
from apply_new_mapping import get_mapping_dict

# Function to transform a given nwk-tree string to a ranked-tree in cluster notation

dellen = re.compile('\:[0-9]+\.[0-9]+(e-[0-9]+)?')


def sort_tree(tree):
    # Returns unique string representation of tree by sorting taxa and converting sets to trees.
    output = []
    for i in tree:
        j = sorted(list(i))
        output.append(set(j))
    return output


def add_leaves(tree):
    # Add sets containing just one leaf to make ultrametric trees become non-ultrametric RNNI trees
    # needed to use RNNI moves as they are implemented for RNNI and not for uRNNI
    tree_with_leaves = []
    n = len(tree) + 1
    for j in range(1,n+1):
        # tree_with_leaves.append(set([str(j)]))
        tree_with_leaves.append(set([j]))
    for k in tree:
        tree_with_leaves.append(set(k))
    return tree_with_leaves


def expand(tree):
    # Takes a tree in list of lists format, adds leave clusters and returns it in list of sets format
    if (len(tree[0]) == 1):
        return (sort_tree(tree))
    tree = add_leaves(tree)
    tree = sort_tree(tree)
    return tree


def name_clade(clade, tree):
    # Function that names the clade in a given tree, part of a recursion

    """
    clade is a given node in the tree
    clade.caldes[0] and clade.clades[1] are the two children of the clade
    tree.distance() returns the distance of two given nodes, if only on is given it returns distance to the root

    Description of the if statements
    l.27 & l.58 : which child of the clade is further away from the root --> will get the smaller rank
    l.30 & l.63 : case that both children are leaves
    l.35 & l.67 : the child node clade[0] is a leaf
    l.42 & l.75 : the child node clade[1] is a leaf
    l.50 & l.83 : both child nodes are subtrees with an already set clade notation
    """

    # if (clade.clades[0].branch_length > clade.clades[1].branch_length):
    if tree.distance(clade.clades[0]) > tree.distance(clade.clades[1]):
        # if len(clade.clades[0].name) == 1 and len(clade.clades[1].name) == 1:
        if not '{' in clade.clades[0].name and not '{' in clade.clades[1].name:
            # clade.branch_length += clade.clades[0].branch_length
            clade.name = '{' + clade.clades[0].name + ',' + clade.clades[1].name + '}:' + str(tree.distance(clade))

        elif not '{' in clade.clades[0].name:
            big1 = clade.clades[1].name.split(sep=',{')
            big1 = big1[len(big1) - 1].replace('}', '').replace('{', '')
            big1 = re.sub(dellen, '', big1)
            # clade.branch_length += clade.clades[0].branch_length
            clade.name = clade.clades[1].name + ',{' + clade.clades[0].name + ',' + big1 + '}:' + str(
                tree.distance(clade))

        elif not '{' in clade.clades[1].name:
            big0 = clade.clades[0].name.split(sep=',{')
            big0 = big0[len(big0) - 1].replace('}', '').replace('{', '')
            big0 = re.sub(dellen, '', big0)
            # clade.branch_length += clade.clades[0].branch_length
            clade.name = clade.clades[0].name + ',{' + big0 + ',' + clade.clades[1].name + '}:' + str(
                tree.distance(clade))

        else:
            big0 = clade.clades[0].name.split(sep=',{')
            big0 = big0[len(big0) - 1].replace('}', '').replace('{', '')
            big0 = re.sub(dellen, '', big0)
            big1 = clade.clades[1].name.split(sep=',{')
            big1 = big1[len(big1) - 1].replace('}', '').replace('{', '')
            big1 = re.sub(dellen, '', big1)
            # clade.branch_length += clade.clades[0].branch_length
            clade.name = clade.clades[0].name + ',' + clade.clades[1].name + ',{' + big0 + ',' + big1 + '}:' + str(
                tree.distance(clade))
            # print(big0, '--2--', big1)

    else:
        if not '{' in clade.clades[0].name and not '{' in clade.clades[1].name:
            # clade.branch_length += clade.clades[1].branch_length
            clade.name = '{' + clade.clades[0].name + ',' + clade.clades[1].name + '}:' + str(tree.distance(clade))

        elif not '{' in clade.clades[0].name:
            # redundant ?
            big1 = clade.clades[1].name.split(sep=',{')
            big1 = big1[len(big1) - 1].replace('}', '').replace('{', '')
            big1 = re.sub(dellen, '', big1)
            # clade.branch_length += clade.clades[1].branch_length
            clade.name = clade.clades[1].name + ',{' + clade.clades[0].name + ',' + big1 + '}:' + str(
                tree.distance(clade))

        elif not '{' in clade.clades[1].name:
            # redundant ?
            big0 = clade.clades[0].name.split(sep=',{')
            big0 = big0[len(big0) - 1].replace('}', '').replace('{', '')
            big0 = re.sub(dellen, '', big0)
            # clade.branch_length += clade.clades[1].branch_length
            clade.name = clade.clades[0].name + ',{' + big0 + ',' + clade.clades[1].name + '}:' + str(
                tree.distance(clade))

        else:
            big0 = clade.clades[0].name.split(sep=',{')
            big0 = big0[len(big0) - 1].replace('}', '').replace('{', '').replace(':.*$', '')
            big0 = re.sub(dellen, '', big0)
            big1 = clade.clades[1].name.split(sep=',{')
            big1 = big1[len(big1) - 1].replace('}', '').replace('{', '')
            big1 = re.sub(dellen, '', big1)
            # clade.branch_length += clade.clades[1].branch_length
            clade.name = clade.clades[1].name + ',' + clade.clades[0].name + ',{' + big0 + ',' + big1 + '}:' + str(
                tree.distance(clade))
            # print(big0,'--1--', big1)


def rec_nameing(clade, tree):
    # Function that recursively names the inner nodes of a tree with the cluster notation, starting top down from the given clade
    if not clade.clades[0].name:
        if not clade.clades[1].name:
            rec_nameing(clade.clades[0], tree)
            rec_nameing(clade.clades[1], tree)
            name_clade(clade, tree)
        else:
            rec_nameing(clade.clades[0], tree)
            name_clade(clade, tree)
    else:
        if not clade.clades[1].name:
            rec_nameing(clade.clades[1], tree)
            name_clade(clade, tree)
        else:
            # both clades have a name
            name_clade(clade, tree)


def to_cluster(raw_cluster):
    # Function that converts the output from the rec_naming function and sorts the clusters accordingly
    out = {}
    next = raw_cluster.split(sep=',{')
    for i in next:
        key, value = i.split(sep=':')
        if not key[0] == '{':
            key = '{' + key
        out[key] = float(value)
    sort = sorted(out.items(), key=lambda x: x[1], reverse=True)
    clust = ''
    for key in sort:
        clust += key[0] + ','
    # Need to remove the last ,
    clust = clust[:-1]
    clust = '[' + clust + ']'
    return (clust)


def nwk_to_cluster(nwk):
    # Function that converts a tree nwk string to its cluster representation
    tree_str = io.StringIO(nwk)
    tree = Phylo.read(tree_str, 'newick')
    tree.root.branch_length = float(0)
    rec_nameing(tree.root, tree)
    clus = to_cluster(tree.root.name)
    return (clus)


# Function only for external input, like BEAST output
def nwk_parser(in_file, out_file, re_tree=re.compile('.')):
    # Function that converts in_file of nwk trees to out_file cluster trees
    # re_tree is an optional tree identification string, i.e. tree lines from BEAST output start with 'tree ...'
    # re_tree Default takes each line as a tree!

    # ntaxa = re.compile('\tDimensions ntax=[0-9]+\;$')
    trees = 0
    out_f = open(out_file, 'w+')
    # out_f.write('taxa\ntrees\n') # First two lines for later usage of number of trees and taxa
    with open(in_file) as in_f:
        for cnt, line in enumerate(in_f):
            if re_tree.match(line):
                out_f.write(nwk_to_cluster(line) + '\n')
                trees += 1
            # elif ntaxa.match(line):
            # n = re.compile('[0-9]+').search(line).group()
            # else:
            # print(cnt)
    # lines = out_f.readline()
    # lines[0] = n
    # lines[1] = str(trees)
    out_f.close()


class Node:
    def __init__(self, clade, rank):
        self.left = None
        self.right = None
        self.rank = rank
        self.clade = clade


class Leaf:
    def __init__(self, taxa):
        self.taxa = taxa


def cluster_to_tree(tree):
    working_treelist = []
    rank = 1
    for cluster in tree:
        if len(cluster) == 2:
            # Internal Node that will be formed from two leafs
            new_cherry = Node(cluster.copy(), rank)
            new_cherry.left = Leaf(cluster.pop())
            new_cherry.right = Leaf(cluster.pop())
            working_treelist.append(new_cherry)
        else:
            # Internal Node that will be formed from two subtrees, one may be a leaf
            to_merge = []
            for index, node in enumerate(working_treelist):
                if not cluster.isdisjoint(node.clade):
                    to_merge.append(index)
                    if len(cluster) - len(node.clade) == 1:
                        break
            if len(to_merge) == 1:
                # Add leaf
                new_tree = Node(cluster.copy(), rank)
                new_tree.left = working_treelist.pop(to_merge[0])
                new_tree.right = Leaf(cluster.difference(new_tree.left.clade).pop())
                working_treelist.append(new_tree)
            elif len(to_merge) == 2:
                # Merge two trees
                new_tree = Node(cluster.copy(), rank)
                new_tree.left = working_treelist.pop(to_merge[0])
                new_tree.right = working_treelist.pop(to_merge[1] - 1)
                working_treelist.append(new_tree)
            else:
                # ERROR
                print('ERROR!')
        rank += 1
    return working_treelist[0]


# Ways to check if either of the children are leafs
# isinstance(node.left, set)
# isinstance(node.right, set)


def print_tree_from_root(root):
    """Recursive Function to get the nwk representation of a node based tree"""
    if isinstance(root.left, Leaf) and isinstance(root.right, Leaf):
        # Merge two leafs
        return f'({root.left.taxa}:{root.rank},{root.right.taxa}:{root.rank})'
    elif isinstance(root.left, Leaf) and not isinstance(root.right, Leaf):
        # Merge left leaf and right subtree
        return f'({root.left.taxa}:{root.rank},{print_tree_from_root(root.right)}:{root.rank - root.right.rank})'
    elif not isinstance(root.left, Leaf) and isinstance(root.right, Leaf):
        # Merge left subtree and right leaf
        return f'({print_tree_from_root(root.left)}:{root.rank - root.left.rank},{root.right.taxa}:{root.rank})'
    else:
        # Merge two subtrees
        return f'({print_tree_from_root(root.left)}:{root.rank - root.left.rank},{print_tree_from_root(root.right)}:{root.rank - root.right.rank})'


def write_tree_to_nwk(tree, out_file, map_file):
    mapping = get_mapping_dict(map_file)
    cur_t = cluster_to_tree(tree)
    with open(out_file, 'w+') as f:
        # Nexus Header
        f.write("#NEXUS\n\nBEGIN TREES;\n\tTRANSLATE\n")
        for i in sorted(mapping.keys()):
            f.write(f"\t\t{i} {mapping[i]}")
            if i is not len(mapping.keys()):
                f.write(',')
            f.write('\n')
        f.write(";\n")

        f.write('TREE tree1 = ')
        f.write(print_tree_from_root(cur_t))

        f.write(';\nEND;\n')
    return 0
