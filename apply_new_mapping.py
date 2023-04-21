__author__ = 'Lars Berling'

import re


# This code will overwrite the file

def get_mapping_dict(file):
    # Extract a mapping dict from a file dict: int --> taxa
    # Begin trees;
    #   Translate

    # ;
    # trees
    # End;
    begin_map = re.compile('\tTranslate\n', re.I)
    end = re.compile('\t*?;\n?')

    mapping = {}

    begin = False
    with open(file) as f:
        for line in f:
            if begin:
                if end.match(line):
                    break
                split = line.split()

                mapping[int(split[0])] = split[1][:-1] if split[1][-1] == "," else split[1]

            if begin_map.match(line):
                begin = True
    return mapping


def apply_new_mapping(file_to_apply, file_extract_mapping, output):
    """
    Remaps the taxa according to another mapping

    :param file_to_apply: containing tree that will be remapped
    :param file_extract_mapping: containing the mapping that will be used
    :param output: output file that will contain the tree with new mapping
    :return: Writes a new .tree file with different mapping
    """
    # Will apply a given mapping to a file
    new_mapping = get_mapping_dict(file_extract_mapping)
    old_mapping = get_mapping_dict(file_to_apply)
    new_mapping_reversed = {v: k for (k, v) in new_mapping.items()}

    re_tree = re.compile('\t?tree .*$', re.I)
    re_taxa = re.compile('([0-9]+)(\[|:)')

    # writing the Nexus header with the new mapping before the trees
    out = open(output, "w+")
    out.write("#Nexus\n\nBegin trees;\n\tTranslate\n")
    for i in sorted(new_mapping.keys()):
        out.write(f"\t\t{i} {new_mapping[i]}")
        if i is not len(new_mapping.keys()):
            out.write(',')
        out.write('\n')
    out.write(";\n")
    # The following will remap the trees
    with open(file_to_apply) as f:
        for line in f:
            if re_tree.match(line):
                out.write(re_taxa.sub(lambda m: m.group().replace(m.group(1),
                                                                  str(new_mapping_reversed[
                                                                          old_mapping[int(m.group(1))]]))
                                      , line))
                # replaces the m.group(1) which is the taxon number with the new maping taxon

    out.write("End;\n")
    out.close()
    return 0
