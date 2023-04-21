from simulationGenerator import work_folder, run_treeannotator
import error_measures as em

from tree import read_trees_from_file

from tetres.summary.compute_sos import compute_sos_mt
from tetres.trees.time_trees import TimeTreeSet

from Bio import Phylo
from scipy import stats
import sys, re, os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

def collapse(tree):
    # Takes a tree and returns list of list format without the leaf cluster
    if not len(tree[0]) == 1:
        return(tree)
    output = []
    s = len(tree[-1])
    for i in range(s,len(tree)):
        output.append(set(tree[i]))
    return output


def is_proper_simulation(folder, dir):
    if os.path.exists(f"{folder}/{dir}"):
        if os.path.exists(f"{folder}/{dir}/{dir}.trees"):
            if os.path.exists(f"{folder}/{dir}/{dir}.fas"):
                if os.path.exists(f"{folder}/{dir}/conv_centroid.tree"):
                    if os.path.exists(f"{folder}/{dir}/MCC_ca.tree"):
                        return True
    print(folder, dir, "PROBLEM!")
    return False


def list_simulations(n_greater=5, sim_greater=5):
    # Returns the n_pairs for all simulations with n > n_greater
    list = os.listdir(work_folder)
    count = {}
    for item in list:
        if os.path.isdir(f'{work_folder}/{item}'):
            n, _, _, sim = item.split('_')
            if int(n) > n_greater:
                if is_proper_simulation(work_folder, item):
                    if int(n) in count:
                        count[int(n)] += 1
                    else:
                        count[int(n)] = 1
    n_pairs = [(k, v) for k, v in count.items() if v>sim_greater]
    n_pairs.sort(key=lambda x: x[0])
    return n_pairs


def error_measures_commonancestor_summary(n_pairs, l=800, r=0.005, x_labels=["LogLik"]):
    data = []

    for n, n_simulations in n_pairs:
        sim = 0
        worked_sim = n_simulations
        while worked_sim > 0:
            working_dir = f"{n}_{l}_{str(r - int(r)).split('.')[1]}_{sim}"

            # Checking if MCC tree exists already
            if not os.path.exists(
                    f"{work_folder}/{working_dir}/MCC_ca.tree"):
                run_treeannotator(n, l, r, sim, ca=True)

            # Everything depending on the Rank topology

            centroid = collapse(read_trees_from_file(f"{work_folder}/{working_dir}/conv_centroid.tree")[0])
            true = collapse(read_trees_from_file(f"{work_folder}/{working_dir}/conv_true_{working_dir}.tree")[0])

            # if not os.path.exists(f"{work_folder}/{working_dir}/conv_MCC_ca.tree"):
            #     run_treeannotator(n, l, r, sim, ca=True)
            mcc = collapse(read_trees_from_file(f"{work_folder}/{working_dir}/conv_MCC_ca.tree")[0])

            if "CME" in x_labels:
                cen_cme = int(em.clades_missed_error(centroid, true))
                mcc_cme = int(em.clades_missed_error(mcc, true))
                data.append(["CME", f"{'Cen' if cen_cme < mcc_cme else 'MCC'}"])

            if "CCE" in x_labels:
                cen_cce = int(em.clades_called_error(centroid, true))
                mcc_cce = int(em.clades_called_error(mcc, true))
                data.append(["CCE", f"{'Cen' if cen_cce > mcc_cce else 'MCC'}"])

            if "CRE" in x_labels:
                cen_cre = int(em.rank_error(true, centroid))
                mcc_cre = int(em.rank_error(true, mcc))
                data.append(["CRE", f"{'Cen' if cen_cre < mcc_cre else 'MCC'}"])

            if "RNNI" in x_labels:
                cen_tts = TimeTreeSet(f"{work_folder}/{working_dir}/centroid.tree")
                mcc_tts = TimeTreeSet(f"{work_folder}/{working_dir}/MCC_ca.tree")
                tru_tts = TimeTreeSet(f"{work_folder}/{working_dir}/true_{working_dir}.tree")

                d_cen_tru = tru_tts[0].fp_distance(cen_tts[0])
                d_mcc_tru = tru_tts[0].fp_distance(mcc_tts[0])
                data.append(["RNNI", f"{'Cen' if d_cen_tru < d_mcc_tru else 'MCC'}"])

            if "CAE" in x_labels:
                cen_cae = em.ca_error(f"{work_folder}/{working_dir}/true_{working_dir}.tree",
                                      f"{work_folder}/{working_dir}/centroid.tree")
                mcc_cae = em.ca_error(f"{work_folder}/{working_dir}/true_{working_dir}.tree",
                                      f"{work_folder}/{working_dir}/MCC_ca.tree")
                data.append(["CAE", f"{'Cen' if cen_cae < mcc_cae else 'MCC'}"])

            if "LogLik" in x_labels:
                mcc_ll = em.log_likelihood(
                    f"{work_folder}/{working_dir}/MCC_ca.tree",
                    f"{work_folder}/{working_dir}/{working_dir}.fas")
                cen_ll = em.log_likelihood(f"{work_folder}/{working_dir}/centroid.tree",
                                               f"{work_folder}/{working_dir}/{working_dir}.fas")
                data.append(["LogLik", f"{'Cen' if cen_ll > mcc_ll else 'MCC'}"])


            mcc_dist = em.tree_distances(f"{work_folder}/{working_dir}/MCC_ca.tree",
                                         f"{work_folder}/{working_dir}/true_{working_dir}.tree")
            cen_dist = em.tree_distances(f"{work_folder}/{working_dir}/centroid.tree",
                                         f"{work_folder}/{working_dir}/true_{working_dir}.tree")
            for metric in mcc_dist.keys():
                if metric in x_labels:
                    data.append([metric,
                                 f"{'Cen' if cen_dist[metric] < mcc_dist[metric] else 'MCC'}"])

            sim += 1
            worked_sim -= 1

    data = pd.DataFrame(data, columns=["Error", "Winner"])

    sns.countplot(data=data, x="Error", hue="Winner", hue_order=['Cen', 'MCC'], order=x_labels)

    plt.xlabel("Error", fontsize=15)
    plt.ylabel("Count", fontsize=15)

    plt.tick_params(labelsize=15)
    plt.tight_layout()
    # plt.show()
    plt.savefig(f"Simulations/ErrorMeasures_summary_{n_pairs}.pdf", format="pdf", bbox_inches='tight')
    plt.clf()


def error_measures_commonancestor(n_pairs, l=800, r=0.005, x_labels=["LogLik"]):
    data = []

    for n, n_simulations in n_pairs:
        sim = 0
        worked_sim = n_simulations
        while worked_sim > 0:
            working_dir = f"{n}_{l}_{str(r - int(r)).split('.')[1]}_{sim}"

            # Checking if MCC tree exists already
            if not os.path.exists(
                    f"{work_folder}/{working_dir}/MCC_ca.tree"):
                run_treeannotator(n, l, r, sim, ca=True)

            # Everything depending on the Rank topology

            centroid = collapse(read_trees_from_file(f"{work_folder}/{working_dir}/conv_centroid.tree")[0])
            true = collapse(read_trees_from_file(f"{work_folder}/{working_dir}/conv_true_{working_dir}.tree")[0])

            if not os.path.exists(f"{work_folder}/{working_dir}/conv_MCC_ca.tree"):
                run_treeannotator(n, l, r, sim, ca=True)
            mcc = collapse(read_trees_from_file(f"{work_folder}/{working_dir}/conv_MCC_ca.tree")[0])

            if "CME" in x_labels:
                cen_cme = int(em.clades_missed_error(centroid, true))
                mcc_cme = int(em.clades_missed_error(mcc, true))
                data.append(["CME", cen_cme,  mcc_cme])


            if "CCE" in x_labels:
                cen_cce = int(em.clades_called_error(centroid, true))
                mcc_cce = int(em.clades_called_error(mcc, true))
                data.append(["CCE", cen_cce, mcc_cce])

            if "CRE" in x_labels:
                cen_cre = int(em.rank_error(true, centroid))
                mcc_cre = int(em.rank_error(true, mcc))
                data.append(["CRE", cen_cre, mcc_cre])

            if "RNNI" in x_labels:
                cen_tts = TimeTreeSet(f"{work_folder}/{working_dir}/centroid.tree")
                mcc_tts = TimeTreeSet(f"{work_folder}/{working_dir}/MCC_ca.tree")
                tru_tts = TimeTreeSet(f"{work_folder}/{working_dir}/true_{working_dir}.tree")

                d_cen_tru = tru_tts[0].fp_distance(cen_tts[0])
                d_mcc_tru = tru_tts[0].fp_distance(mcc_tts[0])
                data.append(["RNNI", d_cen_tru, d_mcc_tru])

            if "CAE" in x_labels:
                cen_cae = em.ca_error(f"{work_folder}/{working_dir}/true_{working_dir}.tree",
                                      f"{work_folder}/{working_dir}/centroid.tree")
                mcc_cae = em.ca_error(f"{work_folder}/{working_dir}/true_{working_dir}.tree",
                                      f"{work_folder}/{working_dir}/MCC_ca.tree")
                data.append(["CAE", cen_cae, mcc_cae])

            if "LogLik" in x_labels:
                mcc_ll = em.log_likelihood(
                    f"{work_folder}/{working_dir}/MCC_ca.tree",
                    f"{work_folder}/{working_dir}/{working_dir}.fas")
                cen_ll = em.log_likelihood(f"{work_folder}/{working_dir}/centroid.tree",
                                               f"{work_folder}/{working_dir}/{working_dir}.fas")
                data.append(["LogLik", f"{'Cen' if cen_ll > mcc_ll else 'MCC'}"])


            mcc_dist = em.tree_distances(f"{work_folder}/{working_dir}/MCC_ca.tree",
                                         f"{work_folder}/{working_dir}/true_{working_dir}.tree")
            cen_dist = em.tree_distances(f"{work_folder}/{working_dir}/centroid.tree",
                                         f"{work_folder}/{working_dir}/true_{working_dir}.tree")
            for metric in mcc_dist.keys():
                if metric in x_labels:
                    data.append([metric, cen_dist[metric], mcc_dist[metric]])

            sim += 1
            worked_sim -= 1

    data = pd.DataFrame(data, columns=["Error", "Centroid", "MCC"])

    fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(15, 12))

    for er, ax in zip(x_labels, axs.ravel()):
        sns.scatterplot(data=data[data["Error"] == er], x="Centroid", y="MCC", ax=ax)

        ax.axline([0, 0], [1, 1], color='red')  # Plot f(x)=x
        ax.set_title(f"{er}")
        ax.set_xlabel("")
        ax.set_ylabel("")
        # plt.axis(ax)  # Resetting the axis limits
        ax.tick_params(labelsize=13)

    fig.supylabel("MCC", fontsize=15)
    fig.supxlabel("Centroid", fontsize=15)

    plt.tight_layout()
    # plt.show()
    plt.savefig("Simulations/ErrorMeasures_xyplot-testing.pdf", format="pdf", bbox_inches='tight')
    plt.clf()


def compare_likelihood_values(n_pairs, l=800, r=0.005, t1_key = 'MCC_ca', t2_key = "centroid"):
    data = []
    for n, n_simulations in n_pairs:
        sim = 0
        worked_sim = n_simulations
        # while worked_sim > 0:
        for i in range(n_simulations):
            working_dir = f"{n}_{l}_{str(r - int(r)).split('.')[1]}_{sim}"

            t1_ll = em.log_likelihood(f"{work_folder}/{working_dir}/{t1_key}.tree",
                                      f"{work_folder}/{working_dir}/{working_dir}.fas")
            data.append([f"{t1_key}", n, t1_ll, f"{work_folder}/{working_dir}/{t1_key}.tree"])

            t2_ll = em.log_likelihood(f"{work_folder}/{working_dir}/{t2_key}.tree",
                                      f"{work_folder}/{working_dir}/{working_dir}.fas")
            data.append([f"{t2_key}", n, t2_ll, f"{work_folder}/{working_dir}/{t2_key}.tree"])

            sim += 1
            worked_sim -= 1

    data = pd.DataFrame(data, columns=['Tree', 'taxa', 'LogLik', 'file'])

    # Histogram plot of the LogLik values
    # sns.boxplot(data=data, x='taxa', y="LogLik", hue='Tree', dodge=True)

    x_list = list(data[(data["Tree"] == t1_key)]['LogLik'])
    y_list = list(data[(data["Tree"] == t2_key)]['LogLik'])

    tl1 = list(data[(data["Tree"] == t1_key)]['file'])
    tl2 = list(data[(data["Tree"] == t2_key)]['file'])

    # RNNI distances only work with a resolved tree, conversion of mean and median with negative branch lengths is pointless
    # from bic_clustering import normalized_distance
    # distances = [normalized_distance(read_trees_from_file(t1)[0], read_trees_from_file(t1)[0]) for t1, t2 in zip(tl1, tl2)]

    # distances = [em.tree_distances(t1, t2)["nwRF"] for t1, t2 in zip(tl1, tl2)]

    sns.scatterplot(x=x_list, y=y_list)  # , hue=distances)

    plt.xlabel(f'LogLik({t1_key})', fontsize=15)
    plt.ylabel(f'LogLik({t2_key})', fontsize=15)

    ax = plt.axis()  # Getting current axis limits
    plt.axline([0, 0], [1, 1], color='red')  # Plot f(x)=x
    plt.axis(ax)  # Resetting the axis limits
    
    plt.tick_params(labelsize=15)
    plt.tight_layout()
    # plt.show()
    # return 1

    plt.savefig(f'Simulations/compare_ll_{t1_key}_{t2_key}_{n_pairs}.pdf',
                dpi=400, bbox_inches='tight', format="pdf")
    plt.clf()


if __name__ == '__main__':
    dirs = list_simulations(n_greater=0, sim_greater=0)  # List of all simulations with more than n_greater taxa and more than sim_greater simulations

    compare_likelihood_values(dirs, t1_key = 'MCC_ca', t2_key = "centroid")  # Create scatter plot of likelihood values for mcc and centroid
    error_measures_commonancestor_summary(dirs, x_labels=["PD", "RF", "wRF", "RNNI", "LogLik"])  # Histogram plot for summary of given error measures
    error_measures_commonancestor(dirs, x_labels=["PD", "RF", "wRF", "RNNI", "LogLik"])  # Scatter plot for the given error measures

