from tetres.trees.time_trees import TimeTreeSet
from tetres.summary.centroid import Centroid
from tetres.summary.annotate_centroid import annotate_centroid
import os


if __name__ == "__main__":
	# Given the path to your treefile, in this case it is in the same directory as this script
	treefile = "DS1.trees"
	# Creating a TimeTreeSet from your trees
	mytts = TimeTreeSet(treefile)
	# Printing the number of trees that were loaded
	print(f"Loaded {len(mytts)} trees from {treefile}.")
	# Printing the number of leafs of the first tree in this set
	print(f"First tree has {len(mytts[0])} leafs.")
	# Creating a Centroid object, restricting it to 2 threads
	mycen = Centroid(n_cores=2)
	# Setting up the number of repetitions, here it is set to 1
	# If you increase the number you want to pick the tree with lowest sos as given in the filename
	repetitions = 1
	for _ in range(repetitions):
		# Compute the centroid of the trees
		cen, sos = mycen.compute_centroid(mytts)
		# If a centroid with this sum of squared distances has not been found yet, save it to a file
		if not os.path.exists(f"centroid_{sos}.tree"):
			# Writing the centroid to a nexus file, file_name defines where this new file will be created
			# name is the name of the centroid tree in the nexus file itself
			cen.write_nexus(mytts.map, file_name=f"centroid_{sos}.tree", name=f'centroid{sos}')
		else:
			print("Already found centroid with same SoS value.")
