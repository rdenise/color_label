# -*- coding: utf-8 -*-


##########################################################################################
##########################################################################################
##
##								Libraries used
##
##########################################################################################
##########################################################################################

from Bio import SeqIO, AlignIO
import sys
import os
import numpy as np
import re
import matplotlib.cm as cm
import matplotlib.colors as colors

##########################################################################################
##########################################################################################
##
##								Functions
##
##########################################################################################
##########################################################################################

def create_folder(mypath):

	"""
	Created the folder that I need to store my result if it doesn't exist
	:param mypath: path where I want the folder (write at the end of the path)
	:type: string
	:return: Nothing
	"""

	try:
		os.makedirs(mypath)
	except OSError:
		pass

##########################################################################################
##########################################################################################

def create_colorstrip_itol_file(info_tab):

	"""
	Function that create a file in itol colorstrip format for the sequence in the info_tab
	with color for each kind of

	:param info_tab: table of the annotation table
	:type: numpy.ndarray
	:return: Nothing
	"""

	print("\n#################")
	print("# COLOR STRIP FILE")
	print("#################\n")

	with open(PREFIX+"_colorstrip.txt", 'w') as writing_file:
		writing_file.write("DATASET_COLORSTRIP\n")
		writing_file.write("SEPARATOR TAB\n")
		writing_file.write("DATASET_LABEL\tT2SS_T4P_Tad_Com\n")
		writing_file.write("COLOR\t#ff0000\n")
		writing_file.write("LEGEND_TITLE\tSecretion_system\n")
		writing_file.write("LEGEND_SHAPES"+"\t1"*len(DICT_COLORSTRIP)+"\n")
		writing_file.write("LEGEND_COLORS\t"+"\t".join(DICT_COLORSTRIP.values())+"\n")
		writing_file.write("LEGEND_LABELS\t"+"\t".join(DICT_COLORSTRIP.keys())+"\n")
		writing_file.write("DATA\n")

		progression=0

		for seq in info_tab :

			progression += 1
			sys.stdout.write("{:.2f}% : {}/{} sequences\r".format(progression/float(info_tab.shape[0])*100, progression,info_tab.shape[0]))
			sys.stdout.flush()

			if "generic" in seq[0] or "choice" in seq[0] or "generique" in seq[0]:
				continue
			else:
				writing_file.write(seq[0]+"\t"+DICT_COLORSTRIP[seq[-1]]+"\n")

	print("Done !")
	return

##########################################################################################
##########################################################################################

def create_binary_itol_file(info_tab):

	"""
	Function that create a file in itol color label format for the sequence in the info_tab
	and set all the verified sequence will have grey square

	:param info_tab: table of the annotation table
	:type: numpy.ndarray
	:return: Nothing
	"""

	print("\n#################")
	print("# LABEL BINARY FILE")
	print("#################\n")

	with open(PREFIX+"_labelbinary.txt", 'w') as writing_file:
		writing_file.write("DATASET_BINARY\n")
		writing_file.write("SEPARATOR TAB\n")
		writing_file.write("COLOR\t#a4a4a4\n")
		writing_file.write("DATASET_LABEL\tVerify_sequence\n")
		writing_file.write("FIELD_SHAPES\t1\n")
		writing_file.write("FIELD_LABELS\tverified\n")
		writing_file.write("LEGEND_TITLE\tExperimentaly_verified_sequences\n")
		writing_file.write("LEGEND_SHAPES\t1\n")
		writing_file.write("LEGEND_COLORS\t#a4a4a4\n")
		writing_file.write("LEGEND_LABELS\tverify\n")
		writing_file.write("DATA\n")

		progression = 0

		for seq in info_tab :

			progression += 1
			sys.stdout.write("{:.2f}% : {}/{} sequences\r".format(progression/float(info_tab.shape[0])*100, progression,info_tab.shape[0]))
			sys.stdout.flush()

			if "_V_" in seq[0] :
				writing_file.write(seq[0]+"\t1\n")
			else :
				writing_file.write(seq[0]+"\t0\n")

	print("Done !")
	return

##########################################################################################
##########################################################################################



def create_colorrange_itol_file(info_tab):

	"""
	Function that create a file in itol color range format for the sequence in the info_tab
	and with color for each kingdom or subkingdom
	:param info_tab: table of the annotation table
	:type: numpy.ndarray
	:return: Nothing
	"""

	print("\n#################")
	print("# COLOR RANGE FILE")
	print("#################\n")

	with open(PREFIX+"_colorrange.txt", 'w') as writing_file:
		writing_file.write("TREE_COLORS\n")
		writing_file.write("SEPARATOR SPACE\n")
		writing_file.write("DATA\n")

		progression = 0

		for seq in info_tab :

			progression += 1
			sys.stdout.write("{:.2f}% : {}/{} sequences\r".format(progression/float(info_tab.shape[0])*100, progression,info_tab.shape[0]))
			sys.stdout.flush()

			lineage = seq[-2]

			if lineage in DICT_COLORRANGE :
				writing_file.write(seq[0]+" range "+DICT_COLORRANGE[lineage]+" "+lineage+"\n")
			elif seq[-3]=="Archaea":
				writing_file.write(seq[0]+" range "+DICT_COLORRANGE["Archaea"]+" Archaea\n")
			else :
				writing_file.write(seq[0]+" range "+DICT_COLORRANGE["Bacteria"]+" Bacteria\n")

	print("Done !")
	return

##########################################################################################
##########################################################################################



def create_labels_itol_file(info_tab):


	"""
	Function that let the possibility to change the name of the leafs' label

	:param info_tab: table of the annotation table
	:type: numpy.ndarray
	:return: Nothing
	"""

	print("\n#################")
	print("# LABELS FILE")
	print("#################\n")

	with open(PREFIX+"_id_label.txt", 'w') as writing_file:
		writing_file.write("LABELS\n")
		writing_file.write("SEPARATOR TAB\n")
		writing_file.write("DATA\n")

		progression=0

		for seq in info_tab :

			progression += 1
			sys.stdout.write("{:.2f}% : {}/{} sequences\r".format(progression/float(info_tab.shape[0])*100, progression,info_tab.shape[0]))
			sys.stdout.flush()

			writing_file.write(seq[0]+"\t"+seq[2]+"\n")

	print("Done !")
	return


##########################################################################################
##########################################################################################



def create_labels_itol_file_reverse(info_tab):


	"""
	Function that let the possibility to change the name of the leafs' label

	:param info_tab: table of the annotation table
	:type: numpy.ndarray
	:return: Nothing
	"""

	print("\n#################")
	print("# LABELS REVERSE FILE")
	print("#################\n")

	with open(PREFIX+"_id_label_reverse.txt", 'w') as writing_file:
		writing_file.write("LABELS\n")
		writing_file.write("SEPARATOR TAB\n")
		writing_file.write("DATA\n")

		progression=0

		for seq in info_tab :

			progression += 1
			sys.stdout.write("{:.2f}% : {}/{} sequences\r".format(progression/float(info_tab.shape[0])*100, progression,info_tab.shape[0]))
			sys.stdout.flush()

			writing_file.write(seq[2]+"\t"+seq[0]+"\n")

	print("Done !")
	return

##########################################################################################
##########################################################################################

def get_color_cmap(name, n_colors=6):

	"""
	Return discrete colors from a matplotlib palette.

	:param name: Name of the palette. This should be a named matplotlib colormap.
	:type: str
    :param n_colors: Number of discrete colors in the palette.
	:type: int
	:return: List-like object of colors as hexadecimal tuples
	:type: list
    """

	brewer_qual_pals = {"Accent": 8, "Dark2": 8, "Paired": 12,
	                    "Pastel1": 9, "Pastel2": 8,
	                    "Set1": 9, "Set2": 8, "Set3": 12}

	cmap = getattr(cm, name)
	if name in brewer_qual_pals:
	    bins = np.linspace(0, 1, brewer_qual_pals[name])[:n_colors]
	else:
	    bins = np.linspace(0, 1, n_colors + 2)[1:-1]
	palette = list(map(tuple, cmap(bins)[:, :3]))
	pal_cycle = cycle(palette)
	palette = [next(pal_cycle) for _ in range(n_colors)]
	return [colors.rgb2hex(rgb) for rgb in palette]


##########################################################################################
##########################################################################################

def conversion_alignment(alignment_file, format):

	"""
	Function that convert the alignement file into fasta format

	:param alignment_file: the alignment file
	:type: str
	:param format: the format of the alignment
	:type: str
	"""

	input_handle = open(alignment_file, "rU")
	output_handle = open(os.path.join(os.path.dirname(alignment_file),"tmp.fasta"), "w")

	alignments = AlignIO.parse(input_handle, format)
	AlignIO.write(alignments, output_handle, "fasta")

	output_handle.close()
	input_handle.close()

	return os.path.join(os.path.dirname(alignment_file),"tmp.fasta")

##########################################################################################
##########################################################################################

def read_color_file(color_file) :

	"""
	Function that read the color file and set the dictionnary with it

	:param color_file: the name of the color file
	:type: str
	:return: the dictionary contructed with name as keys and color as values
	:rtype: dict

	"""


	my_color = np.loadtxt(color_file, delimiter="\t", comments="//", dtype="string")
	return {line[0]:line[1] for line in my_color}

##########################################################################################
##########################################################################################

def create_color_dict(cmap, col_infoTab, name) :

	"""
	Function that create the color file and set the dictionnary with it

	:param cmap: the cmap in which the color will be choose
	:type: str
	:param col_infoTab: the column of the INFO_TAB with the name of the system
	or phylum
	:type: numpy.ndarray
	:param name: name of the file *.color
	:type: str
	:return: the dictionary contructed with name as keys and color as values
	:rtype: dict

	"""

	color_dir = os.path.join(os.path.dirname(PREFIX), "file_color")
	create_folder(color_dir)
	uniq_infoTab = np.unique(col_infoTab)
	uniq_infoTab = uniq_infoTab[uniq_infoTab != "generic"]
	my_color = get_color_cmap(cmap, len(uniq_infoTab))
	uniq_infoTab = np.c_[uniq_infoTab, my_color]
	np.savetxt(os.path.join(color_dir, name), uniq_infoTab, delimiter='\t', fmt='%s')
	return {line[0]:line[1] for line in uniq_infoTab}

##########################################################################################
##########################################################################################

def write_big_new_file(file_tab, file_f, write_file) :

	"""
	Function that take the information about the system in the fasta file and put all the information about taxonomy
	and name

	:param file_tab: Name of the file with the taxonomic information
	:type: str
	:param file_f: name of the alignement file in fasta format
	:type: str
	:param write_file: name of the file where we want to write the information
	:type: str
	"""

	tab_info = np.loadtxt(file_tab, delimiter="\t", dtype="string", comments="##")

	with open(write_file, "w") as w_file :
		line = "#{}\t{}\t{}\t{}\t{}\t{}\n".format("leaf_name", "species_id", "species_name", "kingdom", "phylum", "system")
		w_file.write(line)
		for seq in SeqIO.parse(file_f, format="fasta") :
			if "NC_" in seq.id :
				name_species = "_".join(seq.id.split("_")[:2])
				if 'T4SS' in seq.id :
					system_name = re.search("vir[Bb][0-9][0-9]?", seq.id).group(0)
				else :
					system_name = seq.id.split("_")[seq.id.split("_").index("D")-1]
			else :
				name_species = seq.id.split("_")[0][:4]
				system_name = seq.id.split("_")[seq.id.split("_").index("V")-1]

			index_species = tab_info[:,0].tolist().index(name_species)
			line = "{}\t{}\t{}\t{}\t{}\t{}\n".format(seq.id, tab_info[index_species,0], " ".join(tab_info[index_species,1].split(" ")[:2]), tab_info[index_species,2], tab_info[index_species,3], system_name)
			w_file.write(line)
	return
