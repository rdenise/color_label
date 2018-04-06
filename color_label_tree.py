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
import pandas as pd
import re
import matplotlib.cm as cm
import matplotlib.colors as colors
from itertools import cycle

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
	:type: str
	:return: Nothing
	"""

	try:
		os.makedirs(mypath)
	except OSError:
		pass

	return

##########################################################################################
##########################################################################################

def create_colorstrip_itol_file_systems(info_df, PREFIX, DICT_COLORSTRIP):

	"""
	Function that create a file in itol colorstrip format for the sequence in the info_df
	with color for each kind of systems

	:param info_df: table of the annotation table
	:type: pandas.DataFrame
	:param PREFIX: the path to the forlder of the new file name
	:type: str
	:param DICT_COLORRANGE: dictionnary that contain the name of the systems with a color associate
	:type: dict
	:return: Nothing
	"""

	print("\n#################")
	print("# COLOR STRIP FILE")
	print("#################\n")

	with open(os.path.join(PREFIX,"colorstrip_systems.txt"), 'w') as writing_file:
		writing_file.write("DATASET_COLORSTRIP\n")
		writing_file.write("SEPARATOR TAB\n")
		writing_file.write("DATASET_LABEL\tSystems\n")
		writing_file.write("COLOR\t#ff0000\n")
		writing_file.write("LEGEND_TITLE\tSecretion_system\n")
		writing_file.write("LEGEND_SHAPES{}\n".format("\t1"*len(DICT_COLORSTRIP)))
		writing_file.write("LEGEND_COLORS\t{}\n".format("\t".join(DICT_COLORSTRIP.values())))
		writing_file.write("LEGEND_LABELS\t{}\n".format("\t".join(DICT_COLORSTRIP.keys())))
		writing_file.write("STRIP_WIDTH\t50\n")
		writing_file.write("MARGIN\t25\n")
		writing_file.write("DATA\n")

		#info_df = info_df[~info_df.Predicted_system.isin(["generic", "generique", "choice"])].reset_index(drop=True)
		info_df["color_strip"] = info_df.apply(lambda x : DICT_COLORSTRIP[x.Predicted_system], axis=1)
		info_df[["NewName", "color_strip"]].to_csv(writing_file, sep="\t", index=None, header=None)
		writing_file.write("\n")

	print()
	print("Done !")
	return


##########################################################################################
##########################################################################################

def create_colorstrip_itol_file_phylum(info_df, PREFIX, DICT_COLORSTRIP):

	"""
	Function that create a file in itol colorstrip format for the sequence in the info_df
	with color for each kind of systems

	:param info_df: table of the annotation table
	:type: pandas.DataFrame
	:param PREFIX: the path to the forlder of the new file name
	:type: str
	:param DICT_COLORRANGE: dictionnary that contain the name of the systems with a color associate
	:type: dict
	:return: Nothing
	"""

	print("\n#################")
	print("# COLOR STRIP FILE")
	print("#################\n")

	with open(os.path.join(PREFIX,"colorstrip_phylum.txt"), 'w') as writing_file:
		writing_file.write("DATASET_COLORSTRIP\n")
		writing_file.write("SEPARATOR TAB\n")
		writing_file.write("DATASET_LABEL\tPhylum\n")
		writing_file.write("COLOR\t#ff0000\n")
		writing_file.write("LEGEND_TITLE\tPhylum\n")
		writing_file.write("LEGEND_SHAPES{}\n".format("\t1"*len(DICT_COLORSTRIP)))
		writing_file.write("LEGEND_COLORS\t{}\n".format("\t".join(DICT_COLORSTRIP.values())))
		writing_file.write("LEGEND_LABELS\t{}\n".format("\t".join(DICT_COLORSTRIP.keys())))
		writing_file.write("STRIP_WIDTH\t50\n")
		writing_file.write("MARGIN\t25\n")		
		writing_file.write("DATA\n")

		#info_df = info_df[~info_df.Predicted_system.isin(["generic", "generique", "choice"])].reset_index(drop=True)
		info_df["name_range"] = info_df.apply(lambda x: x.Phylum if x.Phylum in DICT_COLORSTRIP else x.Kingdom, axis=1)
		info_df["color_strip"] = info_df.apply(lambda x: DICT_COLORSTRIP[x.name_range], axis=1)
		info_df[["NewName", "color_strip"]].to_csv(writing_file, sep="\t", index=None, header=None)
		writing_file.write("\n")

	print()
	print("Done !")
	return

##########################################################################################
##########################################################################################

def create_binary_itol_file(info_df, PREFIX):

	"""
	Function that create a file in itol color label format for the sequence in the info_df
	and set all the validated sequence will have grey square

	:param info_df: table of the annotation table
	:type: pandas.DataFrame
	:param PREFIX: the path to the forlder of the new file name
	:type: str
	:return: Nothing
	"""

	print("\n#################")
	print("# LABEL BINARY FILE")
	print("#################\n")

	with open(os.path.join(PREFIX,"labelbinary_validated.txt"), 'w') as writing_file:
		writing_file.write("DATASET_BINARY\n")
		writing_file.write("SEPARATOR TAB\n")
		writing_file.write("COLOR\t#b15928\n")
		writing_file.write("DATASET_LABEL\tValidated sequence\n")
		writing_file.write("FIELD_SHAPES\t1\n")
		writing_file.write("FIELD_LABELS\tvalidated\n")
		writing_file.write("FIELD_COLORS\t#b15928\n")
		writing_file.write("LEGEND_TITLE\tExperimentaly validated sequences\n")
		writing_file.write("LEGEND_SHAPES\t1\n")
		writing_file.write("LEGEND_COLORS\t#b15928\n")
		writing_file.write("LEGEND_LABELS\tverify\n")
		writing_file.write("DATA\n")

		info_df = info_df.replace({"System_status":{"V":1, "D":-1}})
		info_df[["NewName", "System_status"]].to_csv(writing_file, sep="\t", index=None, header=None)
		writing_file.write("\n")

	print()
	print("Done !")
	return

##########################################################################################
##########################################################################################



def create_colorrange_itol_file_phylum(info_df, PREFIX, DICT_COLORRANGE):

	"""
	Function that create a file in itol color range format for the sequence in the info_df
	and with color for each kingdom or subkingdom

	:param info_df: table of the annotation table
	:type: pandas.DataFrame
	:param PREFIX: the path to the forlder of the new file name
	:type: str
	:param DICT_COLORRANGE: dictionnary that contain the name of the phylum with a color associate
	:type: dict
	:return: Nothing
	"""

	print("\n#################")
	print("# COLOR RANGE FILE")
	print("#################\n")

	with open(os.path.join(PREFIX,"colorrange_phylum.txt"), 'w') as writing_file:
		writing_file.write("TREE_COLORS\n")
		writing_file.write("SEPARATOR TAB\n")
		writing_file.write("DATA\n")


		info_df["range_col"] = "range"
		#info_df["name_range"] = info_df.apply(lambda x: x.Phylum if x.Phylum in DICT_COLORRANGE else x.Lineage.split(";")[2] if x.Phylum == "Proteobacteria" else x.Kingdom, axis=1)
		info_df["name_range"] = info_df.apply(lambda x: x.Phylum if x.Phylum in DICT_COLORRANGE else x.Kingdom, axis=1)
		info_df["color_lineage"] = info_df.apply(lambda x: DICT_COLORRANGE[x.name_range], axis=1)

		info_df[["NewName","range_col", "color_lineage", "name_range"]].to_csv(writing_file, sep="\t", index=None, header=None)
		writing_file.write("\n")

	print()
	print("Done !")
	return

##########################################################################################
##########################################################################################



def create_colorrange_itol_file_systems(info_df, PREFIX, DICT_COLORRANGE):

	"""
	Function that create a file in itol color range format for the sequence in the info_df
	and with color for each kingdom or subkingdom

	:param info_df: table of the annotation table
	:type: pandas.DataFrame
	:param PREFIX: the path to the forlder of the new file name
	:type: str
	:param DICT_COLORRANGE: dictionnary that contain the name of the phylum with a color associate
	:type: dict
	:return: Nothing
	"""

	print("\n#################")
	print("# COLOR RANGE FILE")
	print("#################\n")

	with open(os.path.join(PREFIX,"colorrange_systems.txt"), 'w') as writing_file:
		writing_file.write("TREE_COLORS\n")
		writing_file.write("SEPARATOR TAB\n")
		writing_file.write("DATA\n")


		info_df["range_col"] = "range"
		#info_df["name_range"] = info_df.apply(lambda x: x.Phylum if x.Phylum in DICT_COLORRANGE else x.Lineage.split(";")[2] if x.Phylum == "Proteobacteria" else x.Kingdom, axis=1)
		info_df["name_range"] = info_df.apply(lambda x: x.Predicted_system, axis=1)
		info_df["color_lineage"] = info_df.apply(lambda x : DICT_COLORRANGE[x.Predicted_system], axis=1)
		
		info_df[["NewName","range_col", "color_lineage", "name_range"]].to_csv(writing_file, sep="\t", index=None, header=None)
		writing_file.write("\n")

	print()
	print("Done !")
	return

##########################################################################################
##########################################################################################

def create_labels_itol_file(info_df, PREFIX):


	"""
	Function that let the possibility to change the name of the leafs' label

	:param info_df: table of the annotation table
	:type: pandas.DataFrame
	:param PREFIX: the path to the forlder of the new file name
	:type: str
	:return: Nothing
	"""

	print("\n#################")
	print("# LABELS FILE")
	print("#################\n")

	with open(os.path.join(PREFIX,"id_label.txt"), 'w') as writing_file:
		writing_file.write("LABELS\n")
		writing_file.write("SEPARATOR TAB\n")
		writing_file.write("DATA\n")

		progression=0

		for seq in info_df :

			progression += 1
			sys.stdout.write("{:.2f}% : {}/{} sequences\r".format(progression/float(info_df.shape[0])*100, progression,info_df.shape[0]))
			sys.stdout.flush()

			writing_file.write("{}\t{}\n".format(seq[0],seq[2]))

	print()
	print("Done !")
	return


##########################################################################################
##########################################################################################

def create_labels_itol_file_reverse(info_df, PREFIX):


	"""
	Function that let the possibility to change the name of the leafs' label

	:param info_df: table of the annotation table
	:type: pandas.DataFrame
	:param PREFIX: the path to the forlder of the new file name
	:type: str
	:return: Nothing
	"""

	print("\n#################")
	print("# LABELS REVERSE FILE")
	print("#################\n")

	with open(os.path.join(PREFIX,"id_label_reverse.txt"), 'w') as writing_file:
		writing_file.write("LABELS\n")
		writing_file.write("SEPARATOR TAB\n")
		writing_file.write("DATA\n")

		progression=0

		for seq in info_df :

			progression += 1
			sys.stdout.write("{:.2f}% : {}/{} sequences\r".format(progression/float(info_df.shape[0])*100, progression,info_df.shape[0]))
			sys.stdout.flush()

			writing_file.write("{}\t{}\n".format(seq[2],seq[0]))

	print()
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


	my_color = np.genfromtxt(color_file, delimiter="\t", comments="//", dtype="str")
	return {line[0]:line[1] for line in my_color}

##########################################################################################
##########################################################################################

def create_color_dict(cmap, col_infoTab, name, PREFIX) :

	"""
	Function that create the color file and set the dictionnary with it

	:param cmap: the cmap in which the color will be choose
	:type: str
	:param col_infoTab: the column of the INFO_TAB with the name of the system
	or phylum
	:type: pandas.DataFrame
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

def create_binary_itol_file_auto(info_df, PREFIX, columns_names, colors):

	"""
	Function that create a file in itol color label format for the sequence in the info_df
	and set all the information about the evalue and score and multi copy ...

	:param info_df: table of the annotation table
	:type: pandas.DataFrame
	:param PREFIX: the path to the forlder of the new file name
	:type: str
	:param columns_names: Name of the columns we want to set
	:type: list of str
	:param colors: The colors for each params
	:type: list of str
	:return: Nothing
	"""

	print("\n#################")
	print("# LABEL BINARY FILE")
	print("#################\n")

	with open(os.path.join(PREFIX,"labelbinary_params_test.txt"), 'w') as writing_file:
		writing_file.write("DATASET_BINARY\n")
		writing_file.write("SEPARATOR TAB\n")
		writing_file.write("COLOR\t{}\n".format(colors[0]))
		writing_file.write("DATASET_LABEL\tParams test\n")
		writing_file.write("FIELD_SHAPES{}\n".format("\t1"*len(columns_names)))
		writing_file.write("FIELD_LABELS\t{}\n".format("\t".join(columns_names)))
		writing_file.write("FIELD_COLORS\t{}\n".format("\t".join(colors)))
		# XXX Permet d'avoir une legend dans le dessin
		writing_file.write("LEGEND_TITLE\tTest Params Duplicates\n")
		writing_file.write("LEGEND_SHAPES{}\n".format("\t1"*len(columns_names)))
		writing_file.write("LEGEND_COLORS\t{}\n".format("\t".join(colors)))
		writing_file.write("LEGEND_LABELS\t{}\n".format("\t".join(columns_names)))
		writing_file.write("DATA\n")

		if "System_status" in info_df.columns :
			info_df = info_df[info_df.System_status == "D"].reset_index(drop=True)

		for column_name in columns_names :
			if "Multi_copy" in info_df.columns :
				info_df[column_name] = info_df.apply(lambda x : -1 if x.Multi_copy == "No" else 1 if x[column_name] == "Yes" else 0, axis=1)
			else :
				info_df[column_name] = info_df.apply(lambda x : 1 if x[column_name] == "Yes" else 0, axis=1)

		info_df[["NewName"]+columns_names].to_csv(writing_file, sep="\t", index=None, header=None)
		writing_file.write("\n")

	print()
	print("Done !")
	return

##########################################################################################
##########################################################################################

def create_popup_info_itol_file(info_df, PREFIX):


	"""
	Function that let the possibility to change the name of the leafs' label

	:param info_df: table of the annotation table
	:type: pandas.DataFrame
	:param PREFIX: the path to the forlder of the new file name
	:type: str
	:return: Nothing
	"""

	print("\n#################")
	print("# POPUP FILE")
	print("#################\n")

	with open(os.path.join(PREFIX,"popup_info.txt"), 'w') as writing_file:
		writing_file.write("POPUP_INFO\n")
		writing_file.write("SEPARATOR TAB\n")
		writing_file.write("DATA\n")

		progression=0


		##9606,Homo sapiens info popup,<h1>Homo sapiens</h1><p style='color:blue'>More info at <a target='_blank' href='http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9606'>NCBI</a></p>

		info_df["NewName2"] = info_df["NewName"]
		info_df["popup"] = info_df.apply(lambda x : "<h4 style='color:blue'>{}</h4><h4>Info Lineage :</h4><h5>Kingdom</h5><p style='color:blue'>{}</p><h5>Phylum</h5><p style='color:blue'>{}</p><h5>Full lineage</h5><p style='color:blue'>{}</p> ".format(x.Species_name, x.Kingdom, x.Phylum, x.Lineage), axis=1)
		info_df[["NewName","NewName2",  "popup"]].to_csv(writing_file, sep="\t", index=None, header=None)
		writing_file.write("\n")

	print()
	print("Done !")
	return