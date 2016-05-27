"""
USAGE ::: python color_label_tree file_name file_tab myprefix format
 file_name : name of the file used for tree (alignement).
 file_tab : tab with the annotation of the leaves.
 myprefix : name of the start of the color file.
 format : format of the alignement file.
"""
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
import progressbar
import os
import numpy as np
import time
import re

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

def create_colorstrip_itol_file(fileName):

	"""
	Function that create a file in itol colorstrip format for the sequence in the fileName
	with color for each kind of
	:param fileName: base name of the fasta file of the alignement that produce the tree
	:type: string
	:return: Nothing
	"""

	print "\n#################"
	print "# COLOR STRIP FILE"
	print "#################\n"

	create_folder(os.path.join(PATH_ITOL_FILE,NAME_PROTEIN))

	with open(os.path.join(PATH_ITOL_FILE,NAME_PROTEIN,"colorstrip_"+NAME_PROTEIN+".txt"), 'w') as writing_file:
		writing_file.write("DATASET_COLORSTRIP\n")
		writing_file.write("SEPARATOR TAB\n")
		writing_file.write("DATASET_LABEL\tT2SS_T4P_Tad_Com\n")
		writing_file.write("COLOR\t#ff0000\n")
		writing_file.write("LEGEND_TITLE\tSecretion_system\n")
		writing_file.write("LEGEND_SHAPES"+"\t1"*len(DICT_COLORSTRIP)+"\n")
		writing_file.write("LEGEND_COLORS\t"+"\t".join(DICT_COLORSTRIP.values())+"\n")
		writing_file.write("LEGEND_LABELS\t"+"\t".join(DICT_COLORSTRIP.keys())+"\n")
		writing_file.write("DATA\n")

		alignment = SeqIO.parse(os.path.join(PATH_ALIGNEMENT,fileName), "fasta")

		bar = progressbar.ProgressBar(maxval=len(list(alignment)), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		progression=1
		bar.start()

		alignment = SeqIO.parse(os.path.join(PATH_ALIGNEMENT,fileName), "fasta")

		for seq in alignment :

			bar.update(progression)
			progression +=1

			if "generique" in seq.id or "choice" in seq.id:
				continue
			elif "_V_" in seq.id:
				writing_file.write(seq.id+"\t"+DICT_COLORSTRIP[seq.id.split("_")[seq.id.split("_").index('V')-1]]+"\n")
			elif "T4SS" not in seq.id :
				writing_file.write(seq.id+"\t"+DICT_COLORSTRIP[seq.id.split("_")[seq.id.split("_").index('D')-1]]+"\n")
			else :
				try:
					writing_file.write(seq.id+"\t"+DICT_COLORSTRIP[re.search("vir[Bb][0-9][0-9]?", seq.id).group(0)]+"\n")
				except KeyError:
					pass
	bar.finish()

	return

##########################################################################################
##########################################################################################

def create_colorlabel_itol_file(fileName):

	"""
	Function that create a file in itol color label format for the sequence in the fileName
	and set all the verified sequence in red
	:param fileName: base name of the fasta file of the alignement that produce the tree
	:type: string
	:return: Nothing
	"""

	print "\n#################"
	print "# LABEL COLOR FILE"
	print "#################\n"

	create_folder(os.path.join(PATH_ITOL_FILE,NAME_PROTEIN))

	with open(os.path.join(PATH_ITOL_FILE,NAME_PROTEIN,"colorlabel_"+NAME_PROTEIN+".txt"), 'w') as writing_file:
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

		alignment = SeqIO.parse(os.path.join(PATH_ALIGNEMENT,fileName), "fasta")

		bar = progressbar.ProgressBar(maxval=len(list(alignment)), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		progression=1
		bar.start()

		alignment = SeqIO.parse(os.path.join(PATH_ALIGNEMENT,fileName), "fasta")

		for seq in alignment :

			bar.update(progression)
			progression +=1

			if "_V_" in seq.id :
				writing_file.write(seq.id+"\t1\n")
			else :
				writing_file.write(seq.id+"\t0\n")

	bar.finish()

	return

##########################################################################################
##########################################################################################



def create_colorrange_itol_file(fileName, fileTab):

	"""
	Function that create a file in itol color range format for the sequence in the fileName
	and with color for each kingdom or subkingdom
	:param fileName: base name of the fasta file of the alignement that produce the tree
	:type: string
	:param fileTab: Absolute path of the tabular file with the information about the species
	(kingdom, lineage)
	:type: string
	:return: Nothing
	"""

	print "\n#################"
	print "# COLOR RANGE FILE"
	print "#################\n"

	create_folder(os.path.join(PATH_ITOL_FILE,NAME_PROTEIN))

	with open(os.path.join(PATH_ITOL_FILE,NAME_PROTEIN,"colorrange_"+NAME_PROTEIN+".txt"), 'w') as writing_file:
		writing_file.write("TREE_COLORS\n")
		writing_file.write("SEPARATOR SPACE\n")
		writing_file.write("DATA\n")

		alignment = SeqIO.parse(os.path.join(PATH_ALIGNEMENT,fileName), "fasta")
		info_tab = np.loadtxt(fileTab, dtype='string', comments='##', delimiter='\t')

		bar = progressbar.ProgressBar(maxval=len(list(alignment)), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		progression=1
		bar.start()

		alignment = SeqIO.parse(os.path.join(PATH_ALIGNEMENT,fileName), "fasta")

		for seq in alignment :

			bar.update(progression)
			progression +=1

			if "NC_" in seq.id :
				sub_id = "_".join(seq.id.split("_")[:2])
			else :
				sub_id = seq.id.split("_")[0][:4]


			if "Bacteria" in info_tab[:,2][list(info_tab[:,0]).index(sub_id)] :
				lineage = info_tab[:,3][list(info_tab[:,0]).index(sub_id)]
			else :
				lineage = info_tab[:,2][list(info_tab[:,0]).index(sub_id)]


			if lineage in DICT_COLORRANGE :
				writing_file.write(seq.id+" range "+DICT_COLORRANGE[lineage]+" "+lineage+"\n")
			else :
				writing_file.write(seq.id+" range "+DICT_COLORRANGE["Bacteria"]+" Bacteria\n")

	bar.finish()

	return

##########################################################################################
##########################################################################################



def create_labels_itol_file(fileName, fileTab):


	"""
	Function that let the possibility to change the name of the leafs' label

	:param fileName: base name of the fasta file of the alignement that produce the tree
	:type: string
	:param fileTab: Absolute path of the tabular file with the information about the species
	(kingdom, lineage)
	:type: string
	:return: Nothing
	"""

	print "\n#################"
	print "# LABELS FILE"
	print "#################\n"

	create_folder(os.path.join(PATH_ITOL_FILE,NAME_PROTEIN))

	with open(os.path.join(PATH_ITOL_FILE,NAME_PROTEIN,"id_label"+NAME_PROTEIN+".txt"), 'w') as writing_file:
		writing_file.write("LABELS\n")
		writing_file.write("SEPARATOR TAB\n")
		writing_file.write("DATA\n")

		alignment = SeqIO.parse(os.path.join(PATH_ALIGNEMENT,fileName), "fasta")
		info_tab = np.loadtxt(fileTab, dtype='string', comments='##', delimiter='\t')

		bar = progressbar.ProgressBar(maxval=len(list(alignment)), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		progression=1
		bar.start()

		alignment = SeqIO.parse(os.path.join(PATH_ALIGNEMENT,fileName), "fasta")

		for seq in alignment :

			bar.update(progression)
			progression +=1

			if "NC_" in seq.id :
				sub_id = "_".join(seq.id.split("_")[:2])
			else :
				sub_id = seq.id.split("_")[0][:4]

			name = info_tab[:,1][list(info_tab[:,0]).index(sub_id)]

			writing_file.write(seq.id+"\t"+" ".join(name.split(" ")[:2])+"\n")

	bar.finish()

	return


##########################################################################################
##########################################################################################



def create_labels_itol_file_reverse(fileName, fileTab):


	"""
	Function that let the possibility to change the name of the leafs' label

	:param fileName: base name of the fasta file of the alignement that produce the tree
	:type: string
	:param fileTab: Absolute path of the tabular file with the information about the species
	(kingdom, lineage)
	:type: string
	:return: Nothing
	"""

	print "\n#################"
	print "# LABELS REVERSE FILE"
	print "#################\n"

	create_folder(os.path.join(PATH_ITOL_FILE,NAME_PROTEIN))

	with open(os.path.join(PATH_ITOL_FILE,NAME_PROTEIN,"id_label_reverse"+NAME_PROTEIN+".txt"), 'w') as writing_file:
		writing_file.write("LABELS\n")
		writing_file.write("SEPARATOR TAB\n")
		writing_file.write("DATA\n")

		alignment = SeqIO.parse(os.path.join(PATH_ALIGNEMENT,fileName), "fasta")
		info_tab = np.loadtxt(fileTab, dtype='string', comments='##', delimiter='\t')

		bar = progressbar.ProgressBar(maxval=len(list(alignment)), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		progression=1
		bar.start()

		alignment = SeqIO.parse(os.path.join(PATH_ALIGNEMENT,fileName), "fasta")

		for seq in alignment :

			bar.update(progression)
			progression +=1

			if "NC_" in seq.id :
				sub_id = "_".join(seq.id.split("_")[:2])
			else :
				sub_id = seq.id.split("_")[0][:4]

			name = info_tab[:,1][list(info_tab[:,0]).index(sub_id)]

			writing_file.write(" ".join(name.split(" ")[:2])+"\t"+seq.id+"\n")

	bar.finish()

	return

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
##
##								Main
##
##########################################################################################
##########################################################################################

if __name__ == '__main__':
	if len(sys.argv)!=5:
		print "Usage : python %s fichier_alignment INFO_TAB myprefix format_alignment." % sys.argv[0]
		sys.exit(1)

	file_name = os.path.basename(sys.argv[1])
	PATH_ALIGNEMENT = os.path.dirname(os.path.abspath(sys.argv[1]))

	file_tab = os.path.abspath(sys.argv[2])

	NAME_PROTEIN = sys.argv[3]

	FORMAT = sys.argv[4].lower()

	if FORMAT != "fasta" :
		file_name_abspath = conversion_alignment(os.path.join(PATH_ALIGNEMENT,file_name), FORMAT)
		file_name = os.path.basename(file_name_abspath)

	create_colorstrip_itol_file(file_name)

	create_colorlabel_itol_file(file_name)

	create_colorrange_itol_file(file_name, file_tab)

	create_labels_itol_file(file_name, file_tab)

	if FORMAT != "fasta" :
		os.remove(file_name_abspath)

PATH_ITOL_FILE = "/Users/rdenise/Documents/Analysis_tree/virb4/label_itol/%s/" % time.strftime("%d_%m_%y")
