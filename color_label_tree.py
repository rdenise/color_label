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
import seaborn as sns
import argparse
from textwrap import dedent

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

	print "\n#################"
	print "# COLOR STRIP FILE"
	print "#################\n"

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


		bar = progressbar.ProgressBar(maxval=info_tab.shape[0], widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		progression=1
		bar.start()


		for seq in info_tab :

			bar.update(progression)
			progression +=1

			if "generique" in seq[0] or "choice" in seq[0]:
				continue
			else:
				writing_file.write(seq[0]+"\t"+DICT_COLORSTRIP[seq[-1]]+"\n")
	bar.finish()

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

	print "\n#################"
	print "# LABEL BINARY FILE"
	print "#################\n"

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


		bar = progressbar.ProgressBar(maxval=info_tab.shape[0], widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		progression=1
		bar.start()


		for seq in info_tab :

			bar.update(progression)
			progression +=1

			if "_V_" in seq[0] :
				writing_file.write(seq[0]+"\t1\n")
			else :
				writing_file.write(seq[0]+"\t0\n")

	bar.finish()

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

	print "\n#################"
	print "# COLOR RANGE FILE"
	print "#################\n"

	with open(PREFIX+"_colorrange.txt", 'w') as writing_file:
		writing_file.write("TREE_COLORS\n")
		writing_file.write("SEPARATOR SPACE\n")
		writing_file.write("DATA\n")

		bar = progressbar.ProgressBar(maxval=info_tab.shape[0], widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		progression=1
		bar.start()


		for seq in info_tab :

			bar.update(progression)
			progression +=1

			lineage = seq[-2]

			if lineage in DICT_COLORRANGE :
				writing_file.write(seq[0]+" range "+DICT_COLORRANGE[lineage]+" "+lineage+"\n")
			elif seq[-3]=="Archaea":
				writing_file.write(seq[0]+" range "+DICT_COLORRANGE["Archaea"]+" Archaea\n")
			else :
				writing_file.write(seq[0]+" range "+DICT_COLORRANGE["Bacteria"]+" Bacteria\n")

	bar.finish()

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

	print "\n#################"
	print "# LABELS FILE"
	print "#################\n"

	with open(PREFIX+"_id_label.txt", 'w') as writing_file:
		writing_file.write("LABELS\n")
		writing_file.write("SEPARATOR TAB\n")
		writing_file.write("DATA\n")

		bar = progressbar.ProgressBar(maxval=info_tab.shape[0], widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		progression=1
		bar.start()


		for seq in info_tab :

			bar.update(progression)
			progression +=1

			writing_file.write(seq[0]+"\t"+seq[2]+"\n")

	bar.finish()

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

	print "\n#################"
	print "# LABELS REVERSE FILE"
	print "#################\n"

	with open(PREFIX+"_id_label_reverse.txt", 'w') as writing_file:
		writing_file.write("LABELS\n")
		writing_file.write("SEPARATOR TAB\n")
		writing_file.write("DATA\n")

		bar = progressbar.ProgressBar(maxval=info_tab.shape[0], widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		progression=1
		bar.start()


		for seq in info_tab :

			bar.update(progression)
			progression +=1

			writing_file.write(seq[2]+"\t"+seq[0]+"\n")

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
	uniq_infoTab = uniq_infoTab[uniq_infoTab != "generique"]
	my_color = sns.color_palette(cmap, len(uniq_infoTab)).as_hex()
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
		line = "#%s\t%s\t%s\t%s\t%s\t%s\n" % ("leaf_name", "species_id", "species_name", "kingdom", "phylum", "system")
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
			line = "%s\t%s\t%s\t%s\t%s\t%s\n" % (seq.id, tab_info[index_species,0], " ".join(tab_info[index_species,1].split(" ")[:2]), tab_info[index_species,2], tab_info[index_species,3], system_name)
			w_file.write(line)


##########################################################################################
##########################################################################################
##
##								Main
##
##########################################################################################
##########################################################################################


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
     description=dedent("""


---------------------------------------------------------------------------------------

 XXXXXXX             XX                    XXXXXX          XXXXXXXXX           XXXXXX
XX          XXXXX    XX          XXXXX    XX    XX  XXXXXX      XXX  XXXXXXXX XX    XX
XX        XXX   XXX  XX        XXX   XXX  XX     XX   XX       XXX   XX       XX     XX
XX       XX       XX XX       XX       XX XX    XX    XX      XXX    XX       XX    XX
XX       XX       XX XX       XX       XX XX XXX      XX     XXX     XXXX     XX XXX
XX       XX       XX XX       XX       XX XX  XX      XX    XXX      XX       XX  XX
 XXXXXXX  XXX   XXX   XXXXXXX  XXX   XXX  XX   XXXX   XX   XXXXXXXXX XX       XX   XXXX
             XXX                  XXX               XXXXXX           XXXXXXXX

---------------------------------------------------------------------------------------
""") )

general_option = parser.add_argument_group(title = "General input dataset options")
general_option.add_argument("-s",'--seqfile',
 							required=True,
							metavar="<file>",
							dest="seqFile",
							help="File with the sequences used for the alignment")
general_option.add_argument("-f",'--format',
							required=True,
							dest="format",
							help="")
general_option.add_argument("-pre",'--prefix',
 							default=None,
							dest="prefix",
							metavar='<PREFIX>',
							help="Using <PREFIX> for output files (default: seqFile directory)")

annotation_option = parser.add_argument_group(title = "Table annotation options")
annotation_option.add_argument("-old",'--oldannotation',
							metavar="<INFO_TAB>",
							dest="oldInfo",
							default=None,
							help="File with the annotation for each species in old format")
annotation_option.add_argument("-annot",'--annotationtab',
							metavar="<ANNOTATION_TAB>",
							dest="annotFile",
							default=None,
							help="File with the annotation for each leaf in the tree in rigth format")

color_option = parser.add_argument_group(title = "Color options")
color_option.add_argument("-sysCo",'--systemsColor',
							metavar="<Color_system_File>",
							dest="sysColor",
							default=None,
							help="File the name of the of the systems and the color for each systems in hexadecimal")
color_option.add_argument("-phyCo",'--phylumColor',
							metavar="<Color_system_File>",
							dest="phylumColor",
							default=None,
							help="File the name of the of the phylum and the color for each systems in hexadecimal")

args = parser.parse_args()

if not args.prefix :
	PREFIX = os.path.join(os.path.abspath(args.seqFile),"colorize_%s" %(time.strftime("%d_%m_%y")))
else :
	PREFIX = args.prefix

create_folder(os.path.dirname(PREFIX))

FORMAT = args.format.lower()

file_name = os.path.abspath(args.seqFile)

if FORMAT != "fasta" :
	file_name_abspath = conversion_alignment(file_name, FORMAT)
	file_name = file_name_abspath

if not (args.oldInfo or args.annotFile):
    parser.error("you MUST provided annotation table.")
elif args.oldInfo :
	write_big_new_file(args.oldInfo, file_name, os.path.join(os.path.dirname(PREFIX),"ANNOTATION_TAB"))
	file_tab = os.path.abspath(os.path.join(os.path.dirname(PREFIX),"ANNOTATION_TAB"))
elif args.annotFile:
	file_tab = args.annotFile

tab_numpy = np.loadtxt(file_tab, delimiter="\t", dtype="string")

#Paired bon pour les systemes < 12 sinon nipy_spectral
#Set3 pour les phylums < 12 sinon rainbow (mais pas beau et vraiment proche)

if not args.sysColor :
	DICT_COLORSTRIP = create_color_dict("nipy_spectral", tab_numpy[:,-1], "systems.color")
else :
	DICT_COLORSTRIP = read_color_file(color_file)

if not args.phylumColor :
	DICT_COLORRANGE = create_color_dict("Paired", tab_numpy[:,-2], "phylum.color")
else :
	DICT_COLORRANGE = read_color_file(color_file)

# Appel des fonctions
######

create_colorstrip_itol_file(tab_numpy)

create_binary_itol_file(tab_numpy)

create_colorrange_itol_file(tab_numpy)

create_labels_itol_file(tab_numpy)

######
# Fin des fonctions

if FORMAT != "fasta" :
	os.remove(file_name_abspath)
