# -*- coding: utf-8 -*-

import argparse
from textwrap import dedent
import sys, os
import time

from color_label_tree import *

##########################################################################################
##########################################################################################
##
##								Main
##
##########################################################################################
##########################################################################################

# NOTE J'ai trouver le dessin pour le nom sur : http://patorjk.com/software/taag/#p=display&f=Graffiti&t=Type%20Something%20

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
     description=dedent("""


-----------------------------------------

 _____       _            _
/  __ \     | |          (_)
| /  \/ ___ | | ___  _ __ _ _______ _ __
| |    / _ \| |/ _ \| '__| |_  / _ \ '__|
| \__/\ (_) | | (_) | |  | |/ /  __/ |
 \____/\___/|_|\___/|_|  |_/___\___|_|


-----------------------------------------
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
general_option.add_argument("-o",'--output',
 							default=None,
							dest="output",
							metavar='<FOLDER>',
							help="Using <FOLDER> for output files (default: seqFile directory)")

annotation_option = parser.add_argument_group(title = "Table annotation options")
annotation_option.add_argument("-annot",'--annotationtab',
							metavar="<ANNOTATION_TAB>",
							dest="annotFile",
							default=None,
							required=True,
							help="File with the annotation for each leaf in the tree in rigth format")
annotation_option.add_argument("-add_col",'--add_columns',
							metavar="<COLUMN_NAME>",
							dest="add_columns",
							nargs="+",
							default=None,
							help="Add columns to set binary with it, each columns need to be in the annotation table with 'Yes' or 'No' values for each nodes")
annotation_option.add_argument("-sys",'--systemtree',
							action='store_true',
							dest="systemtree",
							help="Option to know if the name of the tree is a systems tree or genes tree")

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

if not args.output :
	OUTPUT = os.path.join(os.path.abspath(args.seqFile),"colorize_{}".format(time.strftime("%d_%m_%y")))
else :
	OUTPUT = args.output

create_folder(OUTPUT)

FORMAT = args.format.lower()

file_name = os.path.abspath(args.seqFile)
all_leafs = list(SeqIO.to_dict(SeqIO.parse(file_name, 'fasta')).keys())

if FORMAT != "fasta" :
	file_name_abspath = conversion_alignment(file_name, FORMAT)
	file_name = file_name_abspath

file_tab = args.annotFile

df_tab = pd.read_table(file_tab)

if args.systemtree :
	df_tab = df_tab[df_tab.System_Id.isin(all_leafs)].reset_index(drop=True)
	df_tab.drop_duplicates(subset='System_Id', keep='first', inplace=True)
	df_tab['NewName'] = df_tab["System_Id"]

else :
	df_tab = df_tab[df_tab.NewName.isin(all_leafs)].reset_index(drop=True)

#Paired bon pour les systemes < 12 sinon nipy_spectral
#Set3 pour les phylums < 12 sinon rainbow (mais pas beau et vraiment proche)

if not args.sysColor :
	DICT_COLORSTRIP = create_color_dict("nipy_spectral", df_tab.Predicted_System, "systems.color")
else :
	DICT_COLORSTRIP = read_color_file(args.sysColor)

if not args.phylumColor :
	DICT_COLORRANGE = create_color_dict("Paired", df_tab.Phylum, "phylum.color")
else :
	DICT_COLORRANGE = read_color_file(args.phylumColor)

# Appel des fonctions
######

create_colorstrip_itol_file(df_tab, OUTPUT, DICT_COLORSTRIP)

create_binary_itol_file(df_tab, OUTPUT)

create_colorrange_itol_file(df_tab, OUTPUT, DICT_COLORRANGE)

if args.add_columns :
	number_of_columns = len(args.add_columns)
	if number_of_columns <= 5 :
		colors = ['#5d3f87', '#1a4c90','#d93230','#0b7d3d', '#d2b7b1'][:number_of_columns]
	else :
		colors = get_color_cmap("Paired", number_of_columns)
	create_binary_itol_file_auto(df_tab, OUTPUT, args.add_columns, colors)

#create_labels_itol_file(df_tab, OUTPUT)

######
# Fin des fonctions

if FORMAT != "fasta" :
	os.remove(file_name_abspath)
