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

tab_numpy = np.genfromtxt(file_tab, delimiter="\t", dtype="str")

#Paired bon pour les systemes < 12 sinon nipy_spectral
#Set3 pour les phylums < 12 sinon rainbow (mais pas beau et vraiment proche)

if not args.sysColor :
	DICT_COLORSTRIP = create_color_dict("nipy_spectral", tab_numpy[:,-1], "systems.color")
else :
	DICT_COLORSTRIP = read_color_file(args.sysColor)

if not args.phylumColor :
	DICT_COLORRANGE = create_color_dict("Paired", tab_numpy[:,-2], "phylum.color")
else :
	DICT_COLORRANGE = read_color_file(args.phylumColor)

# Appel des fonctions
######

create_colorstrip_itol_file(tab_numpy, PREFIX, DICT_COLORSTRIP)

create_binary_itol_file(tab_numpy, PREFIX)

create_colorrange_itol_file(tab_numpy, PREFIX, DICT_COLORRANGE)

create_labels_itol_file(tab_numpy, PREFIX)

######
# Fin des fonctions

if FORMAT != "fasta" :
	os.remove(file_name_abspath)
