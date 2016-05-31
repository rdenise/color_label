# Colorize

Script to construct colorstrip, colorrange or binary file

Installation libraries needed
-----------------------------

1. Biopython

    pip install biopython

Format annotation file
----------------------
The annotation file need to be a tabulate separate file with these columns:

leaf_name **[TAB]** species_id **[TAB]** full_name_of_the_species **[TAB]** kingdom **[TAB]** phyllum **[TAB]** systeme_name

old annotation table format

species_id **[TAB]** full_name_of_the_species **[TAB]** kingdom **[TAB]** phyllum

and comment begin with **##**

Format color files
------------------
The color files need to be a tabulate separate file with these columns:

name(phylum or system) **[TAB]** color(hexadecimal code)

comments line need to begin by **//**
