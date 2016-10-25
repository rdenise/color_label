# Colorizer

Script to construct colorstrip, colorrange or binary file in Itol Format

Dependencies :
--------------

- Python 3.5 (The Python 2.7 version is no longer updated)
   - Numpy 1.11.2
   - Biopython 1.68
   - Matplotlib 1.5.3

Format annotation file
----------------------
- The annotation file need to be a tabulate separate file with these columns:

leaf_name **[TAB]** species_id **[TAB]** full_name_of_the_species **[TAB]** kingdom **[TAB]** phyllum **[TAB]** systeme_name

- old annotation table format :

species_id **[TAB]** full_name_of_the_species **[TAB]** kingdom **[TAB]** phyllum

Comments in all the annotation table begin with **##**

Format color files
------------------
The color files need to be a tabulate separate file with these columns:

name(phylum or system) **[TAB]** color(hexadecimal code)

comments line need to begin by **//**
