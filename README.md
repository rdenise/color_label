# Colorizer

Script to construct colorstrip, colorrange or binary file in Itol Format

Dependencies :
--------------

- Python 3.6 (The Python 2.7 version is no longer updated)
   - Numpy 1.12.0
   - Biopython 1.68
   - Matplotlib 2.0.0
   - Pandas 0.17.0

Format annotation file
----------------------
- The annotation file need to be a tabulate separate file with at least these columns:

NewName **[TAB]** Replicon_name **[TAB]** Predicted_system **[TAB]** System_status **[TAB]** Kingdom **[TAB]** Phylum **[TAB]** [More the columns you want to appear in the tree [-options --binary_symbol]]

Comments in all the annotation table begin with **#**

   - NewName : Name of the leaf_name
   - Replicon_name : Name code of the species
   - Kingdom : Name of the Kingdom of the species
   - Phylum : Name of the phylum of the species
   - System_status : V (verified) or D (detected)

Format color files
------------------
The color files need to be a tabulate separate file with these columns:

name(phylum or system) **[TAB]** color(hexadecimal code)

comments line need to begin by **//**

Index(['NewName', 'Hit_Id', 'Replicon_name', 'Position', 'Sequence_length',
       'Gene', 'Reference_system', 'Predicted_system', 'System_Id',
       'System_status', 'Gene_status', 'i-evalue', 'Score', 'Profile_coverage',
       'Sequence_coverage', 'Begin_match', 'End_match', 'Tandem', 'Loner',
       'Loner_unique', 'Multi_copy', 'Max_Score', 'Min_ievalue', 'Kingdom',
       'Phylum'],
      dtype='object')
