import numpy as np
import pandas as pd
import csv 
import pdb2DF 
from itertools import groupby, islice
import re 
pd.options.mode.chained_assignment = None

# Table of 20 amino acids and the number of exchangeable hydrogens in each
H_ex = pd.read_csv('H_D_ex.txt', sep = '\t',header = 0)
H_ex = H_ex.astype({"H(ex)":int})
H_ex.index.name = 'Index'

# Put amino acids in a list
AA_list = H_ex["Amino_acid"].tolist()

# Use pdb2DF function to create a DF from pdb file
# Also return the input pdb file, and user input for is the file an antibody or protein
DF_from_PDB, file, is_anti = pdb2DF.pdb2DF()

# Read the DataFrame file
pdb_DF = pd.read_csv(DF_from_PDB, header = 0, index_col = False)

# Rename the index column
pdb_DF.index.name = 'Index'

# Convert NaNs to zeroes
pdb_DF = pdb_DF.fillna(0)

# Convert atom numbers and amino acid numbers to integers
pdb_DF = pdb_DF.astype({"ATOM_NR": int, "AA_NR": int})

# Create a list of all values in the amino acid name column in DF
AA_column = pdb_DF["AA_NAME"].tolist()

# Create a list of all values that are NOT names of amino acids, i.e. zeroes/NaNs
SOL_list = []
for AA in AA_column:
    if AA not in AA_list:
        SOL_list.append(AA)


# Create a new clean DF 
clean_DF = pd.DataFrame()
# Only include amino acids
clean_DF = pdb_DF[~pdb_DF["AA_NAME"].isin(SOL_list)]

# Create a amino acid sequence ID for antibody and protein files:
if is_anti == 1:
    # Create a new column which combines amino acid with its sequence number
    clean_DF["AA_seq"] = clean_DF["AA_NAME"] +'_'+clean_DF["CHAIN"]+'_'+ clean_DF["AA_NR"].astype(str)

elif is_anti == 2:
    # Create a new column which combines amino acid with its sequence number
    clean_DF["AA_seq"] = clean_DF["AA_NAME"] +'_'+ clean_DF["AA_NR"].astype(str)

# Create a list of the correct sequence of amino acids
AA_seq = clean_DF["AA_seq"].tolist() 
# Remove duplicating values
AA_seq = list(dict.fromkeys(AA_seq))

# DEBUGGING
# fout = open('AA_seq.txt', 'wt')

# for line in AA_seq:
#     new_line = str(line)+'\n'
#     fout.write(new_line)

# fout.close()

# Create a dataframe containing selected columns
H_DF = clean_DF[["ELEMENT",
                 "AA_NAME",
                 "ATOM_NAME",
                 "AA_seq"]]

# Merge with amino acid DataFrame
H_DF = H_DF.merge(H_ex,
                 left_on = "AA_NAME",
                 right_on = "Amino_acid",
                 how = "left")

# Select required columns
H_DF = H_DF[["ELEMENT",
             "AA_NAME",
             "ATOM_NAME",
             "AA_seq",
             "H(ex)"]]

# Create a dictionary for the amino acid sequence ID
dictionary = {}
for AA in AA_seq:
    dictionary[AA] = H_DF[H_DF["AA_seq"] == AA]

# HYDROGEN-DEUTERIUM EXCHANGE
temp_list = []
# Loop over amino acid sequence IDs
for AA in AA_seq:
    # Create a temporary DataFrame for a given ID
    temp_df = dictionary[AA]

    # Go through each amino acid and exchange the correct number of hydrogens
    if temp_df["AA_NAME"].any() == "LYS":
        for ex in range(4):
            # Try to change all of the required hydrogens
            try:
                H_index = temp_df[temp_df["ATOM_NAME"] == "H"].first_valid_index()
                temp_df.loc[H_index,"ATOM_NAME"] = temp_df.loc[H_index,"ATOM_NAME"].replace("H","D")
            # If there are fewer hydrogens in the sequence than needs to be exchanged
            except KeyError:
                print('There number of exchangeable hydrogens is greater than the number of hydrogens in the amino acid. \nAll available hydrogens will be still exchanged.')
    elif temp_df["AA_NAME"].any() == "GLY" or temp_df["AA_NAME"].any() == "ALA" or temp_df["AA_NAME"].any() == "VAL" or temp_df["AA_NAME"].any() == "LEU" or temp_df["AA_NAME"].any() == "ILE" or temp_df["AA_NAME"].any() == "PHE" or temp_df["AA_NAME"].any() == "TYR" or temp_df["AA_NAME"].any() == "ASP" or temp_df["AA_NAME"].any() == "GLU" or temp_df["AA_NAME"].any() == "HIS" or temp_df["AA_NAME"].any() == "HSD" or temp_df["AA_NAME"].any() == "MET":
        try:
            H_index = temp_df[temp_df["ATOM_NAME"] == "H"].first_valid_index()
            temp_df.loc[H_index,"ATOM_NAME"] = temp_df.loc[H_index,"ATOM_NAME"].replace("H","D")
        except KeyError:
            print('There number of exchangeable hydrogens is greater than the number of hydrogens in the amino acid. \nAll available hydrogens will be still exchanged.')
    elif temp_df["AA_NAME"].any() == "TYR" or temp_df["AA_NAME"].any() == "TRP" or temp_df["AA_NAME"].any() == "SER" or temp_df["AA_NAME"].any() == "THR" or temp_df["AA_NAME"].any() == "CYS":
        for ex in range(2):
            try:
                H_index = temp_df[temp_df["ATOM_NAME"] == "H"].first_valid_index()
                temp_df.loc[H_index,"ATOM_NAME"] = temp_df.loc[H_index,"ATOM_NAME"].replace("H","D")
            except KeyError:
                print('There number of exchangeable hydrogens is greater than the number of hydrogens in the amino acid. \nAll available hydrogens will be still exchanged.')
    elif temp_df["AA_NAME"].any() == "ASN" or temp_df["AA_NAME"].any() == "GLN":
        for ex in range(3):
            try:
                H_index = temp_df[temp_df["ATOM_NAME"] == "H"].first_valid_index()
                temp_df.loc[H_index,"ATOM_NAME"] = temp_df.loc[H_index,"ATOM_NAME"].replace("H","D")
            except KeyError:
                print('There number of exchangeable hydrogens is greater than the number of hydrogens in the amino acid. \nAll available hydrogens will be still exchanged.')
    elif temp_df["AA_NAME"].any() == "ARG":
        for ex in range(6):
            try:
                H_index = temp_df[temp_df["ATOM_NAME"] == "H"].first_valid_index()
                temp_df.loc[H_index,"ATOM_NAME"] = temp_df.loc[H_index,"ATOM_NAME"].replace("H","D")
            except KeyError:
                print('There number of exchangeable hydrogens is greater than the number of hydrogens in the amino acid. \nAll available hydrogens will be still exchanged.')
    # Append all the DataFrames into a temporary list            
    temp_list.append(temp_df)
        
# Create a new DataFrame from the temporary list
D_DF = pd.DataFrame()
D_DF = pd.concat(temp_list)

# Rename the index columns
clean_DF.index.name = "Index"
D_DF.index.name = "Index"

# Merge the exchanged hydrogens with the clean DataFrame
final_DF = clean_DF.merge(D_DF,on = "Index")

# Put the elements into a list
last_col = final_DF["ATOM_NAME_y"].tolist()
third_col = final_DF["ELEMENT_y"].tolist()

# WRITING TO PDB
# Correct space for the last and third column in a PDB file
space = '           '
# Create lists for searching and replacing the entries in last column of PDB file
third_input = H_DF["ELEMENT"].tolist()
# print(third_search_list)
third_search_list = []
for val in third_input:

    if len(val) == 1:
        new_val = str('  ')+str(val)+str('   ')
    elif len(val) == 2:
        new_val = str('  ')+str(val)+str('  ')
    elif len(val) == 3:
        new_val = str('  ')+str(val)+str(' ')
    elif len(val) == 4:
        new_val = str(' ')+str(val)+str(' ')

    third_search_list.append(new_val)

replace_list = []
H_to_list = H_DF["ATOM_NAME"].tolist()
search_list = []

# Add the correct spaces to the last & third column
for i in range(len(last_col)):
    new_replace_string = str(space)+str(last_col[i])
    new_search_string = str(space)+str(H_to_list[i])
    replace_list.append(new_replace_string)
    search_list.append(new_search_string)

third_replace_list = []
for val in last_col:
    new_val = str('  ')+val+str('   ')
    third_replace_list.append(new_val)

# Open the input PDB file
pdb_file = open(file, 'rt')

# Get user input for path and file to be outputted
path = input('Please enter the path to the directory where you want to save the file: ')
path = pdb2DF.is_path(path)
file = input('Please enter the name with which you want to save the PDB file. NOTE: please include the .pdb extension: ')
outfile_and_path = str(path)+str(file)

# Open output file for writing to
pdb_out = open(outfile_and_path, 'wt')

# Read all the lines in the PDB file
lines = pdb_file.readlines()

# Placeholder list for all lines
strings = []
# Cut off string for a PDB file created from a GRO file
cutoff_string = 'MODEL'
number_of_lines = 0
cutoff_line = 0
# Loop over lines in file to get location of cutoff
for line in lines:
    # Get number of lines in file
    number_of_lines += 1
    # Find cutoff string in the file
    if cutoff_string in line:
        cutoff_line = number_of_lines
    # Store the correct preamble for outputting    
    preamble = lines[:cutoff_line]  
    # Create an array containing the data
    data = lines[cutoff_line:number_of_lines]

# Write preamble to out file
for line in preamble:
    pdb_out.write(line)

# Replace last & third column of input PDB file with the H-D exchanged list of elements
for i in range(len(search_list)):
    new_last = data[i].replace(search_list[i],replace_list[i]).replace(third_search_list[i],third_replace_list[i])
    pdb_out.write(new_last)

# Add last two lines to PDB file
end = ["TER\n","ENDMDL"]
for line in end:
    pdb_out.write(line)

# Close input and output files
pdb_file.close()
pdb_out.close()
