import numpy as np
import pandas as pd
import csv 
import os.path


def is_path(filepath):
    is_path = os.path.isdir(filepath)
    while is_path == False:
        filepath = input('Invalid directory %s. Please enter again: ' %(filepath))
        is_path = os.path.isdir(filepath)
    return filepath

def is_file(filepath, filename):
    is_file = os.path.isfile(str(filepath)+str(filename))
    while is_file == False:
        filename = input('Invalid filename %s. Please enter again: ' %(filename))
        is_file = os.path.isfile(str(filepath)+str(filename))
    return filename

def is_int(input_string, input_var):

    input_valid = False

    while input_valid == False:
        try:
            input_var = int(input_var)
            input_valid = True
        except ValueError:
            print('Did not understand the command. Try again.')
            input_var = input(input_string)

    return input_var



def pdb2DF():
    pdb_path = input('Please enter the full path to the folder containing the PDB file: ')
    pdb_path = is_path(pdb_path)

    file = input('Please enter the name of the PDB file: ')
    file = is_file(pdb_path,file)

    fullpath = str(pdb_path)+str(file)

    cutoff_string = 'MODEL'

    # Open file
    f = open(fullpath, 'r')

    # Read in all lines from file to a list
    lines = f.readlines()
    number_of_lines = 0
    cutoff_line = 0

    # Loop over lines in file to get location of cutoff
    for line in lines:
        # Get number of lines in file
        number_of_lines += 1

        if cutoff_string in line:

            cutoff_line = number_of_lines
    
    # Create an array containing only x and y data
    rows = lines[cutoff_line:number_of_lines]

    split_rows = np.char.split(rows)

    output = str(pdb_path)+'LYS_data_clean.csv'

    is_anti = input('Is the file a protein or antibody? \n 1. Antibody \n 2. Protein \nPlease enter: ')
    is_anti = is_int('Is the file a protein or antibody? \n 1. Antibody \n 2. Protein \nPlease enter: ', is_anti)

    if is_anti == 1:

        columns = [['TYPE', 'ATOM_NR', 'ELEMENT', 'AA_NAME','CHAIN','AA_NR', 'X', 'Y', 'Z', 'ONE', 'ZERO', 'ATOM_NAME']]

        with open(output,"w+", newline='') as my_csv:
            csvWriter = csv.writer(my_csv, delimiter=',')
            csvWriter.writerows(columns)
            csvWriter.writerows(split_rows)

    elif is_anti == 2:

        columns = [['TYPE', 'ATOM_NR', 'ELEMENT', 'AA_NAME','AA_NR','X', 'Y', 'Z', 'ONE', 'ZERO', 'ATOM_NAME']]

        with open(output,"w+", newline='') as my_csv:
            csvWriter = csv.writer(my_csv, delimiter=',')
            csvWriter.writerows(columns)
            csvWriter.writerows(split_rows)


    return output, fullpath, is_anti
