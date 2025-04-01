#!/usr/bin/env python3

'''
utilities.py: Contains functions to load the input_data.pkl file and
to check the "Normal termination" line in the output files.
It also contains a function to create a Sumbitscript.sh to run
the calculations simultaneously.
'''

from pathlib import Path
from collections import namedtuple          # queues system
import pickle
import os
import re
from shutil import copy         # copy chk files

from pyssian import GaussianInFile, GaussianOutFile         # pyssian parser
from pyssian.classutils import Geometry         # pyssian xyz


def load_dict():
    """
    Load the binary file input_data.pkl.

    Returns
    -------
    input_data: dict
        Dictionary with electronic and geometrical information
    """
    with open('input_data.pkl', 'rb') as f:
        input_data = pickle.load(f)

    return input_data

def check_file_exists(filepath):
    """
    Check if the the output file exists '.out', '_solvent.out' with Normal termination
    
    Parameters
    ----------
    filepath : str
        Path to the file to check
    
    Returns
    -------
    bool
        True if the file exists, False otherwise
    """
    if os.path.exists(filepath):
        with open(filepath) as F:
            lines = F.readlines()
            if "Normal termination" in lines[-1]:
                print(f'Finished calculation {filepath.stem}')
                return True

    return False

def rename_chk(reactant, product):
    """
    Identify the .chk file and rename it to *solvent.chk

    Parameters
    ----------
    reactant : list
        Output of the ArgumentParser.parse_args function
    product : list
        Output of the ArgumentParser.parse_args function

    Returns
    -------
    Four .chk files *solvent.chk
    """

    for i in range(len(reactant)):
        oldfile_chk = Path(reactant[i].stem + '.chk')
        newfile_chk = Path(reactant[i].stem + '_solvent.chk')
        copy(oldfile_chk, newfile_chk)
    for i in range(len(product)):
        oldfile_chk = Path(product[i].stem + '.chk')
        newfile_chk = Path(product[i].stem + '_solvent.chk')
        copy(oldfile_chk, newfile_chk)
    
    return

def replace_to_chkbasis(input_data, key):
    """
    Identify the user-specified basise set with the keywords ExtraBasis, extrabasis,
    gen, Gen, GenECP, genecp and pseudopotential in the command line file and substitute
    them by the keyword chkbasis or ChkBasis

    Parameters
    ----------
    input_data : dict
        Dictionary with electronic and geometrical information
    key: str
        Key in the dictionary to check the keywords
    """
    keywords_to_replace = ['extrabasis', 'genecp', 'gen' 'pseudopotential', 'pseudo=read', 'pseudo']
    
    # Create a regex pattern to match any of the keywords, case insensitive
    pattern = re.compile('|'.join(keywords_to_replace), re.IGNORECASE)
    
    for i in range(len(input_data[key])):
        matches = pattern.findall(input_data[key][i])
        if matches:
            input_data[key][i] = pattern.sub('ChkBasis', input_data[key][i])
            input_data[key][i] = re.sub(r'(ChkBasis\s*)+', 'ChkBasis ', input_data[key][i]).strip()
            print(f'Keywords found and replaced with ChkBasis in {key}, compound {i + 1}')
        else:
            print(f'No keyword found in {key}, compound {i + 1}, the default basis set will be used')

    return input_data

def warning_basisset(input_data, key):
    """
    Identify the user-specified basis set with the keywords ExtraBasis, extrabasis,
    gen, Gen, GenECP, genecp and pseudopotential in the command line file and 
    substitute them by the keyword chkbasis or ChkBasis

    Parameters
    ----------
    input_data : dict
        Dictionary with electronic and geometrical information
    key: str
        Key in the dictionary to check the keywords
    """
    keywords_to_replace = ['extrabasis', 'genecp', 'gen' 'pseudopotential', 
                           'pseudo=read', 'pseudo']
    
    # Create a regex pattern to match any of the keywords, case insensitive
    pattern = re.compile('|'.join(keywords_to_replace), re.IGNORECASE)
    
    # Check each command line in the list
    for i, command in enumerate(input_data[key]):
        matches = pattern.findall(command)
        if matches:
            raise ValueError(f'Keywords found in command line {i} of {key}: {matches}. \n'
                             'Please, add tail file with the basis set.')
        else:
            print(f'No keyword found in command line {i + 1} of {key}, \n'
                  'The default basis set will be used.')

def wait_condition(reactant=None, product=None):
    """
    Check if '_solvent.chk' files exist and "Normal termination" line in the output files.

    Parameters
    ----------
    reactant : list
        Output of the ArgumentParser.parse_args function
    product : list
        Output of the ArgumentParser.parse_args function
    """
    chkFileReactantS, chkFileProductS = [], []
    OFileReactantS, OFileProductS = [], []

    for i in range(len(reactant)):
        chkFileReactantS.append(Path(Path(reactant[i]).stem + '_solvent.chk'))
        OFileReactantS.append([Path(reactant[i]).with_suffix(ext) for ext in ['.out', '.log']])

    for i in range(len(product)):
        chkFileProductS.append(Path(Path(product[i]).stem + '_solvent.chk'))
        #OFileProductS.append(Path(product[i]).with_suffix('.out'))
        OFileProductS.append([Path(reactant[i]).with_suffix(ext) for ext in ['.out', '.log']])

    # check if chk file exists
    chkFile = chkFileReactantS + chkFileProductS

    for file in chkFile:
        while os.path.exists(file) == False:
            print(f'Waiting for {file} file')
            break
    
    # Check if output file exists with or without the '_solvent.out' ending
    OFile = OFileReactantS + OFileProductS

    for file in OFile:
        while not check_file_exists(file):
            print(f'Waiting for {file} file')
            break

def inputsumbit(files):
    """
    Generate a file with the input files to submit the calculations.

    Parameters
    ----------
    files : Namespace
        Output of the ArgumentParser.parse_args function

    Returns
    -------
    File InputstoSubmit.txt
    """
    path_list, line_list = [], []

    input_calc = files
    path_list = [os.path.abspath(input_calc[i]) for i in range(len(input_calc))]
    
    # Define the text file
    with open('InputstoSubmit.txt', 'w') as F:
        for i in path_list:
            F.write(i + '\n')
    print('InputstoSubmit.txt file created')

def submitline(software, file):
    """
    Create a Sumbitscript.sh to run the calculations simultaneously

    Parameters
    ----------
    file : file
        Parser object

    Returns
    -------
        line: Line with the selected queue to run the calculations.
    """
    Queue = namedtuple('Queue','nprocesors memory'.split())
    QUEUES = {4:Queue(4,8),8:Queue(8, 24),
          12:Queue(12,24),20:Queue(20, 48),
          24:Queue(24,128),28:Queue(28, 128),
          36:Queue(36,192),'q4':Queue('q4', 4)}

    i_path = Path(file)

    with GaussianInFile(file) as GIF:

        GIF.read()
        nprocs = int(GIF.nprocs)
        mem = GIF.mem
    memory = int(mem[:-2])             # Assume 'xxxMB'

    if mem[-2:] == 'GB':
        memory = memory*1024    # Translate to MB in case of G
    if nprocs == 4 and software == 'g16':
        queue = QUEUES['q4']
    elif nprocs == 4 and memory < 4000:
        queue = QUEUES['q4']
    else:
        queue = QUEUES[nprocs]
    line = f'qs {software}.c{queue.nprocesors}m{queue.memory}  {i_path.name};'
    
    return line

def submitcalculation(software, files):
    """
    Check input files in the current directory or the provided folder to generate
    a SubmitScript.sh that properly sends them to their queues.
    It checks in which queue they should go according to the values
    of nprocshared and %mem. To use the generated script
    run in the cluster: 'chmod +x SubmitScript.sh; ./SubmitScript;' or
    'bash SubmitScript.sh'

    Parameters
    ----------
    software : Namespace
        Output of the ArgumentParser.parse_args function

    Returns
    -------
    Script SubmitScript.sh
    """
    path_list, line_list = [], []

    input_calc = files
    path_list = [os.path.abspath(input_calc[i]) for i in range(len(input_calc))]
    
    for i in path_list:
        line = submitline(software, i)
        line_list.append(line)
    
    # Now create the submit script
    SubmitLines = ['#!/bin/bash\n',
                   '# Automated Submit Script\n'
                   f'BASEDIR=$PWD;']
    
    for i in line_list:
        SubmitLines.append(i)
    SubmitLines.append('cd ${BASEDIR};')
    SubmitLines.append('echo "Finished submiting calculations";')

    with open('SubmitScript.sh', 'w') as F:
        F.write('\n'.join(SubmitLines))


