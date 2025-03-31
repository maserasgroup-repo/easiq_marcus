#!/usr/bin/env python3

"""
4points.py script that creates input files required to obtain the value of
the free energy of the outer-sphere electron transfer activation barrier using 
the approximation of the 4-points in two parabolas.
"""

import argparse
import os           # manipulate paths
from shutil import copy         # copy chk files
from pathlib import Path    # to add the names
from collections import namedtuple          # queues system

from .get_info import extract_info
from .utilities import load_dict, wait_condition, submitcalculation, inputsumbit, replace_to_chkbasis

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument('Reactant', nargs=2, type=Path, help='Output files: Reactant1, Reactant2')
parser.add_argument('Product', nargs=2, type=Path, help='Output files: Product1, Product2')
parser.add_argument('solvent',help="""Solvent to use in the calculation
                    written as it would be written in Gaussian.""")
parser.add_argument('smodel',choices=['smd','pcm','cpcm'], help="""
                    Solvent model to use""")
parser.add_argument('-UseBasisSet', '--usebasisset', 
                    choices=['yes', 'no'], 
                    default='no',
                    help="""Maintain the basis set described in the previous input line.
                    By default, the basis set are retrived from the chk file.""")

def noneq_points(reactant, product, input_data, FileStruct_noneq):
    """
    Design input files 'non_eq' with the previous information and copy the
    .chk files derived from the previous script get_info.py

    Parameters
    ----------
    reactant : list
        Output of the ArgumentParser.parse_args function
    product : list
        Output of the ArgumentParser.parse_args function
    dict_args : dict
        Dictionary with electronic and geometrical information
    FileStruct_input2 : string
        Structure of the input file

    Returns
    -------
    nameInReactnoneq, nameInProdnoneq
    """
    nameInReactnoneq, nameInProdnoneq = [], []

    outage = '\nNonEq=write \n\n--link1--\n'

    for i in range(len(reactant)):

        h_input3 =  input_data["commandline_prod"][i] + \
            f' scrf=({smodel},solvent={solvent},read) geom=check guess=read'
        h_input31 =  input_data["commandline_prod"][i] + \
            f' scrf=({smodel},solvent={solvent},read,ExternalIteration) geom=check guess=read'
                        
        oldfile_chk = Path(Path(reactant[i]).stem + '_solvent.chk')
        newfile_chk = Path(Path(reactant[i]).stem + '_noneq.chk')
        copy(oldfile_chk, newfile_chk)
        newfile = Path(reactant[i].stem + '_noneq.in')
        
        with open(newfile,'w') as F:
        
            txt = FileStruct_noneq.format(
                                    nprocs=input_data["nprocs_prod"][i],
                                    mem=input_data["mem_prod"][i],
                                    chk= newfile_chk,
                                    Header=h_input3,
                                    Title=newfile,
                                    charge=input_data["charge_prod"][i],
                                    spin=input_data["spin_prod"][i],
                                    #Tail=dict_args['tail'],
                                    Division=outage,
                                    Header2=h_input31,
                                    Title2=newfile,
                                    charge2=input_data["charge_prod"][i],
                                    spin2=input_data["spin_prod"][i],
                                    Tail_input='\nNonEq=read',
                                    #Tail2=dict_args['tail']
                                    )
            F.write(txt)

            nameInReactnoneq.append(newfile)

    for i in range(len(product)):

        h_input3 =  input_data["commandline_prod"][i] + \
            f' scrf=({smodel},solvent={solvent},read) geom=check guess=read'
        h_input31 =  input_data["commandline_prod"][i] + \
            f' scrf=({smodel},solvent={solvent},read,ExternalIteration) geom=check guess=read'

        oldfile_chk = Path(Path(product[i]).stem + '_solvent.chk')
        newfile_chk = Path(Path(product[i]).stem + '_noneq.chk')
        copy(oldfile_chk, newfile_chk)
        newfile = Path(product[i].stem + '_noneq.in')
        
        with open(newfile,'w') as F:
        
            txt = FileStruct_noneq.format(
                                    nprocs=input_data["nprocs_prod"][i],
                                    mem=input_data["mem_prod"][i],
                                    chk= newfile_chk,
                                    Header=h_input3,
                                    Title=newfile,
                                    charge=input_data["charge_react"][i],
                                    spin=input_data["spin_react"][i],
                                    #Tail=dict_args['tail'],
                                    Division=outage,
                                    Header2=h_input31,
                                    Title2=newfile,
                                    charge2=input_data["charge_react"][i],
                                    spin2=input_data["spin_react"][i],
                                    Tail_input='\nNonEq=read',
                                    #Tail2=dict_args['tail']
                                    )
            F.write(txt)

            nameInProdnoneq.append(newfile)

    return nameInReactnoneq, nameInProdnoneq

if __name__ == "__main__":

    args = parser.parse_args()
    reactant = args.Reactant
    product = args.Product
    solvent = args.solvent
    smodel = args.smodel

    dict_args = vars(args) # convert to dict

    FileStruct = '%nprocshared={nprocs}\n%mem={mem}\n{Header}\n\n{Title}\n\n{charge} {spin}\n{Coords}\n\n{Tail}\n\n\n\n'
    FileStruct_noneq = '%nprocshared={nprocs}\n%mem={mem}\n%chk={chk}\n{Header}\n\n{Title}\n\n{charge} {spin}\n{Division}%nprocshared={nprocs}\n%mem={mem}\n%chk={chk}\n{Header2}\n\n{Title2}\n\n{charge2} {spin2}\n{Tail_input}\n\n\n\n\n'

    wait_condition(reactant, product)

    # change the input_data with tail_basisset if dict_args['tail'] is not None
    if dict_args['usebasisset'] == 'no':
        if os.path.exists("input_data.pkl"):
            input_data = load_dict()
        else:
            input_data = extract_info(reactant, product)

        input_data = replace_to_chkbasis(input_data, "commandline_prod")
        input_data = replace_to_chkbasis(input_data, "commandline_react")
    else:
        if os.path.exists("input_data.pkl"):
            input_data = load_dict()
        else:
            input_data = extract_info(reactant, product)

        print(f"Basis set are retrieved from the input line of the previous calculation. "
        f"Check if the information is provided.")
        input_data = replace_to_chkbasis(input_data, "commandline_prod")
        input_data = replace_to_chkbasis(input_data, "commandline_react")

    nameInReactnoneq, nameInProdnoneq = noneq_points(reactant, product, input_data, FileStruct_noneq)
    files_to_send = nameInReactnoneq + nameInProdnoneq
    #submitcalculation(software, files_to_send)
    inputsumbit(files_to_send)      # generate a list with the input names