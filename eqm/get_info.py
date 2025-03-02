#!/usr/bin/env python3

"""
get_info.py extracts the information from the Gaussian output files of the reactants
and products and save it in a dictionary. This dictionary is used in the next script.
"""

import argparse
from pathlib import Path    # to add the names
import os

from pyssian import GaussianOutFile         # pyssian parser
from pyssian.classutils import Geometry         # pyssian xyz parser
import pickle         # save a dictionary document

from .utilities import inputsumbit, rename_chk, warning_basisset

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument('Reactant', nargs=2, type=Path, help='Output files: Reactant1, Reactant2')
parser.add_argument('Product', nargs=2, type=Path, help='Output files: Product1, Product2')
parser.add_argument('solvent',help="""Solvent to use in the calculation
                    written as it would be written in Gaussian.""")
parser.add_argument('smodel',choices=['smd','pcm','cpcm'], help="""
                    Solvent model to use""")
#parser.add_argument('software', choices=['g09', 'g16'], default='g09',
#                    nargs='?', help="""Gaussian version program""")
parser.add_argument('-Chk', '--chk', 
                    choices=['yes', 'no'], 
                    help=""".chk files to use in the calculation. 
                    By default input files are created.""", 
                    default=None)
parser.add_argument('-Tail', '--tail', help="""Tail File that contains the extra
                    options, such as pseudopotentials, basis sets. Applied to all files""",
                    default=None)
parser.add_argument('-Tail1', '--tail1', help="""Tail1 File that contains the extra
                    options, such as pseudopotentials, basis sets for reactant1 and product1""",
                    default=None)
parser.add_argument('-Tail2', '--tail2', help="""Tail2 File that contains the extra
                    options, such as pseudopotentials, basis sets for reactant2 and product2""",
                    default=None)

def prepare_commandline(link):
    """
    This function remove the keywords: opt, freq and scrf of the commandline

    Parameters
    ----------
    link : object
        Parser object

    Returns
    -------
    raw_commandline : list
        A list with the raw commands of the calculation
    """
    commandline = link.commandline
    split_commandline = commandline.split(" ")

    new_commandline = [x.casefold() for x in split_commandline
                        if 'opt' not in x.casefold()
                        if 'freq' not in x.casefold()
                        if 'scrf' not in x.casefold()]
                        
    raw_commandline = " ".join(str(_) for _ in new_commandline)

    return raw_commandline

def extract_info(reactant, product):
    """
    Extract the output file information (name, charge, spin, xyz, atoms, coords, nprocs, mem).

    Parameters
    ----------
    reactant : list
        Output of the ArgumentParser.parse_args function
    product : list
        Output of the ArgumentParser.parse_args function

    Returns
    -------
    input_data, OFileProduct, OFileReactant
        Dictionary 'input_data'.
        Name of the new input files.
    """
    names_prod, charge_prod = list(), list()
    spin_prod, xyz_prod = list(), list()
    atoms_prod, coords_prod = list(), list()
    nprocs_prod, mem_prod = list(), list()
    names_react, charge_react = list(), list()
    spin_react, xyz_react = list(), list()
    atoms_prod, coords_prod = list(), list()
    commline_prod, commline_react = list(), list()

    for i in range(len(product)):

        with GaussianOutFile(product[i], parselist=[1, 101, 202]) as GOF:
            GOF.read()

        # Get the Link1, 101 and 202 of the Gaussian OutputFile
        l1 = GOF.get_links(1)[0]
        l101 = GOF.get_links(101)[0]
        l202 = GOF.get_links(202)[-1]

        geom = Geometry.from_L202(l202)
        xyz_prod.append(str(geom))
        atoms_prod = geom.atoms
        coords_prod = geom.coordinates
        charge_prod.append(l101.charge)
        spin_prod.append(l101.spin)
        nprocs_prod.append(l1.nprocs)
        mem_prod.append(l1.mem)
        raw_cl_prod = prepare_commandline(l1)
        commline_prod.append(raw_cl_prod)
        names_prod.append(product[i])

    for i in range(len(reactant)):

        with GaussianOutFile(reactant[i], parselist=[1, 101, 202]) as GOF:
            GOF.read()

        # Get the Link101 and 202 of the Gaussian OutputFile
        l101 = GOF.get_links(101)[0]
        l202 = GOF.get_links(202)[-1]

        geom = Geometry.from_L202(l202)
        xyz_react.append(str(geom))
        atoms_react = geom.atoms
        coords_react = geom.coordinates
        charge_react.append(l101.charge)
        spin_react.append(l101.spin)
        raw_cl_react = prepare_commandline(l1)
        commline_react.append(raw_cl_react)
        names_react.append(reactant[i])

    input_data = {}
    
    keys = ("names_prod", "nprocs_prod", "commandline_prod",
            "mem_prod", "charge_prod", "spin_prod",
            "atoms_prod", "xyz_prod", "coords_prod",
            "names_react", "commandline_react",
            "charge_react", "spin_react",
            "atoms_react", "xyz_react", "coords_react")
    values = [names_prod, nprocs_prod, commline_prod,
            mem_prod, charge_prod, spin_prod,
            atoms_prod, xyz_prod, coords_prod,
            names_react, commline_react,
            charge_react, spin_react,
            atoms_react, xyz_react, coords_react]

    for a, b in zip(keys, values):

        input_data[a] = b

    f = open("input_data.pkl","wb")
    pickle.dump(input_data,f)
    f.close()

    return input_data

def compute_chk(reactant, product, dict_args, input_data, FileStruct_chk):
    """
    Get the information from the input_data dictionary and build four input
    files (*solvent.in) with the required spin, charge and chk specified.

    Parameters
    ----------
    reactant : list
        Output of the ArgumentParser.parse_args function
    product : list
        Output of the ArgumentParser.parse_args function
    dict_args : dict
        Dictionary with electronic and geometrical information
    input_data : dict
        Electronic and geometrical information of the reactants and products.
        Computational details of the calculations
    FileStruct_chk : string
        Structure of the new inputfile

    Returns
    -------
    Four input files *solvent.in
    """
    namesIn_prodS, namesIn_reactS = [], []

    for i in range(len(product)):

        h_solvent =input_data["commandline_prod"][i] + f' scrf=({smodel},solvent={solvent})'
        newfile = Path(product[i].stem + '_solvent.in')
        newfile_chk = Path(product[i].stem + '_solvent.chk')

        tail_content = dict_args['tail']
        if i == 0:
            tail_content += dict_args['tail1']
        elif i == 1:
            tail_content += dict_args['tail2']

        with open(newfile, 'w') as F:

            txt = FileStruct_chk.format(
                                    nprocs=input_data["nprocs_prod"][i],
                                    mem=input_data["mem_prod"][i],
                                    chk=newfile_chk,
                                    Header=h_solvent,
                                    charge=input_data["charge_prod"][i],
                                    spin=input_data["spin_prod"][i],
                                    Title=newfile,
                                    Coords=input_data["xyz_prod"][i],
                                    Tail=tail_content #dict_args['tail']
                                    )
            F.write(txt)

        namesIn_prodS.append(newfile)

    for i in range(len(reactant)):

        h_solvent =input_data["commandline_prod"][i] + f' scrf=({smodel},solvent={solvent})'
        newfile = Path(reactant[i].stem + '_solvent.in')
        newfile_chk = Path(reactant[i].stem + '_solvent.chk')
        
        tail_content = dict_args['tail']
        if i == 0:
            tail_content +=  dict_args['tail1']
        elif i == 1:
            tail_content +=  dict_args['tail2']

        with open(newfile, 'w') as F:

            txt = FileStruct_chk.format(
                                    nprocs=input_data["nprocs_prod"][i],
                                    mem=input_data["mem_prod"][i],
                                    chk=newfile_chk,
                                    Header=h_solvent,
                                    charge=input_data["charge_react"][i],
                                    spin=input_data["spin_react"][i],
                                    Title=newfile,
                                    Coords=input_data["xyz_react"][i],
                                    Tail=tail_content #dict_args['tail']
                                    )
            F.write(txt)

        namesIn_reactS.append(newfile)

    return namesIn_prodS, namesIn_reactS

if __name__ == "__main__":

    args = parser.parse_args()
    reactant = args.Reactant
    product = args.Product
    #software = args.software
    solvent = args.solvent
    chk = args.chk
    smodel = args.smodel

    dict_args = vars(args) # convert to dict

    # Process tail, tail1, and tail2 arguments
    for tail_arg in ['tail', 'tail1', 'tail2']:
        if dict_args[tail_arg] is None:
            dict_args[tail_arg] = ''
        else:
            with open(os.path.abspath(dict_args[tail_arg]), 'r') as F:
                Aux = [line.strip() for line in F]
            while len(Aux) >= 1 and not Aux[0]:
                _ = Aux.pop(0)
            while len(Aux) >= 1 and not Aux[-1]:
                _ = Aux.pop(-1)
            dict_args[tail_arg] = '\n'.join(Aux)

    FileStruct = '%nprocshared={nprocs}\n%mem={mem}\n{Header}\n\n{Title}\n\n{charge} {spin}\n{Coords}\n\n{Tail}\n\n\n'
    FileStruct_chk = '%nprocshared={nprocs}\n%mem={mem}\n%chk={chk}\n{Header}\n\n{Title}\n\n{charge} {spin}\n{Coords}\n\n{Tail}\n\n\n'

    input_data = extract_info(reactant, product)

    # check if the .chk file is provided
    if chk == 'yes':
        rename_chk(reactant, product) # rename the .chk file to *solvent.chk
    
    else:
        # check if the tail file is provided and if it contains the keywords
        try: 
            if dict_args['tail'] == '' and dict_args['tail1'] == '' and dict_args['tail2'] == '':
                print("No tail file provided")
                # check if input_data contains the keywords
                warning_basisset(input_data, "commandline_react")
                warning_basisset(input_data, "commandline_prod")
        
                namesIn_prodS, namesIn_reactS = compute_chk(reactant, product, dict_args, input_data, FileStruct_chk)
                filestorun = namesIn_prodS + namesIn_reactS
                #submitcalculation(software, filestorun)
                inputsumbit(filestorun) # generate a list with the input names

            else:
                print(f'Tail file provided: {dict_args[tail_arg]}')
                namesIn_prodS, namesIn_reactS = compute_chk(reactant, product, dict_args, input_data, FileStruct_chk)
                filestorun = namesIn_prodS + namesIn_reactS
                #submitcalculation(software, filestorun)
                inputsumbit(filestorun) # generate a list with the input names

        except ValueError as e:
            print(f'Error: {e}')
            exit(1) # exit the program with an error



