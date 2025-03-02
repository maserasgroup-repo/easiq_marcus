#!/usr/bin/env python3

"""
get_results.py calculates free-energy barriers for Single-Electron Transfer (SET)
and Energy transfer (EnT) events using reported approximations based on Marcus theory.
"""

from pathlib import Path
import os
import argparse

from pyssian import GaussianOutFile
from pyssianutils.functions import thermochemistry, potential_energy, write_2_file

parser = argparse.ArgumentParser(description=__doc__)

requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-r', '--Reactant', help='Reactants', nargs='+', type=Path,
                    required=True)
requiredNamed.add_argument('-p', '--Product', help='Products', nargs='+', type=Path)

parser.add_argument('-v', '--verbose', action='store_true', help="""Verbose mode""")
parser.add_argument('-fp', '--four_points', action='store_true',
                    help="""4-points DFT application of Marcus theory assuming parabolas 
                    of the same width for reactants and products""")
parser.add_argument('-fpa', '--four_points_asymm', action='store_true',
                    help="""4-points DFT application of Marcus theory assuming parabolas 
                    of different width for reactants and products""")
parser.add_argument('-hs', '--hard_sphere', action='store_true', help="""application of 
                    the Hard-Sphere model""")
parser.add_argument('-shs', '--simpl_hard_sphere', type=float,
                    choices=[95, 96, 97, 98, 99],
                    help="""application of the simplified Hard-Sphere model. A value 
                    for the "A" constant is expected.""")
parser.add_argument('-v_d', '--volume_donor', type=float, help="""Volume of the donor
                    molecule in Angstrom""")
parser.add_argument('-v_a', '--volume_acceptor', type=float, help="""Volume of the acceptor
                    molecule in Angstrom""")
parser.add_argument('-r_d', '--radius_donor', type=float, help="""Radius of the donor
                    molecule in Angstrom""")
parser.add_argument('-r_a', '--radius_acceptor', type=float, help="""Radius of the acceptor
                    molecule in Angstrom""")
parser.add_argument('-d', '--dielectric_cte', type=float, help="""Dielectric constant of the media""")
parser.add_argument('-d_opt', '--dielectric_opt', type=float, help="""Dielectric optical
                    constant of the media""")
parser.add_argument('-bde', '--bde_value', type=float, help="""Bond Dissociation
                    Free Energy (kcal/mol) of the bond breaking during SET/EnT event 
                    within the frame of the Savéant's model. Compatible with both 
                    Hard-Sphere and Simplified Hard-Sphere approximations.""")
parser.add_argument('-method', help="""If the Final Potential energy is not the
                    Energy of the 'SCF Done:' the target Potential Energy
                    has to be specified, otherwise defaults to the energy of the
                    'SCF Done:""",
                    choices=['oniom','mp2','mp2scs','mp4','ccsdt'],
                    default='default',type=lambda x: x.lower())
parser.add_argument('-O','--OutFile',help="""File to write the Data. 
                    The data will be appended to the file. 
                    If not specified it will print to the console""",default=None)

def reaction_energy(dict_react, dict_prod):
    """
    Obtain the reaction Gibbs free energy
    Parameters
    ----------
    dict_react : dict
        Dict with all the energy parameters of the reactant output files
    dict_prod : dict
        Dict with all the energy parameters of the product output files
    Returns
    -------
    r_energy : energy of the reaction
    """
    r_energy = sum([float(i) for i in dict_prod["G"]]) - \
                sum([float(i) for i in dict_react["G"]])
    r_energy = r_energy * 627.509

    return r_energy

class Tools:
    """
    Set of Tools to calculate energy parameters relevant to the SET/EnT process

    Attributes
    ----------
    reactant : list
        List of the reactant output files
    product : list
        List of the product output files
    v_donor : float
        Volume of the donor molecule in Angstroms
    v_accept : float
        Volume of the acceptor molecule in Angstroms
    r_donor : float
        Radius of the donor molecule in Angstroms
    r_accep : float
        Radius of the acceptor molecule in Angstroms
    diel : float
        Dielectric constant of the media
    diel_opt : float
        Dielectric optical constant of the media
    _bde_value : float
        Bond Dissociation Free Energy (kcal/mol) of the bond breaking during SET/EnT event 
        within the frame of the Savéant's model
    fourpoints : bool
        4-points DFT application of Marcus theory assuming parabolas of the same width 
        for reactants and products
    fourpointsasymm : bool
        4-points DFT application of Marcus theory assuming parabolas of different width 
        for reactants and products
    hs : bool
        application of the Hard-Sphere model
    shs : float
        application of the simplified Hard-Sphere model. A value for the "A" constant is expected. 
    OutFile : str
        File to write the Data. If it exists, the data will be appended. If
        not specified it will print to the console
    _method : str
        If the Final Potential energy is not the Energy of the 'SCF Done:' the
        target Potential Energy has to be specified, otherwise defaults to the
        energy of the 'SCF Done:"
    _verbose : bool
        Verbose mode
    """
    def __init__(self):

        args = parser.parse_args()

        self.reactant = [i for i in args.Reactant]
        self.product = [i for i in args.Product]
        self.v_donor = args.volume_donor
        self.v_accept = args.volume_acceptor
        self.r_donor = args.radius_donor
        self.r_accep = args.radius_acceptor
        self.diel = args.dielectric_cte
        self.diel_opt = args.dielectric_opt
        self._bde_value = args.bde_value
        self.fourpoints = args.four_points
        self.fourpointsasymm = args.four_points_asymm
        self.hs = args.hard_sphere
        self.shs = args.simpl_hard_sphere
        self.OutFile = args.OutFile
        self._method = args.method
        self._verbose = args.verbose

    def radiu_volume_donor(self):
        """
        Calculate the radius of a sphere from its volume
        """
        self.r_donor = (3*self.v_donor/(4*3.141592653589793))**(1/3)
        return self.r_donor

    def radiu_volume_acceptor(self):
        """
        Calculate the radius of a sphere from its volume
        """
        self.r_accep = (3*self.v_accept/(4*3.141592653589793))**(1/3)
        return self.r_accep

    def get_U(self, InFilepath):
        """
        Parser the output files and get the potential energy

        Parameters
        ----------
        InFilepath : Path
            Path of the output files given as arguments

        Returns
        -------
        U : int
            DFT energy of each calculation
        InFilepath: Path
            Path of the output files
        """
        with GaussianOutFile(InFilepath,[1,120,502,508,716,804,913,9999]) as GOF:
        
            if str(InFilepath).endswith('_noneq.out'):

                with open(InFilepath) as f:
                    for line in f:
                        if 'After PCM corrections, the energy is' in line:
                            U = float(line.split()[6])
                            U = value_fmt.format(number_fmt.format(float(U)))
                            break
            else:
                GOF.read()
                U = potential_energy(GOF,self._method)

        return U

    def parser_energy(self, Files):
        """
        Parser the output files

        Parameters
        ----------
        files : list
            List of all the outputs to parser

        Returns
        -------
        out_dict: dict
            Dictionary with all the energy parameters of the output files
        """
        outfile_list =  []
        U_list, Z_list, H_list, G_list = [], [], [], []

        for outfile in Files:

            if not outfile: #In the case of an empty filename, write an empty line
                WriteOutput('')
                continue
            outfile_path = os.path.abspath(outfile)

            with GaussianOutFile(outfile_path,[1,120,502,508,716,804,913,9999]) as GOF:
                GOF.read()
                U = potential_energy(GOF,self._method)
            if U is not None:
                U = value_fmt.format(number_fmt.format(U))
            elif self._verbose:
                raise RuntimeError(f'Potential Energy not found in file {outfile}')
            else:
                U = value_fmt.format('')
            try:
                Z,H,G = thermochemistry(GOF)
            except IndexError as e:
                if self._verbose:
                    raise e
                Z = value_fmt.format('')
                H = value_fmt.format('')
                G = value_fmt.format('')

            else:
                Z = value_fmt.format(number_fmt.format(Z))
                H = value_fmt.format(number_fmt.format(H))
                G = value_fmt.format(number_fmt.format(G))

            outfile_list.append(outfile.stem)
            U_list.append(U)
            Z_list.append(Z)
            H_list.append(G)
            G_list.append(G)
            out_dict = {}
            keys = ("name", "U", "ZPE", "H", "G")
            values = [outfile_list, U_list, Z_list, H_list, G_list]
            for a, b in zip(keys, values):
                out_dict[a] = b

        return out_dict

    def mh_equation(self, r_energy, reorg):
        """
        Calculate the energy barrier of the SET/EnT process

        Parameters
        ----------
        r_energy : float
            Energy of the reaction (kcal/mol)
        reorg : float
            Reorganization energy (kcal/mol)
        """
        if bool(self._bde_value) is True:
            reorg = self._bde_value + reorg

        intrins_ene = reorg/4
        ene_barrier = intrins_ene * (1 + r_energy/(4*intrins_ene))**2

        return ene_barrier, intrins_ene, reorg


    def four_points(self, react_dict, prod_dict, reactant, product):
        """
        Application of the 4-points Marcus theory approximation (assuming the same with for both parabolas)

        Parameters
        ----------
        react_dict : dict
            Dict with all the energy parameters of the reactant output files
        prod_dict : dict
            Dict with all the energy parameters of the product output files
        reactant : list
            Output of the ArgumentParser.parse_args function
        product : list
            Output of the ArgumentParser.parse_args function
        """
        r_energy = reaction_energy(react_dict, prod_dict)

        def reorg_energy(react_dict, prod_dict):
            """
            Calculation of the reorganization parameters

            Returns
            -------
            reorgT : float
                Total reorganization energy
            reorgN : float
                Nuclear reorganization energy
            reorgT_41 : float
                Total reorganization energy (point 4 - point 1)
            reorgT_32 : float
                Total reorganization energy (point 3 - point 2)
            """
            react_energyU = sum([float(i) for i in react_dict["U"]])
            prod_energyU = sum([float(i) for i in prod_dict["U"]])

            U_Rnoneq_lst, U_Req_lst = list(), list()
            U_Pnoneq_lst, U_Peq_lst = list(), list()

            for outfile in reactant:

                outfileRnoneq_path = Path(outfile.stem + '_noneq.out')
                U_Rnoneq_lst.append(float(self.get_U(outfileRnoneq_path)))

            for outfile in product:

                outfilePnoneq_path = Path(outfile.stem + '_noneq.out')
                U_Pnoneq_lst.append(float(self.get_U(outfilePnoneq_path)))
            
            # reorganization total energy
            reorgT_41 = (U_Pnoneq_lst[0] + U_Pnoneq_lst[1] - react_energyU) * 627.509
            reorgT_32 = (U_Rnoneq_lst[0] + U_Rnoneq_lst[1] - prod_energyU )* 627.509

            reorgT = (reorgT_41 + reorgT_32)/2

            return reorgT, reorgT_41, reorgT_32

        reorgT, reorgT_41, reorgT_32 = reorg_energy(react_dict, prod_dict)

        def activ_energy(reorgT, reorgT_41, reorgT_32, r_energy):
            """
            Provide the activation energy of the SET/EnT process
            """
            act_energy = (reorgT + r_energy)**2/(4*reorgT)
            print('\n')
            print('--'*21,'4-points approximation','--'*21)
            line_fmt2 = "  {:^20} {:^20} {:^20} {:^20} {:^20} "
            line_fmt3 = "  {:^20} {:^20} {:^20} {:^20} {:^20} "
            line_fmt4 = "  {:^20} {:^20} {:^20} {:^20} {:^20} "
            WriteOutput(line_fmt2.format('lambda R',
                        'lambda P', 'lambda tot.', 'Reaction Free',
                        'Free energy'))
            WriteOutput(line_fmt3.format('(kcal/mol)',
                        '(kcal/mol)', '(kcal/mol)', 'energy (kcal/mol)',
                        'barrier (kcal/mol)'))
            print('--'*54)
            WriteOutput(line_fmt4.format(round(reorgT_41,1),
                        round(reorgT_32,1), round(reorgT,1),
                        round(r_energy,1), round(act_energy,1)))
            print('--'*54)
            print('lambda R: reorganization energy measured on the reactants parabola')
            print('lambda P: reorganization energy measured on the products parabola')
            print(f'lambda tot.: total reorganization energy assuming the same width for both reactants' f'and products parabolas; lambda tot. = (lambdaR+lambdaP)/2')
            print('Disclaimer: if lambdaR and lambdaP differ significantly, please'
            f'consider using the asymmetric version of the 4-points approximation (-fpa flag instead of -fp)')

        activ_energy(reorgT, reorgT_41, reorgT_32, r_energy)


    def hard_sphere(self, react_dict, prod_dict):
        """
        Application of the Hard-Sphere model

        Parameters
        ----------
        react_dict : dict
            Dict with all the energy parameters of the reactant output files
        prod_dict : dict
            Dict with all the energy parameters of the product output files
        """
        r_energy = reaction_energy(react_dict, prod_dict)

        def reorg_energy():
            """
            Calculation of the reorganization parameters
            """
            cte_factor = ((6.02214076*10**(23)*((1.602176634*10**(-19))**2))*0.000239\
                        *10**10)/ (4*3.141592653589793*8.854187817620389*10**(-12))

            if bool(self.v_accept) is True:
                self.r_accep = self.radiu_volume_acceptor()
            else:
                self.r_accep = float(self.r_accep)

            if bool(self.v_donor) is True:
                self.r_donor = self.radiu_volume_donor()
            else:
                self.r_donor = float(self.r_donor)

            r_total = self.r_accep + self.r_donor
            radius = ((self.r_donor*2)**-1 + (self.r_accep*2)**-1 - (r_total)**-1)
            dielec = ((self.diel_opt)**-1 - (self.diel)**-1)
            reorg =  cte_factor*radius*dielec

            act_energy, intrins_ene, reorg = self.mh_equation(r_energy, reorg)

            print('\n')
            print('--'*13,'Hard-Sphere model','--'*13)
            line_fmt2 = "  {:^20} {:^25} {:^20} "
            line_fmt3 = "  {:^20} {:^25} {:^20} "
            line_fmt4 = "  {:^20} {:^25} {:^20} "
            WriteOutput(line_fmt2.format('Reorganization', 'Reaction Free', 'Free energy'))
            WriteOutput(line_fmt3.format('energy (kcal/mol)', 'energy (kcal/mol)',
                        'barrier (kcal/mol)'))
            print('--'*36)
            WriteOutput(line_fmt4.format(round(reorg,1), round(r_energy,1),
                        round(act_energy,1)))
            print('--'*36)

        reorg_energy()

    def simpl_hard_sphere(self, react_dict, prod_dict):
        """
        Application of the simplified Hard-Shpere model

        Parameters
        ----------
        react_dict : dict
            Dict with all the energy parameters of the reactant output files
        prod_dict : dict
            Dict with all the energy parameters of the product output files
        """
        r_energy = reaction_energy(react_dict, prod_dict)

        def reorg_energy():
            """
            Calculation of the reorganization parameters
            """
            if bool(self.v_accept) is True:
                self.r_accep = self.radiu_volume_acceptor()
            else:
                self.r_accep = float(self.r_accep)

            if bool(self.v_donor) is True:
                self.r_donor = self.radiu_volume_donor()
            else:
                self.r_donor = float(self.r_donor)

            r_total = self.r_accep + self.r_donor
            radius = ((self.r_donor*2)**-1 + (self.r_accep*2)**-1 - (r_total)**-1)
            reorg = self.shs*radius

            act_energy, intrins_ene, reorg = self.mh_equation(r_energy, reorg)

            print('\n')
            print('--'*9,'Simplified Hard-Sphere model, A =',round(self.shs),'--'*8)
            line_fmt2 = "  {:^20} {:^25} {:^20} "
            line_fmt3 = "  {:^20} {:^25} {:^20} "
            line_fmt4 = "  {:^20} {:^25} {:^20} "
            WriteOutput(line_fmt2.format('Reorganization', 'Reaction Free', 'Free energy'))
            WriteOutput(line_fmt3.format('energy (kcal/mol)', 'energy (kcal/mol)',
                        'barrier (kcal/mol)'))
            print('--'*36)
            WriteOutput(line_fmt4.format(round(reorg,1), round(r_energy,1),
                        round(act_energy,1)))
            print('--'*36)
            print('Note: The values of "A" have been empirically determined for the DMF solvent. \
                  Applying this method to reactions in other solvents may lead to unrealistic results.')
            
        reorg_energy()

    def four_points_asymm(self, react_dict, prod_dict, reactant, product):
        """
        Application of the 4-points Marcus theory approximation considering parabolas of different width

        Parameters
        ----------
        react_dict : dict
            Dict with all the energy parameters of the reactant output files
        prod_dict : dict
            Dict with all the energy parameters of the product output files
        reactant : list
            Output of the ArgumentParser.parse_args function
        product : list
            Output of the ArgumentParser.parse_args function
        """
        r_energy = reaction_energy(react_dict, prod_dict)

        def reorg_energy(react_dict, prod_dict):
            """
            Calculation of the reorganization parameters

            Returns
            -------
            reorgT_41 : float
                Total reorganization energy (point 4 - point 1)
            reorgT_32 : float
                Total reorganization energy (point 3 - point 2)
            """
            react_energyU = sum([float(i) for i in react_dict["U"]])
            prod_energyU = sum([float(i) for i in prod_dict["U"]])

            U_Rnoneq_lst, U_Req_lst = list(), list()
            U_Pnoneq_lst, U_Peq_lst = list(), list()

            for outfile in reactant:

                # reactant equilibrium and non-equilibrium
                outfileRnoneq_path = Path(outfile.stem + '_noneq.out')
                U_Rnoneq_lst.append(float(self.get_U(outfileRnoneq_path)))

            for outfile in product:

                # product equilibrium and non-equilibrium
                outfilePnoneq_path = Path(outfile.stem + '_noneq.out')
                U_Pnoneq_lst.append(float(self.get_U(outfilePnoneq_path)))
            
            # reorganization total energy
            #reorg reactants surf
            reorgT_41 = (U_Pnoneq_lst[0] + U_Pnoneq_lst[1] - react_energyU) * 627.509
            #reorg products surf
            reorgT_32 = (U_Rnoneq_lst[0] + U_Rnoneq_lst[1] - prod_energyU )* 627.509
            
            return reorgT_41, reorgT_32

        reorgT_41, reorgT_32 = reorg_energy(react_dict, prod_dict)

        def activ_energy(reorgT_41, reorgT_32, r_energy):
            """
            Provide the activation energy of the SET/EnT process
            """
            act_energy = reorgT_41*(((reorgT_32*(-1))+((reorgT_41*reorgT_32)+(r_energy*(reorgT_41-reorgT_32)))**(1/2))/(reorgT_41-reorgT_32))**2
            print('\n')
            print('--'*14,'Asymmetric 4-points approximation','--'*11)
            line_fmt2 = "  {:^20} {:^20} {:^20} {:^20} "
            line_fmt3 = "  {:^20} {:^20} {:^20} {:^20} "
            line_fmt4 = "  {:^20} {:^20} {:^20} {:^20} "
            WriteOutput(line_fmt2.format('lambda R',
                        'lambda P', 'Reaction Free',
                        'Free energy'))
            WriteOutput(line_fmt3.format('(kcal/mol)',
                        '(kcal/mol)', 'energy (kcal/mol)',
                        'barrier (kcal/mol)'))
            print('--'*41+'---')
            WriteOutput(line_fmt4.format(round(reorgT_41,1),
                        round(reorgT_32,1),
                        round(r_energy,1), round(act_energy,1)))
            print('--'*41+'---')
            print('lambda R: reorganization energy measured on the reactants parabola')
            print('lambda P: reorganization energy measured on the products parabola')

        activ_energy(reorgT_41, reorgT_32, r_energy)


if __name__ == "__main__":

    print('· Reading output files...')

    obj = Tools()

    # get the files
    files_opt = [obj.reactant, obj.product]
    print('Number of reactants: ', len(obj.reactant), \
        'Number of products:', len(obj.product))

    if obj.OutFile is not None:
        # create the file if it does not exist
        if not os.path.exists(obj.OutFile):
            # create the file
            with open(obj.OutFile, 'w') as f:
                f.write('')

        OutFile = os.path.abspath(obj.OutFile)
        WriteOutput = write_2_file(OutFile)

        # print the content of the file
        with open(OutFile, 'r') as f:
            print(f.read())

    else:
        WriteOutput = print

    for i in range(len(obj.reactant)):

        n = largest_filename = len(str(obj.reactant[i].stem))
        name_format = f'{{: <{n}}}'

    # Values format
    number_fmt = '{: 03.9f}'
    largest_value = len(number_fmt.format(10000))
    value_fmt = f'{{: ^{largest_value}}}'
    spacer = ' '
    line_fmt = spacer.join([name_format,]+[value_fmt,]*4)

    # Methods to run
    reactant_dict = obj.parser_energy(files_opt[0])
    product_dict = obj.parser_energy(files_opt[1])
    #print(reactant_dict, product_dict)

    # 4points
    if obj.fourpoints:
        obj.four_points(reactant_dict, product_dict, obj.reactant, obj.product)

    # 4points asymm
    if obj.fourpointsasymm:
        obj.four_points_asymm(reactant_dict, product_dict, obj.reactant, obj.product)

    # Hard-Sphere
    if obj.hs:
        if obj.r_donor is None and obj.v_donor is None:
            raise ValueError('Neither Donor radius or volume not found')
        elif obj.r_accep is None and obj.v_accept is None:
            raise ValueError('Neither Acceptor radius or volume not found')
        elif obj.diel is None and obj.diel_opt is None:
            raise ValueError('Dielectric constant or dielectric optical constant of the solvent not found')
        else:
            obj.hard_sphere(reactant_dict, product_dict)

    # Simplified Hard-Sphere
    if obj.shs:
        if not isinstance(obj.shs, (int, float)):
            raise ValueError('Factor not found')
        elif obj.r_donor is None and obj.v_donor is None:
            raise ValueError('Neither Donor radius or volume not found')
        elif obj.r_accep is None and obj.v_accept is None:
            raise ValueError('Neither Acceptor radius or volume not found')
        else:
            obj.simpl_hard_sphere(reactant_dict, product_dict)
