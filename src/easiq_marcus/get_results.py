#!/usr/bin/env python3

"""
get_results.py calculates free-energy barriers for Single-Electron Transfer (SET)
and Energy transfer (EnT) events using reported approximations based on Marcus theory.
"""

from __future__ import annotations

import argparse
import math
import os
import re
from pathlib import Path


HARTREE_TO_KCAL_MOL = 627.509
WriteOutput = print
number_fmt = "{: 03.9f}"
largest_value = len(number_fmt.format(10000))
value_fmt = f"{{: ^{largest_value}}}"
line_fmt = ""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)

    requiredNamed = parser.add_argument_group("required named arguments")
    requiredNamed.add_argument(
        "-r", "--Reactant", help="Reactants", nargs="+", type=Path, required=True
    )
    requiredNamed.add_argument("-p", "--Product", help="Products", nargs="+", type=Path)

    #parser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode")
    parser.add_argument(
        "-fp",
        "--four_points",
        action="store_true",
        help="4-points DFT application of Marcus theory assuming parabolas of the same width for reactants and products",
    )
    parser.add_argument(
        "-fpa",
        "--four_points_asymm",
        action="store_true",
        help="4-points DFT application of Marcus theory assuming parabolas of different width for reactants and products",
    )
    parser.add_argument(
        "-hs", "--hard_sphere", action="store_true", help="application of the Hard-Sphere model"
    )
    parser.add_argument(
        "-shs",
        "--simpl_hard_sphere",
        type=float,
        choices=[95, 96, 97, 98, 99],
        help='application of the simplified Hard-Sphere model. A value for the "A" constant is expected.',
    )
    parser.add_argument(
        "-v_d", "--volume_donor", type=float, help="Volume of the donor molecule in Angstrom"
    )
    parser.add_argument(
        "-v_a",
        "--volume_acceptor",
        type=float,
        help="Volume of the acceptor molecule in Angstrom",
    )
    parser.add_argument(
        "-r_d", "--radius_donor", type=float, help="Radius of the donor molecule in Angstrom"
    )
    parser.add_argument(
        "-r_a",
        "--radius_acceptor",
        type=float,
        help="Radius of the acceptor molecule in Angstrom",
    )
    parser.add_argument(
        "-d", "--dielectric_cte", type=float, help="Dielectric constant of the media"
    )
    parser.add_argument(
        "-d_opt",
        "--dielectric_opt",
        type=float,
        help="Dielectric optical constant of the media",
    )
    parser.add_argument(
        "-bde",
        "--bde_value",
        type=float,
        help="Bond Dissociation Free Energy (kcal/mol) of the bond breaking during SET/EnT event within the frame of the Savéant's model. Compatible with both Hard-Sphere and Simplified Hard-Sphere approximations.",
    )
    parser.add_argument(
        "-method",
        help="If the Final Potential energy is not the Energy of the 'SCF Done:' the target Potential Energy has to be specified, otherwise defaults to the energy of the 'SCF Done:'",
        choices=["oniom", "mp2", "mp2scs", "mp4", "ccsdt"],
        default="default",
        type=lambda x: x.lower(),
    )
    parser.add_argument(
        "-O",
        "--OutFile",
        help="File to write the Data. If it exists, the data will be appended. If not specified it will print to the console",
        default=None,
    )
    #parser.add_argument(
    #    "-details",
    #    "--details_lambda",
    #    action="store_true",
    #    help="Gives additional information about the calculated lambdas. Specifically, it splits lambda into inner (solute) and outer (solvent) contributions.",
    #)
    return parser


FLOAT_RE = r"[-+]?\d+(?:\.\d+)?(?:[DdEe][-+]?\d+)?"


def convert_gaussian_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def write_2_file(outfile: str):
    def writer(text=""):
        with open(outfile, "a", encoding="utf-8") as handle:
            handle.write(f"{text}\n")

    return writer


def parse_potential_energy(text: str, method: str) -> float | None:
    patterns: dict[str, list[str]] = {
        "default": [rf"SCF Done:\s+E\([^)]+\)\s*=\s*({FLOAT_RE})"],
        "oniom": [rf"ONIOM:\s+extrapolated energy\s*=\s*({FLOAT_RE})"],
        "mp2": [rf"EUMP2\s*=\s*({FLOAT_RE})", rf"EMP2\s*=\s*({FLOAT_RE})"],
        "mp2scs": [
            rf"SCS-?MP2[^=\n]*=\s*({FLOAT_RE})",
            rf"E\(SCS-?MP2\)[^=\n]*=\s*({FLOAT_RE})",
        ],
        "mp4": [
            rf"UMP4\(SDTQ\)\s*=\s*({FLOAT_RE})",
            rf"UMP4\(SDQ\)\s*=\s*({FLOAT_RE})",
            rf"EMP4\(SDTQ\)\s*=\s*({FLOAT_RE})",
            rf"EMP4\(SDQ\)\s*=\s*({FLOAT_RE})",
        ],
        "ccsdt": [
            rf"CCSD\(T\)\s*=\s*({FLOAT_RE})",
            rf"E\(CCSD\(T\)\)\s*=\s*({FLOAT_RE})",
        ],
    }

    search_order = patterns.get(method, patterns["default"])
    for pattern in search_order:
        matches = re.findall(pattern, text)
        if matches:
            return convert_gaussian_float(matches[-1])
    return None


def parse_thermochemistry(text: str) -> tuple[float, float, float]:
    primary_patterns = {
        "zpe": rf"Sum of electronic and zero-point Energies=\s*({FLOAT_RE})",
        "enthalpy": rf"Sum of electronic and thermal Enthalpies=\s*({FLOAT_RE})",
        "gibbs": rf"Sum of electronic and thermal Free Energies=\s*({FLOAT_RE})",
    }
    fallback_patterns = {
        "zpe": rf"Zero-point correction=\s*({FLOAT_RE})",
        "enthalpy": rf"Thermal correction to Enthalpy=\s*({FLOAT_RE})",
        "gibbs": rf"Thermal correction to Gibbs Free Energy=\s*({FLOAT_RE})",
    }

    values: dict[str, float] = {}
    for key in ("zpe", "enthalpy", "gibbs"):
        matches = re.findall(primary_patterns[key], text)
        if matches:
            values[key] = convert_gaussian_float(matches[-1])
            continue

        matches = re.findall(fallback_patterns[key], text)
        if not matches:
            raise IndexError(f"Could not find {key} thermochemistry value")
        values[key] = convert_gaussian_float(matches[-1])

    return values["zpe"], values["enthalpy"], values["gibbs"]


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
    r_energy = sum(float(i) for i in dict_prod["G"]) - sum(float(i) for i in dict_react["G"])
    return r_energy * HARTREE_TO_KCAL_MOL


class Tools:
    """
    Set of Tools to calculate energy parameters relevant to the SET/EnT process
    """

    def __init__(self, args):
        self.reactant = [i for i in args.Reactant]
        self.product = [i for i in (args.Product or [])]
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
        """Calculate the radius of a sphere from its volume"""
        self.r_donor = (3 * self.v_donor / (4 * math.pi)) ** (1 / 3)
        return self.r_donor

    def radiu_volume_acceptor(self):
        """Calculate the radius of a sphere from its volume"""
        self.r_accep = (3 * self.v_accept / (4 * math.pi)) ** (1 / 3)
        return self.r_accep

    def get_U(self, InFilepath):
        """
        Parse the output files and get the potential energy
        """
        text = Path(InFilepath).read_text(encoding="utf-8", errors="replace")

        if str(InFilepath).endswith("_noneq.out"):
            match = re.findall(
                rf"After PCM corrections, the energy is\s+({FLOAT_RE})", text
            )
            if match:
                U = convert_gaussian_float(match[-1])
                return value_fmt.format(number_fmt.format(U))

        U = parse_potential_energy(text, self._method)
        if U is None:
            raise RuntimeError(f"Potential Energy not found in file {InFilepath}")
        return value_fmt.format(number_fmt.format(U))

    def parser_energy(self, Files):
        """
        Parse the output files
        """
        outfile_list = []
        U_list, Z_list, H_list, G_list = [], [], [], []

        for outfile in Files:
            if not outfile:
                WriteOutput("")
                continue

            outfile_path = os.path.abspath(outfile)
            text = Path(outfile_path).read_text(encoding="utf-8", errors="replace")

            U = parse_potential_energy(text, self._method)
            if U is not None:
                U = value_fmt.format(number_fmt.format(U))
            elif self._verbose:
                raise RuntimeError(f"Potential Energy not found in file {outfile}")
            else:
                U = value_fmt.format("")

            try:
                Z, H, G = parse_thermochemistry(text)
            except IndexError as e:
                if self._verbose:
                    raise e
                Z = value_fmt.format("")
                H = value_fmt.format("")
                G = value_fmt.format("")
            else:
                Z = value_fmt.format(number_fmt.format(Z))
                H = value_fmt.format(number_fmt.format(H))
                G = value_fmt.format(number_fmt.format(G))

            outfile_list.append(Path(outfile).stem)
            U_list.append(U)
            Z_list.append(Z)
            H_list.append(H)
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
        """
        if bool(self._bde_value) is True:
            reorg = self._bde_value + reorg
            print(reorg)

        intrins_ene = reorg / 4
        ene_barrier = intrins_ene * (1 + r_energy / (4 * intrins_ene)) ** 2

        return ene_barrier, intrins_ene, reorg

    def four_points(self, react_dict, prod_dict, reactant, product):
        """
        Application of the 4-points Marcus theory approximation
        """
        r_energy = reaction_energy(react_dict, prod_dict)

        def reorg_energy(react_dict, prod_dict):
            U_Rnoneq_lst, U_Req_lst = [], []
            U_Pnoneq_lst, U_Peq_lst = [], []

            for outfile in reactant:
                outfileRnoneq_path = outfile.with_name(f"{outfile.stem}_noneq.out")
                outfileReq_path = outfile.with_name(f"{outfile.stem}_eq.out")
                U_Rnoneq_lst.append(float(self.get_U(outfileRnoneq_path)))
                U_Req_lst.append(float(self.get_U(outfileReq_path)))

            for outfile in product:
                outfilePnoneq_path = outfile.with_name(f"{outfile.stem}_noneq.out")
                outfilePeq_path = outfile.with_name(f"{outfile.stem}_eq.out")
                U_Pnoneq_lst.append(float(self.get_U(outfilePnoneq_path)))
                U_Peq_lst.append(float(self.get_U(outfilePeq_path)))

            reorgT_41 = (
                sum(U_Rnoneq_lst) - sum(U_Req_lst)
            ) * HARTREE_TO_KCAL_MOL
            reorgT_32 = (
                sum(U_Pnoneq_lst) - sum(U_Peq_lst)
            ) * HARTREE_TO_KCAL_MOL
            reorgT = (reorgT_41 + reorgT_32) / 2

            reorgN_41 = reorgT_41
            reorgN_32 = reorgT_32
            reorgN = (reorgN_41 + reorgN_32) / 2

            return reorgT, reorgN, reorgT_41, reorgT_32

        reorgT, reorgN, reorgT_41, reorgT_32 = reorg_energy(react_dict, prod_dict)

        def activ_energy(reorgN, reorgT, reorgT_41, reorgT_32, r_energy):
            act_energy = (reorgT + r_energy) ** 2 / (4 * reorgT)
            print("\n")
            print("--" * 21, "4-points approximation", "--" * 21)
            line_fmt2 = "  {:^20} {:^20} {:^20} {:^20} {:^20} "
            line_fmt3 = "  {:^20} {:^20} {:^20} {:^20} {:^20} "
            line_fmt4 = "  {:^20} {:^20} {:^20} {:^20} {:^20} "
            WriteOutput(
                line_fmt2.format(
                    "lambda R",
                    "lambda P",
                    "lambda tot.",
                    "Reaction Free",
                    "Free energy",
                )
            )
            WriteOutput(
                line_fmt3.format(
                    "(kcal/mol)",
                    "(kcal/mol)",
                    "(kcal/mol)",
                    "energy (kcal/mol)",
                    "barrier (kcal/mol)",
                )
            )
            print("--" * 54)
            WriteOutput(
                line_fmt4.format(
                    round(reorgT_41, 1),
                    round(reorgT_32, 1),
                    round(reorgT, 1),
                    round(r_energy, 1),
                    round(act_energy, 1),
                )
            )
            print("--" * 54)
            print(
                "lambda R: reorganization energy measured on the reactants parabola"
            )
            print("lambda P: reorganization energy measured on the products parabola")
            print(
                "lambda tot.: total reorganization energy assuming the same width for both reactants and products parabolas; lambda tot. = (lambdaR+lambdaP)/2"
            )
            print(
                "Note: if lambdaR and lambdaP differ significantly, please consider using the asymmetric version of the 4-points approximation (-fpa flag instead of -fp)"
            )

        activ_energy(reorgN, reorgT, reorgT_41, reorgT_32, r_energy)

    def hard_sphere(self, react_dict, prod_dict):
        """
        Application of the Hard-Sphere model
        """
        r_energy = reaction_energy(react_dict, prod_dict)

        def reorg_energy():
            cte_factor = (
                (
                    6.02214076 * 10 ** (23) * ((1.602176634 * 10 ** (-19)) ** 2)
                )
                * 0.000239
                * 10 ** 10
            ) / (4 * math.pi * 8.854187817620389 * 10 ** (-12))

            if bool(self.v_accept) is True:
                self.r_accep = self.radiu_volume_acceptor()
            else:
                self.r_accep = float(self.r_accep)

            if bool(self.v_donor) is True:
                self.r_donor = self.radiu_volume_donor()
            else:
                self.r_donor = float(self.r_donor)

            r_total = self.r_accep + self.r_donor
            radius = ((self.r_donor * 2) ** -1 + (self.r_accep * 2) ** -1 - (r_total) ** -1)
            dielec = (self.diel_opt ** -1) - (self.diel ** -1)
            reorg = cte_factor * radius * dielec

            act_energy, intrins_ene, reorg = self.mh_equation(r_energy, reorg)

            print("\n")
            print("--" * 13, "Hard-Sphere model", "--" * 13)
            line_fmt2 = "  {:^20} {:^25} {:^20} "
            line_fmt3 = "  {:^20} {:^25} {:^20} "
            line_fmt4 = "  {:^20} {:^25} {:^20} "
            WriteOutput(line_fmt2.format("Reorganization", "Reaction Free", "Free energy"))
            WriteOutput(
                line_fmt3.format(
                    "energy (kcal/mol)", "energy (kcal/mol)", "barrier (kcal/mol)"
                )
            )
            print("--" * 36)
            WriteOutput(line_fmt4.format(round(reorg, 1), round(r_energy, 1), round(act_energy, 1)))
            print("--" * 36)

        reorg_energy()

    def simpl_hard_sphere(self, react_dict, prod_dict):
        """
        Application of the simplified Hard-Sphere model
        """
        r_energy = reaction_energy(react_dict, prod_dict)

        def reorg_energy():
            if bool(self.v_accept) is True:
                self.r_accep = self.radiu_volume_acceptor()
            else:
                self.r_accep = float(self.r_accep)

            if bool(self.v_donor) is True:
                self.r_donor = self.radiu_volume_donor()
            else:
                self.r_donor = float(self.r_donor)

            r_total = self.r_accep + self.r_donor
            radius = ((self.r_donor * 2) ** -1 + (self.r_accep * 2) ** -1 - (r_total) ** -1)
            reorg = self.shs * radius

            act_energy, intrins_ene, reorg = self.mh_equation(r_energy, reorg)

            print("\n")
            print("--" * 9, "Simplified Hard-Sphere model, A =", round(self.shs), "--" * 8)
            line_fmt2 = "  {:^20} {:^25} {:^20} "
            line_fmt3 = "  {:^20} {:^25} {:^20} "
            line_fmt4 = "  {:^20} {:^25} {:^20} "
            WriteOutput(line_fmt2.format("Reorganization", "Reaction Free", "Free energy"))
            WriteOutput(
                line_fmt3.format(
                    "energy (kcal/mol)", "energy (kcal/mol)", "barrier (kcal/mol)"
                )
            )
            print("--" * 36)
            WriteOutput(line_fmt4.format(round(reorg, 1), round(r_energy, 1), round(act_energy, 1)))
            print("--" * 36)
            print(
                'Note: The values of "A" have been empirically determined for the DMF solvent. Applying this method to reactions in other solvents may lead to unrealistic results.'
            )

        reorg_energy()

    def four_points_asymm(self, react_dict, prod_dict, reactant, product):
        """
        Application of the 4-points Marcus theory approximation considering parabolas of different width
        """
        r_energy = reaction_energy(react_dict, prod_dict)

        def reorg_energy(react_dict, prod_dict):
            U_Rnoneq_lst, U_Req_lst = [], []
            U_Pnoneq_lst, U_Peq_lst = [], []

            for outfile in reactant:
                outfileRnoneq_path = outfile.with_name(f"{outfile.stem}_noneq.out")
                outfileReq_path = outfile.with_name(f"{outfile.stem}_eq.out")
                U_Rnoneq_lst.append(float(self.get_U(outfileRnoneq_path)))
                U_Req_lst.append(float(self.get_U(outfileReq_path)))

            for outfile in product:
                outfilePnoneq_path = outfile.with_name(f"{outfile.stem}_noneq.out")
                outfilePeq_path = outfile.with_name(f"{outfile.stem}_eq.out")
                U_Pnoneq_lst.append(float(self.get_U(outfilePnoneq_path)))
                U_Peq_lst.append(float(self.get_U(outfilePeq_path)))

            reorgT_41 = (
                sum(U_Rnoneq_lst) - sum(U_Req_lst)
            ) * HARTREE_TO_KCAL_MOL
            reorgT_32 = (
                sum(U_Pnoneq_lst) - sum(U_Peq_lst)
            ) * HARTREE_TO_KCAL_MOL
            return reorgT_41, reorgT_32

        reorgT_41, reorgT_32 = reorg_energy(react_dict, prod_dict)

        def activ_energy(reorgT_41, reorgT_32, r_energy):
            act_energy = reorgT_41 * (
                (
                    (-reorgT_32)
                    + ((reorgT_41 * reorgT_32) + (r_energy * (reorgT_41 - reorgT_32))) ** 0.5
                )
                / (reorgT_41 - reorgT_32)
            ) ** 2
            print("\n")
            print("--" * 14, "Asymmetric 4-points approximation", "--" * 11)
            line_fmt2 = "  {:^20} {:^20} {:^20} {:^20} "
            line_fmt3 = "  {:^20} {:^20} {:^20} {:^20} "
            line_fmt4 = "  {:^20} {:^20} {:^20} {:^20} "
            WriteOutput(
                line_fmt2.format("lambda R", "lambda P", "Reaction Free", "Free energy")
            )
            WriteOutput(
                line_fmt3.format(
                    "(kcal/mol)",
                    "(kcal/mol)",
                    "energy (kcal/mol)",
                    "barrier (kcal/mol)",
                )
            )
            print("--" * 41 + "---")
            WriteOutput(
                line_fmt4.format(
                    round(reorgT_41, 1),
                    round(reorgT_32, 1),
                    round(r_energy, 1),
                    round(act_energy, 1),
                )
            )
            print("--" * 41 + "---")
            print("lambda R: reorganization energy measured on the reactants parabola")
            print("lambda P: reorganization energy measured on the products parabola")

        activ_energy(reorgT_41, reorgT_32, r_energy)


def main(argv: list[str] | None = None) -> int:
    global WriteOutput
    global number_fmt
    global value_fmt
    global line_fmt

    parser = build_parser()
    args = parser.parse_args(argv)
    obj = Tools(args)
    print("· Reading output files...")

    files_opt = [obj.reactant, obj.product]
    print(
        "Number of reactants: ",
        len(obj.reactant),
        "Number of products:",
        len(obj.product),
    )

    if obj.OutFile is not None:
        OutFile = os.path.abspath(obj.OutFile)
        WriteOutput = write_2_file(OutFile)
    else:
        WriteOutput = print

    largest_filename = max((len(str(path.stem)) for path in obj.reactant), default=1)
    name_format = f"{{: <{largest_filename}}}"

    number_fmt = "{: 03.9f}"
    largest_value = len(number_fmt.format(10000))
    value_fmt = f"{{: ^{largest_value}}}"
    spacer = " "
    line_fmt = spacer.join([name_format] + [value_fmt] * 4)

    reactant_dict = obj.parser_energy(files_opt[0])
    product_dict = obj.parser_energy(files_opt[1])

    if obj.fourpoints:
        obj.four_points(reactant_dict, product_dict, obj.reactant, obj.product)

    if obj.fourpointsasymm:
        obj.four_points_asymm(reactant_dict, product_dict, obj.reactant, obj.product)

    if obj.hs:
        if obj.r_donor is None and obj.v_donor is None:
            raise RuntimeError("Neither Donor radius or volume not found")
        if obj.r_accep is None and obj.v_accept is None:
            raise RuntimeError("Neither Acceptor radius or volume not found")
        if obj.diel is None and obj.diel_opt is None:
            raise RuntimeError(
                "Dielectric constant or dielectric optical constant of the solvent not found"
            )
        obj.hard_sphere(reactant_dict, product_dict)

    if obj.shs:
        if not isinstance(obj.shs, (int, float)):
            raise RuntimeError("Factor not found")
        if obj.r_donor is None and obj.v_donor is None:
            raise RuntimeError("Neither Donor radius or volume not found")
        if obj.r_accep is None and obj.v_accept is None:
            raise RuntimeError("Neither Acceptor radius or volume not found")
        obj.simpl_hard_sphere(reactant_dict, product_dict)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
