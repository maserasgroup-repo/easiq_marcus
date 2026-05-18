#!/usr/bin/env python3
"""Generate Gaussian equilibrium and nonequilibrium inputs from four pairs."""

from __future__ import annotations

import shutil
import sys
from dataclasses import dataclass
from pathlib import Path


PERIODIC_TABLE = {
    1: "H",
    2: "He",
    3: "Li",
    4: "Be",
    5: "B",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    10: "Ne",
    11: "Na",
    12: "Mg",
    13: "Al",
    14: "Si",
    15: "P",
    16: "S",
    17: "Cl",
    18: "Ar",
    19: "K",
    20: "Ca",
    21: "Sc",
    22: "Ti",
    23: "V",
    24: "Cr",
    25: "Mn",
    26: "Fe",
    27: "Co",
    28: "Ni",
    29: "Cu",
    30: "Zn",
    31: "Ga",
    32: "Ge",
    33: "As",
    34: "Se",
    35: "Br",
    36: "Kr",
    37: "Rb",
    38: "Sr",
    39: "Y",
    40: "Zr",
    41: "Nb",
    42: "Mo",
    43: "Tc",
    44: "Ru",
    45: "Rh",
    46: "Pd",
    47: "Ag",
    48: "Cd",
    49: "In",
    50: "Sn",
    51: "Sb",
    52: "Te",
    53: "I",
    54: "Xe",
    55: "Cs",
    56: "Ba",
    57: "La",
    58: "Ce",
    59: "Pr",
    60: "Nd",
    61: "Pm",
    62: "Sm",
    63: "Eu",
    64: "Gd",
    65: "Tb",
    66: "Dy",
    67: "Ho",
    68: "Er",
    69: "Tm",
    70: "Yb",
    71: "Lu",
    72: "Hf",
    73: "Ta",
    74: "W",
    75: "Re",
    76: "Os",
    77: "Ir",
    78: "Pt",
    79: "Au",
    80: "Hg",
    81: "Tl",
    82: "Pb",
    83: "Bi",
    84: "Po",
    85: "At",
    86: "Rn",
    87: "Fr",
    88: "Ra",
    89: "Ac",
    90: "Th",
    91: "Pa",
    92: "U",
    93: "Np",
    94: "Pu",
    95: "Am",
    96: "Cm",
    97: "Bk",
    98: "Cf",
    99: "Es",
    100: "Fm",
    101: "Md",
    102: "No",
    103: "Lr",
    104: "Rf",
    105: "Db",
    106: "Sg",
    107: "Bh",
    108: "Hs",
    109: "Mt",
    110: "Ds",
    111: "Rg",
    112: "Cn",
    113: "Nh",
    114: "Fl",
    115: "Mc",
    116: "Lv",
    117: "Ts",
    118: "Og",
}


@dataclass
class GaussianInput:
    link0_lines: list[str]
    route_lines: list[str]
    title_lines: list[str]
    charge_mult_line: str
    tail_lines: list[str]


def split_top_level_tokens(text: str) -> list[str]:
    tokens: list[str] = []
    current: list[str] = []
    depth = 0
    for char in text:
        if char.isspace() and depth == 0:
            if current:
                tokens.append("".join(current))
                current = []
            continue
        if char == "(":
            depth += 1
        elif char == ")" and depth > 0:
            depth -= 1
        current.append(char)
    if current:
        tokens.append("".join(current))
    return tokens


def parse_gaussian_input(path: Path) -> GaussianInput:
    lines = path.read_text(encoding="utf-8").splitlines()
    idx = 0
    link0_lines: list[str] = []
    while idx < len(lines) and lines[idx].lstrip().startswith("%"):
        link0_lines.append(lines[idx])
        idx += 1

    if idx >= len(lines) or not lines[idx].lstrip().startswith("#"):
        raise ValueError(f"{path} does not look like a Gaussian input file.")

    route_lines: list[str] = []
    while idx < len(lines) and lines[idx].strip():
        route_lines.append(lines[idx])
        idx += 1

    while idx < len(lines) and not lines[idx].strip():
        idx += 1

    title_lines: list[str] = []
    while idx < len(lines) and lines[idx].strip():
        title_lines.append(lines[idx])
        idx += 1

    while idx < len(lines) and not lines[idx].strip():
        idx += 1

    if idx >= len(lines):
        raise ValueError(f"{path} is missing the charge/multiplicity line.")
    charge_mult_line = lines[idx]
    idx += 1

    while idx < len(lines) and lines[idx].strip():
        idx += 1
    while idx < len(lines) and not lines[idx].strip():
        idx += 1

    tail_lines = lines[idx:]
    return GaussianInput(
        link0_lines=link0_lines,
        route_lines=route_lines,
        title_lines=title_lines,
        charge_mult_line=charge_mult_line,
        tail_lines=tail_lines,
    )


def update_chk_lines(link0_lines: list[str], chk_name: str) -> list[str]:
    new_lines: list[str] = []
    replaced = False
    for line in link0_lines:
        stripped = line.lstrip()
        if stripped.lower().startswith("%chk="):
            prefix = line[: len(line) - len(stripped)]
            new_lines.append(f"{prefix}%chk={chk_name}")
            replaced = True
        else:
            new_lines.append(line)
    if not replaced:
        new_lines.insert(0, f"%chk={chk_name}")
    return new_lines


def add_scrf_read(token: str) -> str:
    if "=" not in token:
        return token
    key, value = token.split("=", 1)
    if not value.startswith("(") or not value.endswith(")"):
        return token
    inner = value[1:-1].strip()
    options = [part.strip() for part in inner.split(",") if part.strip()]
    if not any(part.lower() == "read" for part in options):
        options.append("read")
    return f"{key}=({','.join(options)})"


def update_route_lines(route_lines: list[str], checkpoint_mode: bool) -> list[str]:
    route_text = " ".join(line.strip() for line in route_lines)
    tokens = split_top_level_tokens(route_text)
    if not tokens:
        raise ValueError("Empty route section.")

    prefix = tokens[0]
    kept_tokens: list[str] = []
    has_geom_checkpoint = False
    has_guess_read = False

    for token in tokens[1:]:
        lowered = token.lower()
        head = lowered.split("=", 1)[0]
        if head in {"opt", "freq"}:
            continue
        if head == "scrf":
            token = add_scrf_read(token)
            lowered = token.lower()
        if lowered == "geom=checkpoint":
            has_geom_checkpoint = True
        if lowered == "guess=read":
            has_guess_read = True
        kept_tokens.append(token)

    if checkpoint_mode and not has_geom_checkpoint:
        kept_tokens.append("geom=checkpoint")
    if checkpoint_mode and not has_guess_read:
        kept_tokens.append("guess=read")

    return [f"{prefix} {' '.join(kept_tokens)}".rstrip()]


def parse_orientation_blocks(path: Path, heading: str) -> list[list[str]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    blocks: list[list[str]] = []
    idx = 0
    while idx < len(lines):
        if heading in lines[idx]:
            idx += 1
            while idx < len(lines) and set(lines[idx].strip()) != {"-"}:
                idx += 1
            if idx >= len(lines):
                break
            idx += 1
            while idx < len(lines) and set(lines[idx].strip()) != {"-"}:
                idx += 1
            if idx >= len(lines):
                break
            idx += 1
            current_block: list[str] = []
            while idx < len(lines):
                stripped = lines[idx].strip()
                if set(stripped) == {"-"}:
                    break
                parts = lines[idx].split()
                if len(parts) >= 6:
                    atomic_number = int(parts[1])
                    symbol = PERIODIC_TABLE.get(atomic_number)
                    if symbol is None:
                        raise ValueError(
                            f"Unsupported atomic number {atomic_number} in {path}."
                        )
                    x, y, z = parts[-3:]
                    current_block.append(
                        f"{symbol:<2} {float(x):>15.6f} {float(y):>15.6f} {float(z):>15.6f}"
                    )
                idx += 1
            if current_block:
                blocks.append(current_block)
        idx += 1
    return blocks


def extract_last_geometry(path: Path) -> list[str]:
    standard_blocks = parse_orientation_blocks(path, "Standard orientation:")
    if standard_blocks:
        return standard_blocks[-1]
    input_blocks = parse_orientation_blocks(path, "Input orientation:")
    if input_blocks:
        return input_blocks[-1]
    raise ValueError(
        f"Could not find 'Standard orientation' or 'Input orientation' in {path}."
    )


def trim_trailing_blank_lines(lines: list[str]) -> list[str]:
    trimmed = list(lines)
    while trimmed and not trimmed[-1].strip():
        trimmed.pop()
    return trimmed


def finalize_tail(lines: list[str], noneq_value: str) -> list[str]:
    trimmed = trim_trailing_blank_lines(lines)
    if trimmed and set(trimmed[-1].strip()) == {"*"}:
        blank_count = 2
    else:
        blank_count = 1
    return trimmed + [""] * blank_count + [f"NonEq={noneq_value}"] + [""] * 4


def render_gaussian_input(
    parsed: GaussianInput,
    route_lines: list[str],
    chk_name: str,
    geometry_lines: list[str] | None,
    noneq_value: str,
) -> str:
    link0_lines = update_chk_lines(parsed.link0_lines, chk_name)
    final_tail = finalize_tail(parsed.tail_lines, noneq_value)

    lines: list[str] = []
    lines.extend(link0_lines)
    lines.extend(route_lines)
    lines.append("")
    lines.extend(parsed.title_lines)
    lines.append("")
    lines.append(parsed.charge_mult_line)
    if geometry_lines:
        lines.extend(geometry_lines)
    lines.append("")
    lines.extend(final_tail)
    return "\n".join(lines) + "\n"


def paired_names(paths: list[Path]) -> list[tuple[Path, Path]]:
    return [(paths[i], paths[i + 1]) for i in range(0, len(paths), 2)]


def ensure_same_extension_output(input_path: Path, suffix: str) -> Path:
    return input_path.with_name(f"{input_path.stem}{suffix}{input_path.suffix}")


def print_next_steps(eq_outputs: list[tuple[Path, str]], noneq_outputs: list[Path]) -> None:
    eq_names = ", ".join(path.name for path, _ in eq_outputs)
    noneq_names = ", ".join(path.name for path in noneq_outputs)
    print("Generated Gaussian inputs successfully.")
    print(f"Submit these equilibrium inputs first: {eq_names}")
    print("Once those four jobs finish without errors, submit the nonequilibrium inputs:")
    print(noneq_names)
    print(
        "After all 8 calculations finish without errors, run get_results.py to analyze the results."
    )


def main(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv[1:]

    if len(argv) != 8:
        print(
            "Usage: fourpoints input1 aux1 input2 aux2 input3 aux3 input4 aux4",
            file=sys.stderr,
        )
        return 1

    paths = [Path(arg) for arg in argv]
    for path in paths:
        if not path.exists():
            print(f"Missing file: {path}", file=sys.stderr)
            return 1

    pairs = paired_names(paths)
    parsed_inputs = [parse_gaussian_input(input_path) for input_path, _ in pairs]

    eq_outputs: list[tuple[Path, str]] = []
    eq_chk_names: list[str] = []
    checkpoint_modes: list[bool] = []

    for parsed, (input_path, aux_path) in zip(parsed_inputs, pairs):
        eq_input_path = ensure_same_extension_output(input_path, "_eq")
        eq_chk_name = f"{input_path.stem}_eq.chk"
        eq_chk_path = input_path.with_name(eq_chk_name)
        checkpoint_mode = aux_path.suffix.lower() == ".chk"
        checkpoint_modes.append(checkpoint_mode)
        eq_chk_names.append(eq_chk_name)

        route_lines = update_route_lines(parsed.route_lines, checkpoint_mode)
        geometry_lines = None if checkpoint_mode else extract_last_geometry(aux_path)
        eq_text = render_gaussian_input(
            parsed=parsed,
            route_lines=route_lines,
            chk_name=eq_chk_name,
            geometry_lines=geometry_lines,
            noneq_value="write",
        )
        eq_outputs.append((eq_input_path, eq_text))

        if checkpoint_mode:
            shutil.copy2(aux_path, eq_chk_path)

    for output_path, content in eq_outputs:
        output_path.write_text(content, encoding="utf-8")

    partner_indices = {0: 1, 1: 0, 2: 3, 3: 2}
    noneq_outputs: list[Path] = []
    for idx, ((input_path, _), parsed) in enumerate(zip(pairs, parsed_inputs)):
        partner_idx = partner_indices[idx]
        noneq_input_path = ensure_same_extension_output(input_path, "_noneq")
        swapped_chk_name = eq_chk_names[partner_idx]
        route_lines = update_route_lines(parsed.route_lines, checkpoint_mode=True)
        geometry_lines = None
        noneq_text = render_gaussian_input(
            parsed=parsed,
            route_lines=route_lines,
            chk_name=swapped_chk_name,
            geometry_lines=geometry_lines,
            noneq_value="read",
        )
        noneq_input_path.write_text(noneq_text, encoding="utf-8")
        noneq_outputs.append(noneq_input_path)

    print_next_steps(eq_outputs, noneq_outputs)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
