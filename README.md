# easiq-marcus

`easiq-marcus` is a small Python package to analyze single-electron transfer (SET) and related processes using Marcus theory from Gaussian 16 derived files.

The current workflow is intentionally simple and includes two commands:

- `fourpoints` prepares the equilibrium and nonequilibrium Gaussian input files required for a four-point Marcus analysis.
- `get_results` parses Gaussian output files and reports reaction free energies, reorganization energies, and free-energy barriers. Note that beyond the four-points approach, this code also enables the application of empirical hard-sphere models.

## Installation

Install from a local clone:

```bash
git clone https://github.com/maserasgroup-repo/easiq_marcus.git
cd easiq-marcus
pip install .
```

Install directly from GitHub:

```bash
pip install git+https://github.com/maserasgroup-repo/easiq_marcus.git
```

For development:

```bash
pip install -e .
```

## Commands

The package installs two commands:

- `fourpoints`
- `get_results`

You can verify the installation with:

```bash
fourpoints
get_results -h
```

## What the code does

`easiq-marcus` supports two main approaches to apply Marcus theory:


### 1. Four-point Marcus analysis

This workflow calculates reorganization energies in the package. It uses:

- equilibrium Gaussian calculations for the reactant and product charge/spin states
- nonequilibrium single-point calculations built from those states

From those files, `get_results` extracts:

- the reaction free energy from Gaussian total free energies
- `lambda_R` and `lambda_P` from electronic energies
- the barrier from either the symmetric (`-fp`) or asymmetric (`-fpa`) Marcus expression

For the standard symmetric Marcus picture, the free-energy barrier is

```text
ΔG‡ = (λ + ΔGr°)^2 / (4λ)
```

where:

- `ΔGr°` is the reaction free energy
- `λ` is the reorganization energy

In the four-point route, the code computes reactant- and product-side reorganization energies:

```text
λR = VR(qP) - VR(qR)
λP = VP(qR) - VP(qP)
```

and, in the symmetric approximation, uses

```text
λ = (λR + λP) / 2
```

If the two sides are significantly different, the asymmetric Marcus expression can be used:

```text
ΔG‡ = λR * [(-λP + sqrt(λR λP + (λR - λP) ΔGr°)) / (λR - λP)]^2
```

The code exposes these two routes through:

- `-fp` for the symmetric four-point treatment
- `-fpa` for the asymmetric four-point treatment

  
### 2. Hard-sphere style estimates

`get_results` also supports:

- the hard-sphere model with `-hs`
- the simplified hard-sphere model with `-shs`
- the Savéant correction for dissociative transfer events with `-bde`

These routes are useful when you want a fast estimate without running the full four-point cycle. They require solvent and size information such as radii or volumes, and they are available as options in the example workflow below. For the present README, the worked example is focused on the four-point route because it is the most fully documented and reproducible one.

For the hard-sphere model, the solvent contribution is estimated as

```text
λS = (NA e^2 / 4πϵ0) * (1/2rD + 1/2rA - 1/R) * (1/ϵopt - 1/ϵ)
```

with `R = rD + rA`.

For the simplified hard-sphere model, the code uses

```text
λS = A * (1/2rD + 1/2rA - 1/R)
```

where `A` is one of `95, 96, 97, 98, 99`.

For dissociative electron- or energy-transfer events, the Savéant correction is available as

```text
λ = λS + BDFE
```

where `BDFE` is provided with the `-bde` option in kcal/mol.

## General usage

### `fourpoints`

Generate the equilibrium and nonequilibrium Gaussian inputs from four input/auxiliary pairs:

```bash
fourpoints input1.gjf aux1.out input2.gjf aux2.out input3.gjf aux3.out input4.gjf aux4.out
```

Each auxiliary file can be either:

- a Gaussian output file containing the optimized geometry
- a Gaussian checkpoint file (`.chk`)

This choice matters:

- if you pass an `.out` file, `fourpoints` extracts the final geometry from that output and writes it explicitly into the generated `*_eq.gjf` file
- if you pass a `.chk` file, `fourpoints` copies that checkpoint to the corresponding `*_eq.chk` file and modifies the generated `*_eq.gjf` route section so Gaussian reads geometry and wavefunction from the checkpoint using `geom=checkpoint` and `guess=read`

In other words, checkpoint usage is not automatic just because `.chk` files are present in the folder. The `.chk` file is used only if you explicitly provide it as the auxiliary argument to `fourpoints`.

The command writes:

- `*_eq.gjf` files for the equilibrium jobs
- `*_noneq.gjf` files for the nonequilibrium jobs

### `get_results`

Analyze Gaussian output files and compute Marcus-theory quantities:

```bash
get_results -r reactant1.out reactant2.out -p product1.out product2.out -fp
```

Useful options:

- `-fp` for the symmetric four-point approximation
- `-fpa` for the asymmetric four-point approximation
- `-hs` for the hard-sphere model
- `-shs <95|96|97|98|99>` for the simplified hard-sphere model
- `-bde <value>` to apply the Savéant correction on top of `-hs` or `-shs`
- `-O results.txt` to append output to a file instead of printing

For four-point analysis, the corresponding `*_eq.out` and `*_noneq.out` files must be present in the same directory as the main output files.

## Worked example: SET from a BI(·) radical to a Re-bipyridine complex

The repository includes a complete example in [example](example/). This example describes a single-electron transfer from a BI radical donor to a Re-based catalyst.

The initial equilibrium-state files are:

- `bi_rad.gjf`
- `bi_rad.out`
- `bi_cat.gjf`
- `bi_cat.out`
- `re-bpy-co3-cl.gjf`
- `re-bpy-co3-cl.out`
- `re-bpy-co3-cl-1e.gjf`
- `re-bpy-co3-cl-1e.out`

The example folder also includes the generated `*_eq.gjf`, `*_noneq.gjf`, `*_eq.out`, and `*_noneq.out` files, so you can inspect the full workflow end to end.

### Step 1. Move into the example directory

```bash
cd example
```

### Step 2. Generate the four-point Gaussian inputs

For this example, the four-point inputs were generated with Gaussian output files as the auxiliary geometry source:

```bash
fourpoints bi_rad.gjf bi_rad.out \
           bi_cat.gjf bi_cat.out \
           re-bpy-co3-cl.gjf re-bpy-co3-cl.out \
           re-bpy-co3-cl-1e.gjf re-bpy-co3-cl-1e.out
```

This creates:

- `bi_rad_eq.gjf`
- `bi_rad_noneq.gjf`
- `bi_cat_eq.gjf`
- `bi_cat_noneq.gjf`
- `re-bpy-co3-cl_eq.gjf`
- `re-bpy-co3-cl_noneq.gjf`
- `re-bpy-co3-cl-1e_eq.gjf`
- `re-bpy-co3-cl-1e_noneq.gjf`

The `*_eq.gjf` files should be submitted first. Once those calculations finish successfully, submit the `*_noneq.gjf` files.

If you prefer to drive the workflow from checkpoint files instead, you can do:

```bash
fourpoints bi_rad.gjf bi_rad.chk \
           bi_cat.gjf bi_cat.chk \
           re-bpy-co3-cl.gjf re-bpy-co3-cl.chk \
           re-bpy-co3-cl-1e.gjf re-bpy-co3-cl-1e.chk
```

In that case, the generated `*_eq.gjf` files will read geometry and wavefunction from the copied checkpoint files rather than from coordinates written into the input file.

### Step 3. Collect the Gaussian outputs

After the Gaussian calculations finish, the four-point analysis requires the following output files:

- `bi_rad.out`
- `bi_cat.out`
- `re-bpy-co3-cl.out`
- `re-bpy-co3-cl-1e.out`
- `bi_rad_eq.out`
- `bi_rad_noneq.out`
- `bi_cat_eq.out`
- `bi_cat_noneq.out`
- `re-bpy-co3-cl_eq.out`
- `re-bpy-co3-cl_noneq.out`
- `re-bpy-co3-cl-1e_eq.out`
- `re-bpy-co3-cl-1e_noneq.out`

All of these files are included in the example folder included in the repository.

### Step 4. Compute the symmetric four-point barrier

Run:

```bash
get_results -r bi_rad.out re-bpy-co3-cl.out \
            -p bi_cat.out re-bpy-co3-cl-1e.out -fp
```

This computes:

- the reaction free energy from the total Gaussian free energies
- `lambda_R` and `lambda_P`
- the symmetric Marcus barrier using the average reorganization energy

For this example, the calculation gives approximately:

- `ΔG_r° = -6.0 kcal/mol`
- `lambda_R = 48.4 kcal/mol`
- `lambda_P = 47.1 kcal/mol`
- `lambda = 47.7 kcal/mol`
- `ΔG‡ ≈ 9.1 kcal/mol`

### Step 5. Compute the asymmetric four-point barrier

Run:

```bash
get_results -r bi_rad.out re-bpy-co3-cl.out \
            -p bi_cat.out re-bpy-co3-cl-1e.out -fpa
```

This uses `lambda_R` and `lambda_P` separately in the asymmetric Marcus expression.

For this example, the asymmetric result is also about:

- `ΔG‡ ≈ 9.1 kcal/mol`

The symmetric and asymmetric barriers are almost identical here because `lambda_R` and `lambda_P` are very similar.

### Step 6. Optional: hard-sphere models

If you want a simpler Marcus estimate without the explicit four-point cycle, you can instead use:

```bash
get_results -r bi_rad.out re-bpy-co3-cl.out \
            -p bi_cat.out re-bpy-co3-cl-1e.out -hs ...
```

or

```bash
get_results -r bi_rad.out re-bpy-co3-cl.out \
            -p bi_cat.out re-bpy-co3-cl-1e.out -shs 95 ...
```

These modes require additional solvent and molecular size information:

- `-r_d` and `-r_a`, or alternatively `-v_d` and `-v_a`
- `-d` and `-d_opt` for `-hs`

We do not give a single recommended hard-sphere result in this README because those values depend directly on the radius/volume and dielectric choices supplied by the user. The option is nonetheless available in the code and can be useful as a faster or more approximate screening route.

### Step 7. Optional: Savéant correction for dissociative events

If the transfer event is concerted with bond cleavage, you can add a bond dissociation free energy correction with `-bde`.

With the hard-sphere model:

```bash
get_results -r reactant1.out reactant2.out \
            -p product1.out product2.out product3.out\
            -hs ... -bde 25.0
```

With the simplified hard-sphere model:

```bash
get_results -r reactant1.out reactant2.out \
            -p product1.out product2.out product3.out\
            -shs 95 ... -bde 25.0
```

Here:

- `25.0` is an example `BDFE` value in kcal/mol
- the final reorganization energy used in the barrier expression becomes `λ = λS + BDFE`

This option is intended for dissociative SET/EnT situations where a pure solvent-only reorganization estimate would be unrealistically small.

## What `get_results` extracts from Gaussian outputs

For the current code:

- reaction free energies are taken from lines like:

```text
Sum of electronic and thermal Free Energies=
```

- electronic energies are taken from lines like:

```text
SCF Done:
```

For some nonequilibrium calculations, the code will also use:

```text
After PCM corrections, the energy is
```

if that line is present.

## Notes and limitations

- The package is currently designed around Gaussian-style output parsing.
- Four-point analysis expects the `*_eq.out` and `*_noneq.out` files to follow the naming convention generated by `fourpoints`.
- The example in this repository is a SET from BI(·) to a Re-based catalyst; it is intended as both a regression test and a tutorial.

## Authors

Lucía Morán-González
Albert Solé-Daura
Feliu Maseras

## Citation

If you use `easiq-marcus` in academic work, please cite the associated software paper once available, and cite Gaussian as appropriate for the underlying quantum-chemical calculations.
