# Rotated N-Layered Systems

Scripts to generate and analyze twisted N-layer van der Waals systems (graphene and TMDs) for LAMMPS molecular dynamics simulations. These scripts were used in the following papers:

- [Flat bands and gaps in twisted double bilayer graphene](https://doi.org/10.1039/C9NR10830K)
- [Flat bands and gaps in twisted double trilayer graphene](https://doi.org/10.1103/PhysRevB.111.075111)
- [Multi scale calculation of light-induced interlayer changes in low-angle twisted bilayer WSe₂](https://doi.org/10.48550/arXiv.2604.28143)



## Files Description

| File | Purpose |
|------|---------|
| `rotated_n_layer.py` | Generate twisted double-trilayer graphene (t3+3LG) supercells for LAMMPS |
| `generate_twisted_2L_WSe2.py` | Generate twisted bilayer WSe₂ using integer-based moiré construction (recommended) |
| `generate_twisted_2L_WSe2_new.py` | Newer brute-force version for twisted bilayer WSe₂ |
| `generate_twisted_2L_MoS2.py` | Generate twisted bilayer MoS₂ supercells (SW+KC relaxed parameters, AA registry at origin) |
| `generate_twisted_MoSe2WSe2.py` | Generate twisted MoSe₂/WSe₂ heterobilayer supercells (averaged lattice constant, AB_3R stacking) |
| `build_heterobilayer.py` | General heterobilayer builder: brute-force search for commensurate (N, M) pairs satisfying N·a₁ = M·a₂ across two different lattice constants |
| `plot_cos_theta_L_twisted_2L.py` | Find (N, M) pairs for a target angle and plot supercell size vs. twist angle |
| `energia_por_atomo_t6LG.py` | Post-process LAMMPS results: plot energy per atom vs. twist angle for t6LG |
| `run_code.bash` | Batch script to generate t6LG structures for multiple (N, M) pairs |


## Theory

For a hexagonal lattice, commensurate twist angles are defined by integer pairs (N, M). The supercell vectors are:

```
A1 = N·a1 + M·a2
A2 = -M·a1 + (N+M)·a2
```

The twist angle θ and supercell lattice parameter L are:

```
θ = arccos((N² + 4NM + M²) / (2(N² + NM + M²)))

L = a · √(N² + NM + M²)
```

The supercell contains exactly `N² + NM + M²` primitive cells. Each layer is rotated by ±θ/2 about the z-axis. The Moiré pattern periodicity is `L / |M - N|`.


## Usage

### Twisted bilayer WSe₂ (`generate_twisted_2L_WSe2.py`)

Takes the lattice parameter and (N, M) indices as command-line arguments:

```bash
python generate_twisted_2L_WSe2.py <a_angstrom> <N> <M>
```

Example (magic angle ~1.1°, a = 3.29 Å):

```bash
python generate_twisted_2L_WSe2.py 3.29 34 33
```

Outputs:
- `N_M_twisted.top` — LAMMPS topology file (triclinic box, atom types: W-layer1=1, Se-layer1=2, W-layer2=3, Se-layer2=4)
- `N_M_twisted.xyz` — XYZ file for visualization
- `1st_layer_twisted.png` — diagnostic plot of layer 1 atomic positions

The script uses an exact integer-based enumeration of primitive-cell origins inside the moiré supercell, guaranteeing exactly `3(N² + NM + M²)` atoms per layer (3 atoms per primitive cell: W, Se, Se).

### Twisted 6-layer graphene (`rotated_n_layer.py`)

Edit `N` and `M` directly in the script (default: N=61, M=59) and run:

```bash
python rotated_n_layer.py
```

Outputs:
- `t6LG_N_M.xyz` — XYZ positions for the supercell
- `visual_t6LG_N_M.xyz` — XYZ for visualization (includes atoms outside supercell)
- `N_M/t6LG_N_M.top` — LAMMPS topology file (6 atom types, ABA stacking per trilayer)
- `cell_basis` — supercell lattice vectors

### Batch generation (`run_code.bash`)

```bash
bash run_code.bash
```

Loops over a predefined list of (N, M) pairs and calls the structure generator for each.

### Finding (N, M) for a target angle (`plot_cos_theta_L_twisted_2L.py`)

Edit the target angle in the script (default: 1.1°) and run:

```bash
python plot_cos_theta_L_twisted_2L.py
```

Prints all (N, M) pairs within 0.05° of the target and plots L vs. θ for the full (N, M) space explored.

### Twisted bilayer MoS₂ (`generate_twisted_2L_MoS2.py`)

Same interface as the WSe₂ generator:

```bash
python generate_twisted_2L_MoS2.py <N> <M>
```

Uses SW+KC relaxed lattice parameters for AB_3R stacking (a = 3.117 Å, d_Mo-Mo = 6.087 Å). Outputs:
- `N_M_twisted.top` — LAMMPS topology file (atom types: Mo-L1=1, S-L1=2,3; Mo-L2=4, S-L2=5,6)
- `t2LMoS2_N_M.xyz` — XYZ file for visualization

### Twisted MoSe₂/WSe₂ heterobilayer (`generate_twisted_MoSe2WSe2.py`)

```bash
python generate_twisted_MoSe2WSe2.py <N> <M>
```

Uses AB_3R bilayer relaxed parameters with the averaged lattice constant (a = 3.2977 Å). Outputs to `N_M/` subdirectory:
- `N_M/twisted_MoSe2WSe2_initial.top` — LAMMPS topology file (6 atom types: Mo1, Se2, Se1 for MoSe₂ layer; W2, Se4, Se3 for WSe₂ layer)
- `N_M/t2L_MoSe2WSe2_N_M.xyz` — XYZ file for visualization

### General heterobilayer builder (`build_heterobilayer.py`)

For systems with different lattice constants in each layer (e.g. MoS₂/WSe₂, Δa ≈ 5.35%). Uses a brute-force search to find integers (N, M) such that N·a₁ ≈ M·a₂ within a tolerance, producing a commensurate supercell.

```bash
python build_heterobilayer.py [--a1 A1] [--a2 A2] [--theta THETA] \
    [--h1 H1] [--h2 H2] [--d_int D_INT] [--tol TOL] [--N_max N_MAX] \
    [--output FILENAME] [--m1_metal M] [--m1_chalc M] [--m2_metal M] [--m2_chalc M]
```

Key arguments:

| Argument | Default | Description |
|----------|---------|-------------|
| `--a1` | 3.31 | Lattice constant of layer 1 (Å) |
| `--a2` | 3.11 | Lattice constant of layer 2 (Å) |
| `--theta` | 1.0 | Twist angle (degrees) |
| `--h1` | 1.6 | Chalcogen–metal height in layer 1 (Å) |
| `--h2` | 1.5 | Chalcogen–metal height in layer 2 (Å) |
| `--d_int` | 7.5 | Metal–metal interlayer distance (Å) |
| `--tol` | 0.05 | Commensurability tolerance (Å) |
| `--N_max` | 100 | Maximum integer search range |
| `--output` | auto | Output `.dat` filename |

Example (MoS₂/WSe₂, θ = 1°):

```bash
python build_heterobilayer.py --a1 3.288 --a2 3.117 --theta 1.0 \
    --h1 1.644 --h2 1.559 --d_int 6.376 \
    --m1_metal 183.84 --m1_chalc 78.97 --m2_metal 95.96 --m2_chalc 32.06
```

Outputs a LAMMPS `.dat` file (atom types: W=1, Se=2,3 for layer 1; Mo=4, S=5,6 for layer 2).

### Post-processing energy (`energia_por_atomo_t6LG.py`)

Reads LAMMPS log files from subdirectories named `N_M/` and plots energy per atom vs. twist angle:

```bash
python energia_por_atomo_t6LG.py
```

Requires completed LAMMPS runs with log files in each `N_M/` directory.


## Dependencies

- Python 3
- NumPy
- Matplotlib


## Citation

```bibtex
@article{Culchac2020,
    title = {Flat bands and gaps in twisted double bilayer graphene},
    volume = {12},
    ISSN = {2040-3372},
    url = {http://dx.doi.org/10.1039/C9NR10830K},
    DOI = {10.1039/c9nr10830k},
    number = {8},
    journal = {Nanoscale},
    publisher = {Royal Society of Chemistry (RSC)},
    author = {Culchac, F. J. and Del Grande, R. R. and Capaz, Rodrigo B. and Chico, Leonor and Morell, E. Suárez},
    year = {2020},
    pages = {5014–5020}
}
```

```bibtex
@article{Culchac2025,
    title = {Flat bands and gaps in twisted double trilayer graphene},
    author = {Culchac, F. J. and Del Grande, R. R. and Menezes, Marcos G. and Capaz, Rodrigo B.},
    journal = {Phys. Rev. B},
    volume = {111},
    issue = {7},
    pages = {075111},
    numpages = {8},
    year = {2025},
    month = {Feb},
    publisher = {American Physical Society},
    doi = {10.1103/PhysRevB.111.075111},
    url = {https://link.aps.org/doi/10.1103/PhysRevB.111.075111}
}
```

```bibtex
@misc{delgrande2026multiscalecalculationlightinducedstructural,
      title={Multi-scale calculation of light-induced structural changes in low-angle twisted bilayer WSe$_2$}, 
      author={Rafael R. Del Grande and David A. Strubbe},
      year={2026},
      eprint={2604.28143},
      archivePrefix={arXiv},
      primaryClass={cond-mat.mtrl-sci},
      url={https://arxiv.org/abs/2604.28143}, 
}
```
