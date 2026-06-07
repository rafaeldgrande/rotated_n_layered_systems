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
