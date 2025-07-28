# Rotated N-Layered Systems

Set of scripts to study twisted N-layer systems (graphene and TMDs). Those scripts were used in the following papers (see bibtex entries in the end of this page):
- [Flat bands and gaps in twisted double bilayer graphene](https://doi.org/10.1039/C9NR10830K)
- [Flat bands and gaps in twisted double trilayer graphene](https://doi.org/10.1103/PhysRevB.111.075111)



## Files Description

| File | Purpose |
|------|---------|
| `rotated_n_layer.py` | Main script for generating twisted 6-layer graphene structures |
| `generate_twisted_2L_WSe2.py` | Generate twisted bilayer WSe₂ systems |
| `energia_por_atomo_t6LG.py` | Calculate and analyze energy per atom for different configurations |
| `plot_cos_theta_L_twisted_2L.py` | Plot relationships between twist angle and supercell size |
| `run_code.bash` | Batch processing script for multiple configurations |



## Theory

The twist angle θ between layers is calculated using:

```
θ = arccos((N² + 4NM + M²)/(2(N² + NM + M²)))
```

The supercell size is given by:

```
L = a√(3(N² + NM + M²))
```

Where:
- `N`, `M` are integer indices defining the supercell
- `a` is the lattice parameter
- The Moiré pattern periodicity is `L/|M-N|`




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
    author = {Culchac,  F. J. and Del Grande,  R. R. and Capaz,  Rodrigo B. and Chico,  Leonor and Morell,  E. Suárez},
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


