# CompareCBS
Structural comparison of binding sites based on atom types and distances

## Prerequisites

To use CompareCBS, you need to have R installed along with the following packages:
- `igraph`
- `bio3d`
- `cowplot`
- `fastcluster`
- `openxlsx`
- `vanddraabe`

### Recommended Installation using Mamba

To set up the required environment, follow these steps:

1. Create a new environment and install the necessary R packages:
   ```bash
   mamba create -n legacy_pepit r-base r-essentials r-bio3d r-igraph r-cowplot r-fastcluster r-openxlsx
   mamba activate legacy_pepit
   ```

2. Install the `vanddraabe` package manually, as it has been removed from the CRAN repository:
   ```bash
   wget https://cran.r-project.org/src/contrib/Archive/vanddraabe/vanddraabe_1.0.0.tar.gz
   R
   ```

3. In the R console, install the downloaded package:
   ```R
   install.packages('vanddraabe_1.0.0.tar.gz')
   ```

### Pepit Installation

To install the legacy `pepit` package, open R in the folder containing the package file and run:
```R
install.packages('pepit_1.0.tar.gz')
```

## File Format

Binding site files should be formatted as classical PDB files with one modification: atom typing should be specified in column 77 (1-based index).

### DIONYSUS Atom Typing

Protein:
- Cα: A
- Aromatic carbon: a
- Cβ: b
- Sidechain oxygen: o
- Backbone carbon: C
- Backbone oxygen: O
- Backbone nitrogen: N
- Sidechain carbon (non-aromatic): c

Carbohydrate:
- C1: s
- Other carbons: g
- Oxygens: u

Different chains can be selected during the comparison process. It is advisable to rename chains according to the precise selections you intend to compare.

## Binding Site Comparison

To compare binding sites, use the following command:
```bash
Rscript sugar.R [target_pdb] [chains] [query_pdb] [prefix]
```

### Example
```bash
Rscript sugar.R data/5EYX_1_MAN_1.pdb A,B data/1DGL_1_MAN_3.pdb examples/jacalin_vs_l-type
```

### Output Files

The comparison script outputs four types of files:
- `.al`: Alignment file containing the correspondence between each mapped atom.
- `.index`: Summary of the `.score` file.
- `.score`: Contains alignment metrics, notably the coverage (number of atoms mapped), alignment mean distortion (mean dist), and its score (coverage weighted by atomic pair distortion).
- `.pdb`: Contains the aligned query to the target.

By default, if the same prefix is used multiple times, `sugar.R` will append the `.al`, `.index`, and `.score` files and create new `.pdb` files.

### References
[1] Rasolohery, Inès, Gautier Moroy, and Frédéric Guyon. "PatchSearch: a fast computational method for off-target detection." J Chem. Info. Model. 2017.
[2] Rey, Julien, Inès Rasolohery, Pierre Tufféry, Frédéric Guyon, and Gautier Moroy. "PatchSearch: a web server for off-target protein identification." Nucleic Acids Res, 2019.
[3] Guyon, Frédéric, and Gautier Moroy. "Non-sequential alignment of binding sites for fast peptide screening." bioRxiv, 2023.
[4] More recent version of pepit https://github.com/DSIMB/pepit
