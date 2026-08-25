# Household resistome permutation tests

Permutation (randomization) tests to assess whether antibiotic resistance gene
profiles (resistomes) are more similar **within households** than expected by
chance, across four sample compartments: **Patient, Animal, Environment,
Food**.

## Input data

A csv file (semicolon-separated, read with `read.csv2`) with the following
columns, in this order:

| sample_ID | patient_ID | sample_type | gene1 | gene2 | ... | geneX |
|-----------|------------|-------------|-------|-------|-----|-------|

- **sample_ID**: unique identifier of the sample.
- **patient_ID**: household identifier (all samples collected from the same
  household — patient, animals, environment, food — share the same
  `patient_ID`).
- **sample_type**: one of `"Patient"`, `"Animal"`, `"Environment"`, `"Food"`.
- **gene1...geneX**: presence (1) / absence (0) of each antibiotic resistance
  gene.

## Method

For each script:

1. Compute pairwise **Jaccard distances** between all samples' gene profiles
   (`vegan::vegdist`).
2. Keep only intra-household comparisons (`patient_ID1 == patient_ID2`),
   removing self-comparisons and duplicate pairs (A–B kept, B–A dropped).
3. For each household, average the distance within each of the 10 pairwise
   comparison categories: `HuHu, HuAni, HuEnv, HuFood, AniAni, AniEnv,
   AniFood, EnvEnv, EnvFood, FoodFood` (Hu = Patient/Human).
4. Average across households to get one **observed value (Mtot)** per
   category.
5. Build a **null distribution** by repeating steps 1–4 on 10,000 randomized
   versions of the dataset (see randomization strategies below).
6. Compute an **empirical p-value** per category: the proportion of
   randomized Mtot values ≤ the observed Mtot.
7. Export diagnostic plots: null distribution histograms (observed value in
   blue, 5%/95% thresholds in red, 2.5%/97.5% in orange) and p-value
   stabilization curves across the 10,000 iterations (to check that
   `n_iter` is large enough for the p-value to have converged).

## The 8 scripts

The 8 scripts differ only in **which compartment is randomized** and **how**
it is randomized. Within each group of 4, scripts are identical except for
the randomization block (~lines 100–135).

### Group 1 — Individual-level randomization

Shuffles gene profiles **between samples of the same compartment**, across
all households, while keeping the number of samples per household per
compartment unchanged. This breaks any household-specific signal in the
targeted compartment while preserving household composition.

- `permutation_test_randomize_patient.R`
- `permutation_test_randomize_animal.R`
- `permutation_test_randomize_environment.R`
- `permutation_test_randomize_food.R`

### Group 2 — Block-level randomization

Reassigns the **entire compartment block of a household to a different,
randomly drawn household** (e.g., all animal samples from household A are
moved together to household B), instead of shuffling individual gene
profiles. This tests whether households as a whole (rather than individual
gene profiles) drive the observed similarity.

- `permutation_test_block_randomize_patient.R`
- `permutation_test_block_randomize_animal.R`
- `permutation_test_block_randomize_environment.R`
- `permutation_test_block_randomize_food.R`

> **Note:** unlike individual-level randomization, block-level randomization
> does not guarantee that each household keeps the same number of samples
> per compartment after randomization, if household sizes are unequal in the
> original dataset (a household inherits another household's block "as is").

## Requirements

```r
install.packages(c("vegan", "dplyr", "ggplot2", "tidyr", "viridis"))
```

## Usage

1. Place your data file in the working directory and update the path in:
   ```r
   table <- as.data.frame(read.csv2("data_file"))
   ```
2. Run the script for the compartment/strategy you want to test (e.g. in R
   or via `Rscript permutation_test_randomize_animal.R`).
3. Each run produces (in the working directory):
   - `Mtot_jaccard.csv` — observed mean intra-household distance per
     comparison category.
   - `table_pvalues.csv` — empirical p-value per comparison category.
   - `jac_<Category>.tiff` (×10) — null distribution histogram per category.
   - `p-values_evolution_jac.tiff` and `p-values_evolution_jac_faceting.tiff`
     — p-value stabilization plots across permutation iterations.

Since all 8 scripts write to the same output filenames, run them in separate
working directories (or rename outputs) if you want to keep results from
several scripts side by side.

## Runtime

Each iteration recomputes a full pairwise distance matrix and re-aggregates
it by household, so runtime scales with `n_iter × (number of samples)²`.
With `n_iter = 10,000`, this can take a while on large datasets — consider
lowering `n_iter` for quick tests before a full run.
