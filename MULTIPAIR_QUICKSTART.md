# Multi-protein-pair workflow — quick start (FHNW Exercise 4)

This is an **additional branch** of the original EIIA-Glc/HPr Snakemake project. The already verified HPr NMR ensemble and its 600 poses are **not modified**.

## What has been added

A complete `{pair}`-wildcard pipeline:

```text
RCSB PDB download (only if missing)
     | 
     +--> select/clean receptor chain -> optional PDBFixer heavy-atom repair --+
     +--> select/clean ligand chain   -> optional PDBFixer heavy-atom repair --+ 
                                                                              |
                                                            LightDock setup -> run
                                                                              |
                                                         score + rank candidate poses
                                                                              |
                                                         reconstruct best chains A/B
                                                                              |
                                                         automated PyMOL PNG
```

The existing rules `docking_setup` and `docking_run` are reused, rather than copied. The extra rules are in `workflow/rules/multipair.smk`. Their logs and benchmarks are under `logs/multipair/` and `benchmarks/multipair/`.

## Configured protein pairs and verification status

- **`two_uuy` — COMPLETED:** The student executed a real Snakemake run using `data/raw/2UUY/2UUY_rec.pdb` (chain A) and `data/raw/2UUY/2UUY_lig.pdb` (chain B) and verified **7/7 completed jobs**, **20 ranked poses**, **top score +4.83015167 (swarm 0, pose index 7)**, a reconstructed complex with chains A/B and an automated PyMOL figure. These are workflow outputs, not independent structural validation.

- `e2a_hpr`: EIIA-Glc/HPr (the existing input proteins). The original complete **30-model ensemble** is the main, tested workflow under `results/e2a_hpr/ensemble`. The optional single-pair pipeline under `results/e2a_hpr/automated` is separate.
- `barnase_barstar` (**configured but NOT executed**): barnase and barstar from RCSB [PDB 1BRS](https://www.rcsb.org/structure/1BRS). Receptor barnase is **chain A** and ligand barstar is **chain D**. These are two *different protein species*, unlike the HPr NMR conformations. They are extracted from a **bound experimental complex**, so a new docking test is a *workflow demonstration*, **not independent prediction validation**. Source: Buckle et al. (1994), PDB DOI `10.2210/pdb1BRS/pdb`.

The `1BRS.pdb` file is **not included** in this package. **2UUY has already been executed and its local raw input files and outputs are included.** When first needed, Snakemake downloads it from `https://files.rcsb.org/download/1BRS.pdb` using the new `download_pair_pdb` rule. This requires an internet connection in WSL; if the download fails, obtain this file manually and put it in `data/raw/barnase_barstar/1BRS.pdb`.

## Run from VS Code / WSL Ubuntu

```bash
cd ~/antibody-docking-snakemake
conda activate snakemake

# 1. Safe re-check of the ALREADY COMPLETED second pair, 2UUY:
snakemake --cores 1 --use-conda -n \
  results/two_uuy/automated/best_complex_automated.png \
  --rerun-triggers mtime

# 2. Target ONLY 2UUY (reuses its current results if up to date):
snakemake --cores 1 --use-conda \
  results/two_uuy/automated/best_complex_automated.png \
  --rerun-triggers mtime

# 3. Optional: dry-run ALL configured pairs. This includes the NOT YET RUN
#    barnase_barstar example and may schedule its PDB download / LightDock.
snakemake --cores 1 --use-conda -n all_pairs --rerun-triggers mtime

# 4. Optional: preview the unexecuted barnase–barstar pair without running it:
snakemake --cores 1 --use-conda -n \
  results/barnase_barstar/automated/best_complex_automated.png \
  --rerun-triggers mtime

# 5. Visualize the 2UUY DAG:
snakemake --dag results/two_uuy/automated/best_complex_automated.png \
  | dot -Tpng -o two_uuy_dag.png
```

**DO NOT use `--forceall`, `--forcerun` or `--touch` during routine verification.** The project's original `snakemake all` still targets the already verified 30-model EIIA/HPr ensemble. A dry-run (`-n`) never performs docking. Without `-n`, the new example will perform real computational work and can take substantial computing time depending on the machine.

## Add another protein pair later

Edit `config/config.yaml` and add:

```yaml
protein_pairs:
  my_pair:
    receptor: data/prepared/my_pair/receptor.pdb
    ligand: data/prepared/my_pair/ligand.pdb
    preparation:
      receptor:
        source: data/raw/my_pair/1ABC.pdb
        chain: A
        model: 1
        repair_missing_atoms: false
        add_hydrogens: false
      ligand:
        source: data/raw/my_pair/2XYZ.pdb
        chain: B
        model: 1
        repair_missing_atoms: false
        add_hydrogens: false

multi_pair_targets:
  - my_pair
```

`1ABC` and `2XYZ` above are **illustrative placeholders**, not validated PDB choices. Replace them with actual RCSB accessions, choose the correct chains, and check the structures are protein inputs that LightDock can process. For a PDB file not hosted at RCSB, place the file at the configured path manually. Then dry-run `results/my_pair/automated/best_complex_automated.png`.

The preprocessing rule uses `ATOM` records, a chosen chain/MODEL, removes HETATM and alternate locations other than blank/A, and can optionally call PDBFixer. It does **not** promise that arbitrary PDB structures will work without additional biological cleaning, missing-residue rebuilding, or handling of ligands/cofactors. Hydrogens are optional in prepared PDBs and are removed by the existing `--noh` docking setup.

## Final output files for each pair

```text
results/<pair>/automated/docking_scores.csv
results/<pair>/automated/best_complex_chainsAB.pdb
results/<pair>/automated/best_complex_automated.png
logs/multipair/{download,prepare,score,build,render}_*.log
benchmarks/multipair/*.tsv
```

Scores are from LightDock's *computational scoring function*, not binding energies. The rebuilt PDB assumes LightDock's generated complex uses a receptor–ligand atom-serial reset (as demonstrated by the existing EIIA–HPr output). Other LightDock versions / formats may require adapting this parser.

## What has actually been verified

- The **original** EIIA/HPr NMR ensemble was executed (30 conformers, 60 swarm files and 600 ranked poses); its final PDB and PyMOL PNG are present.
- **2UUY was executed end to end on WSL on 23 November 2026**, using the generic `{pair}` workflow: **7/7 Snakemake jobs**, two swarms, 20 ranked poses, best score **+4.83015167** (swarm 0, pose index 7), best-complex PDB with chains A/B and automated PyMOL PNG. The input structures were supplied as local files.
- Seven local Python unit/regression tests passed; these alone would not prove a real docking run, so the completed Snakemake execution is recorded separately in `MULTIPAIR_VALIDATION_2UUY.md`.
- **Barnase–barstar is configured but has not been run**. The project has not been independently reproduced from scratch on a different machine, and the docking geometries are not biologically validated.
