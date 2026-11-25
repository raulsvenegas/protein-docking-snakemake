# Automated Protein–Protein Docking with Snakemake and LightDock
**Authors:** Raúl Sánchez Venegas and Jovana Arsenov

**FHNW · MSc Life Sciences – Data Science · Computer and Software Architecture · Exercise 4**

## Project goal

This project implements a Snakemake workflow for preparing two interacting proteins, performing exploratory protein–protein docking with LightDock, collecting and ranking the computed poses, reconstructing the highest-scoring complex, and automatically producing a molecular figure with PyMOL.

The objective is **workflow automation and software architecture**, rather than establishing an experimentally validated protein–protein interaction geometry.

## Input structures and preprocessing

| Role | Protein | Source structure | Preparation |
| --- | --- | --- | --- |
| Receptor | EIIA-Glc | RCSB PDB `1F3G` | Keep protein `ATOM` records, then use PDBFixer to repair 23 missing heavy atoms |
| Ligand | HPr | RCSB PDB `1HDN` | Extract all 30 structural models from the NMR ensemble |

Original source files are expected at `data/raw/e2a_hpr/1F3G.pdb` and `data/raw/e2a_hpr/1HDN.pdb`. The prepared receptor contains one protein chain, A. The extracted HPr model 1 has one chain, A; the extracted file contains hydrogens, but LightDock setup uses `--noh` to remove them for docking.

**Preprocessing limitation:** EIIA-Glc's unresolved complete N-terminal residues were deliberately **not** reconstructed; the script calls `findMissingResidues()` but then clears `missingResidues`. The receptor repair adds missing heavy atoms, not hydrogens. Preparation uses the experimentally resolved protein coordinates, so neither a complete rebuilding of missing residues nor a universally applicable protonation workflow has been demonstrated.

## How the workflow works

```text
EIIA-Glc (1F3G) -> clean -> repair ------------------+
                                                     |
HPr NMR ensemble (1HDN) -> extract 30 models -------+-> LightDock setup
                                                              |
                                                              v
                                               LightDock run (2 swarms/model)
                                                              |
                                                              v
                                                   Parse and rank scores
                                                              |
                                                              v
                                                  Reconstruct best PDB
                                                              |
                                                              v
                                                    Render PyMOL PNG
                                                              |
                                                              v
                                                      Final outputs
```

Snakemake builds a directed acyclic graph (DAG) from file dependencies and avoids rerunning existing outputs unless they are missing or considered out of date.

The original 30-model ensemble remains **specific to EIIA-Glc/HPr**, using the `{model}` wildcard. A separate **end-to-end `{pair}`-wildcard workflow** is now available for independent protein pairs configured in `config/config.yaml`. It was **executed successfully in WSL with 2UUY**, including preparation, LightDock, scoring, PDB reconstruction and automatic PyMOL rendering (7 of 7 jobs completed). **The 30 HPr NMR conformations are not 30 different protein pairs.**

## Configuration

The configuration file is `config/config.yaml`. Tested settings:

| Parameter | Value |
| --- | ---: |
| HPr NMR models | 30 (models 1–30) |
| LightDock swarms per model | 2 |
| Glowworms per swarm | 10 |
| Optimization steps | 10 |
| Random seed | 12345 |
| CPU cores | 1 |

The docking rules read the model range (for HPr), number of swarms, glowworms, steps, seed and CPU allocation from the YAML configuration. These are **small exploratory runs**, not exhaustive docking benchmarks. The source structures are referenced by PDB accession, and the 2UUY input chains were supplied as separate local PDB files.

## Installation and execution

Tested on Ubuntu under Windows Subsystem for Linux (WSL2) with Snakemake 9.27.1 and LightDock 1.0.0rc2. Python scripts, shell commands, Conda environment definitions and Snakemake rules are included.

Conda environment definitions:

- `workflow/envs/pdbfixer.yaml` — structure repair
- `workflow/envs/lightdock_snapshot.yaml` — LightDock simulation and best-pose reconstruction
- `workflow/envs/pymol_snapshot.yaml` — automated PyMOL figure rendering

The Conda-managed LightDock and PyMOL rules were executed successfully for the 2UUY demonstration on the student’s WSL machine. The original EIIA/HPr ensemble was previously run and its final PyMOL rendering was also executed through Snakemake. **Reproduction from scratch on a different computer has not been demonstrated.** Some exported environment snapshots may require platform-specific adjustment.

From the repository root, with Snakemake activated:

```bash
# Check dependencies without executing any simulations
snakemake --cores 1 --use-conda -n all --rerun-triggers mtime

# Run the main target when inputs or outputs require work
snakemake --cores 1 --use-conda all

# Run / verify the separately tested 2UUY pair
snakemake --cores 1 --use-conda -n \
  results/two_uuy/automated/best_complex_automated.png \
  --rerun-triggers mtime

# Run only 2UUY when a real execution is needed (not a dry run)
snakemake --cores 1 --use-conda \
  results/two_uuy/automated/best_complex_automated.png \
  --rerun-triggers mtime

# Produce the rule graph for the original EIIA/HPr main target
snakemake --rulegraph dot all | dot -Tpng \
  -o results/e2a_hpr/ensemble/workflow_rulegraph_final.png
```

**Important:** An ordinary Snakemake run may schedule fresh docking if input, code or software-environment provenance has changed. The `--rerun-triggers mtime` option was used during final verification to avoid unnecessary reruns of previously generated results; using this option does not replace an independent reproduction test. To recreate results from scratch, adequate computational resources and functioning environments are required.

## Results

| Quantity | Observed result |
| --- | ---: |
| Extracted HPr structures | 30 |
| LightDock swarm outputs | 60 |
| Ranked candidate poses | 600 |
| Highest-scoring HPr model | 9 |
| Winning swarm | 0 |
| Winning pose index | 5 |
| Highest LightDock score | **-4.13045687** |
| Atoms in final receptor chain A | 1,137 |
| Atoms in final ligand chain B | 639 |

### Second protein-pair test: 2UUY

A **second independent protein pair**, supplied as `data/raw/2UUY/2UUY_rec.pdb` (chain A; 1,628 input protein ATOM records) and `data/raw/2UUY/2UUY_lig.pdb` (chain B; 415 input protein ATOM records), was processed using the **generic `{pair}` workflow**, not the HPr 30-model ensemble.

| Quantity | Observed 2UUY result |
| --- | ---: |
| Snakemake jobs completed | **7 / 7** |
| LightDock swarms | 2 |
| Candidate poses ranked | **20** |
| Top-ranked swarm | 0 |
| Top-ranked pose index | 7 |
| Highest LightDock score | **+4.83015167** |
| Final reconstructed complex | chains A and B |
| Automatic PyMOL PNG | produced and verified |

The full 2UUY run completed successfully on WSL on **23 November 2026**. The run demonstrates end-to-end **software integration and reuse** with a second pair, not validation of its biological binding geometry. Its 20 candidate poses should not be compared as an equivalent sampling experiment with the 600-pose HPr ensemble.

**Scientific interpretation:** Scores are **LightDock computational scores**, not measured binding affinities; +4.8302 for 2UUY must **not** be interpreted as stronger binding than −4.1305 for EIIA/HPr. An exploratory EIIA/HPr structural inspection identified **54 inter-protein heavy-atom pairs below 2 Å**, suggesting severe steric overlaps. The 2UUY predicted interface also requires independent structural validation. Neither selected pose is an experimentally validated biological prediction.

## Outputs and evidence

- `results/e2a_hpr/ensemble/docking_scores.csv` — all 600 poses ranked by docking score
- `results/e2a_hpr/ensemble/best_scores_by_model.csv` — best candidate from each NMR model
- `results/e2a_hpr/ensemble/model_*/swarm_*/gso_10.out` — 60 docking swarm output files
- `results/e2a_hpr/ensemble/best_complex_chainsAB.pdb` — best-ranked complex, chains A and B
- `results/e2a_hpr/ensemble/best_complex_automated.png` — **automatically rendered by PyMOL through Snakemake**
- `results/e2a_hpr/ensemble/workflow_rulegraph_final.png` — current simplified rule graph
- `results/e2a_hpr/ensemble/model_9/swarm_0/complex_interface.png` — earlier manual interface visualization (when available)
- `results/two_uuy/automated/docking_scores.csv` — **20** 2UUY candidate scores
- `results/two_uuy/automated/best_complex_chainsAB.pdb` — 2UUY best-ranked complex
- `results/two_uuy/automated/best_complex_automated.png` — **automatically generated** 2UUY molecular figure
- `MULTIPAIR_VALIDATION_2UUY.md` — concise record of the completed second-pair run
- `logs/` — command output and error logs
- `benchmarks/` — Snakemake runtime/resource benchmarks where recorded

Snakemake successfully ran `render_best_complex` in its managed PyMOL Conda environment, producing its output, log and benchmark, and then completed the `all` target. The original docking runs predate the docking rules' current benchmark directives, so **do not infer 30 measured setup/run runtimes from those directives**.

## Key source files

| File | Purpose |
| --- | --- |
| `Snakefile` | Workflow entry point; EIIA preparation and optional `{pair}` docking path |
| `config/config.yaml` | Protein pair definitions, 30-model range and docking resources |
| `workflow/rules/hpr_batch.smk` | `{model}` wildcard, extraction, 30 LightDock setups and runs |
| `workflow/rules/score_analysis.smk` | Aggregate scoring rule |
| `workflow/rules/best_complex.smk` | Best-scoring complex reconstruction rule |
| `workflow/rules/render_complex.smk` | Conda-managed PyMOL figure rule |
| `workflow/rules/final_project.smk` | Named aggregate target (`final_project`) |
| `workflow/rules/multipair.smk` | Generic `{pair}` rules for preparation, score parsing, reconstruction and PyMOL |
| `workflow/scripts/prepare_pair_protein.py` | Select and clean configured PDB chain; optional repair/protonation |
| `workflow/scripts/score_pair.py` | Rank candidates for one configured protein pair |
| `workflow/scripts/build_pair_complex.py` | Rebuild a generic pair’s highest-scoring complex |
| `tests/test_multipair.py` | Seven Python unit and regression tests |
| `workflow/scripts/collect_docking_scores.py` | Parse and rank LightDock scores |
| `workflow/scripts/build_best_complex.py` | Generate/convert the best pose to a two-chain PDB |
| `workflow/scripts/render_best_complex.py` | Render the selected complex in PyMOL without the GUI |
| `workflow/scripts/repair_eiia.py` | Repair missing receptor heavy atoms using PDBFixer |

## What was verified, and what remains limited

**Verified:** 30 HPr NMR models processed, 60 swarm files, 600 ranked poses, best PDB and automated PyMOL PNG; a **second real end-to-end execution with 2UUY** using the generic `{pair}` rules (7/7 Snakemake jobs, 20 ranked poses, final PDB chains A/B and PNG); seven unit/regression tests passed; Conda-based LightDock/PyMOL execution on WSL; results and available logs/benchmarks.

**Not established:** reproducibility from scratch on another machine; guaranteed support for **arbitrary** PDB pairs without structure-specific checks; success of the configured but **unexecuted barnase–barstar example**; complete reconstruction of unresolved residue stretches; comprehensive protonation; exhaustive conformational sampling; or scientific validity of either predicted protein interface.

## Configuration-driven multi-pair workflow

The optional `all_pairs` target uses the `multi_pair_targets` list from `config/config.yaml` and produces scores, a two-chain best-complex PDB and a PyMOL PNG for each named protein pair. It does not replace the original 30-model HPr ensemble under `results/e2a_hpr/ensemble`.

Configured pairs at submission:

- **`e2a_hpr`** — the original EIIA/HPr proteins; the main **30-conformer ensemble** is validated separately. The optional `results/e2a_hpr/automated/` single-pair route is not the same 30-model experiment.
- **`two_uuy`** — **executed successfully end to end** using the raw local `2UUY_rec.pdb` and `2UUY_lig.pdb` inputs. The seven-job run and its CSV/PDB/PNG outputs are included.
- **`barnase_barstar`** — optional example from RCSB **1BRS** chains A and D, **configured but not executed**. Its download and docking would require additional runtime work and input checks.

```bash
# Inspect the DAG for 2UUY without recomputing anything
snakemake --cores 1 --use-conda -n \
  results/two_uuy/automated/best_complex_automated.png \
  --rerun-triggers mtime

# Preview all configured pairs (this includes the unexecuted barnase–barstar pair)
snakemake --cores 1 --use-conda -n all_pairs --rerun-triggers mtime
```

**Caution:** Running `snakemake all_pairs` without `-n` will also attempt the **barnase–barstar** experiment and any missing single-pair outputs for `e2a_hpr`. It does **not** just verify the two completed cases. For a safe check of the already executed second example, target `results/two_uuy/automated/best_complex_automated.png` directly. Running `all` by itself still targets the original HPr 30-model ensemble.

The generic preprocessing can select a configured MODEL and chain, retain protein `ATOM` records and optionally call PDBFixer. However, chain assignment, missing residues, cofactors and interface validity require biological judgment. An experimentally bound complex used as docking input is not independent validation of a predicted interface. See [`MULTIPAIR_QUICKSTART.md`](MULTIPAIR_QUICKSTART.md) for setup and [`MULTIPAIR_VALIDATION_2UUY.md`](MULTIPAIR_VALIDATION_2UUY.md) for the second-pair result.

## Recommendations for future students

1. Inspect the contents, chains, missing atoms and missing residues of input PDB structures before simulation.
2. Start with a small pilot, check the DAG and dry run, and distinguish workflow correctness from biological correctness.
3. Define clear filenames, wildcards and parameters in YAML early; do not mistake ensemble conformations for independent protein pairs.
4. Test environment creation and keep logging/benchmarking active from the first run.
5. Evaluate the molecular geometry of the selected pose; the highest numerical score is not proof of a correct complex.
