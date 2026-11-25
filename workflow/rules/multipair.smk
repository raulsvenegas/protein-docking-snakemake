# FHNW Exercise 4 — optional, end-to-end, configuration-driven docking of
# independent protein pairs. This branch is separate from the verified HPr
# NMR ensemble, so existing results/e2a_hpr/ensemble outputs stay unchanged.
#
# To process one pair without affecting the others:
#   snakemake --cores 1 --use-conda \
#     results/barnase_barstar/automated/best_complex_automated.png \
#     --rerun-triggers mtime
#
# To process every listed pair:
#   snakemake --cores 1 --use-conda all_pairs --rerun-triggers mtime

MULTIPAIR_NAMES = list(config.get("multi_pair_targets", []))

# Explicit final targets make the optional batch easy to invoke and show in DAG.
rule all_pairs:
    input:
        expand(
            f"{RESULTS_DIR}/{{pair}}/automated/docking_scores.csv",
            pair=MULTIPAIR_NAMES,
        ),
        expand(
            f"{RESULTS_DIR}/{{pair}}/automated/best_complex_chainsAB.pdb",
            pair=MULTIPAIR_NAMES,
        ),
        expand(
            f"{RESULTS_DIR}/{{pair}}/automated/best_complex_automated.png",
            pair=MULTIPAIR_NAMES,
        )


# Download a missing RCSB PDB file from its four-character accession.
# Files that are already present are reused without a download.
rule download_pair_pdb:
    output:
        "data/raw/{pair}/{pdb_id}.pdb"
    wildcard_constraints:
        pdb_id="[0-9][A-Za-z0-9]{3}"
    log:
        "logs/multipair/download_{pair}_{pdb_id}.log"
    benchmark:
        "benchmarks/multipair/download_{pair}_{pdb_id}.tsv"
    shell:
        """
        mkdir -p logs/multipair benchmarks/multipair
        python workflow/scripts/download_rcsb_pdb.py \
            --pdb-id {wildcards.pdb_id:q} --output {output:q} > {log:q} 2>&1
        """


# Extract one configured protein chain from an experimental PDB and clean it.
# The PDBFixer repair and protonation options are opt-in through config.yaml.
rule prepare_pair_protein:
    input:
        lambda wc: config["protein_pairs"][wc.pair]["preparation"][wc.role]["source"]
    output:
        "data/prepared/{pair}/{role}.pdb"
    params:
        chain=lambda wc: config["protein_pairs"][wc.pair]["preparation"][wc.role]["chain"],
        model=lambda wc: int(config["protein_pairs"][wc.pair]["preparation"][wc.role].get("model", 1)),
        repair=lambda wc: int(config["protein_pairs"][wc.pair]["preparation"][wc.role].get("repair_missing_atoms", False)),
        hydrogens=lambda wc: int(config["protein_pairs"][wc.pair]["preparation"][wc.role].get("add_hydrogens", False))
    wildcard_constraints:
        role="receptor|ligand"
    conda:
        "../envs/pdbfixer.yaml"
    log:
        "logs/multipair/prepare_{pair}_{role}.log"
    benchmark:
        "benchmarks/multipair/prepare_{pair}_{role}.tsv"
    shell:
        """
        mkdir -p logs/multipair benchmarks/multipair
        python workflow/scripts/prepare_pair_protein.py \
            --source {input:q} --output {output:q} \
            --chain {params.chain:q} --model {params.model} \
            --repair {params.repair} --hydrogens {params.hydrogens} \
            > {log:q} 2>&1
        test -s {output:q}
        """


# Keep the existing configuration-driven docking_setup / docking_run rules
# from Snakefile. Their outputs are results/{pair}/automated/swarm_*/gso_*.out.


rule score_pair:
    input:
        DOCK_GSO_PATTERNS
    output:
        f"{RESULTS_DIR}/{{pair}}/automated/docking_scores.csv"
    params:
        basedir=lambda wc: f"{RESULTS_DIR}/{wc.pair}/automated",
        steps=STEPS
    log:
        "logs/multipair/score_{pair}.log"
    benchmark:
        "benchmarks/multipair/score_{pair}.tsv"
    shell:
        """
        mkdir -p logs/multipair benchmarks/multipair
        python workflow/scripts/score_pair.py \
            --base {params.basedir:q} --steps {params.steps} \
            --output {output:q} > {log:q} 2>&1
        test -s {output:q}
        """


rule build_pair_complex:
    input:
        scores=f"{RESULTS_DIR}/{{pair}}/automated/docking_scores.csv"
    output:
        f"{RESULTS_DIR}/{{pair}}/automated/best_complex_chainsAB.pdb"
    params:
        basedir=lambda wc: f"{RESULTS_DIR}/{wc.pair}/automated",
        steps=STEPS,
        glowworms=int(config["docking"]["glowworms"])
    conda:
        "../envs/lightdock_snapshot.yaml"
    log:
        "logs/multipair/build_{pair}.log"
    benchmark:
        "benchmarks/multipair/build_{pair}.tsv"
    shell:
        """
        mkdir -p logs/multipair benchmarks/multipair
        python workflow/scripts/build_pair_complex.py \
            --scores {input.scores:q} --base {params.basedir:q} \
            --steps {params.steps} --glowworms {params.glowworms} \
            --output {output:q} > {log:q} 2>&1
        test -s {output:q}
        """


rule render_pair_complex:
    input:
        f"{RESULTS_DIR}/{{pair}}/automated/best_complex_chainsAB.pdb"
    output:
        f"{RESULTS_DIR}/{{pair}}/automated/best_complex_automated.png"
    conda:
        "../envs/pymol_snapshot.yaml"
    log:
        "logs/multipair/render_{pair}.log"
    benchmark:
        "benchmarks/multipair/render_{pair}.tsv"
    shell:
        """
        mkdir -p logs/multipair benchmarks/multipair
        python workflow/scripts/render_best_complex.py \
            {input:q} {output:q} > {log:q} 2>&1
        test -s {output:q}
        """
