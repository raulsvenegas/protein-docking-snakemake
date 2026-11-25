
# FHNW Exercise 4 - Protein-Protein Docking Workflow
# Automated workflow using Snakemake and LightDock

configfile: "config/config.yaml"


RESULTS_DIR = config.get("results_dir", "results")

SWARMS = int(config["docking"]["swarms"])
STEPS = int(config["docking"]["steps"])

# Output patterns used by the docking simulation rule
DOCK_GSO_PATTERNS = [
    f"{RESULTS_DIR}/{{pair}}/automated/swarm_{swarm}/gso_{STEPS}.out"
    for swarm in range(SWARMS)
]


# --------------------------------------------------
# RULE 1 - Final workflow targets
# --------------------------------------------------

rule all:
    input:
        "results/e2a_hpr/ensemble/docking_scores.csv",
        "results/e2a_hpr/ensemble/best_scores_by_model.csv",
        "results/e2a_hpr/ensemble/best_complex_chainsAB.pdb",
        "results/e2a_hpr/ensemble/best_complex_automated.png"


# --------------------------------------------------
# RULE 4 - Extract HPr NMR model
# --------------------------------------------------

rule extract_hpr_model:
    input:
        "data/raw/e2a_hpr/1HDN.pdb"
    output:
        "data/prepared/e2a_hpr/hpr_model1.pdb"
    params:
        model=1
    script:
        "workflow/scripts/extract_nmr_model.py"


# --------------------------------------------------
# RULE 5 - Clean EIIA-Glc
# --------------------------------------------------

rule prepare_eiia:
    input:
        "data/raw/e2a_hpr/1F3G.pdb"
    output:
        "data/prepared/e2a_hpr/eiia_clean.pdb"
    shell:
        """
        mkdir -p data/prepared/e2a_hpr
        grep '^ATOM  ' {input} > {output}
        echo END >> {output}
        test $(grep -c '^ATOM  ' {output}) -gt 0
        """


# --------------------------------------------------
# RULE 6 - Repair missing EIIA-Glc atoms
# --------------------------------------------------

rule repair_eiia:
    input:
        "data/prepared/e2a_hpr/eiia_clean.pdb"
    output:
        "data/prepared/e2a_hpr/eiia_repaired.pdb"
    conda:
        "workflow/envs/pdbfixer.yaml"
    script:
        "workflow/scripts/repair_eiia.py"


# --------------------------------------------------
# RULE 7 - Automated LightDock setup
# --------------------------------------------------

rule docking_setup:
    input:
        receptor=lambda wc: config["protein_pairs"][wc.pair]["receptor"],
        ligand=lambda wc: config["protein_pairs"][wc.pair]["ligand"]
    output:
        f"{RESULTS_DIR}/{{pair}}/automated/setup.json"
    params:
        outdir=lambda wc: f"{RESULTS_DIR}/{wc.pair}/automated",
        swarms=SWARMS,
        glowworms=config["docking"]["glowworms"],
        seed=config["docking"]["seed"]
    log:
        "logs/docking_setup/{pair}.log"
    conda:
        "workflow/envs/lightdock_snapshot.yaml"
    benchmark:
        "benchmarks/docking_setup/{pair}.tsv"
    shell:
        """
        mkdir -p {params.outdir} logs/docking_setup

        cp {input.receptor} {params.outdir}/receptor.pdb
        cp {input.ligand} {params.outdir}/ligand.pdb

        (
            cd {params.outdir}
            lgd_setup.py \
                receptor.pdb ligand.pdb \
                -s {params.swarms} \
                -g {params.glowworms} \
                --noh --noxt \
                --seed_points {params.seed}
        ) > {log} 2>&1

        test -s {output}
        """


# --------------------------------------------------
# RULE 8 - Automated LightDock optimization
# --------------------------------------------------

rule docking_run:
    input:
        f"{RESULTS_DIR}/{{pair}}/automated/setup.json"
    output:
        DOCK_GSO_PATTERNS
    params:
        outdir=lambda wc: f"{RESULTS_DIR}/{wc.pair}/automated",
        steps=STEPS,
        seed=config["docking"]["seed"]
    threads: int(config["docking"]["cores"])
    log:
        "logs/docking_run/{pair}.log"
    conda:
        "workflow/envs/lightdock_snapshot.yaml"
    benchmark:
        "benchmarks/docking_run/{pair}.tsv"
    shell:
        """
        mkdir -p logs/docking_run

        (
            cd {params.outdir}
            lgd_run.py \
                setup.json {params.steps} \
                -c {threads} \
                -sg {params.seed}
        ) > {log} 2>&1

        for result in {output}; do
            test -s "$result" || exit 1
        done
        """

include: "workflow/rules/hpr_batch.smk"

include: "workflow/rules/score_analysis.smk"

include: "workflow/rules/best_complex.smk"

include: "workflow/rules/final_project.smk"

include: "workflow/rules/render_complex.smk"

# Optional full workflow for multiple independent protein pairs.
include: "workflow/rules/multipair.smk"
