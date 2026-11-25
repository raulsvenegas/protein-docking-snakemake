
# FHNW Exercise 4
# Automated batch docking of HPr NMR models against EIIA-Glc

# --------------------------------------------------
# CONFIGURATION
# --------------------------------------------------

HPR_FIRST = int(config["hpr_models"]["first"])
HPR_LAST = int(config["hpr_models"]["last"])

HPR_MODELS = list(range(HPR_FIRST, HPR_LAST + 1))

HPR_SOURCE = config["hpr_models"]["source"]
HPR_RECEPTOR = config["protein_pairs"]["e2a_hpr"]["receptor"]

HPR_SWARMS = int(config["docking"]["swarms"])
HPR_GLOWWORMS = int(config["docking"]["glowworms"])
HPR_STEPS = int(config["docking"]["steps"])
HPR_SEED = int(config["docking"]["seed"])
HPR_CORES = int(config["docking"]["cores"])

HPR_BASE = f'{config.get("results_dir", "results")}/e2a_hpr/ensemble'

# Expected outputs for all selected models and swarms
HPR_BATCH_TARGETS = [
    f"{HPR_BASE}/model_{model}/swarm_{swarm}/gso_{HPR_STEPS}.out"
    for model in HPR_MODELS
    for swarm in range(HPR_SWARMS)
]

HPR_SWARM_OUTPUTS = [
    f"{HPR_BASE}/model_{{model}}/swarm_{swarm}/gso_{HPR_STEPS}.out"
    for swarm in range(HPR_SWARMS)
]


# --------------------------------------------------
# RULE 1 - Complete the batch
# --------------------------------------------------

rule batch_hpr:
    input:
        HPR_BATCH_TARGETS


# --------------------------------------------------
# RULE 2 - Extract each HPr NMR model
# --------------------------------------------------

rule extract_hpr_ensemble:
    input:
        HPR_SOURCE
    output:
        "data/prepared/e2a_hpr/ensemble/hpr_model{model}.pdb"
    params:
        model=lambda wc: int(wc.model)
    script:
        "../scripts/extract_nmr_model.py"


# --------------------------------------------------
# RULE 3 - LightDock setup for each model
# --------------------------------------------------

rule setup_hpr_ensemble:
    input:
        receptor=HPR_RECEPTOR,
        ligand="data/prepared/e2a_hpr/ensemble/hpr_model{model}.pdb"
    output:
        f"{HPR_BASE}/model_{{model}}/setup.json"
    params:
        outdir=lambda wc: f"{HPR_BASE}/model_{wc.model}",
        swarms=HPR_SWARMS,
        glowworms=HPR_GLOWWORMS,
        seed=HPR_SEED
    log:
        "logs/hpr_batch/setup_model_{model}.log"
    conda: "../envs/lightdock_snapshot.yaml"
    benchmark: "benchmarks/hpr_batch/setup_model_{model}.tsv"
    shell:
        """
        mkdir -p {params.outdir} logs/hpr_batch

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
# RULE 4 - Run LightDock for each model
# --------------------------------------------------

rule run_hpr_ensemble:
    input:
        f"{HPR_BASE}/model_{{model}}/setup.json"
    output:
        HPR_SWARM_OUTPUTS
    params:
        outdir=lambda wc: f"{HPR_BASE}/model_{wc.model}",
        steps=HPR_STEPS,
        seed=HPR_SEED
    threads: HPR_CORES
    log:
        "logs/hpr_batch/run_model_{model}.log"
    conda: "../envs/lightdock_snapshot.yaml"
    benchmark: "benchmarks/hpr_batch/run_model_{model}.tsv"
    shell:
        """
        mkdir -p logs/hpr_batch

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
