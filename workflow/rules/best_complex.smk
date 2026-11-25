# FHNW Exercise 4
# Reconstruct the best-ranked docking pose using the Snakemake-managed LightDock environment.

rule build_best_complex:
    input:
        scores="results/e2a_hpr/ensemble/docking_scores.csv",
        docking_files=HPR_BATCH_TARGETS
    output:
        "results/e2a_hpr/ensemble/best_complex_chainsAB.pdb"
    params:
        base=HPR_BASE,
        steps=HPR_STEPS,
        glowworms=HPR_GLOWWORMS
    conda:
        "../envs/lightdock_snapshot.yaml"
    log:
        "logs/best_complex.log"
    benchmark:
        "benchmarks/best_complex.tsv"
    script:
        "../scripts/build_best_complex.py"
