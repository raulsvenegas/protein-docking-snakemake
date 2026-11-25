# FHNW Exercise 4
# Automated visualization of the best docking complex

rule render_best_complex:
    input:
        "results/e2a_hpr/ensemble/best_complex_chainsAB.pdb"
    output:
        "results/e2a_hpr/ensemble/best_complex_automated.png"
    conda:
        "../envs/pymol_snapshot.yaml"
    log:
        "logs/render_best_complex.log"
    benchmark:
        "benchmarks/render_best_complex.tsv"
    shell:
        """
        mkdir -p logs benchmarks

        python workflow/scripts/render_best_complex.py \
            {input} {output} > {log} 2>&1

        test -s {output}
        """
