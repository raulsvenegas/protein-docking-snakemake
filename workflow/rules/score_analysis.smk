
# FHNW Exercise 4
# Automated docking score analysis

rule collect_hpr_scores:
    input:
        HPR_BATCH_TARGETS

    output:
        all_scores="results/e2a_hpr/ensemble/docking_scores.csv",
        best_scores="results/e2a_hpr/ensemble/best_scores_by_model.csv"

    params:
        steps=HPR_STEPS

    log:
        "logs/hpr_score_analysis.log"

    benchmark:
        "benchmarks/hpr_score_analysis.tsv"

    shell:
        """
        mkdir -p logs benchmarks

        python workflow/scripts/collect_docking_scores.py --steps {params.steps} \
            > {log} 2>&1

        test -s {output.all_scores}
        test -s {output.best_scores}
        """
