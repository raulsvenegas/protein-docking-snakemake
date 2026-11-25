
# FHNW Exercise 4
# Final protein-protein docking workflow target

rule final_project:
    input:
        scores="results/e2a_hpr/ensemble/docking_scores.csv",
        best_scores="results/e2a_hpr/ensemble/best_scores_by_model.csv",
        best_complex="results/e2a_hpr/ensemble/best_complex_chainsAB.pdb",
        figure="results/e2a_hpr/ensemble/best_complex_automated.png"
