
from pathlib import Path

# Input and output paths provided by Snakemake
input_file = Path(snakemake.input[0])
output_file = Path(snakemake.output[0])
model_number = int(snakemake.params.model)

selected = False
atoms = []

# Read the original NMR structure
with input_file.open() as pdb:
    for line in pdb:
        if line.startswith("MODEL "):
            selected = int(line.split()[1]) == model_number

        elif line.startswith("ENDMDL") and selected:
            break

        elif selected and line.startswith(("ATOM  ", "HETATM", "TER   ")):
            atoms.append(line)

# Check that the selected model contains protein atoms
if not any(line.startswith("ATOM  ") for line in atoms):
    raise ValueError(f"Model {model_number} contains no protein atoms")

# Save the selected model as a separate PDB file
output_file.parent.mkdir(parents=True, exist_ok=True)
output_file.write_text("".join(atoms) + "END\n")
