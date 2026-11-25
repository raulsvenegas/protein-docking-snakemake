
from pathlib import Path
from pdbfixer import PDBFixer
from openmm.app import PDBFile

# --------------------------------------------------
# STEP 1: Define input and output paths
# --------------------------------------------------

# Snakemake provides the paths automatically.
# Manual execution uses our EIIA-Glc example.

smk = globals().get("snakemake")

if smk is not None:
    input_pdb = Path(smk.input[0])
    output_pdb = Path(smk.output[0])
else:
    input_pdb = Path("data/prepared/e2a_hpr/eiia_clean.pdb")
    output_pdb = Path("data/prepared/e2a_hpr/eiia_repaired.pdb")

# --------------------------------------------------
# STEP 2: Load the protein
# --------------------------------------------------

print(f"Loading protein: {input_pdb}")

fixer = PDBFixer(filename=str(input_pdb))

# --------------------------------------------------
# STEP 3: Identify missing residues
# --------------------------------------------------

fixer.findMissingResidues()

# Do not reconstruct unresolved whole residues.
# We will document the missing N-terminal segment.
fixer.missingResidues = {}

# --------------------------------------------------
# STEP 4: Identify missing heavy atoms
# --------------------------------------------------

fixer.findMissingAtoms()

print("Missing heavy atoms:")

for residue, atoms in fixer.missingAtoms.items():
    atom_names = ", ".join(atom.name for atom in atoms)
    print(f"{residue.name} {residue.id}: {atom_names}")

# --------------------------------------------------
# STEP 5: Reconstruct missing heavy atoms
# --------------------------------------------------

fixer.addMissingAtoms()

# --------------------------------------------------
# STEP 6: Save the repaired protein structure
# --------------------------------------------------

output_pdb.parent.mkdir(parents=True, exist_ok=True)

with output_pdb.open("w") as file:
    PDBFile.writeFile(
        fixer.topology,
        fixer.positions,
        file,
        keepIds=True
    )

print(f"Repaired structure saved: {output_pdb}")
print("Protein repair completed successfully!")
