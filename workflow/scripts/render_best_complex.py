# FHNW Exercise 4
# Automated visualization of the best docking complex

import sys
from pathlib import Path
import pymol
from pymol import cmd

input_pdb = Path(sys.argv[1]).resolve()
output_png = Path(sys.argv[2]).resolve()

if not input_pdb.is_file():
    raise FileNotFoundError(input_pdb)

output_png.parent.mkdir(parents=True, exist_ok=True)

# Start PyMOL without opening a graphical window
pymol.finish_launching(["pymol", "-cq"])

# Load the selected protein complex
cmd.load(str(input_pdb), "complex")

# Verify that both protein chains are present
receptor_atoms = cmd.count_atoms("complex and chain A")
ligand_atoms = cmd.count_atoms("complex and chain B")

if receptor_atoms == 0 or ligand_atoms == 0:
    raise ValueError("Both protein chains must be present.")

# Create a scientific protein visualization
cmd.hide("everything", "all")
cmd.show("cartoon", "complex")

# Receptor chain A: cyan
cmd.color("cyan", "complex and chain A")

# Ligand chain B: orange
cmd.color("orange", "complex and chain B")

# Improve presentation quality
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("antialias", 2)

cmd.orient("complex")
cmd.zoom("complex", buffer=5)

# Export high-resolution PNG
cmd.png(
    str(output_png),
    width=1800,
    height=1200,
    dpi=180,
    ray=1
)

if not output_png.is_file() or output_png.stat().st_size == 0:
    raise RuntimeError("PyMOL did not generate the PNG.")

print("SUCCESS: PyMOL visualization generated")
print(f"Receptor atoms: {receptor_atoms}")
print(f"Ligand atoms: {ligand_atoms}")
print(f"Output: {output_png}")
