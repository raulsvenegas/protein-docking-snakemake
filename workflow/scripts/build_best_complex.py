# FHNW Exercise 4 — automatic reconstruction of the highest-scoring docking pose.
# This script is executed by Snakemake in the LightDock Conda environment.
import csv
import subprocess
from pathlib import Path

with open(snakemake.input.scores, newline="") as handle:
    best = next(csv.DictReader(handle))

model = int(best["model"])
swarm = int(best["swarm"])
pose = int(best["pose_index"])
steps = int(snakemake.params.steps)
glowworms = int(snakemake.params.glowworms)

swarm_dir = Path(snakemake.params.base) / f"model_{model}" / f"swarm_{swarm}"
source = swarm_dir / f"lightdock_{pose}.pdb"
log_path = Path(str(snakemake.log[0]))
log_path.parent.mkdir(parents=True, exist_ok=True)

if not source.is_file():
    # Direct executable call: Snakemake activates our rule's Conda environment.
    with log_path.open("w") as log_handle:
        subprocess.run(
            [
                "lgd_generate_conformations.py",
                "--setup", "../setup.json",
                "../receptor.pdb", "../ligand.pdb",
                f"gso_{steps}.out",
                str(max(glowworms, pose + 1)),
            ],
            cwd=swarm_dir,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            check=True,
        )

if not source.is_file():
    raise FileNotFoundError(source)

atoms = [line for line in source.read_text().splitlines() if line.startswith("ATOM  ")]
boundaries = [
    i for i in range(1, len(atoms))
    if int(atoms[i][6:11]) <= int(atoms[i - 1][6:11])
]
if len(boundaries) != 1:
    raise ValueError(f"Expected one receptor–ligand boundary; found {len(boundaries)}")

receptor, ligand = atoms[:boundaries[0]], atoms[boundaries[0]:]
lines = []
serial = 1
for chain_id, group in [("A", receptor), ("B", ligand)]:
    for line in group:
        lines.append(line[:6] + f"{serial:5d}" + line[11:21] + chain_id + line[22:])
        serial += 1
    if chain_id == "A":
        lines.append(f"TER   {serial:5d}      {group[-1][17:26]}")
        serial += 1
lines.append("END")

destination = Path(str(snakemake.output[0]))
destination.parent.mkdir(parents=True, exist_ok=True)
destination.write_text("\n".join(lines) + "\n")

message = (
    f"Best docking pose: model={model}, swarm={swarm}, pose={pose}\n"
    f"Score: {best['score']}\n"
    f"Receptor atoms: {len(receptor)}\n"
    f"Ligand atoms: {len(ligand)}\n"
    f"Output: {destination}\n"
)
with log_path.open("a") as handle:
    handle.write(message)
print(message)
