"""Reconstruct the highest-ranked docked complex for one configured pair.

Reuses the proven LightDock conformation command from build_best_complex.py.
Requires the LightDock Conda environment for lgd_generate_conformations.py.
"""
import argparse
import csv
import subprocess
from pathlib import Path


def rewrite_two_chain_pdb(source: Path, destination: Path):
    atoms = [l for l in source.read_text().splitlines() if l.startswith("ATOM  ")]
    if not atoms:
        raise ValueError(f"No protein atoms in {source}")
    # LightDock writes receptor atoms followed by ligand atoms, restarting the
    # serial number at the boundary (verified in our EIIA-HPr output).
    boundaries = [i for i in range(1, len(atoms))
                  if int(atoms[i][6:11]) <= int(atoms[i - 1][6:11])]
    if len(boundaries) != 1:
        raise ValueError("Expected exactly one receptor/ligand atom-serial reset "
                         f"in {source}, found {len(boundaries)}")
    receptor = atoms[:boundaries[0]]
    ligand = atoms[boundaries[0]:]
    if not receptor or not ligand:
        raise ValueError("Either receptor or ligand has zero atoms")
    out_lines = []
    serial = 1
    for chain, group in (("A", receptor), ("B", ligand)):
        for line in group:
            out_lines.append(line[:6] + f"{serial:5d}" + line[11:21]
                             + chain + line[22:])
            serial += 1
        if chain == "A":
            out_lines.append(f"TER   {serial:5d}      {group[-1][17:26]}")
            serial += 1
    out_lines.append("END")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(out_lines) + "\n")
    print(f"Receptor (chain A): {len(receptor)} atoms; "
          f"ligand (chain B): {len(ligand)} atoms")


def build(scores: Path, base: Path, destination: Path,
          steps: int, glowworms: int):
    with scores.open(newline="") as handle:
        best = next(csv.DictReader(handle))
    swarm = int(best["swarm"])
    pose = int(best["pose_index"])
    swarm_dir = base / f"swarm_{swarm}"
    source = swarm_dir / f"lightdock_{pose}.pdb"
    # Reconstruct again if needed, so that it always matches the current GSO.
    args = ["lgd_generate_conformations.py", "--setup", "../setup.json",
            "../receptor.pdb", "../ligand.pdb", f"gso_{steps}.out",
            str(max(glowworms, pose + 1))]
    subprocess.run(args, cwd=swarm_dir, check=True)
    if not source.is_file():
        raise FileNotFoundError(source)
    rewrite_two_chain_pdb(source, destination)
    print(f"Best pose: swarm={swarm}, pose={pose}, score={best['score']}")
    print(f"Output: {destination}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--scores", type=Path, required=True)
    p.add_argument("--base", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--steps", type=int, required=True)
    p.add_argument("--glowworms", type=int, required=True)
    a = p.parse_args()
    build(a.scores, a.base, a.output, a.steps, a.glowworms)
