"""Parse all LightDock GSO files for ONE independently configured protein pair.

The ranking convention is the same as in the verified HPr ensemble pipeline:
more positive (numerically higher) LightDock scores are ranked first.
"""
import argparse
import csv
from pathlib import Path


COLUMNS = [
    "rank", "swarm", "pose_index", "score", "rec_id", "lig_id", "luciferin",
    "neighbors", "vision_range", "coordinates", "source_file",
]


def collect(base: Path, steps: int, output: Path):
    files = sorted(base.glob(f"swarm_*/gso_{steps}.out"),
                   key=lambda f: int(f.parent.name.split("_")[-1]))
    if not files:
        raise ValueError(f"No GSO files found under {base}")
    rows = []
    for file in files:
        swarm = int(file.parent.name.split("_")[-1])
        pose_idx = 0
        for line in file.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if not line.startswith("(") or ")" not in line:
                raise ValueError(f"Unexpected GSO line in {file}: {line[:100]}")
            coords, values = line.split(")", 1)
            fields = values.split()
            if len(fields) != 6:
                raise ValueError(f"Expected six GSO metadata columns in {file}")
            rows.append(dict(
                swarm=swarm, pose_index=pose_idx, score=float(fields[5]),
                rec_id=int(fields[0]), lig_id=int(fields[1]),
                luciferin=float(fields[2]), neighbors=int(fields[3]),
                vision_range=float(fields[4]), coordinates=coords + ")",
                source_file=str(file),
            ))
            pose_idx += 1
    if not rows:
        raise ValueError("No docking poses found")
    rows.sort(key=lambda row: row["score"], reverse=True)
    for rank, row in enumerate(rows, 1):
        row["rank"] = rank
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Scored {len(rows)} poses in {len(files)} swarms for {base}")
    print(f"Top pose: swarm {rows[0]['swarm']}, index {rows[0]['pose_index']}, "
          f"score {rows[0]['score']:.8f}")
    return rows


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--base", type=Path, required=True)
    p.add_argument("--steps", type=int, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    collect(a.base, a.steps, a.output)
