
# FHNW Exercise 4
# Collect and rank LightDock docking scores
# Standard Python library only

import argparse
import csv
from pathlib import Path

parser = argparse.ArgumentParser(description="Collect and rank LightDock results")
parser.add_argument("--steps", required=True, type=int, help="Docking steps configured in config.yaml")
args = parser.parse_args()

BASE_DIR = Path("results/e2a_hpr/ensemble")
OUTPUT_ALL = BASE_DIR / "docking_scores.csv"
OUTPUT_BEST = BASE_DIR / "best_scores_by_model.csv"

files = sorted(BASE_DIR.glob(f"model_*/swarm_*/gso_{args.steps}.out"))

if not files:
    raise SystemExit("ERROR: No LightDock output files found.")

rows = []

for file_path in files:
    model = int(file_path.parent.parent.name.split("_")[-1])
    swarm = int(file_path.parent.name.split("_")[-1])

    pose_index = 0

    with file_path.open() as handle:
        for line in handle:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            if not line.startswith("(") or ")" not in line:
                raise ValueError(f"Invalid pose in {file_path}: {line}")

            coordinates, values = line.split(")", 1)
            coordinates = coordinates + ")"

            fields = values.split()

            if len(fields) != 6:
                raise ValueError(
                    f"Unexpected number of columns in {file_path}"
                )

            rec_id = int(fields[0])
            lig_id = int(fields[1])
            luciferin = float(fields[2])
            neighbors = int(fields[3])
            vision_range = float(fields[4])
            score = float(fields[5])

            rows.append({
                "model": model,
                "swarm": swarm,
                "pose_index": pose_index,
                "score": score,
                "rec_id": rec_id,
                "lig_id": lig_id,
                "luciferin": luciferin,
                "neighbors": neighbors,
                "vision_range": vision_range,
                "coordinates": coordinates,
                "source_file": str(file_path),
            })

            pose_index += 1

if not rows:
    raise SystemExit("ERROR: No docking poses found.")

# Higher LightDock scores are ranked first.
rows.sort(key=lambda row: row["score"], reverse=True)

for rank, row in enumerate(rows, start=1):
    row["rank"] = rank

columns = [
    "rank",
    "model",
    "swarm",
    "pose_index",
    "score",
    "rec_id",
    "lig_id",
    "luciferin",
    "neighbors",
    "vision_range",
    "coordinates",
    "source_file",
]

with OUTPUT_ALL.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=columns)
    writer.writeheader()
    writer.writerows(rows)

# Select the highest-scoring pose for each HPr model.
best_by_model = {}

for row in rows:
    model = row["model"]
    if model not in best_by_model:
        best_by_model[model] = row

with OUTPUT_BEST.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=columns)
    writer.writeheader()

    for model in sorted(best_by_model):
        writer.writerow(best_by_model[model])

print("Docking score collection completed!")
print(f"Input files: {len(files)}")
print(f"Total poses: {len(rows)}")
print(f"HPr models: {len(best_by_model)}")
print(f"All scores: {OUTPUT_ALL}")
print(f"Best scores: {OUTPUT_BEST}")

print("\nTop 5 scoring poses:")

for row in rows[:5]:
    print(
        f"Rank {row['rank']}: "
        f"Model {row['model']}, "
        f"Swarm {row['swarm']}, "
        f"Pose {row['pose_index']}, "
        f"Score {row['score']:.4f}"
    )
