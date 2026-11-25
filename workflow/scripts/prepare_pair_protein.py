"""Select and clean a protein chain from a PDB; optionally repair/protonate.

- Works with X-ray structures and multi-model NMR PDBs.
- Selects a chain and an optional MODEL, retaining only protein ATOM records.
- Excludes heteroatoms/waters and non-selected alternate locations.
- Renumbers ATOM serials and normalizes the selected chain to A, because the
  receptor and ligand are supplied as separate files to LightDock.
- Optional PDBFixer processing repairs missing *heavy atoms*; whole missing
  residues are deliberately not rebuilt. Hydrogens are optional, and LightDock
  setup uses --noh so they will not contribute to docking.
"""
import argparse
from pathlib import Path


def extract_chain(source: Path, output: Path, chain: str, model: int = 1):
    if len(chain) != 1:
        raise ValueError("A PDB chain ID must be exactly one character")
    lines = source.read_text().splitlines()
    contains_models = any(x.startswith("MODEL ") for x in lines)
    in_model = not contains_models
    selected_model_seen = not contains_models
    atoms = []
    for line in lines:
        if line.startswith("MODEL "):
            in_model = int(line[10:14].strip()) == model
            selected_model_seen |= in_model
            continue
        if line.startswith("ENDMDL"):
            in_model = False
            continue
        if not in_model or not line.startswith("ATOM  "):
            continue
        if len(line) < 54 or line[21] != chain:
            continue
        altloc = line[16]
        if altloc not in (" ", "A"):
            continue
        atoms.append(line)
    if not selected_model_seen:
        raise ValueError(f"MODEL {model} not found in {source}")
    if not atoms:
        raise ValueError(f"No protein ATOM records for chain {chain} in {source}")
    if len(atoms) >= 99999:
        raise ValueError("PDB fixed-width atom count limit exceeded")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as handle:
        for serial, line in enumerate(atoms, 1):
            # ATOM serial [6:11], alternate location [16], chain [21]
            handle.write(line[:6] + f"{serial:5d}" + line[11:16]
                         + " " + line[17:21] + "A" + line[22:] + "\n")
        handle.write("END\n")
    return len(atoms)


def process(source: Path, output: Path, chain: str, model: int,
            repair: bool = False, hydrogens: bool = False):
    if not repair and not hydrogens:
        count = extract_chain(source, output, chain, model)
        print(f"Extracted chain {chain}, model {model}: {count} protein atoms -> {output}")
        return

    import tempfile
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    with tempfile.TemporaryDirectory() as temp:
        clean = Path(temp) / "clean.pdb"
        count = extract_chain(source, clean, chain, model)
        fixer = PDBFixer(filename=str(clean))
        fixer.findMissingResidues()
        fixer.missingResidues = {}  # avoid guessing unresolved entire segments
        fixer.findMissingAtoms()
        if repair:
            fixer.addMissingAtoms()
        elif fixer.missingAtoms:
            raise ValueError("Cannot add hydrogens with heavy atoms missing; enable repair")
        if hydrogens:
            fixer.addMissingHydrogens(7.0)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w") as handle:
            PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
    print(f"Extracted {count} records, repaired={repair}, hydrogens={hydrogens}: {output}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--chain", required=True)
    p.add_argument("--model", type=int, default=1)
    p.add_argument("--repair", type=int, choices=[0, 1], default=0)
    p.add_argument("--hydrogens", type=int, choices=[0, 1], default=0)
    a = p.parse_args()
    process(a.source, a.output, a.chain, a.model, bool(a.repair), bool(a.hydrogens))
