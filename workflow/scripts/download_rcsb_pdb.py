"""Download a structure from the RCSB PDB when Snakemake needs it.

Examples: python workflow/scripts/download_rcsb_pdb.py --pdb-id 1BRS \
            --output data/raw/barnase_barstar/1BRS.pdb
"""
import argparse
import re
from pathlib import Path
from urllib.request import urlopen


def download(pdb_id: str, output: Path) -> None:
    if not re.fullmatch(r"[0-9][A-Za-z0-9]{3}", pdb_id):
        raise ValueError("Invalid four-character PDB accession")
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    with urlopen(url, timeout=45) as response:
        content = response.read()
    if not any(line.startswith(b"ATOM  ") for line in content.splitlines()):
        raise ValueError(f"No protein atoms in downloaded PDB: {pdb_id}")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(".pdb.partial")
    temporary.write_bytes(content)
    temporary.replace(output)
    print(f"Source: {url}\nSaved PDB: {output} ({len(content)} bytes)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb-id", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    download(args.pdb_id, args.output)
