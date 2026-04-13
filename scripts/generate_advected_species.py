#!/usr/bin/env python3
"""Generate advected_species.txt from a MICM mechanism configuration.

Reads species from either:
  - v1 format: config.json with top-level "species" array
  - v0 format: micm/species.json with "camp-data" array

Writes advected_species.txt containing one lowercase species name per line
for every species with "__is_advected": true.

Usage:
    python3 generate_advected_species.py <chemistry_data_dir>

Example:
    python3 generate_advected_species.py chemistry_data/chapman
"""
import json
import sys
from pathlib import Path


def load_species(config_dir: Path) -> list[dict]:
    """Load species list from v1 or v0 config format."""
    v1 = config_dir / "config.json"
    v0 = config_dir / "micm" / "species.json"
    if v1.exists():
        with open(v1) as f:
            data = json.load(f)
        # v1 may have top-level "species", or nested under "camp-data"
        return data.get("species", data.get("camp-data", []))
    elif v0.exists():
        with open(v0) as f:
            data = json.load(f)
        return data.get("camp-data", data.get("species", []))
    else:
        print(f"ERROR: No config.json or micm/species.json in {config_dir}", file=sys.stderr)
        sys.exit(1)


def main():
    if len(sys.argv) != 2:
        print(__doc__, file=sys.stderr)
        sys.exit(1)

    config_dir = Path(sys.argv[1])
    species = load_species(config_dir)
    advected = sorted(
        s["name"].lower()
        for s in species
        if s.get("__is_advected")
    )

    out_path = config_dir / "advected_species.txt"
    with open(out_path, "w") as f:
        f.write("# Auto-generated from mechanism config — do not edit\n")
        for name in advected:
            f.write(name + "\n")

    print(f"Wrote {len(advected)} advected species to {out_path}")


if __name__ == "__main__":
    main()
