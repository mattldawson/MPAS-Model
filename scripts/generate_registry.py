#!/usr/bin/env python3
"""Generate MPAS Registry.xml and core_interface snippets from a MICM mechanism config.

Reads a mechanism config folder (auto-detects v0 vs v1 format), identifies
advected species (those with "__is_advected": true), and generates:

  1. registry_snippet.xml  — <package> + <var> entries for Registry.xml
  2. core_interface_snippet.F — package activation code for mpas_atm_core_interface.F

Usage:
    python scripts/generate_registry.py <mechanism_folder> [--name <mechanism_name>]

Example:
    python scripts/generate_registry.py chemistry_data/chapman_emis_dep
"""
import argparse
import json
import os
import sys


def load_species(mechanism_path):
    """Load species list from a mechanism config (v0 or v1 format).

    Returns a list of dicts, each with at least 'name' and properties.
    """
    # Try v1 format first: <path>/config.json with top-level "species" key
    v1_config = os.path.join(mechanism_path, "config.json")
    if os.path.isfile(v1_config):
        with open(v1_config) as f:
            data = json.load(f)
        if "species" in data:
            return data["species"]

    # Try v0 format: <path>/micm/species.json with "camp-data" key
    v0_config = os.path.join(mechanism_path, "micm", "species.json")
    if os.path.isfile(v0_config):
        with open(v0_config) as f:
            data = json.load(f)
        if "camp-data" in data:
            return data["camp-data"]

    print(f"ERROR: Cannot find species config in {mechanism_path}", file=sys.stderr)
    sys.exit(1)


def get_advected_species(species_list):
    """Filter species that have '__is_advected': true."""
    advected = []
    for sp in species_list:
        if sp.get("__is_advected", False):
            advected.append(sp["name"])
    return advected


def to_mpas_name(name):
    """Convert a MICM species name to an MPAS-safe variable name.

    Lowercases the name and replaces dots with underscores, matching
    the Fortran to_mpas_name function in mpas_chemistry_utils.F90.
    """
    return name.lower().replace(".", "_")


def generate_registry_snippet(mechanism_name, advected_species):
    """Generate Registry.xml snippet text."""
    pkg_name = f"chem_{mechanism_name}_in"
    lines = [
        f'<!-- Auto-generated for mechanism: {mechanism_name} -->',
        f'<package name="{pkg_name}"',
        f'         description="Chemistry species for {mechanism_name} mechanism"/>',
        '',
        '<!-- Add these inside the scalars var_array in Registry.xml -->',
    ]
    for sp_name in advected_species:
        var_name = to_mpas_name(sp_name)
        lines.append(
            f'<var name="{var_name}" '
            f'array_group="chem_{mechanism_name}" '
            f'units="kg kg^{{-1}}" '
            f'description="{sp_name} mass mixing ratio" '
            f'packages="{pkg_name}"/>'
        )
    return "\n".join(lines) + "\n"


def generate_core_interface_snippet(mechanism_name):
    """Generate mpas_atm_core_interface.F snippet text."""
    pkg_name = f"chem_{mechanism_name}_in"
    lines = [
        f"! Auto-generated package activation for: {mechanism_name}",
        f"nullify({pkg_name}Active)",
        f"call mpas_pool_get_package(packages, '{pkg_name}Active', &",
        f"                           {pkg_name}Active)",
        f"if (associated({pkg_name}Active)) then",
        f"   {pkg_name}Active = .false.",
        f"end if",
    ]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(
        description="Generate MPAS Registry snippets from MICM mechanism config"
    )
    parser.add_argument(
        "mechanism_path",
        help="Path to mechanism config folder (e.g., chemistry_data/chapman_emis_dep)",
    )
    parser.add_argument(
        "--name",
        default=None,
        help="Override mechanism name (default: folder basename)",
    )
    args = parser.parse_args()

    mechanism_path = args.mechanism_path.rstrip("/")
    mechanism_name = args.name or os.path.basename(mechanism_path)

    species_list = load_species(mechanism_path)
    advected = get_advected_species(species_list)

    if not advected:
        print(f"WARNING: No advected species found in {mechanism_path}", file=sys.stderr)

    # Create output directory
    out_dir = os.path.join(mechanism_path, "generated")
    os.makedirs(out_dir, exist_ok=True)

    # Write registry snippet
    reg_path = os.path.join(out_dir, "registry_snippet.xml")
    with open(reg_path, "w") as f:
        f.write(generate_registry_snippet(mechanism_name, advected))
    print(f"  Registry snippet: {reg_path}")

    # Write core interface snippet
    ci_path = os.path.join(out_dir, "core_interface_snippet.F")
    with open(ci_path, "w") as f:
        f.write(generate_core_interface_snippet(mechanism_name))
    print(f"  Core interface snippet: {ci_path}")

    print(f"\nAdvected species for '{mechanism_name}': {', '.join(advected)}")
    print(f"Generated {len(advected)} var entries.")


if __name__ == "__main__":
    main()
