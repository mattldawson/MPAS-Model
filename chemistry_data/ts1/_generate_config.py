#!/usr/bin/env python3
"""Generate MPAS-specific TS1 config files from the musica TS1 mechanism.

Creates:
  chemistry_data/ts1/config.json          - TS1 mechanism with MPAS properties
  chemistry_data/ts1/tuvx_micm_mapping.json - TUV-x → MICM rate parameter mapping
  chemistry_data/ts1/initial_conditions.csv - Height-profile ICs for MPAS init
  chemistry_data/ts1/static_rate_params.csv - Default rates for non-TUV-x reactions
  chemistry_data/ts1/aerosol_profile.csv   - Aerosol stub (radius, number density vs height)
  chemistry_data/ts1/generated/registry_snippet.xml
  chemistry_data/ts1/generated/core_interface_snippet.F
"""
import json
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# This file lives in MPAS-Model/chemistry_data/ts1/, so MPAS root is two
# levels up.
MPAS_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
MUSICA_DIR = os.path.dirname(MPAS_DIR)

TS1_JSON = os.path.join(MUSICA_DIR, "configs", "v1", "ts1", "ts1.json")
TS1_IC = os.path.join(MUSICA_DIR, "configs", "v1", "ts1", "initial_conditions.csv")
TUVX_JSON = os.path.join(MUSICA_DIR, "configs", "tuvx", "ts1_tsmlt.json")
OUT_DIR = SCRIPT_DIR

# ──────────── Species classification ────────────

# Third body
THIRD_BODY = {"M"}

# Constant VMR species (not advected, set from VMR)
CONSTANT_VMR = {
    "N2":  0.7808,
    "O2":  0.2095,
    "H2O": 0.01,       # ~1% tropospheric average; JW test has no moist physics
    "CO2": 4.15e-4,
    "H2":  5.5e-7,
}

# Molecular weights for species that have MW=0 in ts1.json but need it for MPAS
MW_FIXES = {
    "N2":    0.0280134,
    "O2":    0.0319988,
    "H2O":   0.0180153,
    "CO2":   0.0440095,
    "H2":    0.0020159,
    "M":     0.0289644,  # air
    "COF2":  0.0660071,
    "COFCL": 0.0824621,
    "HF":    0.0200063,
    "F":     0.0189984,
    "OCS":   0.0600751,
    "S":     0.0320650,
    "SF6":   0.1460554,
    "SO":    0.0480640,
    "SO3":   0.0800630,
    "CH4":   0.0160425,  # listed without MW in some versions
}

# Short-lived radicals: diagnosed locally, NOT advected
SHORT_LIVED = {
    # Atoms
    "O", "O1D", "H", "N", "CL", "BR", "F", "S",
    # Simple radicals
    "OH", "HO2", "NO3", "SO", "SO3",
    # Halogen radicals
    "CLO", "BRO", "OCLO", "CL2O2",
    # Peroxy radicals (RO2 family)
    "CH3O2", "C2H5O2", "C3H7O2", "CH3CO3", "MCO3",
    "ISOPAO2", "ISOPBO2", "ISOPNO3",
    "MACRO2", "NTERPO2", "TERPO2", "TERP2O2",
    "TOLO2", "XYLENO2", "XYLOLO2",
    "PO2", "EO2", "ENEO2", "MEKO2", "XO2",
    "RO2", "ALKO2", "MALO2", "DICARBO2", "MDIALO2",
    "C6H5O2", "ACBZO2", "BENZO2", "BZOO", "PHENO2",
    "HOCH2OO",
    # Unstable intermediates
    "EO", "PHENO",
}

# Species that are purely internal to MICM (VBS intermediates, passive tracers)
# These keep MW=0 and are not tracked by MPAS
MICM_INTERNAL = {
    "BCARYO2VBS", "BENZO2VBS", "ISOPO2VBS", "IVOCO2VBS",
    "MTERPO2VBS", "TOLUO2VBS", "XYLEO2VBS",
    "E90", "NH_5", "NH_50", "ST80_25", "sink",
    "NH4",
    "soa1_a1", "soa1_a2", "soa2_a1", "soa2_a2",
    "soa3_a1", "soa3_a2", "soa4_a1", "soa4_a2",
    "soa5_a1", "soa5_a2",
    "SOAG0", "SOAG1", "SOAG2", "SOAG3", "SOAG4",
    "IVOC", "SVOC", "MTERP", "BCARY", "SF6",
}

# TUV-x profile species (provide gas profiles to TUV-x)
TUVX_PROFILES = {
    "O3": {"profile_name": "O3", "scale_height_km": 4.5, "column_density_method": "arithmetic"},
    "O2": {"profile_name": "O2", "scale_height_km": 8.01, "column_density_method": "geometric"},
}

# Surface emissions (molec cm⁻² s⁻¹) — order-of-magnitude values
EMISSIONS = {
    "NO":     1.0e11,
    "CO":     5.0e11,
    "CH2O":   1.0e10,
    "ISOP":   5.0e11,
    "C2H4":   5.0e9,
    "C2H6":   5.0e9,
    "C3H6":   2.0e9,
    "C3H8":   2.0e9,
    "BIGALK":  1.0e10,
    "BIGENE":  2.0e9,
    "TOLUENE": 5.0e9,
    "XYLENES": 3.0e9,
    "BENZENE": 2.0e9,
    "CH3CHO":  5.0e9,
    "CH3COCH3": 3.0e9,
    "MEK":     2.0e9,
    "DMS":     5.0e10,
    "SO2":     5.0e10,
    "NH3":     5.0e10,
    "HCN":     1.0e9,
    "CH3OH":   5.0e9,
    "C2H5OH":  3.0e9,
    "HCOOH":   2.0e9,
    "CH3COOH": 2.0e9,
}

# Dry deposition velocities (cm s⁻¹)
DEPOSITION = {
    "O3":     0.4,
    "NO2":    0.3,
    "HNO3":   2.0,
    "H2O2":   1.0,
    "SO2":    0.8,
    "CH2O":   0.5,
    "HCOOH":  0.3,
    "CH3COOH": 0.2,
    "HCL":    2.0,
    "NH3":    1.0,
    "N2O5":   2.0,
    "PAN":    0.3,
    "MPAN":   0.3,
    "HNO3":   2.0,
    "CH3OOH": 0.3,
    "H2SO4":  2.0,
    "CH3CHO": 0.2,
    "GLYALD":  0.3,
    "GLYOXAL": 0.3,
    "HYAC":   0.3,
    "MACR":   0.2,
    "MVK":    0.2,
}


def load_ts1():
    with open(TS1_JSON) as f:
        return json.load(f)


def load_tuvx_reactions():
    with open(TUVX_JSON) as f:
        data = json.load(f)
    rxns = data.get("photolysis", {}).get("reactions", [])
    return {r["name"] for r in rxns}


def load_ts1_initial_conditions():
    """Parse musica-format initial_conditions.csv into dicts."""
    conc = {}
    user_rates = {}
    surf_params = {}
    env = {}
    with open(TS1_IC) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            key = parts[0]
            if key.startswith("CONC."):
                species = key[5:]
                conc[species] = float(parts[1])
            elif key.startswith("USER."):
                name = key[5:]
                user_rates[name] = float(parts[1])
            elif key.startswith("SURF."):
                name = key[5:]
                surf_params[name] = (float(parts[1]), float(parts[2]))
            elif key.startswith("ENV."):
                env[key[4:]] = float(parts[1])
    return conc, user_rates, surf_params, env


def classify_species(ts1_data):
    """Classify each species and determine MPAS properties."""
    results = []
    for sp in ts1_data["species"]:
        name = sp["name"]
        mw = sp.get("molecular weight [kg mol-1]", 0)

        # Fix MW if needed
        if mw == 0 and name in MW_FIXES:
            mw = MW_FIXES[name]

        props = {
            "name": name,
            "molecular weight [kg mol-1]": mw,
        }

        # Copy description if present
        if "__description" in sp:
            props["__description"] = sp["__description"]

        # Absolute tolerance
        props["__absolute tolerance"] = sp.get("__absolute tolerance", 1e-12)

        if name in THIRD_BODY:
            props["is third body"] = True
        elif name in CONSTANT_VMR:
            props["__mpas_constant_vmr"] = CONSTANT_VMR[name]
            if name in TUVX_PROFILES:
                tp = TUVX_PROFILES[name]
                props["__is_tuvx_profile"] = tp["profile_name"]
                props["__tuvx_scale_height [km]"] = tp["scale_height_km"]
                props["__tuvx_column_density_method"] = tp["column_density_method"]
        elif name in MICM_INTERNAL:
            pass  # Keep MW=0 (or whatever it is), won't be tracked by MPAS
        elif name in SHORT_LIVED:
            # Short-lived radicals: give proper MW but don't mark as advected
            if mw == 0 and name in MW_FIXES:
                props["molecular weight [kg mol-1]"] = MW_FIXES[name]
            # Not advected - diagnosed locally
        else:
            # Long-lived species: advected
            if mw > 0:
                props["__is_advected"] = True
                if name in TUVX_PROFILES:
                    tp = TUVX_PROFILES[name]
                    props["__is_tuvx_profile"] = tp["profile_name"]
                    props["__tuvx_scale_height [km]"] = tp["scale_height_km"]
                    props["__tuvx_column_density_method"] = tp["column_density_method"]
                if name in EMISSIONS:
                    props["__mpas_surface_emission_flux [molec cm-2 s-1]"] = EMISSIONS[name]
                if name in DEPOSITION:
                    props["__mpas_surface_deposition_velocity [cm s-1]"] = DEPOSITION[name]

        results.append(props)
    return results


def build_config(ts1_data, species_list):
    """Build MPAS-specific config.json in v1 format.

    Appends EMISSION reactions for species in EMISSIONS and
    FIRST_ORDER_LOSS reactions for species in DEPOSITION.
    These reactions provide the MICM rate parameters (EMIS.<name>,
    LOSS.<name>) that the MPAS emissions/deposition modules set.
    """
    reactions = list(ts1_data["reactions"])

    # Add EMISSION reactions
    for name in sorted(EMISSIONS):
        reactions.append({
            "type": "EMISSION",
            "name": name,
            "gas phase": "gas",
            "products": [{"species name": name, "coefficient": 1}],
        })

    # Add FIRST_ORDER_LOSS reactions
    for name in sorted(DEPOSITION):
        reactions.append({
            "type": "FIRST_ORDER_LOSS",
            "name": name,
            "gas phase": "gas",
            "reactants": [{"species name": name, "coefficient": 1}],
        })

    config = {
        "version": ts1_data["version"],
        "name": ts1_data["name"] + " (MPAS-adapted)",
        "species": species_list,
        "phases": ts1_data["phases"],
        "reactions": reactions,
    }
    return config


def build_tuvx_mapping(tuvx_reactions, ts1_data):
    """Build TUV-x → MICM rate parameter mapping."""
    # Get MICM USER_DEFINED reaction names
    micm_user = set()
    for rxn in ts1_data["reactions"]:
        if rxn["type"] == "USER_DEFINED":
            micm_user.add(rxn["name"])

    mapping = []
    for tuvx_name in sorted(tuvx_reactions):
        # Direct match
        if tuvx_name in micm_user:
            mapping.append({"source": tuvx_name, "target": f"USER.{tuvx_name}"})
        elif tuvx_name == "jno_i" and "jno" in micm_user:
            mapping.append({"source": "jno_i", "target": "USER.jno"})
        else:
            print(f"  WARNING: TUV-x reaction '{tuvx_name}' has no MICM counterpart", file=sys.stderr)

    return mapping


def build_initial_conditions(species_list, conc_data):
    """Build height-profile initial_conditions.csv.

    For TS1, we use a simplified 2-height profile: surface + stratosphere.
    Species with known tropospheric/stratospheric values get a basic profile;
    others get a uniform value from the musica IC data (VMR).
    """
    # Get advected species
    advected = [sp for sp in species_list if sp.get("__is_advected")]
    if not advected:
        return ""

    # Heights (km): surface to ~80 km with key levels
    heights = [0, 2, 5, 8, 12, 16, 20, 25, 30, 40, 50, 60, 80]

    # Species that have meaningful vertical structure
    # For most TS1 species, use uniform profile from IC data
    # Special profiles for O3 (stratospheric peak) and a few others

    # O3 profile (VMR) — realistic Chapman-like
    o3_profile = [
        5.0e-8, 6.1e-8, 7.9e-8, 9.6e-8, 1.2e-7,  # troposphere
        1.7e-7, 1.9e-6, 4.5e-6,                      # tropopause-lower strat
        8.9e-6, 7.1e-6, 2.0e-6, 5.8e-7, 4.8e-8       # mid-upper strat
    ]

    # CO profile — decreases with altitude
    co_vmr = conc_data.get("CO", 3.25e-6)
    co_profile = [co_vmr * f for f in [1.0, 1.0, 0.9, 0.8, 0.7, 0.5, 0.3, 0.15, 0.08, 0.03, 0.01, 0.005, 0.002]]

    # NO/NO2 — surface peak + stratospheric source
    no_vmr = conc_data.get("NO", 1.5e-8)
    no_profile = [no_vmr * f for f in [1.0, 0.8, 0.5, 0.3, 0.1, 0.05, 0.1, 0.3, 0.5, 0.3, 0.1, 0.05, 0.01]]
    no2_vmr = conc_data.get("NO2", 6.4e-8)
    no2_profile = [no2_vmr * f for f in [1.0, 0.8, 0.5, 0.3, 0.1, 0.05, 0.1, 0.2, 0.3, 0.15, 0.05, 0.02, 0.005]]

    # CH4 — well-mixed in troposphere, decreases in stratosphere
    ch4_vmr = conc_data.get("CH4", 7.0e-5)
    ch4_profile = [ch4_vmr * f for f in [1.0, 1.0, 1.0, 1.0, 1.0, 0.95, 0.8, 0.5, 0.3, 0.1, 0.05, 0.02, 0.01]]

    special_profiles = {
        "O3": o3_profile,
        "CO": co_profile,
        "NO": no_profile,
        "NO2": no2_profile,
        "CH4": ch4_profile,
    }

    # Build CSV
    names = [sp["name"] for sp in advected]
    header = "height_km," + ",".join(names)
    rows = [header]

    for hi, h in enumerate(heights):
        vals = []
        for sp_name in names:
            if sp_name in special_profiles:
                vals.append(f"{special_profiles[sp_name][hi]:.6e}")
            else:
                # Uniform profile from musica IC data
                vmr = conc_data.get(sp_name, 1.0e-15)
                # Decrease to near-zero above ~50 km for tropospheric species
                if h > 50:
                    factor = 0.001
                elif h > 30:
                    factor = 0.01
                elif h > 20:
                    factor = 0.1
                else:
                    factor = 1.0
                vals.append(f"{vmr * factor:.6e}")
        rows.append(f"{h}," + ",".join(vals))

    return "\n".join(rows) + "\n"


def build_static_rate_params(tuvx_reactions, ts1_data, user_rates):
    """Build static_rate_params.csv for non-TUV-x USER_DEFINED reactions."""
    # Get MICM USER_DEFINED reactions not covered by TUV-x
    tuvx_names = set(tuvx_reactions)
    # Account for jno_i → jno mapping
    tuvx_names.add("jno")  # mapped from jno_i

    lines = ["# Static rate parameters for USER_DEFINED reactions not computed by TUV-x"]
    lines.append("# Format: USER.<name>,<rate_value>")
    lines.append("# These are set once at init and remain constant (TUV-x rates override dynamically)")

    for rxn in sorted(ts1_data["reactions"], key=lambda r: r.get("name", "")):
        if rxn["type"] != "USER_DEFINED":
            continue
        name = rxn["name"]
        if name in tuvx_names:
            continue  # Will be computed by TUV-x
        rate = user_rates.get(name, 0.0)
        lines.append(f"USER.{name},{rate:.6e}")

    return "\n".join(lines) + "\n"


def build_aerosol_profile(surf_params):
    """Build aerosol_profile.csv for SURFACE reaction stub."""
    lines = [
        "# Aerosol stub profile for SURFACE reactions",
        "# Format: height_km,effective_radius_m,number_density_m-3",
        "# Provides a simple background aerosol for heterogeneous chemistry",
        "# Default: 100 nm radius, 1e9 particles/m³ in boundary layer, decreasing with altitude",
    ]

    # Get default values from SURF params (they're all the same in the IC file)
    default_radius = 2.82e-8  # ~28 nm (from IC file)
    default_num = 1.0e11      # particles/m³ (from IC file)

    heights_factors = [
        (0, 1.0), (1, 1.0), (2, 0.8), (3, 0.5),
        (5, 0.3), (8, 0.1), (12, 0.05), (16, 0.02),
        (20, 0.01), (30, 0.005), (50, 0.001), (80, 0.0001),
    ]

    lines.append("height_km,effective_radius_m,number_density_m-3")
    for h, f in heights_factors:
        lines.append(f"{h},{default_radius:.6e},{default_num * f:.6e}")

    return "\n".join(lines) + "\n"


def build_registry_snippet(species_list):
    """Generate Registry.xml snippet for TS1 species."""
    lines = ['<!-- Auto-generated for mechanism: ts1 -->',
             '<package name="chem_ts1_in"',
             '         description="TS1 gas-phase chemistry species"/>',
             '',
             '<!-- Add these inside the scalars var_array in Registry.xml -->']

    advected = [sp for sp in species_list if sp.get("__is_advected")]
    for sp in sorted(advected, key=lambda s: s["name"].lower()):
        name_lower = sp["name"].lower()
        desc = sp.get("__description", f"{sp['name']} mass mixing ratio")
        if "mass mixing ratio" not in desc.lower():
            desc = f"{desc} mass mixing ratio"
        lines.append(f'<var name="{name_lower}" array_group="ts1" units="kg kg^{{-1}}"'
                      f' description="{desc}"'
                      f' packages="chem_ts1_in"/>')

    return "\n".join(lines) + "\n"


def build_core_interface_snippet():
    """Generate core_interface_snippet.F."""
    return """! Auto-generated package activation for: ts1
nullify(chem_ts1_inActive)
call mpas_pool_get_package(packages, 'chem_ts1_inActive', &
                           chem_ts1_inActive)
if (associated(chem_ts1_inActive)) then
   chem_ts1_inActive = .false.
end if
"""


def main():
    print("Loading TS1 mechanism...")
    ts1_data = load_ts1()
    print(f"  {len(ts1_data['species'])} species, {len(ts1_data['reactions'])} reactions")

    print("Loading TUV-x config...")
    tuvx_reactions = load_tuvx_reactions()
    print(f"  {len(tuvx_reactions)} photolysis reactions")

    print("Loading TS1 initial conditions...")
    conc_data, user_rates, surf_params, env = load_ts1_initial_conditions()
    print(f"  {len(conc_data)} concentrations, {len(user_rates)} rates, {len(surf_params)} surface params")

    print("Classifying species...")
    species_list = classify_species(ts1_data)
    advected = [sp for sp in species_list if sp.get("__is_advected")]
    constant = [sp for sp in species_list if "__mpas_constant_vmr" in sp]
    short = [sp for sp in species_list
             if sp["name"] in SHORT_LIVED and sp.get("molecular weight [kg mol-1]", 0) > 0]
    internal = [sp for sp in species_list if sp["name"] in MICM_INTERNAL]
    print(f"  Advected: {len(advected)}")
    print(f"  Constant VMR: {len(constant)}")
    print(f"  Short-lived (not advected): {len(short)}")
    print(f"  MICM-internal: {len(internal)}")

    # Create output directories
    os.makedirs(os.path.join(OUT_DIR, "generated"), exist_ok=True)

    # 1. config.json
    print("Writing config.json...")
    config = build_config(ts1_data, species_list)
    with open(os.path.join(OUT_DIR, "config.json"), "w") as f:
        json.dump(config, f, indent=2)

    # 2. TUV-x mapping
    print("Writing tuvx_micm_mapping.json...")
    mapping = build_tuvx_mapping(tuvx_reactions, ts1_data)
    print(f"  {len(mapping)} mapped reactions")
    with open(os.path.join(OUT_DIR, "tuvx_micm_mapping.json"), "w") as f:
        json.dump(mapping, f, indent=2)

    # 3. Initial conditions
    print("Writing initial_conditions.csv...")
    ic_csv = build_initial_conditions(species_list, conc_data)
    with open(os.path.join(OUT_DIR, "initial_conditions.csv"), "w") as f:
        f.write(ic_csv)

    # 4. Static rate parameters
    print("Writing static_rate_params.csv...")
    static_csv = build_static_rate_params(tuvx_reactions, ts1_data, user_rates)
    with open(os.path.join(OUT_DIR, "static_rate_params.csv"), "w") as f:
        f.write(static_csv)

    # 5. Aerosol profile
    print("Writing aerosol_profile.csv...")
    aerosol_csv = build_aerosol_profile(surf_params)
    with open(os.path.join(OUT_DIR, "aerosol_profile.csv"), "w") as f:
        f.write(aerosol_csv)

    # 6. Registry snippet
    print("Writing generated/registry_snippet.xml...")
    snippet = build_registry_snippet(species_list)
    with open(os.path.join(OUT_DIR, "generated", "registry_snippet.xml"), "w") as f:
        f.write(snippet)

    # 7. Core interface snippet
    print("Writing generated/core_interface_snippet.F...")
    snippet_f = build_core_interface_snippet()
    with open(os.path.join(OUT_DIR, "generated", "core_interface_snippet.F"), "w") as f:
        f.write(snippet_f)

    print(f"\nDone! Files written to {OUT_DIR}/")
    print(f"\nAdvected species for Registry.xml ({len(advected)}):")
    for sp in sorted(advected, key=lambda s: s["name"]):
        flags = []
        if sp["name"] in EMISSIONS:
            flags.append("E")
        if sp["name"] in DEPOSITION:
            flags.append("D")
        if sp["name"] in TUVX_PROFILES:
            flags.append("T")
        flag_str = f" [{','.join(flags)}]" if flags else ""
        print(f"  {sp['name']}{flag_str}")


if __name__ == "__main__":
    main()
