#!/usr/bin/env python3
"""
convert_all_simple.py

A simplified Python script that:
  - Loops over a list of ROOT files
  - Infers beam energies from the presence of 'on41', 'on100', or 'on275'
  - Chooses 'pi_n' or 'k_lambda' from filename
  - Records mapping (original file -> output prefix) in a .info.txt
  - Calls root_hepmc_converter.py with minimal options

Requires the 'sh' library (pip install sh).
"""

import os
import sys
import sh

FILES = [
    "/w/eic-scshelf2104/users/singhav/JLEIC/USERS/trottar/OUTPUTS/raty_eic/k_lambda_crossing_0_10.0on100.0_x0.0001-0.9000_q1.0-500.0.root",
    "/w/eic-scshelf2104/users/singhav/JLEIC/USERS/trottar/OUTPUTS/raty_eic/k_lambda_crossing_0_18.0on275.0_x0.0001-0.9000_q1.0-500.0.root",
    "/w/eic-scshelf2104/users/singhav/JLEIC/USERS/trottar/OUTPUTS/raty_eic/k_lambda_crossing_0_5.0on41.0_x0.0001-0.9000_q1.0-500.0.root",
    #"/w/eic-scshelf2104/users/singhav/JLEIC/USERS/trottar/OUTPUTS/raty_eic/pi_n_crossing_0_10.0on100.0_x0.0001-0.9000_q1.0-500.0.root",
    #"/w/eic-scshelf2104/users/singhav/JLEIC/USERS/trottar/OUTPUTS/raty_eic/pi_n_crossing_0_18.0on275.0_x0.0001-0.9000_q1.0-500.0.root",
    #"/w/eic-scshelf2104/users/singhav/JLEIC/USERS/trottar/OUTPUTS/raty_eic/pi_n_crossing_0_5.0on41.0_x0.0001-0.9000_q1.0-500.0.root",
]

OUTDIR = "/volatile/eic/romanov/meson-structure-2025-03/eg-hepmc"
EVENTS_PER_FILE = 5000

def infer_energies(filename: str) -> str:
    """Return a short beam-energy string based on 'on41', 'on100', or 'on275'."""
    if "on41" in filename:
        return "5x41"
    elif "on100" in filename:
        return "10x100"
    elif "on275" in filename:
        return "18x275"
    else:
        return "??x??"

def infer_particle(filename: str) -> str:
    """Return 'pi_n' if present, otherwise 'k_lambda' (or ??? if not found)."""
    if "pi_n" in filename:
        return "pi_n"
    elif "k_lambda" in filename:
        return "k_lambda"
    else:
        return "???"

def main():
    info_file = os.path.join(OUTDIR, "hepmc_files_simple.info.txt")
    # Overwrite .info.txt at start
    sh.echo("Mapping of ROOT files -> output prefix", _out=info_file, _append=False)

    for f in FILES:
        base = os.path.basename(f)
        
        # Choose energies
        energies = infer_energies(base)
        # Choose particle
        particle = infer_particle(base)

        # Construct output prefix, e.g. "/volatile/eic/.../pi_n_5x41"
        prefix = os.path.join(OUTDIR, f"{particle}_{energies}_{EVENTS_PER_FILE}evt")

        # Append info line
        sh.echo(f"FILE={base} --> PREFIX={prefix}", _out=info_file, _append=True)

        # Print status
        print(f"Converting: {base}")
        print(f"  Using energies='{energies}' and particle='{particle}'")
        print(f"  Output prefix:\n    {prefix}\n")

        # Run your converter with 'sh'
        sh.python(
            "root_hepmc_converter.py",
            "--input-files", f,
            "--output-prefix", prefix,
            "--events-per-file", str(EVENTS_PER_FILE)
        , _out=sys.stdout)

        print(f"Done with {base}")
        print("----------------------------------")

    print("All conversions complete!")

if __name__ == "__main__":
    main()
