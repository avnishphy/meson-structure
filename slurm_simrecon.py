#!/usr/bin/env python3
"""
create_eic_workflow.py

Generates a set of Bash scripts for each provided HepMC file to run a typical EIC workflow:
  1) afterburner -> produces *.afterburner.hepmc
  2) npsim       -> produces *.edm4hep.root
  3) eicrecon    -> produces *.edm4eic.root

For each input file, two scripts are created in the output directory (which must be given
as a full path and mounted inside the container):

  1. <basename>.container.sh
     - Intended to run inside the container.
     - Sources the environment, then runs:
           afterburner, npsim, and eicrecon
     - Produces output files:
           *.afterburner.hepmc, *.edm4hep.root, *.edm4eic.root

  2. <basename>.slurm.sh
     - A Slurm job script that:
           loads singularity,
           executes the container-run script using its full path.

Finally, a master script "submit_all_slurm_jobs.sh" is created to submit all Slurm jobs.

Usage:
  python create_eic_workflow.py --outdir /full/path/to/workflow_scripts \
      --time 12:00:00 --cpus-per-task 4 file1.hepmc file2.hepmc

Ensure that the output directory is specified with a full path and is properly mounted
inside your Singularity container.
"""

import argparse
import os
import textwrap

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate container-run and Slurm job scripts for EIC workflow."
    )
    parser.add_argument("input_files", nargs="+", help="List of input HepMC files (e.g., 'file1.hepmc file2.hepmc').")
    parser.add_argument("-o", "--outdir", required=True, help="Full path to the directory where scripts will be generated. "
                             "This directory must be mounted inside the container.")
    parser.add_argument("--time", default="24:00:00", help="Slurm time limit for each job (default: 24:00:00).")
    parser.add_argument("--cpus-per-task", type=int, default=1, help="Slurm --cpus-per-task (default: 1).")
    parser.add_argument("--mem-per-cpu", default="2G", help="Slurm memory per CPU (default: 2G).")
    parser.add_argument("--account", default="eic", help="Slurm account name (default: eic).")
    parser.add_argument("--container",
                        default="/cvmfs/singularity.opensciencegrid.org/eicweb/eic_xl:nightly",
                        help="Singularity container image (default: eicweb/eic_xl:nightly).")
    return parser.parse_args()

def make_container_script(input_file, output_dir):
    """
    Generates the container-run script for one HepMC file.
    All paths are absolute.
    """
    basename = os.path.splitext(os.path.basename(input_file))[0]
    # Use full paths for input and outputs.
    infile_full = os.path.abspath(input_file)
    afterburn_file = os.path.join(output_dir, f"{basename}.afterburner.hepmc")
    edm4hep_file   = os.path.join(output_dir, f"{basename}.edm4hep.root")
    edm4eic_file   = os.path.join(output_dir, f"{basename}.edm4eic.root")

    container_script = textwrap.dedent(f"""\
    #!/bin/bash
    set -e
    # This script is intended to run INSIDE the Singularity container.
    
    echo "Sourcing EIC environment..."
    # Adjust the path if needed:
    source /opt/detector/epic-main/bin/thisepic.sh

    echo "Running afterburner on {infile_full} -> {afterburn_file}"
    afterburner --input {infile_full} --output {afterburn_file}

    echo "Running npsim on {afterburn_file} -> {edm4hep_file}"
    npsim --runType run --inputFiles {afterburn_file} --outputFile {edm4hep_file} --numberOfEvents 10000

    echo "Running eicrecon on {edm4hep_file} -> {edm4eic_file}"
    eicrecon -Ppodio:output_file={edm4eic_file} -Pplugins=janadot {edm4hep_file}

    echo "All steps completed for {basename}!"
    """)
    return container_script

def make_slurm_script(infile, outdir, container_image, time_limit, cpus_per_task, mem_per_cpu, account):
    """
    Generates the Slurm job script for one HepMC file.
    Uses full paths to ensure the container can find the scripts.
    """

    # Use full path for the container-run script

    basename = os.path.splitext(os.path.basename(infile))[0]
    container_script_full = os.path.join(outdir, f"{basename}.container.sh")
    slurm_script = textwrap.dedent(f"""\
    #!/bin/bash
    #SBATCH --account={account}
    #SBATCH --job-name={basename}
    #SBATCH --time={time_limit}
    #SBATCH --cpus-per-task={cpus_per_task}
    #SBATCH --mem-per-cpu={mem_per_cpu}
    #SBATCH --output={os.path.join(outdir, f"{basename}.slurm.log")}

    set -e

    # Ensure singularity is available
    if ! command -v singularity &> /dev/null; then
      echo "singularity not found. Please load the module or install singularity."
      exit 1
    fi

    echo "Running job {basename} on $(hostname)"
    echo "Using container image: {container_image}"
    echo "Executing container-run script: {container_script_full}"

    # Execute the container-run script inside the container.
    singularity exec -B {outdir}:{outdir} {container_image} {container_script_full}

    echo "Slurm job finished for {basename}!"
    """)
    return slurm_script

def main():
    args = parse_args()
    # Ensure outdir is an absolute path.
    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    slurm_scripts = []
    for infile in args.input_files:
        # Create container-run script with full paths.
        container_text = make_container_script(
            input_file=infile,
            output_dir=outdir
        )
        basename = os.path.splitext(os.path.basename(infile))[0]
        container_filename = os.path.join(outdir, f"{basename}.container.sh")
        with open(container_filename, "w") as f:
            f.write(container_text)
        os.chmod(container_filename, 0o755)

        # Create Slurm job script.
        slurm_text = make_slurm_script(
            infile=infile,
            outdir=outdir,
            container_image=args.container,
            time_limit=args.time,
            cpus_per_task=args.cpus_per_task,
            mem_per_cpu=args.mem_per_cpu,
            account=args.account
        )
        slurm_filename = os.path.join(outdir, f"{basename}.slurm.sh")
        with open(slurm_filename, "w") as f:
            f.write(slurm_text)
        os.chmod(slurm_filename, 0o755)
        slurm_scripts.append(os.path.basename(slurm_filename))

    # Create a master submit script to submit all slurm jobs.
    master_script = os.path.join(outdir, "submit_all_slurm_jobs.sh")
    with open(master_script, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("set -e\n\n")
        f.write("# This script submits all generated slurm jobs.\n\n")
        for s in slurm_scripts:
            f.write(f"sbatch {s}\n")
        f.write("\necho \"All slurm jobs submitted!\"\n")
    os.chmod(master_script, 0o755)

    print(f"Done! Created scripts in {outdir}:")
    print(f"  Container-run scripts: {', '.join([os.path.splitext(s)[0] + '.container.sh' for s in slurm_scripts])}")
    print(f"  Slurm scripts: {', '.join(slurm_scripts)}")
    print(f"  Master submit script: {master_script}")

if __name__ == "__main__":
    main()
