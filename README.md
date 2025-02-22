# meson-structure
Meson structure analysis


## Processing data files

### Overview

The general EIC processing chain look like this: 

<img src="chain.svg" >

In order to run it on multiple farm jobs it is convenient to 
split it in multiple files. So the chain looks like this

<img src="ifarm-processing.svg" >

Scripts are called in this order: 

1. **root_hepmc_converter.py** - converts original root files and split to small *.hepmc chunks
2. **ifarm_process_eg.py** - is a wrapper of `root_hepmc_converter` to do many EG files at once, handling energies and names. 
2. **create_jobs.py** - for each hepmc create scripts how we process it in container and a slurm job
3. **collect_job_stats.py** - Check the actual status of simulation (how many events now generated)
4. Done! Should be ready to analyze

**Running with special EICrecon branch**

As example one can see zdc_lambda/ folder for 
example of Singularity image creation that compiles
EICrecon particular branch. Then `create_jobs.py` has
`--container` flag where you can specify what image to use.  

### Data chain (full)
... and split ...

Our files of interest located at

```
/w/eic-scshelf2104/users/singhav/JLEIC/USERS/trottar/OUTPUTS/raty_eic
```

```bash
python convert_to_hepmc3.py \
      --input-files file_5x41.root file_10x100.root \
      --chunk-size 50000 \
      --events 100000 \
      --output-prefix out_hepmc \
      --events-per-file 20000
```

`ifarm_process_eg.py` is a handy wrapper with hard
coded direcotires for the ifarm. It processes
all MC event generators files at once. 

```bash
python ifarm_process_eg.py
```

```bash

# Generate slurm scripts for 1 file
# We need to bind root directory /volatile/eic/romanov/meson-structure-2025-02
python create_jobs.py \
       -b /volatile/eic/romanov/meson-structure-2025-02 \
       -o /volatile/eic/romanov/meson-structure-2025-02/reco \
       -e 5000 \
       /volatile/eic/romanov/meson-structure-2025-02/eg-hepmc/k_lambda_10x100_5000evt_001.hepmc

# Generate slurm scripts for all files
python create_jobs.py \
       -b /volatile/eic/romanov/meson-structure-2025-02 \
       -o /volatile/eic/romanov/meson-structure-2025-02/reco \
       -e 5000 \
       /volatile/eic/romanov/meson-structure-2025-02/eg-hepmc/*.hepmc

. /volatile/eic/romanov/meson-structure-2025-02/reco/submit_all_slurm_jobs.sh
```

Now for custom Singularity image (for specific EICrecon branch) we
run `create_jobs.py` with --container flag

```bash
python create_jobs.py \
       --container /work/eic/users/romanov/meson-structure-work/eicrecon_lambda_zdc.sif \
       -b /volatile/eic/romanov/meson-structure-2025-02 \
       -o /volatile/eic/romanov/meson-structure-2025-02/reco-zdc-lambda \
       -e 5000 \
       /volatile/eic/romanov/meson-structure-2025-02/eg-hepmc/*.hepmc

. /volatile/eic/romanov/meson-structure-2025-02/reco-zdc-lambda/submit_all_slurm_jobs.sh
```

There are scripts to show how many events exactly were processed by jobs: 

```bash 
python collect_job_stats.py /volatile/eic/romanov/meson-structure-2025-02/reco-zdc-lambda
```

```
python root_hepmc_converter.py -o /w/eic-scshelf2104/users/romanov/meson_structure_2025/temp /w/eic-scshelf2104/users/romanov/meson_structure_2025/temp/k-lambda_10x100.root
```

generate eic-shell and slurm scripts
```
python create_jobs.py -o /w/eic-scshelf2104/users/romanov/meson_structure_2025/temp /w/eic-scshelf2104/users/romanov/meson_structure_2025/temp/temp_001.hepmc
```