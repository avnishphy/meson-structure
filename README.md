# meson-structure
Meson structure analysis


## Processing data files

### Overview

The general EIC processing chain look like this: 

<img src="chain.svg" >

Scripts are called in this order: 

1. **root_hepmc_converter.py** - converts original root files and split to small *.hepmc chunks
2. **create_jobs.py** - for each hepmc file create what to do
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

# Generate slurm scripts for 1 file
# We need to bind root directory /volatile/eic/romanov/meson-structure-2025-02
python slurm_simrecon.py \
       -b /volatile/eic/romanov/meson-structure-2025-02 \
       -o /volatile/eic/romanov/meson-structure-2025-02/reco \
       -e 5000 \
       /volatile/eic/romanov/meson-structure-2025-02/eg-hepmc/k_lambda_10x100_5000evt_001.hepmc

# Generate slurm scripts for all files
python slurm_simrecon.py \
       -b /volatile/eic/romanov/meson-structure-2025-02 \
       -o /volatile/eic/romanov/meson-structure-2025-02/reco \
       -e 5000 \
       /volatile/eic/romanov/meson-structure-2025-02/eg-hepmc/*.hepmc
```

```bash
python convert_to_hepmc3.py \
      --input-files file_5x41.root file_10x100.root \
      --chunk-size 50000 \
      --events 100000 \
      --output-prefix out_hepmc \
      --events-per-file 20000
```


```
python root_hepmc_converter.py -o /w/eic-scshelf2104/users/romanov/meson_structure_2025/temp /w/eic-scshelf2104/users/romanov/meson_structure_2025/temp/k-lambda_10x100.root
```

generate eic-shell and slurm scripts
```
python slurm_simrecon.py -o /w/eic-scshelf2104/users/romanov/meson_structure_2025/temp /w/eic-scshelf2104/users/romanov/meson_structure_2025/temp/temp_001.hepmc
```