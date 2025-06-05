# Simulation campaigns

This page documents the 2025-03 meson structure simulation campaign.

> (!) For the list of files go to [DATA PAGE](data.md) 


## Campaign 2025-05

### Overview

The Campaign 2025-05 run to make data current with EIC EPIC software updates. 

The campaign includes simulations with three beam energy configurations:

1. 5x41 GeV
2. 10x100 GeV
3. 18x275 GeV

Each configuration has multiple files (indexed 001-200) with 5000 events per file.

```yaml
timestamp: '2025-06-04T12:05:06.867483'
input_file: /volatile/eic/romanov/meson-structure-2025-06/eg-hepmc/*.hepmc
container_image: /cvmfs/singularity.opensciencegrid.org/eicweb/eic_xl:nightly
```

### Data Location

The campaign data is stored in the following locations (on JLab farm):

- HEPMC files:   
  `/volatile/eic/romanov/meson-structure-2025-03/eg-hepmc`
- reco info: 
  `/volatile/eic/romanov/meson-structure-2025-06/reco`



### Processing Commands

The exact commands used in this campaign:

```bash
mkdir /volatile/eic/romanov/meson-structure-2025-06
mkdir /volatile/eic/romanov/meson-structure-2025-06/reco

# Using previous converted hepmc files
cp -r /volatile/eic/romanov/meson-structure-2025-03/eg-hepmc /volatile/eic/romanov/meson-structure-2025-06

# Creating jobs (using latest eic_xl container)
cd /home/romanov/meson-structure-work/meson-structure/full-sim-pipeline
python create_jobs.py \
       -b /volatile/eic/romanov/meson-structure-2025-06 \
       -o /volatile/eic/romanov/meson-structure-2025-06/reco \
       -e 5000 \
       /volatile/eic/romanov/meson-structure-2025-06/eg-hepmc/*.hepmc

# Submit jobs
cd /volatile/eic/romanov/meson-structure-2025-06/reco/
submit_all_slurm_jobs.sh
```



## Campaign 2025-03


### Overview

The Campaign 2025-03 is focused on testing the new ZDC lambda reconstruction 
algorithm using the latest ePIC software. 
This campaign reuses meson-structure-2025-02, 
reusing some of the existing `hepmc` files while implementing improved reconstruction techniques.

### Processing Details

The campaign uses a processing pipeline that converts Monte Carlo event generator files 
to a format suitable for full detector simulation and reconstruction:

1. MCEG files are converted to HEPMC format (splitting large files into manageable chunks)
2. The HEPMC files are processed through the latest ePIC reconstruction software
3. Output files include both EDM4EIC format and histogram files

### Data Location

The campaign data is stored in the following locations:

- HEPMC files:   
  `/volatile/eic/romanov/meson-structure-2025-03/eg-hepmc`
  
   Note: These are linked from the previous campaign:
   
   `/volatile/eic/romanov/meson-structure-2025-02/eg-hepmc`
- Reconstruction output:  
  `/volatile/eic/romanov/meson-structure-2025-03/reco`

### Processing Commands

The exact commands used in this campaign:

```bash
# Original MCEG files location
# Note: Using the same hepmc files as the previous campaign
cd /volatile/eic/romanov/meson-structure-2025-03
ln -s /volatile/eic/romanov/meson-structure-2025-02/eg-hepmc eg-hepmc

cd /home/romanov/meson-structure-work/meson-structure/full-sim-pipeline

# Generate job scripts (using latest eic_xl container)
python create_jobs.py \
       -b /volatile/eic/romanov/meson-structure-2025-03 \
       -o /volatile/eic/romanov/meson-structure-2025-03/reco \
       -e 5000 \
       /volatile/eic/romanov/meson-structure-2025-03/eg-hepmc/*.hepmc

# Submit jobs
cd /volatile/eic/romanov/meson-structure-2025-03/reco/
submit_all_slurm_jobs.sh
```


## Accessing the Data

The data can be accessed using XRootD:

```bash
# XRootD root URL
root://dtn-eic.jlab.org

# Browse available files
xrdfs root://dtn-eic.jlab.org ls /volatile/eic/romanov/meson-structure-2025-03/reco

# Example download
xrdcp root://dtn-eic.jlab.org//volatile/eic/romanov/meson-structure-2025-03/reco/k_lambda_5x41_5000evt_001.edm4eic.root ./
```

For more details on accessing the data, see the [Data](./data) page.