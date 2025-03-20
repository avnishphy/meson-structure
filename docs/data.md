# Data

This page contains information about accessing and working with the meson structure data.

## Data Location

The meson structure data is available from the following locations:

- Original MCEG files: `/work/eic/users/singhav/JLEIC/USERS/trottar/OUTPUTS/raty_eic/` (last update of October 2024)
- Processed reconstruction files: `/volatile/eic/romanov/meson-structure-2025-03/reco`

## Accessing Data with XRootD

The data is available through XRootD for remote access. Here's how to access it:

### Exploring Available Files

To browse the available files:

```bash
xrdfs root://dtn-eic.jlab.org
ls /volatile/eic/romanov/meson-structure-2025-03/reco
```

### Downloading Files

To download files using XRootD:

```bash
xrdcp root://dtn-eic.jlab.org//volatile/eic/romanov/meson-structure-2025-03/reco/k_lambda_5x41_5000evt_200.edm4eic.root ./
```

## File Naming Convention

Files follow this naming pattern:

```
k_lambda_{beam_E}_5000evt_{index}.edm4eic.root
```

Where:
- `beam_E`: Beam energy configuration [5x41, 10x100, 18x275]
- `index`: File index [001-200]

Example filenames:
- `k_lambda_5x41_5000evt_001.edm4eic.root`
- `k_lambda_10x100_5000evt_001.edm4eic.root`
- `k_lambda_18x275_5000evt_001.edm4eic.root`

## File Types

The reconstruction produces several file types:
- `.edm4eic.root`: Main reconstruction output with full event data
- `.afterburner.hist.root`: Histograms from afterburner analysis

## Example Usage

Here's an example of downloading files using XRootD:

```bash
# Download a small file (~12MB)
xrdcp root://dtn-eic.jlab.org//volatile/eic/romanov/meson-structure-2025-03/reco/k_lambda_5x41_5000evt_200.afterburner.hist.root ./

# Download a medium file (~357MB)
xrdcp root://dtn-eic.jlab.org//volatile/eic/romanov/meson-structure-2025-03/reco/k_lambda_5x41_5000evt_200.edm4eic.root ./

# Download a larger file (~626MB)
xrdcp root://dtn-eic.jlab.org//volatile/eic/romanov/meson-structure-2025-03/reco/k_lambda_10x100_5000evt_001.edm4eic.root ./

# Download a very large file (~1.1GB)
xrdcp root://dtn-eic.jlab.org//volatile/eic/romanov/meson-structure-2025-03/reco/k_lambda_18x275_5000evt_001.edm4eic.root ./
```