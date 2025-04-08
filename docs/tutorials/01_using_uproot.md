# Using uproot

Corresponding python example: 
- [tutorials/00_uproot.py](https://github.com/JeffersonLab/meson-structure/tree/main/tutorials/00_uproot.py)
- [tutorials/01_plot_mcparticles.py](https://github.com/JeffersonLab/meson-structure/tree/main/tutorials/01_plot_mcparticles.py)


## Prerequisites

Install these packages:

```bash
pip install uproot awkward numpy matplotlib hist
```

::: info
We use [uproot] to process reconstruction files saved in CERN ROOT format `.root`.
While python is a slow language if compared to C++, uproot can achieve comparable 
performance in event processing. It is based on [awkward-arrays][awkward] arrays which core
is written in C and utilizes vectorized data processing the same way as numpy. What is even more
important, uproot can be easily installed via `pip install`, runs on all operating systems, 
very compatible with main python data science and AI tools. 
:::

Related documentation: 

- [uproot]
- [awkward]


[uproot-github]: https://github.com/scikit-hep/uproot5
[uproot]: https://uproot.readthedocs.io/en/latest/basic.html
[awkward]: https://github.com/scikit-hep/awkward


## Plot MCParticles


Here is the description of 1st tutorial:

1. Reading particle data from ROOT files
2. Processing data in manageable chunks
3. Creating histograms directly with the `hist` library
4. Creating visualizations with minimal code


### Reading root file

`uproot` provides several ways to open and read data from files.
Uproot tutorials starts with `array` or `arrays` method which reads all required data at once. 
But it could easily take too much time and memory on e.g. a laptop if files are large. 

The most efficient way to develop analysis scripts and process large number of files/events is
to use `iterate` method which reads data in chunks, which could be
processed in a vectorized way (using numpy or better suited awkward array library)
So we e.g. read 1000 events at once, process them, add data to histos, process
next 1000 events, etc. 


[uproot iterate method](https://uproot.readthedocs.io/en/latest/uproot.behaviors.TBranch.iterate.html):

The minimal you need to iterate

```python
# The simplest way to process a file
for chunk in uproot.iterate(
        {file_name: "events"},      # File Name : TTree name (EDM4EIC ttree is "events")
        branches,                   # List of branches you want to process
        step_size=1000,             # How many events to process per chunk
        entry_stop=10_000           # On what event to stop (there is also etry_start) variable
    ):
    # process chunk by chunk here
```

See [tutorials/00_uproot.py](https://github.com/JeffersonLab/meson-structure/tree/main/tutorials/00_uproot.py)

Full code stub: 

```python
import uproot

# What branches to process
branches = [
    "MCParticles.PDG",
    "MCParticles.momentum.z",
    "MCParticles.endpoint.z",
]

# Read and process file in chunks
for chunk_i, chunk in enumerate(uproot.iterate(
        {"my_file.root": "events"}, 
        branches, 
        step_size=100, 
        entry_stop=200)):

    print(f"Сhunk {chunk_i} read")

    # Print data shape. It is going to be
    # [n-events-in-chunk]x{branch:[n-particles]}
    chunk.type.show()

    # Show a value of a single particle
    particle_pdg = chunk[0]["MCParticles.PDG"][2]
    print(f"  PDG of the 3d particle of the 1st event in this chunk: {particle_pdg}")
```

How data is organized:

```
chunk.type.show() output is: 

100 * {
"MCParticles.PDG": var * int32,
"MCParticles.momentum.z": var * float32,
"MCParticles.endpoint.z": var * float64
}
```

This means that data can be accessed as:

```python
chunk[event_index][branch_name][particle_index]

# Example: 3d particle of the 1st event: 
particle_pdg = chunk[0]["MCParticles.PDG"][2]
```

What is more important, that one can use `chunk[branch_name]` to get `[events]x[particle data]` 
awkward array that can be processed in vectorized way: 

```python
events_pdgs = chunk["MCParticles.PDG"]

# [event_0, event_1, ...] where event_0=[particle_0_pdg, particle_1_pdg, ... etc]
# e.g. [[2212, 11, 11, 321, 3122, 22, 22, 22], ..., [2212, 11, 11, ..., 11, 11, 11]]

# Vectorized way of processing the data
lambda_filter = chunk["MCParticles.PDG"] == 3122
lam_pz = ak.mean(chunk[lambda_filter]["MCParticles.momentum.z"])
print(f"  Lambdas pz in this chunk: {lam_pz}")
```

The key features are:
- Using uproot iterate to process data in chunks
- Using awkward arrays for vectorized operations
- Using boolean masks to select specific particle types

### Histograms

Look [tutorials/01_plot_mcparticles.py](https://github.com/JeffersonLab/meson-structure/tree/main/tutorials/01_plot_mcparticles.py)
for full details. 

We use [hist](https://hist.readthedocs.io/en/latest/user-guide/notebooks/Plots.html) library, 
which uses [boost-histogram](https://boost-histogram.readthedocs.io/en/latest/index.html)
under the hood and provides familiar for HENP way to create and fill histograms when iterating the 
file:

```python
# Create histograms
pz_hist = hist.Hist(hist.axis.Regular(100, -50, 50, name="momentum_z"))

# Fill histograms
pz_hist.fill(branch_pz[lambda_mask])
```

Hist package can plot its own histograms in jupyter. Here is an example how to save histogrms to file:

```python
import matplotlib.pyplot as plt

# Figure and Axes in matplotlib are analog of Canvas and Pad in ROOT
fig, ax = plt.subplots(figsize=(10, 6))
pz_hist.plot(ax=ax)
fig.savefig("lambda_pz.png")
```

## Tips for Efficiency

1. **Adjust chunk size**: Find the right balance between processing speed and memory usage.

2. **Use vectorized operations**: The awkward array library enables fast, NumPy-like operations on jagged arrays.

3. **Focus on what you need**: Only request the branches you actually need from the ROOT files.

4. **Use boolean masks**: They're much faster than loops for filtering particles.
