# Using uproot

::: info
Corresponding python file with example: 
- [01_plot_mcparticles.py](https://github.com/JeffersonLab/meson-structure/tree/main/tutorials/01_plot_mcparticles.py)
:::

## Prerequisites

Install these packages:

```bash
pip install uproot awkward numpy matplotlib hist
```

## 01 Plot MCParticles


Here is the description of 1st tutorial:

1. Reading particle data from ROOT files
2. Processing data in manageable chunks
3. Creating histograms directly with the `hist` library
4. Creating visualizations with minimal code

### Histograms


We use global histogram variables:

```python
# Global histograms
pz_hist = hist.Hist(hist.axis.Regular(100, -50, 50, name="momentum_z"))
lambda_decay_z = hist.Hist(hist.axis.Regular(100, 0, 40000, name="lambda_decay_z"))
lambda_pz = hist.Hist(hist.axis.Regular(100, -50, 50, name="lambda_momentum_z"))
proton_pz = hist.Hist(hist.axis.Regular(100, -50, 50, name="proton_momentum_z"))
```

### Efficient Chunk Processing

The most efficient way to process large number of files/events is
to use `iterate` method which reads data in chunks, which could be
processed in a vectorized way (using numpy or better suited awkward array library)
So we e.g. read 1000 events at once, process them, add data to histos, process
next 1000 events, etc.

First we create a function that processes such chunks:

```python
def process_chunk(chunk):
    # Extract arrays from the chunk
    pdg = chunk["MCParticles.PDG"]
    pz = chunk["MCParticles.momentum.z"]
    decay_z = chunk["MCParticles.endpoint.z"]
    
    # Flatten arrays for all particles
    flat_pdg = ak.flatten(pdg)
    flat_pz = ak.flatten(pz)
    flat_decay_z = ak.flatten(decay_z)
    
    # Fill histograms directly
    pz_hist.fill(flat_pz)
    
    # Use masks to select particle types
    lambda_mask = (flat_pdg == PDG_LAMBDA)
    if ak.sum(lambda_mask) > 0:
        lambda_pz.fill(flat_pz[lambda_mask])
        lambda_decay_z.fill(flat_decay_z[lambda_mask])
```

The key features are:
- Using awkward arrays for vectorized operations
- Flattening the arrays to process all particles at once
- Using boolean masks to select specific particle types

### Plotting with hist

`hist` package produce pretty histograms out of the box, but
we can enhance how they look configuring underlying figure and ax-es.

The `hist` library provides built-in plotting functionality:

```python
def create_plots(outdir):
    # Plot momentum in z-direction
    fig, ax = plt.subplots(figsize=(10, 6))
    pz_hist.plot(ax=ax, color='blue', alpha=0.7)
    ax.set_xlabel("Momentum in z-direction [GeV/c]")
    ax.set_ylabel("Count")
    ax.set_title("Particle Momentum (pz)")
    ax.grid(True, alpha=0.3)
```


### Reading Files in Chunks

`uproot` has several methods reading file.
`array` and `arrays`, read whole data from file, which might be fine in some
cases but takes too much time and memory in others.

We use `uproot.iterate` to read ROOT files in manageable chunks:

```python
for chunk in uproot.iterate(
        file_dict,
        expressions=[
            "MCParticles.PDG",
            "MCParticles.momentum.z",
            "MCParticles.endpoint.z",
        ],
        step_size=step_size,
        library="ak"):
    
    process_chunk(chunk)
```

### Processing Particles

We focus on processing all particles in a chunk at once rather than event-by-event:

In Uproot, each branch in a TTree corresponds to a column, and “events” are rows.
Awkward leverages this to provide one list‐of‐events dimension on top,
with potentially variable‐length sublists for each event (the “jagged” part).

So if each event has a different number of MCParticles,
you still read them as a single Awkward array—each event’s sublist is just a different length.

For example, if you have
MCParticles.momentum.x, MCParticles.momentum.y, MCParticles.momentum.z,
and MCParticles.pdg stored in Uproot,
you end up with arrays of shape [n_events, counts_per_event].
One event may have 5 particles, the next 12, etc.,
and each coordinate or PDG code can be accessed with the same event alignment: momentum.x[i_event][i_particle].

```python
# Flatten arrays for all particles in all events
flat_pdg = ak.flatten(pdg)
flat_pz = ak.flatten(pz)
flat_decay_z = ak.flatten(decay_z)

# Fill histograms for all particles at once
pz_hist.fill(flat_pz)
```


## Running the Script

Run the script with your ROOT files:

```bash
python 01_plot_mcparticles.py file1.root file2.root --step-size 2000 -o output_plots -n 10000
```

## Tips for Efficiency

1. **Adjust chunk size**: Find the right balance between processing speed and memory usage.

2. **Use vectorized operations**: The awkward array library enables fast, NumPy-like operations on jagged arrays.

3. **Focus on what you need**: Only request the branches you actually need from the ROOT files.

4. **Use boolean masks**: They're much faster than loops for filtering particles.
