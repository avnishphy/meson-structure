# Simplified Particle Physics Analysis: A Practical Guide

This guide introduces a straightforward approach to analyzing particle physics data stored in ROOT files. We focus on simplicity and efficiency, showing how to process large datasets with minimal code.

## Prerequisites

Install these packages:

```bash
pip install uproot awkward numpy matplotlib hist
```

## The Tutorial Script

Our simplified script demonstrates:

1. Reading particle data from ROOT files
2. Processing data in manageable chunks
3. Creating histograms directly with the `hist` library
4. Creating visualizations with minimal code

## Main Concepts

### Global Histograms

We use global histogram variables for clarity and simplicity:

```python
# Global histograms
pz_hist = hist.Hist(hist.axis.Regular(100, -50, 50, name="momentum_z"))
lambda_decay_z = hist.Hist(hist.axis.Regular(100, 0, 40000, name="lambda_decay_z"))
lambda_pz = hist.Hist(hist.axis.Regular(100, -50, 50, name="lambda_momentum_z"))
proton_pz = hist.Hist(hist.axis.Regular(100, -50, 50, name="proton_momentum_z"))
```

This approach makes the code more readable and avoids passing histograms between functions.

### Efficient Chunk Processing

We process data in chunks to manage memory usage:

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

### Simple Plotting with hist

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

This approach is simpler than manually extracting bin edges and values.

## Step-by-Step Guide

### 1. Reading Files in Chunks

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

### 2. Processing Particles

We focus on processing all particles at once rather than event-by-event:

```python
# Flatten arrays for all particles in all events
flat_pdg = ak.flatten(pdg)
flat_pz = ak.flatten(pz)
flat_decay_z = ak.flatten(decay_z)

# Fill histograms for all particles at once
pz_hist.fill(flat_pz)
```

### 3. Creating Visualizations

Using `hist`'s built-in plotting capabilities makes visualization simple:

```python
fig, ax = plt.subplots(figsize=(10, 6))
lambda_pz.plot(ax=ax, color='green', alpha=0.7)
ax.set_xlabel("Momentum in z-direction [GeV/c]")
ax.set_ylabel("Count")
ax.set_title("Lambda Momentum (pz)")
ax.grid(True, alpha=0.3)
```

## Running the Script

Run the script with your ROOT files:

```bash
python simplified_particle_analysis.py file1.root file2.root --step-size 2000 -o output_plots -n 10000
```

## Tips for Efficiency

1. **Adjust chunk size**: Find the right balance between processing speed and memory usage.

2. **Use vectorized operations**: The awkward array library enables fast, NumPy-like operations on jagged arrays.

3. **Focus on what you need**: Only request the branches you actually need from the ROOT files.

4. **Use boolean masks**: They're much faster than loops for filtering particles.

## What's Next?

As you become comfortable with this basic analysis:

1. Try analyzing different particle properties
2. Create more complex histograms (2D, profile plots)
3. Implement more sophisticated event selection
4. Compare data with Monte Carlo simulations

## Conclusion

This simplified approach demonstrates that particle physics analysis doesn't need to be complicated. By focusing on vectorized operations and leveraging modern libraries like `hist` and `awkward`, you can analyze large datasets efficiently with minimal code.

Happy analyzing!