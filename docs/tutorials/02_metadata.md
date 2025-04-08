# Accessing Event Metadata 

... in EDM4eic Files

::: info
Corresponding python file with example:
- [02_metadata.py](https://github.com/JeffersonLab/meson-structure/tree/main/tutorials/02_metadata.py)
:::


## Understanding the Metadata Flow

Metadata in EDM4eic files contains important physics quantities from the event generation stage that propagate through the simulation chain. This includes true values like Q², Bjorken x, and other physics quantities crucial for analysis.

```mermaid
flowchart LR
    gen[Event Generator] -->|true values| hepmc[HepMC Converter]
    hepmc -->|attributes| after[Afterburner]
    after --> dd4hep[DD4Hep]
    dd4hep -->|event metadata| eicrecon[EICRecon]
    
    style gen fill:#f9f,stroke:#333,stroke-width:2px
    style eicrecon fill:#bbf,stroke:#333,stroke-width:2px
```

The diagram above illustrates how metadata travels through the simulation chain:
1. The event generator produces the true physics values
2. These values are converted to HepMC format with attributes
3. After passing through afterburner and DD4Hep simulation
4. The metadata is preserved in the final EICRecon output

## Metadata in EDM4eic Files

Metadata in EDM4eic files exists at both file and event levels:

- **File-level metadata**: Global information about the dataset
- **Event-level metadata**: Physics quantities for each event

The event-level metadata is stored in special branches of the 'events' tree:
- `GPStringKeys`: Contains the metadata field names
- `GPStringValues`: Contains the corresponding values as strings

## Tutorial: Accessing Event Metadata

### Prerequisites

```python
# Required packages
import uproot
import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
from hist import Hist  # scikit-hep/hist package
import os
```

### Reading Metadata from a File

```python
# Open the file and get the events tree
file_path = "your_edm4eic_file.root"
file = uproot.open(file_path)
events = file["events"]

# Read the metadata branches
metadata = events.arrays(["GPStringKeys", "GPStringValues"])

# The metadata is now in awkward arrays format
print(f"Number of events: {len(metadata['GPStringKeys'])}")
```

### Converting String Values to Numeric Data

```python
# Convert string values to float64 (numeric values)
keys_ak = metadata["GPStringKeys"]
values_ak_str = metadata["GPStringValues"]
values_ak = ak.strings_astype(values_ak_str, "float64")

# Get all unique metadata field names from the first event
metadata_fields = ak.to_list(ak.ravel(keys_ak[0]))
print(f"Available metadata fields: {metadata_fields}")
```

### Extracting Values for a Specific Metadata Field

```python
# Example: Extract 'dis_Q2' values for all events
q2_values = values_ak[keys_ak == "dis_Q2"]

# Flatten the array to get a simple list of values
q2_flat = ak.ravel(q2_values)

# Convert to numpy array for easier plotting/analysis
q2_numpy = ak.to_numpy(q2_flat)

print(f"Q² statistics: min={np.min(q2_numpy)}, max={np.max(q2_numpy)}, mean={np.mean(q2_numpy)}")
```

### Creating a Histogram for a Metadata Field

```python
# Create a simple histogram using matplotlib
plt.figure(figsize=(10, 6))
plt.hist(q2_numpy, bins=50)
plt.xlabel("Q² (GeV²)")
plt.ylabel("Counts")
plt.title("Distribution of Q² Values")
plt.grid(True, alpha=0.3)
plt.savefig("q2_distribution.png")
plt.close()
```

### Processing Multiple Files and Large Datasets

For large datasets or multiple files, it's more efficient to process data in chunks:

```python
def process_metadata_chunks(file_paths, field_name="dis_Q2", max_events=None):
    """Process metadata from multiple files in chunks."""
    values_list = []
    total_processed = 0
    
    # Create a dictionary of file paths to tree names
    files_dict = {path: "events" for path in file_paths}
    
    # Iterate through chunks
    for chunk in uproot.iterate(
        files=files_dict,
        expressions=["GPStringKeys", "GPStringValues"],
        step_size="100MB"
    ):
        # Count events in this chunk
        chunk_size = len(chunk["GPStringKeys"])
        
        # Extract values for the requested field
        keys = chunk["GPStringKeys"]
        values_str = chunk["GPStringValues"]
        values = ak.strings_astype(values_str, "float64")
        field_values = values[keys == field_name]
        
        # Add to our collection
        values_list.append(ak.ravel(field_values))
        
        # Update counter and check if we've reached max_events
        total_processed += chunk_size
        if max_events and total_processed >= max_events:
            break
    
    # Combine all chunks
    combined_values = ak.concatenate(values_list)
    return ak.to_numpy(combined_values)

# Example usage
file_paths = ["file1.root", "file2.root", "file3.root"]
q2_values = process_metadata_chunks(file_paths, field_name="dis_Q2", max_events=10000)
```

### Creating Histograms for All Metadata Fields

This example shows how to automatically create histograms for all metadata fields:

```python
def create_metadata_histograms(file_paths, output_dir="metadata_plots", max_events=None):
    """Create histograms for all metadata fields."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Dictionary to hold histograms for each metadata field
    histograms = {}
    
    # Process files in chunks
    files_dict = {path: "events" for path in file_paths}
    total_processed = 0
    
    for chunk in uproot.iterate(
        files=files_dict,
        expressions=["GPStringKeys", "GPStringValues"],
        step_size="100MB"
    ):
        # Extract keys and values
        keys = chunk["GPStringKeys"]
        values_str = chunk["GPStringValues"]
        values = ak.strings_astype(values_str, "float64")
        
        # Get all unique field names if we haven't seen them yet
        if not histograms:
            if len(keys) > 0:
                field_names = ak.to_list(ak.ravel(keys[0]))
                
                # Create histograms for each field
                for name in field_names:
                    field_values = values[keys == name]
                    if len(field_values) > 0:
                        min_val = ak.min(field_values)
                        max_val = ak.max(field_values)
                        # Add padding to range
                        padding = (max_val - min_val) * 0.1
                        min_val -= padding
                        max_val += padding
                        
                        # Create histogram using scikit-hep/hist
                        histograms[name] = Hist.new.Reg(
                            100, min_val, max_val, name=name, label=name
                        ).Double()
                        
                        print(f"Created histogram for {name}: range [{min_val}, {max_val}]")
        
        # Fill histograms with data from this chunk
        for name in histograms:
            field_values = values[keys == name]
            if len(field_values) > 0:
                histograms[name].fill(ak.ravel(field_values))
        
        # Update counter and check if we've reached max_events
        total_processed += len(keys)
        if max_events and total_processed >= max_events:
            break
    
    # Save histograms as images
    for name, hist in histograms.items():
        fig, ax = plt.subplots(figsize=(10, 6))
        hist.plot(ax=ax)
        ax.set_title(f"Distribution of {name}")
        ax.set_xlabel(f"{name} Value")
        ax.set_ylabel("Counts")
        ax.grid(True, alpha=0.3)
        
        output_path = os.path.join(output_dir, f"{name}.png")
        plt.savefig(output_path)
        plt.close(fig)
        
        print(f"Saved histogram for {name} to {output_path}")
    
    return histograms

# Example usage
file_paths = ["your_edm4eic_file.root"]
histograms = create_metadata_histograms(file_paths, output_dir="metadata_plots")
```

## Common Metadata Fields

Here are some common metadata fields you might find in EDM4eic files:

| Field Name | Description |
|------------|-------------|
| dis_Q2     | Four-momentum transfer squared (GeV²) |
| dis_xbj    | Bjorken x, the momentum fraction of the struck parton |
| dis_y      | Inelasticity, the fraction of energy transferred to the hadronic system |
| dis_W2     | Squared invariant mass of the hadronic system (GeV²) |
| dis_nu     | Energy transferred to the hadronic system (GeV) |
| projectile_energy | Energy of the projectile beam (GeV) |
| target_energy | Energy of the target beam (GeV) |

## Advanced Usage: Correlations Between Metadata Fields

You can also study correlations between different metadata fields:

```python
def plot_metadata_correlation(file_paths, field_x="dis_xbj", field_y="dis_Q2", max_events=None):
    """Create a 2D correlation plot between two metadata fields."""
    x_values = []
    y_values = []
    
    # Process files in chunks
    files_dict = {path: "events" for path in file_paths}
    total_processed = 0
    
    for chunk in uproot.iterate(
        files=files_dict,
        expressions=["GPStringKeys", "GPStringValues"],
        step_size="100MB"
    ):
        # Extract keys and values
        keys = chunk["GPStringKeys"]
        values_str = chunk["GPStringValues"]
        values = ak.strings_astype(values_str, "float64")
        
        # Extract values for each field
        x_field_values = values[keys == field_x]
        y_field_values = values[keys == field_y]
        
        # Flatten and append to our lists
        x_values.append(ak.ravel(x_field_values))
        y_values.append(ak.ravel(y_field_values))
        
        # Update counter and check if we've reached max_events
        total_processed += len(keys)
        if max_events and total_processed >= max_events:
            break
    
    # Combine all chunks
    x_combined = ak.to_numpy(ak.concatenate(x_values))
    y_combined = ak.to_numpy(ak.concatenate(y_values))
    
    # Create the correlation plot
    plt.figure(figsize=(10, 8))
    plt.hist2d(x_combined, y_combined, bins=50, cmap='viridis')
    plt.colorbar(label='Counts')
    plt.xlabel(field_x)
    plt.ylabel(field_y)
    plt.title(f"Correlation: {field_x} vs {field_y}")
    plt.grid(True, alpha=0.3)
    
    output_path = f"correlation_{field_x}_vs_{field_y}.png"
    plt.savefig(output_path)
    plt.close()
    
    print(f"Saved correlation plot to {output_path}")
    
    return x_combined, y_combined

# Example: Plot correlation between Bjorken x and Q²
x_values, q2_values = plot_metadata_correlation(
    ["your_edm4eic_file.root"], 
    field_x="dis_xbj", 
    field_y="dis_Q2"
)
```

## Conclusion

Accessing metadata in EDM4eic files allows you to retrieve important physics quantities from the event generation stage. This can be valuable for:

- Understanding the characteristics of your dataset
- Performing truth-level analyses
- Evaluating reconstruction performance by comparing with reconstructed quantities

The metadata values provide the "ground truth" for your analysis, allowing you to better understand and interpret your results.