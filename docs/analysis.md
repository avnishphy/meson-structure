# Analysis


# Understanding Daughter Indices in EDM4hep

We will use MCParticles as example

## 1. Overview of `MCParticles` Data

In EDM4hep, each event has an array of `MCParticles`. Each particle typically stores:

- **PDG code**
- **4-momentum** (momentum.x, y, z, energy, etc.)
- **Endpoint** coordinates (x, y, z) – indicating where it decayed, absorbed, or exiting the world
- **`parents_begin` / `parents_end`** – range of references to the `_MCParticles_parents` array
- **`daughters_begin` / `daughters_end`** – range of references to the `_MCParticles_daughters` array

The **`_MCParticles_{daughters, parents}.index`** arrays store the *actual* integer indexes of the
parent or daughter particles in the `MCParticles` collection. For example, if a particle has
`daughters_begin = 10` and `daughters_end = 12`, that means in the `_MCParticles_daughters.index`
array **from** index 10 **up to** (but **not** including) 12, you’ll find the references (e.g.
`(0, 22)`) telling you the daughter is at index `22` in the main `MCParticles` array.

## 2. Why Daughter Indices Might Be `None` (Empty)

It’s possible to see a particle with **non-zero endpoint** \(z\) (e.g. `endpoint.z = 8700 mm`) but *
*no daughters** in the final `MCParticles` collection. Several common reasons:

1. **Particle Decays Outside the Instrumented Region**  
   Geant4 (or dd4hep) can track a particle until it leaves the “world volume.” If it decays *after*
   crossing the boundary, or if the simulation is configured to discard secondaries outside a region
   of interest, the daughters never appear in the final data.
    - The parent’s endpoint might reflect where it left the geometry (thus “endpoint.z”), but no
      actual decaying is recorded in `MCParticles`.

2. **Particle Treated as “Stable”**  
   Some simulation setups treat \(\Lambda\) or other hyperons as stable. They set an endpoint if it
   leaves the volume, but no official decay is recorded.

3. **Filtering / Cuts Remove the Daughters**  
   Many dd4hep or Geant4 jobs apply logic like “don’t store secondaries that never hit a sensitive
   volume” (tracker, calorimeter). If the \(\Lambda\) decays in empty space and the daughters are
   never registered in a sub-detector, they can be omitted.

4. **Post-Processing or Conversion**  
   Sometimes the mother-daughter link is lost during intermediate steps (e.g. digitization,
   reconstruction, or conversion to EDM4hep). The data ends up with a \(\Lambda\) that has an
   endpoint but no recognized offspring in the final record.

## 3. dd4hep / ddsim Flags to Preserve All Particles

To ensure you keep *as many decays as possible* in the final EDM4hep output, consider using the
following flags or configurations in **dd4hep** or **ddsim**:

1. **`--keepAllParticles`**  
   Forces the simulator to save *all* Geant4-tracked particles, even if they don’t produce hits.

2. **Disable Filters**
    - `--filter.tracker false`
    - `--filter.calo false`  
      Disables built-in filters that discard particles not depositing energy in the tracker or
      calorimeter.

3. **`--enableDetailedShowerMode`**  
   In some iLCSoft-based versions, this turns on more verbose saving of secondaries.

4. **Steering File**  
   Alternatively, if you use a Python steering file for dd4hep, you might set:
   ```python
   from DDSim.DDSimConfiguration import DDSimConfig
   config = DDSimConfig()
   config.keepAllParticles = True
   config.filterCalo = False
   config.filterTracker = False
   config.enableDetailedShowerMode = True
   # ...
   ```

**Example** command line:

```bash
ddsim \
  --compactFile MyDetector.xml \
  --inputFile MyGen.hepmc \
  --outputFile MySim.edm4hep.root \
  --keepAllParticles \
  --filter.tracker false \
  --filter.calo false \
  --enableDetailedShowerMode
```

> **Note**: Not all flags exist in every dd4hep release. Check `ddsim --help` or your local docs.

## 4. Verifying the Configuration

- **Print the MCParticles**: Use a debug script (e.g.
  with [Rich tables](https://github.com/Textualize/rich)) to see each particle’s PDG, endpoint, and
  daughters.
- **Check for $\(\Lambda\) (PDG=3122) or \(\bar{\Lambda}\)$ (PDG=-3122)**. Are they truly decaying
  into `p (2212) + \pi^- (-211)` or `\bar{p} (-2212) + \pi^+ (211)`?
- **Look for mother-daughter slices** in `_MCParticles_daughters.index` and
  `_MCParticles_parents.index`. If they are empty, the simulation never stored those secondaries.

## 5. Summary

Seeing **a $\(\Lambda\)$ with “no daughters”** is typically not a bug in the code, but a reflection of
**how** Geant4 and dd4hep are configured to store or discard certain secondaries. By enabling the
right flags (or removing filters), you can preserve the full decay chain in the `MCParticles`
collection—even if it occurs far downstream or in regions without sensitive detectors.