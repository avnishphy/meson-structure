# Navigating *all* EDM4eic relations with **uproot**, **awkward‑array**, and **Rich**

> **Target audience** – analysts who already skim ROOT files produced by ePIC/eAST/Geant4 and now need to *follow* every link: **MCParticle ⇄ MCParticle**, **Track ⇄ Hit**, **RecoParticle ⇄ MCParticle**, …

---

## 0  Why does everything look indirect?

EDM4eic follows the \[PODIO] flattened‑table philosophy:

1. **Each collection** (e.g. `MCParticles`, `TrackerHits`) appears as a branch of fixed‑length record arrays – one row per object.
2. **Each *relation*** (parents, daughters, hits, clusters, …) is *not* stored inline.
   Instead we get

   * two *offset* arrays per object: `relation_begin` and `relation_end`;
   * **one** supplementary flattened vector – a separate branch whose name starts with an underscore, e.g. `_MCParticles_daughters`.
3. The object’s relations therefore live in the half‑open slice

```text
relation_indices = _SUPPLEMENTAL[b:e]  # where b = relation_begin[i], e = relation_end[i]
```

Exactly the same for **all** links in EDM4eic:

| Example relation                 | Offset fields on object            | Supplemental branch                   |
| -------------------------------- | ---------------------------------- | ------------------------------------- |
| MCParticle → *parents*           | `parents_begin`, `parents_end`     | `_MCParticles_parents.*`              |
| MCParticle → *daughters*         | `daughters_begin`, `daughters_end` | `_MCParticles_daughters.*`            |
| Track → *TrackerHits*            | `hits_begin`, `hits_end`           | `_Track_hits.*`                       |
| RecoParticle → *MCParticle* link | `particles_begin`, `particles_end` | `_ReconstructedParticles_particles.*` |

Once you grok one case, you can traverse them all.

---

## 1  Prerequisites

```bash
pip install uproot awkward rich
```

We assume a single‑file EDM4hep/EDM4eic output containing the `events` TTree.

---

## 2  Annotated walk‑through code

Below is a **minimal but complete** script.  Each line is commented; variable names favour clarity over brevity.

```python
#!/usr/bin/env python3
"""Explore EDM4eic relation tables – print an event‑by‑event particle tree.

* Shows how begin/end offsets map into the supplemental *_MCParticles_daughters
  and *_MCParticles_parents vectors.
* The very same pattern works for *any* relation in EDM4eic (tracks↔hits, etc.).
"""

import argparse
import uproot
import awkward as ak
from rich.table import Table
from rich.console import Console

# ---------------------------------------------------------------------------
# 1. CLI – choose file, tree, and how many events to print
# ---------------------------------------------------------------------------

def parse_cli() -> argparse.Namespace:
    """Return parsed command‑line arguments."""

    cli = argparse.ArgumentParser(
        description="Print PDG, pz, decay‑z and parent/daughter indices for the first N events"
    )
    cli.add_argument("-i", "--input-file", required=True,
                     help="Path to an EDM4hep ROOT file")
    cli.add_argument("-t", "--tree-name", default="events",
                     help="Name of the TTree (default: events)")
    cli.add_argument("--max-events", type=int, default=50,
                     help="How many events to dump (default: 50)")
    return cli.parse_args()

# ---------------------------------------------------------------------------
# 2. Helper: slice supplemental table given begin/end offsets
# ---------------------------------------------------------------------------

def indices_for(obj_begin: int, obj_end: int, flat_vector: ak.Array) -> list[int]:
    """Return a Python list of indices (or empty list) for a single object."""
    if obj_begin >= 0 and obj_end > obj_begin:
        return list(flat_vector[obj_begin:obj_end])
    return []  # no relation recorded

# ---------------------------------------------------------------------------
# 3. Main analysis driver
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_cli()
    print(f"Opening {args.input_file}")

    with uproot.open(args.input_file) as root_file:
        tree = root_file[args.tree_name]

        # --- 3.1 pick only the branches we need – uproot is lazy! ---------
        branches = [
            # core kinematics
            "MCParticles.PDG",                # int
            "MCParticles.momentum.z",         # float
            "MCParticles.endpoint.z",         # float (decay‑vertex z)
            # per‑object relation *offsets*
            "MCParticles.daughters_begin",    # int[nevt][npart]
            "MCParticles.daughters_end",
            "MCParticles.parents_begin",
            "MCParticles.parents_end",
            # flattened relation tables – we only need the .index field
            "_MCParticles_daughters.index",   # int[n_total_relations]
            "_MCParticles_parents.index",
        ]

        arrays = tree.arrays(branches, entry_stop=args.max_events, library="ak")

    # --- 3.2  unpack into convenience variables --------------------------------
    pdg_per_event          = arrays["MCParticles.PDG"]
    pz_per_event           = arrays["MCParticles.momentum.z"]
    decayZ_per_event       = arrays["MCParticles.endpoint.z"]
    dobeg_per_event        = arrays["MCParticles.daughters_begin"]
    doend_per_event        = arrays["MCParticles.daughters_end"]
    pabeg_per_event        = arrays["MCParticles.parents_begin"]
    paend_per_event        = arrays["MCParticles.parents_end"]
    flat_daughter_indices  = arrays["_MCParticles_daughters.index"]
    flat_parent_indices    = arrays["_MCParticles_parents.index"]

    console = Console()

    # --- 3.3  iterate over events ---------------------------------------------
    for ievt in range(len(pdg_per_event)):
        n_particles = len(pdg_per_event[ievt])
        table = Table(title=f"Event {ievt}", expand=True)
        for col in ("Idx", "PDG", "pz", "endZ", "Parents", "Daughters"):
            table.add_column(col, no_wrap=True)

        # loop over particles in this event
        for ip in range(n_particles):
            parents   = indices_for(pabeg_per_event[ievt][ip],
                                    paend_per_event[ievt][ip],
                                    flat_parent_indices[ievt])
            daughters = indices_for(dobeg_per_event[ievt][ip],
                                    doend_per_event[ievt][ip],
                                    flat_daughter_indices[ievt])

            table.add_row(
                str(ip),
                str(pdg_per_event[ievt][ip]),
                f"{pz_per_event[ievt][ip]:.2f}",
                f"{decayZ_per_event[ievt][ip]:.2f}",
                ", ".join(map(str, parents)) or "–",
                ", ".join(map(str, daughters)) or "–",
            )

        console.print(table)

        # Optional: dump the *raw* flattened tables so users can see
        # how begin/end offsets map into rows --------------------------------
        def dump_flat(title: str, flat: ak.Array):
            t = Table(title=title, expand=True)
            t.add_column("Flat row", no_wrap=True)
            t.add_column("Object index", no_wrap=True)
            for i, idx in enumerate(flat):
                t.add_row(str(i), str(idx))
            console.print(t)

        dump_flat(f"Event {ievt} – daughters flat table", flat_daughter_indices[ievt])
        dump_flat(f"Event {ievt} – parents   flat table", flat_parent_indices[ievt])
        console.print()  # blank line between events

if __name__ == "__main__":
    main()
```

---

## 3  Key take‑aways

* **begin/end are *offsets*, not indices** – they tell you where to slice the *supplemental* flattened array.
* The supplemental array holds `podio::ObjectID`s; inside one file the `collectionID` is usually constant, so we often only need `.index`.
* **Every relation in EDM4eic uses the *same* three‑vector pattern**.  Once you implement `indices_for()` you can navigate any link graph:

  * `Track.hits_begin/end` ↔ `_Track_hits.*`  → slice to get hit indices.
  * `ReconstructedParticles.particles_begin/end` ↔ `_ReconstructedParticles_particles.*` to match reco objects back to MC truth.

---

## 4  Next steps

1. **Generalise** the helper so it works for arbitrary collection names – great for small validation notebooks.
2. **Visualise** whole event graphs (e.g. with `networkx`) by feeding edges from `indices_for()`.
3. **Performance tip** – if you only need a few events, add `entry_start` / `entry_stop` in uproot to avoid loading the full tree.

Happy tracing! 🙌
