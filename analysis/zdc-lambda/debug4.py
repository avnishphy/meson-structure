import argparse
import uproot
import awkward as ak
from rich.table import Table
from rich.console import Console


def parse_args():
    parser = argparse.ArgumentParser( description="Minimal narrow debug: Prints PDG, pz, endpoint.z, and partial parent/daughter indices.")
    parser.add_argument("-i", "--input-file", required=True, help="Path to a single EDM4hep ROOT file with MCParticles.")
    parser.add_argument("-t", "--tree-name", default="events", help="Name of the TTree (default 'events').")
    parser.add_argument("--max-events", type=int, default=50, help="Max number of events to read/print (default=50).")
    return parser.parse_args()


def main():
    args = parse_args()

    print(f"Opening file: {args.input_file}")
    with uproot.open(args.input_file) as f:
        tree = f[args.tree_name]

        # Read up to 'max_events' events worth of data
        arrays = tree.arrays(
            expressions=[
                # Basic MC info
                "MCParticles.PDG",
                "MCParticles.momentum.z",  # only pz
                "MCParticles.endpoint.z",  # decay z
                # Daughter indices
                "MCParticles.daughters_begin",
                "MCParticles.daughters_end",
                # Parent indices
                "MCParticles.parents_begin",
                "MCParticles.parents_end",
                # The actual Podio::ObjectID vectors for daughters & parents
                # We skip the .collectionID, focusing only on .index for narrowness
                "_MCParticles_daughters.index",
                "_MCParticles_parents.index",
            ],
            entry_stop=args.max_events,
            library="ak"
        )

    all_pdg = arrays["MCParticles.PDG"]
    all_pz = arrays["MCParticles.momentum.z"]
    all_endz = arrays["MCParticles.endpoint.z"]

    all_daughters_begin = arrays["MCParticles.daughters_begin"]
    all_daughters_end = arrays["MCParticles.daughters_end"]
    all_parents_begin = arrays["MCParticles.parents_begin"]
    all_parents_end = arrays["MCParticles.parents_end"]

    all_daughters_idx = arrays["_MCParticles_daughters.index"]
    all_parents_idx = arrays["_MCParticles_parents.index"]

    n_events = len(all_pdg)
    console = Console()

    # Loop over each event
    for i_evt in range(n_events):
        event_pdg = all_pdg[i_evt]
        event_pz = all_pz[i_evt]
        event_endz = all_endz[i_evt]

        event_daughters_begin = all_daughters_begin[i_evt]
        event_daughters_end = all_daughters_end[i_evt]
        event_parents_begin = all_parents_begin[i_evt]
        event_parents_end = all_parents_end[i_evt]

        event_daughters_idx = all_daughters_idx[i_evt]
        event_parents_idx = all_parents_idx[i_evt]

        # trying uproot to navigate
        n_parts = len(event_pdg)

        # Create a narrow Rich table
        table = Table(title=f"Event {i_evt}")
        table.expand = True  # allow full width
        table.add_column("Idx", no_wrap=True)
        table.add_column("PDG", no_wrap=True)
        table.add_column("pz", no_wrap=True)
        table.add_column("endZ", no_wrap=True)
        table.add_column("Parents (idx)", no_wrap=True)
        table.add_column("Daughters (idx)", no_wrap=True)

        for i_part in range(n_parts):
            dbeg = event_daughters_begin[i_part]
            dend = event_daughters_end[i_part]
            pbeg = event_parents_begin[i_part]
            pend = event_parents_end[i_part]

            # Build a short string for the parent's indexes
            if (pbeg >= 0) and (pend > pbeg):
                these_par_idxs = event_parents_idx[pbeg:pend]
                par_str = ", ".join(str(x) for x in these_par_idxs)

            else:
                par_str = "None"

            # Build a short string for the daughter's indexes
            if (dbeg >= 0) and (dend > dbeg):
                these_dau_idxs = event_daughters_idx[dbeg:dend]
                dau_str = ", ".join(str(x) for x in these_dau_idxs)

            else:
                dau_str = "None"

            table.add_row(
                str(i_part),
                str(event_pdg[i_part]),
                f"{event_pz[i_part]:.2f}",
                f"{event_endz[i_part]:.2f}",
                par_str,
                dau_str
            )

        console.print(table)

        # Add a table for daughters_idx_all
        daughters_table = Table(title=f"Event {i_evt} - daughters_idx_all")
        daughters_table.expand = True
        daughters_table.add_column("Row Index", no_wrap=True)
        daughters_table.add_column("Daughter Index", no_wrap=True)

        for i_idx, daughter_idx in enumerate(event_daughters_idx):
            daughters_table.add_row(
                str(i_idx),
                str(daughter_idx)
            )

        console.print(daughters_table)

        # Add a table for parents_idx_all
        parents_table = Table(title=f"Event {i_evt} - parents_idx_all")
        parents_table.expand = True
        parents_table.add_column("Row Index", no_wrap=True)
        parents_table.add_column("Parent Index", no_wrap=True)

        for i_idx, parent_idx in enumerate(event_parents_idx):
            parents_table.add_row(
                str(i_idx),
                str(parent_idx)
            )

        console.print(parents_table)
        console.print("\n")  # blank line for spacing


if __name__ == "__main__":
    main()