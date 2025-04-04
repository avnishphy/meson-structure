import argparse
import uproot
import awkward as ak
from rich.table import Table
from rich.console import Console


def parse_args():
    parser = argparse.ArgumentParser(description="Minimal narrow debug: Prints PDG, pz, endpoint.z, and partial parent/daughter indices.")
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
                "MCParticles.momentum.z",   # only pz
                "MCParticles.endpoint.z",   # decay z
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

    pdg_all    = arrays["MCParticles.PDG"]
    pz_all     = arrays["MCParticles.momentum.z"]
    endz_all   = arrays["MCParticles.endpoint.z"]

    d_beg_all  = arrays["MCParticles.daughters_begin"]
    d_end_all  = arrays["MCParticles.daughters_end"]
    p_beg_all  = arrays["MCParticles.parents_begin"]
    p_end_all  = arrays["MCParticles.parents_end"]

    daughters_idx_all = arrays["_MCParticles_daughters.index"]
    parents_idx_all   = arrays["_MCParticles_parents.index"]

    n_events = len(pdg_all)
    console = Console()

    # Loop over each event
    for i_evt in range(n_events):
        evt_pdg  = pdg_all[i_evt]
        evt_pz   = pz_all[i_evt]
        evt_endz = endz_all[i_evt]

        evt_db   = d_beg_all[i_evt]
        evt_de   = d_end_all[i_evt]
        evt_pb   = p_beg_all[i_evt]
        evt_pe   = p_end_all[i_evt]

        evt_d_idx= daughters_idx_all[i_evt]
        evt_p_idx= parents_idx_all[i_evt]

        n_parts = len(evt_pdg)

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
            dbeg = evt_db[i_part]
            dend = evt_de[i_part]
            pbeg = evt_pb[i_part]
            pend = evt_pe[i_part]

            # Build a short string for the parent's indexes
            if (pbeg >= 0) and (pend > pbeg):
                these_par_idxs = evt_p_idx[pbeg:pend]
                # show up to 2
                show_par = list(these_par_idxs[:2])
                par_str = ", ".join(str(x) for x in show_par)
                if len(these_par_idxs) > 2:
                    par_str += ", ..."
            else:
                par_str = "None"

            # Build a short string for the daughter's indexes
            if (dbeg >= 0) and (dend > dbeg):
                these_dau_idxs = evt_d_idx[dbeg:dend]
                # show up to 2
                show_dau = list(these_dau_idxs[:2])
                dau_str = ", ".join(str(x) for x in show_dau)
                if len(these_dau_idxs) > 2:
                    dau_str += ", ..."
            else:
                dau_str = "None"

            table.add_row(
                str(i_part),
                str(evt_pdg[i_part]),
                f"{evt_pz[i_part]:.2f}",
                f"{evt_endz[i_part]:.2f}",
                par_str,
                dau_str
            )

        console.print(table)
        console.print("\n")  # blank line for spacing

if __name__ == "__main__":
    main()
