#!/usr/bin/env python3

import argparse
import uproot
import awkward as ak

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-file", "-i", required=True, help="Path to ROOT file.")
    parser.add_argument("--chunk-size", type=int, default=100_000, help="Size of each chunk.")
    parser.add_argument("--max-events", type=int, default=None, help="Stop after N events total.")
    return parser.parse_args()

def main():
    args = parse_args()

    # We’ll read the entire `invts` struct as one branch:
    # plus the TLorentzVectors: P_Inc, e_Inc, e_Scat, etc.
    # Notice we do NOT use "invts.Q2" here.
    branches = [
        "invts",        # the whole struct
        "P_Inc.",        # a TLorentzVector
        "e_Inc.",
        "e_Scat.",
        "k.",
        "lamb_scat."
    ]

    # We'll iterate chunk by chunk
    events_processed = 0
    for data_chunk in uproot.iterate(
            f"{args.input_file}:Evnts",
            expressions=branches,
            step_size=args.chunk_size,
            library="ak"  # get awkward arrays
    ):
        # data_chunk["invts"] is now an awkward array of records,
        # each record having fields like Q2, xBj, etc.
        invts = data_chunk["invts"]
        p_inc = data_chunk["P_Inc."]
        lamb = data_chunk["lamb_scat."]

        # For example, let's extract Q2 and xBj
        q2_array = invts["Q2"]
        xbj_array = invts["xBj"]

        # Show a few sample values
        print("Q2 sample:", q2_array[:5].tolist())
        print("xBj sample:", xbj_array[:5].tolist())

        # Or for the TLorentzVector lamb_scat, we can do:
        print("lamb: ", lamb)
        #lamb_px = lamb["px"]
        #lamb_py = lamb["py"]
        #lamb_pz = lamb["pz"]
        #lamb_e  = lamb["E"]

        # Do something with them...
        print(len(q2_array))
        print("win")
        exit(0)

        events_processed += len(q2_array)
        if args.max_events and events_processed >= args.max_events:
            break

    print(f"Done! Processed {events_processed} events.")

if __name__ == "__main__":
    main()
