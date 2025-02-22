#!/usr/bin/env python3
import uproot
import json

def read_generic_parameters(tree):
    """
    Reads the GP* branches from a TTree (like 'metadata' or 'runs'), which are:
      - GPDoubleKeys, GPDoubleValues
      - GPFloatKeys, GPFloatValues
      - GPIntKeys, GPIntValues
      - GPStringKeys, GPStringValues
    Returns a dictionary with key -> list of values for each type.
    """
    out_dict = {}

    # We assume there's only 1 entry in these metadata TTrees, so we take [0].
    double_keys = tree["GPDoubleKeys"].array(library="python")[0]
    double_vals = tree["GPDoubleValues"].array(library="python")[0]
    for k, v in zip(double_keys, double_vals):
        out_dict[k] = v  # v is list of doubles

    float_keys = tree["GPFloatKeys"].array(library="python")[0]
    float_vals = tree["GPFloatValues"].array(library="python")[0]
    for k, v in zip(float_keys, float_vals):
        out_dict[k] = v  # v is list of floats

    int_keys = tree["GPIntKeys"].array(library="python")[0]
    int_vals = tree["GPIntValues"].array(library="python")[0]
    for k, v in zip(int_keys, int_vals):
        out_dict[k] = v  # v is list of ints

    str_keys = tree["GPStringKeys"].array(library="python")[0]
    str_vals = tree["GPStringValues"].array(library="python")[0]
    for k, v in zip(str_keys, str_vals):
        out_dict[k] = v  # v is list of strings

    return out_dict

def main():
    # Open your EDM4hep or PODIO-based ROOT file
    filename = "example.root"  # Change to your actual file
    f = uproot.open(filename)

    # Print the available keys
    print("ROOT file keys:", f.keys())
    # e.g. you'll see ["metadata;1", "podio_metadata;1", "runs;1", ...]

    # 1) Read 'metadata' TTree
    if "metadata" in f.keys():
        meta_tree = f["metadata"]
        meta_dict = read_generic_parameters(meta_tree)
        print("\n=== metadata (GP parameters) ===")
        for k, v in meta_dict.items():
            print(f"  {k} => {v}")
    else:
        meta_dict = {}
        print("\nNo TTree named 'metadata' found.")

    # 2) Read 'runs' TTree (same structure as 'metadata')
    if "runs" in f.keys():
        runs_tree = f["runs"]
        runs_dict = read_generic_parameters(runs_tree)
        print("\n=== runs (GP parameters) ===")
        for k, v in runs_dict.items():
            print(f"  {k} => {v}")
    else:
        runs_dict = {}
        print("\nNo TTree named 'runs' found.")

    # 3) Read 'podio_metadata' TTree. This typically has arrays describing
    #    collection IDs, versions, etc. We'll just dump the arrays.
    if "podio_metadata" in f.keys():
        podio_tree = f["podio_metadata"]
        print("\n=== podio_metadata contents ===")

        # You can list branches with:
        print("Branches in podio_metadata:", podio_tree.keys())
        # For example: "EDMDefinitions._0", "EDMDefinitions._1",
        #              "PodioBuildVersion.major", etc.

        # We'll read them one by one, printing the first entry if present.
        podio_dict = {}
        for branch_name in podio_tree.keys():
            arr = podio_tree[branch_name].array(library="python")
            # Usually there's 1 entry. We'll just grab arr[0].
            if len(arr) > 0:
                val = arr[0]
            else:
                val = None
            podio_dict[branch_name] = val
            print(f"  {branch_name} => {val}")
    else:
        podio_dict = {}
        print("\nNo TTree named 'podio_metadata' found.")

    # Combine them into one big dictionary if desired
    all_metadata = {
        "metadata": meta_dict,
        "runs": runs_dict,
        "podio_metadata": podio_dict
    }

    # 4) Write to JSON file
    out_json = "metadata_dump.json"
    with open(out_json, "w") as jf:
        json.dump(all_metadata, jf, indent=2)
    print(f"\nWrote combined metadata to {out_json}")

if __name__ == "__main__":
    main()
