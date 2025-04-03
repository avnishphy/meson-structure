import uproot
import argparse
import awkward as ak
import matplotlib.pyplot as plt


def main(input_file, n_events=10):

    # Open root file
    with uproot.open(input_file) as file:

        # Get the events tree from file
        tree = file["events"]

        # Process first 10 events
        n_events = min(n_events, tree.num_entries)
        metadata = tree.arrays(('GPStringKeys', 'GPStringValues'), entry_stop=n_events)

        # Convert values to float
        values = ak.strings_astype(metadata.GPStringValues, "float64")

        # If we look at dimensions of GPStringKeys and GPStringValues we see they mismatch
        # Somehow after all the copying of metadta vlues end up in a single dimension array, like:
        # "dis_q2: [1.2345]" while it should be "dis_q2: 1.2345"
        # Here is the debug output:
        metadata.GPStringKeys.type.show()
        metadata.GPStringValues.type.show()

        # so we remove this last layer by flattening with axis=-1 (remove innermost nesting)
        values = ak.flatten(values, axis=-1)

        # Now values are the same as for keys
        values.type.show()
        metadata.GPStringKeys.type.show()

        # zip keys and values together
        res = ak.zip(
            {
                "key": metadata.GPStringKeys,
                "value": values,
            })

        # Should be nice now! Here is debug print
        res.show(all=True)

        # Now we can get e.g. x vs q2
        xbj = res[res.key == "dis_xbj"].value
        q2  = res[res.key == "dis_q2"].value

        # Convert them to numpy as we use matplotlib which knows how to work with numpy
        xbj_flat = ak.to_numpy(ak.flatten(xbj))
        q2_flat  = ak.to_numpy(ak.flatten(q2))

        # build plot
        plt.hist2d(xbj_flat, q2_flat, bins=50)
        plt.xlabel("xBj")
        plt.ylabel("Q^2")
        plt.title("dis_xbj vs dis_q2")
        plt.savefig("dis_xbj_vs_dis_q2.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help="Shows event level metadata")
    parser.add_argument("input_file", help="The Input file")
    parser.add_argument("--events", "-e", help="Number of events")
    args = parser.parse_args()
    print("Input file: ", args.input_file)
    main(args.input_file, int(args.events))
