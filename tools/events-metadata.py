import uproot
from pprint import pprint
import argparse

def process_gp_branches(tree, entry):
    """Process GP*Keys and GP*Values branches for a specific entry/event."""
    data = {}
    for dtype in ['Double', 'Float', 'Int', 'String']:
        keys_branch = f'GP{dtype}Keys'
        values_branch = f'GP{dtype}Values'
        
        if keys_branch in tree and values_branch in tree:
            # Get data for this specific event
            keys = tree[keys_branch].array()[entry].tolist()
            values = tree[values_branch].array()[entry].tolist()
            data[dtype.lower()] = dict(zip(keys, values))
    return data


def main(input_file):
    with uproot.open(input_file) as file:
        # Get the events tree
        events_tree = file["events"]
        
        # Process first 10 events
        n_events = min(10, events_tree.num_entries)
        events_data = [process_gp_branches(events_tree, i) for i in range(n_events)]

        # Print results
        print("="*50 + "\nFirst 10 Events Metadata:\n" + "="*50)
        for i, event_data in enumerate(events_data):
            print(f"\nEvent {i}:")
            pprint(event_data)
            print("-"*50)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help="Shows event level metadata")
    parser.add_argument("input_file", help="The Input file")
    args = parser.parse_args()
    print("Input file: ", args.input_file)
    main(args.input_file)
