import argparse

import uproot
import json
from pprint import pprint

def process_gp_branches(tree, entry=0):
    """Process GP*Keys and GP*Values branches into a dictionary."""
    data = {}
    
    for dtype in ['Double', 'Float', 'Int', 'String']:
        keys_branch = f'GP{dtype}Keys'
        values_branch = f'GP{dtype}Values'
        if keys_branch in tree and values_branch in tree:
            keys = tree[keys_branch].array()[entry].tolist()
            values = tree[values_branch].array()[entry].tolist()
            data[dtype.lower()] = dict(zip(keys, values))
    return data

def process_podio_metadata(tree, entry=0):
    """Process complex structures in podio_metadata tree."""
    data = {}
    
    # Handle EDMDefinitions (vector of tuples)
    # if "EDMDefinitions._0" in tree:
    #     edm_0 = tree["EDMDefinitions._0"].array()[entry].tolist()
    #     edm_1 = tree["EDMDefinitions._1"].array()[entry].tolist()
    #     data["EDMDefinitions"] = list(zip(edm_0, edm_1))
    
    # Handle PodioBuildVersion
    if "PodioBuildVersion/major" in tree:
        data["PodioBuildVersion"] = {
            'major': int(tree["PodioBuildVersion/major"].array()[entry]),
            'minor': int(tree["PodioBuildVersion/minor"].array()[entry]),
            'patch': int(tree["PodioBuildVersion/patch"].array()[entry])
        }
    
    # Process CollectionTypeInfo for events, metadata, runs
    for col in ['events', 'metadata', 'runs']:
        branch = f"{col}___CollectionTypeInfo"
        if f"{branch}._0" in tree:
            c0 = tree[f"{branch}._0"].array()[entry].tolist()
            c1 = tree[f"{branch}._1"].array()[entry].tolist()
            c2 = tree[f"{branch}._2"].array()[entry].tolist()
            c3 = tree[f"{branch}._3"].array()[entry].tolist()
            data[f"{col}_CollectionTypeInfo"] = list(zip(c0, c1, c2, c3))
    
    # Process idTables for events, metadata, runs
    for col in ['events', 'metadata', 'runs']:
        branch = f"{col}___idTable"
        ids_branch = f"{branch}/m_collectionIDs"
        names_branch = f"{branch}/m_names"
        if ids_branch in tree and names_branch in tree:
            ids = tree[ids_branch].array()[entry].tolist()
            names = tree[names_branch].array()[entry].tolist()
            data[f"{col}_idTable"] = dict(zip(names, ids))
    
    return data

def main(input_file, output_json):
    """Main function to process ROOT file and generate JSON output."""
    with uproot.open(input_file) as file:

        all_data = {}   # Combined metadata

        # Process metadata tree
        if "metadata" in file:
            metadata_tree = file["metadata"]
            metadata_data = process_gp_branches(metadata_tree)
            print("="*50 + "\nMetadata:\n" + "="*50)
            pprint(metadata_data)
            all_data['metadata'] = metadata_data
        else:
            print("No 'metadata' tree found in file")
        
        # Process podio_metadata tree
        podio_metadata_tree = file["podio_metadata"]
        podio_data = process_podio_metadata(podio_metadata_tree)
        
        # Process runs tree
        if "runs" in file:
            runs_tree = file["runs"]
            runs_data = process_gp_branches(runs_tree)
            print("\n" + "="*50 + "\nRuns:\n" + "="*50)
            pprint(runs_data)
        else:
            print("No 'runs' tree found in file")
        


        
        # Print to screen

        print("\n" + "="*50 + "\nPodio Metadata:\n" + "="*50)
        pprint(podio_data)

        
        # Write to JSON
        #with open(output_json, 'w') as f:
        #    json.dump(all_data, f, indent=2, ensure_ascii=False)
        #print(f"\nData written to {output_json}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help="Shows event level metadata")
    parser.add_argument("input_file", help="The Input file")
    parser.add_argument("-o", "--output", dest="output_file", default="", required=False, help="Output JSON file")
    args = parser.parse_args()

    main(args.input_file, args.output_file)
