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
    if "EDMDefinitions._0" in tree:
        edm_0 = tree["EDMDefinitions._0"].array()[entry].tolist()
        edm_1 = tree["EDMDefinitions._1"].array()[entry].tolist()
        data["EDMDefinitions"] = list(zip(edm_0, edm_1))
    
    # Handle PodioBuildVersion
    if "PodioBuildVersion/major" in tree:
        data["PodioBuildVersion"] = {
            'major': tree["PodioBuildVersion/major"].array()[entry],
            'minor': tree["PodioBuildVersion/minor"].array()[entry],
            'patch': tree["PodioBuildVersion/patch"].array()[entry]
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
        # Process metadata tree
        metadata_tree = file["metadata"]
        metadata_data = process_gp_branches(metadata_tree)
        
        # Process podio_metadata tree
        podio_metadata_tree = file["podio_metadata"]
        podio_data = process_podio_metadata(podio_metadata_tree)
        
        # Process runs tree
        runs_tree = file["runs"]
        runs_data = process_gp_branches(runs_tree)
        
        # Combine all data
        all_data = {
            'metadata': metadata_data,
            'podio_metadata': podio_data,
            'runs': runs_data
        }
        
        # Print to screen
        print("="*50 + "\nMetadata:\n" + "="*50)
        pprint(metadata_data)
        print("\n" + "="*50 + "\nPodio Metadata:\n" + "="*50)
        pprint(podio_data)
        print("\n" + "="*50 + "\nRuns:\n" + "="*50)
        pprint(runs_data)
        
        # Write to JSON
        #with open(output_json, 'w') as f:
        #    json.dump(all_data, f, indent=2, ensure_ascii=False)
        #print(f"\nData written to {output_json}")

if __name__ == "__main__":
    input_file = "/volatile/eic/romanov/meson-structure-2025-02/reco/k_lambda_5x41_5000evt_200.edm4hep.root"  # Replace with your ROOT file path
    output_json = "k_lambda_5x41_5000evt_200.edm4hep.json"
    main(input_file, output_json)