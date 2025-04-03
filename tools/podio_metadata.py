import uproot
import argparse
import json
import numpy as np
import pandas as pd
from collections import defaultdict


def inspect_podio_metadata(file):
    """
    Extract and analyze podio_metadata from a ROOT file.

    PODIO metadata typically includes:
    - Class definitions
    - Collection relationships
    - ID information and references
    """
    metadata = {}

    # Check if podio_metadata exists
    if "podio_metadata" not in file:
        print("No podio_metadata tree found in this file.")
        return None

    print("\n=== PODIO Metadata ===")
    podio_meta = file["podio_metadata"]

    # Basic info about the metadata tree
    print(f"Type: {podio_meta.classname}")
    if podio_meta.classname != "TTree":
        print("Warning: podio_metadata is not a TTree as expected")
        return None

    print(f"Entries: {podio_meta.num_entries}")
    print(f"Branches: {len(podio_meta.keys())}")

    # Collect all branch names
    branch_names = list(podio_meta.keys())
    print("\nAvailable metadata branches:")
    for name in branch_names:
        print(f" - {name}")

    metadata["branches"] = branch_names

    # Process known metadata branches
    metadata_contents = {}

    # Try to extract common PODIO metadata branches
    common_branches = [
        "classIDs", "classVersions", "collectionIDs", "collections",
        "schemaVersion", "podioVersion"
    ]

    for branch in branch_names:
        try:
            # Get array from branch
            data = podio_meta[branch].array()

            if len(data) > 0:
                # Special handling for string data
                if isinstance(data[0], (np.ndarray, list)) and len(data[0]) > 0:
                    if isinstance(data[0][0], bytes):
                        # Convert bytes to strings
                        string_data = [[s.decode('utf-8', errors='replace') for s in item] for item in data]
                        print(f"\nBranch '{branch}':")
                        for item in string_data:
                            print(f" - {item}")
                        metadata_contents[branch] = string_data
                    else:
                        print(f"\nBranch '{branch}':")
                        print(data)
                        metadata_contents[branch] = data.tolist()
                else:
                    print(f"\nBranch '{branch}':")
                    print(data)
                    metadata_contents[branch] = data.tolist()
        except Exception as e:
            print(f"Error reading branch '{branch}': {str(e)}")

    metadata["contents"] = metadata_contents

    # Try to build collection-ID mapping
    try:
        if "collections" in metadata_contents and "collectionIDs" in metadata_contents:
            collections = metadata_contents["collections"]
            collection_ids = metadata_contents["collectionIDs"]

            if collections and collection_ids and len(collections) == len(collection_ids):
                collection_map = {}
                for i, collection_list in enumerate(collections):
                    if i < len(collection_ids):
                        for j, coll_name in enumerate(collection_list):
                            if j < len(collection_ids[i]):
                                collection_map[coll_name] = collection_ids[i][j]

                print("\nCollection to ID Mapping:")
                for coll, cid in collection_map.items():
                    print(f" - {coll}: {cid}")

                metadata["collection_map"] = collection_map
    except Exception as e:
        print(f"Error creating collection mapping: {str(e)}")

    # Try to build class-ID mapping
    try:
        if "classIDs" in metadata_contents:
            class_ids = metadata_contents["classIDs"]

            if class_ids:
                print("\nClass IDs:")
                for class_id_list in class_ids:
                    for class_id in class_id_list:
                        print(f" - {class_id}")

                metadata["class_ids"] = class_ids
    except Exception as e:
        print(f"Error processing class IDs: {str(e)}")

    # Check for PODIO schema and version information
    try:
        if "schemaVersion" in metadata_contents:
            schema_version = metadata_contents["schemaVersion"]
            print(f"\nSchema Version: {schema_version}")
            metadata["schema_version"] = schema_version

        if "podioVersion" in metadata_contents:
            podio_version = metadata_contents["podioVersion"]
            print(f"PODIO Version: {podio_version}")
            metadata["podio_version"] = podio_version
    except Exception as e:
        print(f"Error reading version information: {str(e)}")

    return metadata


def analyze_class_relationships(metadata):
    """Analyze the relationships between classes if possible."""
    if not metadata or "contents" not in metadata:
        return

    contents = metadata["contents"]

    # Look for relationship information which may be in different branches
    # depending on the PODIO version
    relationship_branches = [
        "relationShips", "relationships", "relations"
    ]

    for branch in relationship_branches:
        if branch in contents:
            print(f"\nClass Relationships (from {branch}):")
            relationships = contents[branch]

            for rel_list in relationships:
                for rel in rel_list:
                    print(f" - {rel}")


def check_collection_presence(file, metadata):
    """Check if the collections listed in metadata are present in the file."""
    if not metadata or "collection_map" not in metadata:
        return

    collection_map = metadata["collection_map"]

    print("\nChecking for Collections in File:")
    for collection_name in collection_map.keys():
        found = False
        # Try direct access
        try:
            if collection_name in file:
                print(f" ✓ {collection_name} (direct)")
                found = True
            # Try looking in events tree
            elif "events" in file:
                events = file["events"]
                if collection_name in events:
                    print(f" ✓ {collection_name} (in events tree)")
                    found = True
        except:
            pass

        if not found:
            print(f" ✗ {collection_name} (not found)")


def export_metadata_json(metadata, output_file):
    """Export the metadata to a JSON file."""
    if metadata:
        with open(output_file, 'w') as f:
            json.dump(metadata, f, indent=2, default=str)
        print(f"\nMetadata exported to {output_file}")


def main(input_file, export_json=False):
    """Main function to inspect PODIO metadata in a ROOT file."""
    print(f"Analyzing PODIO metadata in: {input_file}")

    try:
        with uproot.open(input_file) as file:
            print(f"Successfully opened file: {file.file_path}")

            # Extract PODIO metadata
            metadata = inspect_podio_metadata(file)

            if metadata:
                # Analyze relationships
                analyze_class_relationships(metadata)

                # Check if collections exist
                check_collection_presence(file, metadata)

                # Export to JSON if requested
                if export_json:
                    output_file = input_file.replace(".root", "_podio_metadata.json")
                    export_metadata_json(metadata, output_file)

    except Exception as e:
        print(f"Error analyzing file: {str(e)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Inspect PODIO metadata in ROOT files")
    parser.add_argument("input_file", help="The input ROOT file to analyze")
    parser.add_argument("-j", "--json", action="store_true",
                        help="Export metadata to JSON file")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Show more detailed information")

    args = parser.parse_args()
    print("Input file:", args.input_file)
    main(args.input_file, export_json=args.json)