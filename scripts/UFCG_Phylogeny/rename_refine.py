import json
import re
import argparse
# Version 31-03-2025 - Mikael Terp
## README ##
#

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Replace accession numbers with organism name and strain in a Newick file.")
parser.add_argument("--newick_path", required=True, help="Path to the Newick file to be updated.")
parser.add_argument("--jsonl_path", required=True, help="Path to the JSONL file containing accession-to-organism mapping.")
args = parser.parse_args()

# Step 1: Create a mapping of accession numbers to organism name and strain
accession_to_name_strain = {}

with open(args.jsonl_path, 'r') as jsonl_file:
    for line in jsonl_file:
        data = json.loads(line)
        accession = data['accession']
        organism_name = data['organism']['organismName']
        strain = None

        # Extract strain from attributes or infraspecificNames
        if 'attributes' in data.get('biosample', {}):
            for attribute in data['biosample']['attributes']:
                if attribute['name'] == 'strain':
                    strain = attribute['value']
                    break
        elif 'infraspecificNames' in data['organism']:
            strain = data['organism']['infraspecificNames'].get('strain')

        # Combine organism name and strain
        name_strain = f"{organism_name} {strain}" if strain else organism_name
        accession_to_name_strain[accession] = name_strain

# Step 2: Read the Newick file
with open(args.newick_path, 'r') as newick_file:
    newick_content = newick_file.read()

# Step 3: Replace accession numbers with organism name and strain
for accession, name_strain in accession_to_name_strain.items():
    newick_content = re.sub(rf'{accession}_[^:]+', name_strain, newick_content)

# Step 4: Write the updated Newick content back to the file
with open(args.newick_path, 'w') as newick_file:
    newick_file.write(newick_content)

print("Accession numbers have been replaced with organism name and strain in the Newick file.")