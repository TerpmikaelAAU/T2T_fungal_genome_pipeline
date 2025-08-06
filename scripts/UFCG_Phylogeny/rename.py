import json
import re
# Version 31-03-2025 - Mikael Terp
## README ##
#



# Step 1: Create a mapping of accession numbers to taxId
accession_to_taxid = {}

with open('path/to/.jsonl', 'r') as jsonl_file:
    for line in jsonl_file:
        data = json.loads(line)
        accession = data['accession']
        taxid = data['organism']['taxId']
        accession_to_taxid[accession] = taxid

# Step 2: Read the Newick file and replace accession numbers with taxId
newick_file_path = 'Path/to/.nwk'
with open(newick_file_path, 'r') as newick_file:
    newick_content = newick_file.read()

# Step 3: Replace accession numbers with taxId and remove leftover parts
for accession, taxid in accession_to_taxid.items():
    newick_content = re.sub(rf'{accession}_[^:]+', str(taxid), newick_content)

# Step 4: Write the updated Newick content back to the file
with open(newick_file_path, 'w') as newick_file:
    newick_file.write(newick_content)

print("Accession numbers have been replaced with taxId in the Newick file.")