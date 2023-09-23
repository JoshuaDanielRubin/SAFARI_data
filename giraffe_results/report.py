import json
import csv
import sys

# Load nodes data
nodes_data = {}
with open('nodes.txt', 'r') as f:
    reader = csv.reader(f, delimiter=' ')
    for row in reader:
        nodes_data[row[0]] = (int(row[1]), int(row[2]))

# Special clade mapping
special_clade_mapping = {'Ovis': 'Bovidae', 'Capra': 'Bovidae', 'Aquatica': 'Lampyridae', 'Luciola': 'Lampyridae', \
                         'Ursus': 'Ursidae', 'Tremarctos' : 'Ursidae'
                        }

# Initialize counters
true_positive = 0
false_positive = 0
true_negative = 0
false_negative = 0

# Details of the reads
read_details = []

# Read and process JSON lines
with open(sys.argv[1], 'r') as f:
    for line in f:
        read = json.loads(line.strip())
        read_name = read.get("name", "")
        path = read.get("path", None)
        
        # Check if it is a contaminant
        if not read_name.startswith("Mito_"):
            if path is None:
                true_negative += 1
                #read_details.append({'read_name': read_name, 'status': 'True Negative'})
            else:
                false_positive += 1
                #read_details.append({'read_name': read_name, 'status': 'False Positive (Contaminant aligned)'})
            continue

        # Get clade from read name
        clade_from = read_name.split('_')[1]
        
        # Check if clade is special and map it to another clade
        clade_from = special_clade_mapping.get(clade_from, clade_from)
        
        # If path does not exist, it is a false negative
        if path is None:
            false_negative += 1
            #read_details.append({'read_name': read_name, 'status': 'False Negative', 'clade_from': clade_from})
            continue
        
        # Get the first node id from the path
        first_node_id = int(path["mapping"][0]["position"]["node_id"])
        
        # Find to which clade it mapped
        clade_to = None
        for clade, (start, end) in nodes_data.items():
            if start <= first_node_id <= end:
                clade_to = clade
                break

        # If clade_to is found, check if it is correctly mapped
        if clade_to:
            if clade_from == clade_to:
                true_positive += 1
                #read_details.append({'read_name': read_name, 'status': 'True Positive', 'clade_from': clade_from, 'clade_to': clade_to})
            else:
                false_positive += 1
                #read_details.append({'read_name': read_name, 'status': 'False Positive (Incorrect clade)', 'clade_from': clade_from, 'clade_to': clade_to})
        else:
            false_positive += 1
            read_details.append({'read_name': read_name, 'status': 'False Positive (Clade not found)', 'clade_from': clade_from})

# Print detailed report
for detail in read_details:
    print(detail)

# Print summary
print(f"True Positives: {true_positive}")
print(f"False Positives: {false_positive}")
print(f"True Negatives: {true_negative}")
print(f"False Negatives: {false_negative}")

