import csv

# Open the input file for reading and the output file for writing
with open('../../euka_dir/euka_db.bins', 'r') as infile, open('nodes.txt', 'w', newline='') as outfile:
    writer = csv.writer(outfile, delimiter=' ')
    # Loop through each line in the input file
    for line in infile:
        # Split the line into parts based on spaces
        parts = line.strip().split()
        
        # The first part is the taxon name
        taxon = parts[0]
        
        # The remaining parts are the triplets of numbers (entropy, start_node, end_node)
        # We are interested in the start and end nodes (positions 1 and 2 in each triplet)
        for i in range(1, len(parts), 3):
            start_node = int(float(parts[i]))
            end_node = int(float(parts[i+1]))
            
            # Write the start and end nodes along with the taxon name to the output file
            writer.writerow([taxon, start_node, end_node])

