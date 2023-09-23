import os
import re

def parse_stat_file(filename):
    # Dictionary to hold the data
    data = {}
    try:
        with open(filename, 'r') as file:
            for line in file:
                # Split the line on ':'
                parts = line.strip().split(':')
                if len(parts) == 2:
                    key, value = parts
                    data[key.strip()] = value.strip()
        return data
    except Exception as e:
        raise ValueError(f"Failed to parse {filename}: {e}")

def parse_directory(directory_path):
    # List to hold all the data
    all_data = []

    for filename in os.listdir(directory_path):
        if filename.endswith('.stat'):
            # Full path to the file
            file_path = os.path.join(directory_path, filename)

            # Parse the file
            file_data = parse_stat_file(file_path)

            # Extract additional info from filename
            match = re.match(r'numtS_and_gen_0_n(\d+)_l(\d+)_d([^_]+)_s([^_]+)_((?:[^_]+_)?[^_]+).stat', filename)

            if match:
                file_data.update({
                    'endogenous_fragments': match.group(1),
                    'fragment_length': match.group(2),
                    'damage_pattern': match.group(3),
                    'fraction_of_reads': match.group(4),
                    'aligner': match.group(5)
                })
                all_data.append(file_data)
            else:
                raise ValueError(f"Filename {filename} does not match the expected format")

    return all_data

def generate_summary(all_data):
    # Check if data is empty
    if not all_data:
        raise ValueError("No data to summarize")

    # Dictionary to hold summary data per aligner
    summary_data = {}

    # Constants for the total number of mappable reads and NUMTs
    TOTAL_MAPPABLE_READS = 10000
    TOTAL_NUMTS = 100

    # Iterate through all data
    for data in all_data:
        aligner = data.get('aligner')
        if aligner not in summary_data:
            summary_data[aligner] = {
                'total_files': 0,
                'mapped_percent_sum': 0,
                'correctmap_percent_sum': 0,
                'unmapped_percent_sum': 0,
                'numt_percent_sum': 0
            }

        # Increment total_files counter
        summary_data[aligner]['total_files'] += 1

        summary_data[aligner]['mapped_percent_sum'] += float(data.get('mapped%', 0))
        summary_data[aligner]['correctmap_percent_sum'] += float(data.get('correctmap%', 0))
        summary_data[aligner]['unmapped_percent_sum'] += float(data.get('unmapped%', 0))
        summary_data[aligner]['numt_percent_sum'] += float(data.get('numt%', 0))

    # Calculate averages, totals, and print summary
    for aligner, data in summary_data.items():
        total_files = data['total_files']
        
        avg_mapped_percent = data['mapped_percent_sum'] / total_files
        avg_correctmap_percent = data['correctmap_percent_sum'] / total_files
        avg_unmapped_percent = data['unmapped_percent_sum'] / total_files
        avg_numt_percent = data['numt_percent_sum'] / total_files

        total_mapped = (avg_mapped_percent / 100) * TOTAL_MAPPABLE_READS
        total_correctmap = (avg_correctmap_percent / 100) * TOTAL_MAPPABLE_READS
        total_unmapped = (avg_unmapped_percent / 100) * TOTAL_MAPPABLE_READS
        total_numt = (avg_numt_percent / 100) * TOTAL_NUMTS

        print(f"Aligner: {aligner}")
        print(f"Total Files: {total_files}")
        #print(f"Total Mappable Reads: {TOTAL_MAPPABLE_READS}")
        #print(f"Average Mapped Percentage: {avg_mapped_percent:.2f}%")
        print(f"Total Mapped: {total_mapped}")
        print(f"Average Correctly Mapped Percentage: {avg_correctmap_percent:.2f}%")
        print(f"Total Correctly Mapped: {total_correctmap}")
        #print(f"Average Unmapped Percentage: {avg_unmapped_percent:.2f}%")
        #print(f"Total Unmapped: {total_unmapped}")
        print(f"Average NUMT Percentage: {avg_numt_percent:.2f}%")
        print(f"Total NUMT: {total_numt}")
        print("-" * 50)


def main():
    directory_path = 'alignments'
    if not os.path.exists(directory_path):
        raise ValueError(f"Directory {directory_path} does not exist")

    all_data = parse_directory(directory_path)
    generate_summary(all_data)

if __name__ == '__main__':
    main()

