import csv
import argparse
import re
import numpy as np

# Define the mapping for modifications
unimod_mapping = {
    "57.0214": "Carbamidomethyl",
    "0.984": "Deamidated",
    "15.9949": "Oxidation"
}

def parse_modifications(modifications_str):
    modifications = []
    if modifications_str:
        for mod in modifications_str.split(','):
            try:
                # Extract position and mass using regular expressions
                match = re.search(r'(\d+)[A-Z]\(([\d\.]+)\)', mod.strip())
                if match:
                    position = int(match.group(1))
                    mass = match.group(2)
                    modification_name = unimod_mapping.get(mass, 'Unknown')

                    if modification_name == 'Unknown':
                        print(f"Warning: Unrecognized mass {mass} in modification '{mod}'")

                    modification = f"{position}|{modification_name}"
                    modifications.append(modification)
                else:
                    print(f"Warning: Unexpected modification format '{mod}'")
            except ValueError as e:
                print(f"Error parsing modification '{mod}': {e}")
                continue
    return modifications

def reformat_csv(input_file, output_file):
    with open(input_file, mode='r', newline='') as infile, open(output_file, mode='w', newline='') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')  # Use tab delimiter for TSV
        headers = reader.fieldnames
        print("Headers in the input file:", headers)  # Debugging: Print headers

        outfile.write('seq,modifications,charge,ion_mobility,observed_mass,Retention\n')  # Write header manually

        for row in reader:
            try:
                seq = row['Peptide']
                charge = int(row['Charge'])
                assigned_modifications = row['Assigned Modifications']
                ion_mobility = float(row['Ion Mobility'])
                observed_mass = float(row['Observed Mass'])  # Include observed mass
                Retention = row['Retention']
                modifications = parse_modifications(assigned_modifications)
                modifications_str = '|'.join(modifications)  # Join modifications with pipe symbol

                outfile.write(f'{seq},{modifications_str},{charge},{ion_mobility},{observed_mass}, {Retention}\n')
            except KeyError as e:
                print(f"Error processing row: Missing expected column {e}")
                continue

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Reformat TSV file.')
    parser.add_argument('-i', '--input', required=True, help='Input TSV file')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    
    args = parser.parse_args()
    
    reformat_csv(args.input, args.output)
