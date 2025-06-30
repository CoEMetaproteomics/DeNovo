import csv
import argparse
import numpy as np

# Constants
MASS_H = 1.007276466812  # Mass of a proton in Daltons

def ook0_to_ccs(ook0, observed_mass, charge):
    if ook0 == 0:
        return 0
    mol_weight_gas = 28.0134
    room_temp = 305
    conversion_factor = 18509.863216340458
    reduced_mass = observed_mass * mol_weight_gas / (observed_mass + mol_weight_gas)
    return conversion_factor * np.abs(charge) / (np.sqrt(reduced_mass * room_temp) * (1/ook0))  # Use 1/ook0 here

def calculate_ccs(input_file, output_file):
    with open(input_file, mode='r', newline='') as infile, open(output_file, mode='w', newline='') as outfile:
        reader = csv.DictReader(infile, delimiter=',')  # Use comma delimiter for CSV
        headers = reader.fieldnames
        print("Headers in the input file:", headers)  # Debugging: Print headers

        # Add CCS column to the output headers
        output_headers = headers + ['CCS']
        writer = csv.DictWriter(outfile, fieldnames=output_headers, delimiter=',')  # Use comma delimiter for CSV
        writer.writeheader()

        for row in reader:
            try:
                print("Processing row:", row)  # Debugging: Print the row being processed
                charge = int(row['charge'])  # Ensure column name matches exactly
                ion_mobility = float(row['ion_mobility'])  # Ensure column name matches exactly
                observed_mass = float(row['observed_mass'])  # Ensure column name matches exactly

                # Calculate CCS
                ccs = ook0_to_ccs(ion_mobility, observed_mass, charge)

                # Add CCS to the row
                row['CCS'] = ccs

                # Write the updated row to the output file
                writer.writerow(row)
            except KeyError as e:
                print(f"Error processing row: Missing expected column {e}")
            except ValueError as e:
                print(f"Error processing row: Invalid data {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate CCS from Ion Mobility and other parameters in a CSV file.')
    parser.add_argument('-i', '--input', required=True, help='Input CSV file')
    parser.add_argument('-o', '--output', required=True, help='Output CSV file')
    
    args = parser.parse_args()
    
    calculate_ccs(args.input, args.output)
