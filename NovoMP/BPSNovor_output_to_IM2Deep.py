import csv
import argparse

# Define the mapping for modifications
unimod_mapping = {
    "57.0215": "Carbamidomethyl",
    "0.984": "Deamidated",
    "15.9949": "Oxidation"
}

def reformat_csv(input_file, output_file):
    with open(input_file, mode='r', newline='') as infile, open(output_file, mode='w', newline='') as outfile:
        reader = csv.DictReader(infile)
        outfile.write('seq,modifications,charge,ion_mobility,observed_mass\n')  # Write header manually

        for row in reader:
            seq = row['stripped_peptide']
            ptms = row['ptms'].split(',')
            ptm_locations = row['ptm_locations'].split(',')

            charge = row['charge']
            ion_mobility = row['ook0']
            observed_mass = row['calc_mh']
            
            modifications = []
            for ptm, loc in zip(ptms, ptm_locations):
                if loc != 'NA':
                    modification = f"{int(loc) + 1}|{unimod_mapping.get(ptm, 'Unknown')}"
                    modifications.append(modification)
            
            modifications_str = '|'.join(modifications)  # Join modifications with pipe symbol
            outfile.write(f'{seq},{modifications_str},{charge},{ion_mobility},{observed_mass}\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Reformat CSV file.')
    parser.add_argument('-i', '--input', required=True, help='Input CSV file')
    parser.add_argument('-o', '--output', required=True, help='Output CSV file')
    
    args = parser.parse_args()
    
    reformat_csv(args.input, args.output)
