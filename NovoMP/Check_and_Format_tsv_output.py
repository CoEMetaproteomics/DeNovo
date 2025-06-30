import pandas as pd

# Load the corrected output file
file_path = "D:/Feng/DDA_PASEF/20250623_HumanFeces_VolcanicEnzyme_Denovo/Volcanic_enzyme_denovo_20250623/novoMP_Output/BlastP/novoMP_peptides_blastp_output.tsv"  # Replace with actual path

# Define the column headers
column_headers = ['qseqid', 'full_qseq', 'sseqid', 'pident', 'nident', 'length', 'mismatch', 'gapopen',
                  'qstart', 'qend', 'sstart', 'send', 'score', 'evalue', 'bitscore', 'staxids',
                  'sscinames', 'skingdoms', 'salltitles']

# Function to check and fix row integrity
def check_and_fix_file(file_path, expected_columns=19):
    # Open the file and process line by line
    corrected_lines = []
    corrected_lines.append('\t'.join(column_headers))  # Add the specified column headers as the first line

    with open(file_path, 'r', encoding='utf-8') as file:
        header = file.readline().strip()  # Skip the original header (or adjust as needed)
        for line in file:
            columns = line.strip().split('\t')
            # Check if the number of columns matches the expected number
            if len(columns) == expected_columns:
                corrected_lines.append(line.strip())  # Keep valid rows
            else:
                # Attempt to fix the row by merging split content into a single line
                # If the previous line is incomplete, append the current line to it
                if len(corrected_lines) > 0 and len(corrected_lines[-1].split('\t')) < expected_columns:
                    corrected_lines[-1] += " " + line.strip()  # Merge with previous row
                else:
                    # If no previous row to merge, just skip this line (or log it)
                    print(f"Skipping malformed line: {line.strip()}")

    # Write the corrected content back to a new file
    corrected_file_path = 'D:/Feng/DDA_PASEF/20250623_HumanFeces_VolcanicEnzyme_Denovo/Volcanic_enzyme_denovo_20250623/novoMP_Output/BlastP/novoMP_peptides_blastp_output_Validated.tsv'  # Replace with actual path
    with open(corrected_file_path, 'w', encoding='utf-8') as output_file:
        for corrected_line in corrected_lines:
            output_file.write(corrected_line + '\n')
    
    print(f"File validation and correction completed. Corrected file saved at {corrected_file_path}")

# Run the check and fix function
check_and_fix_file(file_path)

