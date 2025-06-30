from Bio import Entrez, SeqIO
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
import os

# Specify your email for NCBI
Entrez.email = "feng.xian@univie.ac.at"  # Always provide your email for NCBI Entrez

# Load the CSV file containing the sseqid list
file_path = "D:/Feng/DDA_PASEF/20250623_HumanFeces_VolcanicEnzyme_Denovo/Volcanic_enzyme_denovo_20250623/novoMP_Output/BlastP/novoMP_BlastP_LCA_filtered_output_80Pident.csv"  # Replace with your actual file path
df = pd.read_csv(file_path)

# Clean up sseqid column: remove leading/trailing spaces and drop empty or missing entries
df['sseqid'] = df['sseqid'].str.strip()  # Strip leading/trailing spaces
df = df.dropna(subset=['sseqid'])  # Drop rows where sseqid is missing

# Print the number of total and unique sseqid entries
print(f"Total sseqid entries: {len(df)}")
seq_ids = df['sseqid'].unique()
print(f"Unique sseqid entries after de-duplication: {len(seq_ids)}")

# Paths
output_fasta = "D:/Feng/DDA_PASEF/20250623_HumanFeces_VolcanicEnzyme_Denovo/Volcanic_enzyme_denovo_20250623/novoMP_Output/BlastP/NCBI_retrived_sequences_80pident.fasta"  # Replace with the desired output file path
failed_ids_file = "D:/Feng/DDA_PASEF/20250623_HumanFeces_VolcanicEnzyme_Denovo/Volcanic_enzyme_denovo_20250623/novoMP_Output/BlastP/failed_sseqid.txt"

# Create the output directory if it doesn't exist
output_dir = os.path.dirname(output_fasta)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(f"Created directory: {output_dir}")

# Check if output file exists
if os.path.exists(output_fasta):
    print(f"File already exists: {output_fasta}")
else:
    # Explicitly create an empty FASTA file before anything else happens
    with open(output_fasta, 'w') as f:
        pass  # Create the file and immediately close it
    print(f"Empty FASTA file created: {output_fasta}")

batch_size = 56  # NCBI limits the number of sequences per request
max_workers = 5  # Number of threads
delay_between_batches = 0.5  # Delay (in seconds) between requests to prevent hitting rate limit

# Initialize a list to store failed sseqids
failed_sseqids = []

# Function to fetch sequences for a batch
def fetch_sequences(batch):
    retries = 3
    for attempt in range(retries):
        try:
            with Entrez.efetch(db="protein", id=",".join(batch), rettype="fasta", retmode="text") as handle:
                records = list(SeqIO.parse(handle, "fasta"))
                return records
        except Exception as e:
            print(f"Error fetching batch {batch}, attempt {attempt + 1}: {e}")
            time.sleep(5)  # Wait 5 seconds before retrying
    return None

# Use a ThreadPoolExecutor to fetch sequences in parallel
with ThreadPoolExecutor(max_workers=max_workers) as executor:
    futures = []

    # Submit fetching tasks
    for i in tqdm(range(0, len(seq_ids), batch_size), desc="Submitting batches"):
        batch = seq_ids[i:i + batch_size]
        futures.append(executor.submit(fetch_sequences, batch))
        time.sleep(delay_between_batches)  # Add a delay between batch requests

    # Process the results as they complete
    with open(output_fasta, 'a') as fasta_output:
        for future in tqdm(as_completed(futures), desc="Fetching sequences", total=len(futures)):
            try:
                result = future.result()
                if result:  # Write successful fetches to the FASTA file
                    SeqIO.write(result, fasta_output, "fasta")
                    fasta_output.flush()  # Force immediate writing to the file
                    print(f"Written {len(result)} sequences to {output_fasta}")
                else:
                    # Log the failed batch
                    failed_batch = futures[futures.index(future)].args[0]
                    print(f"Batch failed: {failed_batch}")
                    failed_sseqids.extend(failed_batch)
            except Exception as e:
                # Handle any other exceptions, log the batch as failed
                failed_batch = futures[futures.index(future)].args[0]
                print(f"Error processing batch {failed_batch}: {e}")
                failed_sseqids.extend(failed_batch)

# Save failed sseqid list to a file
if failed_sseqids:
    with open(failed_ids_file, 'w') as failed_file:
        for sseqid in failed_sseqids:
            failed_file.write(sseqid + '\n')
    print(f"Failed sseqid entries saved to {failed_ids_file}")
else:
    print("All sequences fetched successfully, no failed sseqids.")

print(f"FASTA sequences saved to {output_fasta}")
