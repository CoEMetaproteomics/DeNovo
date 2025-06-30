import pandas as pd
from tqdm import tqdm
import re

# Enable tqdm for pandas apply operations
tqdm.pandas()

# Step 1: Load input files
input_file = 'D:/Feng/DDA_PASEF/20250623_HumanFeces_VolcanicEnzyme_Denovo/Volcanic_enzyme_denovo_20250623/novoMP_Output/BlastP/novoMP_peptides_blastp_output_Validated.tsv'
diamond_df = pd.read_csv(input_file, sep='\t', on_bad_lines='skip')
lca_df = pd.read_csv('D:/Feng/DDA_PASEF/20250623_HumanFeces_VolcanicEnzyme_Denovo/Volcanic_enzyme_denovo_20250623/novoMP_Output/BlastP/novoMP_peptides_blastp_output_LCA.tsv', sep='\t')
taxa_lineage_df = pd.read_csv("D:/Feng/DDA_PASEF/20240605_Ultra_mouse_fractionas_Novor_Fasta_mapping_Mgnify/Ultra_Mouse_Feces_Fractions_Novor/Re-Process_Rep1_20240625/20240910_Filtering/BlastP_peptide/NCBI_Taxa_Info/full_taxon_lineage.csv")

# Step 1.1: Remove entries with invalid or missing taxa info
print("Excluding TaxaID=0 entries from the LCA file...")
lca_df = lca_df[lca_df['TaxaID'] != 0]

print("Excluding empty 'superkingdom' entries from the lineage file...")
taxa_lineage_df = taxa_lineage_df[taxa_lineage_df['superkingdom'].notna()]

# Step 1.2: Filter for relevant branches in the taxonomy tree
valid_taxa_ids = {2, 4751, 2157, 10239}
taxa_lineage_df = taxa_lineage_df[
    (taxa_lineage_df['superkingdom'].isin(valid_taxa_ids)) |
    (taxa_lineage_df['kingdom'].isin(valid_taxa_ids)) |
    (taxa_lineage_df['phylum'].isin(valid_taxa_ids)) |
    (taxa_lineage_df['class'].isin(valid_taxa_ids)) |
    (taxa_lineage_df['order'].isin(valid_taxa_ids)) |
    (taxa_lineage_df['family'].isin(valid_taxa_ids)) |
    (taxa_lineage_df['genus'].isin(valid_taxa_ids)) |
    (taxa_lineage_df['species'].isin(valid_taxa_ids))
]

# Step 1.3: Convert staxids and TaxaID to string
diamond_df['staxids'] = diamond_df['staxids'].astype(str)
lca_df['TaxaID'] = lca_df['TaxaID'].astype(str)

# Step 2: Merge Diamond output with LCA results
print("Merging Diamond output with LCA data...")
merged_df = pd.merge(diamond_df, lca_df, on='qseqid', how='left')

# Step 3: Define taxonomic match function
def match_staxids_vectorized(staxids, lca_taxa_id):
    if pd.isna(staxids) or pd.isna(lca_taxa_id):
        return False
    staxid_list = [s.strip() for s in str(staxids).split(';')]
    return str(lca_taxa_id).strip() in staxid_list

print("Finding LCA matches...")
merged_df['lca_matches'] = merged_df.progress_apply(
    lambda row: match_staxids_vectorized(row['staxids'], row['TaxaID']),
    axis=1
)

# Step 3.5: Debug mismatches
def debug_lca_match(row):
    staxid_list = [s.strip() for s in str(row['staxids']).split(';')]
    expected = str(row['TaxaID']).strip()
    return {
        'qseqid': row['qseqid'],
        'staxids': row['staxids'],
        'TaxaID': row['TaxaID'],
        'parsed_staxids': staxid_list,
        'TaxaID_in_list': expected in staxid_list,
        'lca_matches': row['lca_matches']
    }

debug_df = merged_df[merged_df['lca_matches'] == False].apply(debug_lca_match, axis=1, result_type='expand')
debug_df = debug_df[debug_df['TaxaID_in_list'] == True]

if not debug_df.empty:
    print("⚠️ Debugging: These rows should have matched but didn't:")
    print(debug_df[['qseqid', 'staxids', 'TaxaID', 'parsed_staxids', 'lca_matches']])

# Step 4: Higher taxonomic match logic
taxa_set = set(taxa_lineage_df['tax_id'])

def belongs_to_higher_taxon_vectorized(staxids, lca_taxa_id, lineage_df):
    if pd.isna(staxids) or pd.isna(lca_taxa_id):
        return False
    staxid_list = [s.strip() for s in str(staxids).split(';')]
    for staxid in staxid_list:
        try:
            staxid_int = int(staxid)
            if staxid_int in taxa_set:
                lineage_row = lineage_df[lineage_df['tax_id'] == staxid_int]
                if not lineage_row.empty and str(lca_taxa_id) in lineage_row.astype(str).values:
                    return True
        except ValueError:
            continue
    return False

print("Finding higher taxon matches...")
merged_df['higher_taxon_matches'] = merged_df.progress_apply(
    lambda row: belongs_to_higher_taxon_vectorized(row['staxids'], row['TaxaID'], taxa_lineage_df),
    axis=1
)

# Step 5: Prioritize best match
def prioritize_matches(group):
    lca_matches = group[group['lca_matches']]
    if not lca_matches.empty:
        return lca_matches.iloc[0:1]
    higher_taxon_matches = group[group['higher_taxon_matches']]
    if not higher_taxon_matches.empty:
        return higher_taxon_matches.iloc[0:1]
    return group.sort_values(by=['bitscore', 'pident', 'evalue'], ascending=[False, False, True]).iloc[0:1]

print("Applying prioritization to each peptide group...")
filtered_df = merged_df.groupby('qseqid').progress_apply(prioritize_matches).reset_index(drop=True)

# Step 6: Save final output
output_file = 'novoMP_BlastP_LCA_filtered_output_New.tsv'
print("Saving final output...")
filtered_df.to_csv(output_file, sep='\t', index=False)

print(f"✅ Filtering and final output saved at: {output_file}")
