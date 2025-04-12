# %%
import pandas as pd
lookup_table="/home/jake/Projects/Ectopic_recombination/Lookup_table.txt"
# Read the TSV file
df = pd.read_csv('/home/jake/Projects/Ectopic_recombination/results/coverage_results.tsv', sep='\t')

# Read the filter file and create a set of IDs to match
with open('/home/jake/Projects/Ectopic_recombination/results/Gypsys_domain_complete_LTR.txt', 'r') as f:
    filter_ids = {line.strip() for line in f}

# Create a function to extract the base ID from the Sequence ID column
def extract_base_id(sequence_id):
    # Extract the part before "_De"
    return sequence_id.split('_De')[0]

# Filter the dataframe
filtered_df = df[df['Sequence ID'].apply(extract_base_id).isin(filter_ids)]

# Save the filtered results
output_file = '/home/jake/Projects/Ectopic_recombination/results/Gypsy_coverage_results.tsv'
filtered_df.to_csv(output_file, sep='\t', index=False)

# Print the number of matches found
print(f"Found {len(filtered_df)} matching entries")


# %% [markdown]
# Add Species column

# %%
look_up = pd.read_csv(lookup_table, sep="\t")
for index, row in filtered_df.iterrows():  # Use iterrows() to iterate over rows
    gca = row["CSV file"].split('_')[1]  # Access the 'CSV file' column from the current row
    GCA = f"GCA_{gca}"
    
    # Corrected assignment using boolean indexing
    filtered_df.loc[index, "Species"] = look_up.loc[look_up["Genome_Accension"] == GCA, "Species"].values[0] if not look_up.loc[look_up["Genome_Accension"] == GCA, "Species"].empty else None

look_up = pd.read_csv(lookup_table, sep="\t")
display(filtered_df) 
out_file = '/home/jake/Projects/Ectopic_recombination/results/Gypsy_coverage_results_with_species.tsv'
filtered_df.to_csv(out_file, sep='\t', index=False)

# %%
filtered_df

# %%

print(filtered_df['Species'].value_counts())


