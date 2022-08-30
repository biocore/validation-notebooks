## Spike-in - Absolute abundances
A validation notebook for quantifying absolute abundances for shotgun metagenomics data, using plasmid spike-ins with synthetic sequences.

### Inputs
- A feature-table with per-sample counts of synthetic sequences (table_plasmid_sequence_hits.txt)
- A feature-table with your community data (table_community_hits.txt)
- A sample metadata file indicating which plasmid pools and at which dilutions were spike-in to each sample (metadata_samples_plasmid_sequences.txt)
- A metadata file indicating the plasmids present in each pool and respective dilutions (metadata_dilutions_pool.tsv)
- A feature metadata file indicating species names and genome sizes (bp) for each feature (metadata_features.tsv)

### References
This code was adapted from Dr. Livia Zaramela (zaramela at gmail dot com):
https://github.com/lzaramela/SynDNA
