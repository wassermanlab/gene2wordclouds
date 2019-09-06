This file, for each TAXON in JASPAR (i.e. vertebrates, insects, namatodes, fungi, plants), it fetches all TFs, then transforms this TFs to ENTREZ IDs, then identifies the PMIDs associated to individual ENTREZ IDs, then call the R script to get the word counts.

For this to work, you need to:
1) iterate the taxons.
2) create a folder of the taxons.
3) iterate through the individual TFs of each taxon.
4) execute the R script for each PMID associated to each TF.

Files that you need to download and store as JSONs for repeated use:

Since JASPAR stores UniProt IDs, a UniProt ID to Entrez ID mapping.
Since PMIDs are not going to change, download from NCBI the Entrez ID to PMID mapping.
