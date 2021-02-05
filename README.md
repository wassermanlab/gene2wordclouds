<strong>Gene2WordClouds</strong>
=====

<strong>Availability:</strong>

<strong>How to cite:</strong>

Table of Contents
-----
* [Pipeline](#pipeline-overview)
  + [Directed Acyclic Graph](#dag)
* [Dependencies](#dependencies)
  + [Environment](#conda-environment-installation)
* [Usage](#usage)

Pipeline Overview
-----
1. `uniacc2entrezid.py` - convert Uniprot accession to Entrez gene ID (EntrezIDs).
2. `entrezid2pmids.py` - extract all PubMed IDs (PMIDs) associated with query EntrezIDs.
3. `gene2pmids_stats.py` - compute statistics regarding PMIDs and EntrezIDs.
4. `pmid2abstract.py` - download abstracts of the PMID papers.
5. `abstract2words.py` - tally words within each abstract, excluding punctuation and stop words. Stems are obtained for each word.
6. Compute the [Inverse Document Frequency (IDF)](https://en.wikipedia.org/wiki/Tf%E2%80%93idf#Inverse_document_frequency) for each stem.
  - Formula = total number of PMIDs in the analysis over the number of PMIDs containing the specific stem.
7. Compute the [Term Frequency-IDF](https://en.wikipedia.org/wiki/Tf%E2%80%93idf#Term_frequency%E2%80%93Inverse_document_frequency) of each stem for each transcription factor (TF).
  - Formula = number of TF-associated PMIDs containing the specific stem multiplied by the IDF of the stem.
8. `words2cloud.py` - generate a word cloud for each transcription factor with the top 50 words of unique stems.
  - Stems are used to remove redundancy in words such as "ctcf" and "ctcf-binding".



Dependencies
-----
Gene2WordClouds require the following main dependencies.

| Package | Version |
|:-:|:-:|
| biopython | 1.78 |
| click | 7.1.2 |
| click-option-group | 0.5.1 |
| distance | 0.1.3 |
| fuzzywuzzy | 0.18.0 |
| matplotlib | 3.3.3 |
| nltk | 3.4.4 |
| nltk_data | 2019.07.04 |
| numpy | 1.19.5 |
| pandas | 1.2.1 |
| python | 3.8.5 |
| scikit-learn | 0.24.1 |
| tqdm | 4.35.0 |
| wordcloud | 1.8.1 |

#### Conda Environment Installation
All packages are provided within the YML environment file and a conda environment can be created using the following command.
```bash
conda env create -f ./conda/environment.yml
```

Usage
-----
The `wordcloud.py`
