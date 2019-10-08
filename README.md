# JASPAR-wordcloud
Word clouds for individual JASPAR transcription factors are generated as follows:
1) [Gene to PubMed ID (PMID) associations](ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz) were downloaded from NCBI;
2) For each PMID, the number of occurrences of each word in the abstract of that PMID were retrieved using the R package [PubMedWordcloud](http://felixfan.github.io/PubMedWordcloud/) (punctuation and “stop words” were not considered);
3) For each JASPAR taxon (except for urochordates), the [inverse document frequency (IDF)](https://en.wikipedia.org/wiki/Tf%E2%80%93idf#Inverse_document_frequency) of each word was computed as the number of transcription factors in that taxon over the number of transcription factors (in that taxon) associated with an abstract comprising that word;
4) For each transcription factor, the [term frequency-IDF (TF-IDF)](https://en.wikipedia.org/wiki/Tf%E2%80%93idf#Term_frequency%E2%80%93Inverse_document_frequency) of each word was calculated as the number of abstracts associated with that transcription factor comprising that word over the number of words associated with that transcription factor and multiplied by the taxon-specific IDF of that word; and
6) Transcription factor word clouds were generated from the ranks of the top 50 words (in terms of TF-IDF).
Note that words that resembled the gene name of the TF (or one of its synonyms) and redundant words (e.g. “insulator” and “insulator-binding”) were removed using the Python module [fuzzywuzzy](https://github.com/seatgeek/fuzzywuzzy).

## Dependencies
```
install.packages('PubMedWordcloud', repos='http://cran.us.r-project.org')
```
