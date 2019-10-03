library(PubMedWordcloud)
library(tm)
stopwords_regex = paste(stopwords('en'), collapse = '\\b|\\b')
stopwords_regex = paste('\\b', stopwords_regex, '\\b')

args = commandArgs(trailingOnly=TRUE)
pmid = args[1]

word_count_per_pmid<-function(pmid){
    abstracts=getAbstracts(pmid)
    document<-paste(abstracts, collapse=' ')
    document = tolower(document)
    document <- gsub("[^[:alnum:][:space:]-]", "", document)
    document <- gsub("-", "X", document)
    document<-gsub("\\b\\d+\\b", "", document)
    document<-gsub('\\b\\w{2,3}\\b','',document)
    document<-gsub('\\b\\w{1}\\b','',document)
    document <- gsub("X", "-", document)
    document <- Corpus(VectorSource(document))
    document = stringr::str_replace_all(document, stopwords_regex, '')
    one_gene_word_counts <- as.data.frame(table(unlist( strsplit(document, "\ ") )))  # split vector by space
    one_gene_word_counts <- with(one_gene_word_counts, one_gene_word_counts[ Var1 != "", ] )
    name<-paste(pmid, ".rds", sep="")
    saveRDS(one_gene_word_counts, file=name)
}

word_count_per_pmid(pmid)