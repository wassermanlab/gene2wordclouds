library(doParallel)
library(foreach)
library(tm)
library(PubMedWordcloud)
taxons=c("fungi")

for (i in 1:length(taxons)){
  taxon_name = taxons[i]
 filenames = dir(path=taxon_name ,pattern="*.rds")

foreach(i=1:length(filenames)) %dopar%
  {
    library(PubMedWordcloud)
    name=paste(taxon_name,"/",filenames[i], sep="")
    entrez_id<-sapply(strsplit(filenames[i],split= "_"),'[',1)
    cloudname = paste(taxon_name,"/",entrez_id, ".png", sep="")
    png(cloudname, width = 500, height = 500)
    entrez_id<-readRDS(name)
    raw_output=data.frame(entrez_id$Var1,entrez_id$score)
    colnames(raw_output)<-c("word","freq")
    raw_output <- raw_output[order(raw_output$freq, decreasing=TRUE),]
    plotWordCloud(raw_output[2:nrow(raw_output),], min.freq = 0, scale = c(2, .5), max.words = 100)
    dev.off()
    } 
}
