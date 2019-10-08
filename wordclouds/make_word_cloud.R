library(doParallel)
library(foreach)
library(tm)
library(PubMedWordcloud)

args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_file = args[2]
max_words = as.integer(args[3])

wordcloud=paste(output_file, sep="")
png(output_file, width=500, height=500)
rdsdata<-readRDS(input_file)
tfidfs=data.frame(rdsdata$Var1, rdsdata$tfidf)
colnames(tfidfs)<-c("word", "freq")
tfidfs<-tfidfs[order(tfidfs$freq, decreasing=TRUE),]
iranks <- max_words - c(1:nrow(tfidfs))
tfidfs$freq<-iranks
print(tfidfs[14838:15000,])
exit(0)
plotWordCloud(tfidfs[1:nrow(tfidfs),], min.freq=0, scale=c(3, 1), max.words=max_words)
dev.off()