library(doParallel)
library(foreach)
library(tm)
library(PubMedWordcloud)

args = commandArgs(trailingOnly=TRUE)
uniacc = args[1]
max_words = as.integer(args[2])

plot_word_cloud<-function(uniacc){
    rdsfile=paste(uniacc,".rds", sep="")
    wordcloud=paste(uniacc,".png", sep="")
    png(wordcloud, width=500, height=500)
    rdsdata<-readRDS(rdsfile)
    tfidfs=data.frame(rdsdata$Var1, rdsdata$tfidf)
    colnames(tfidfs)<-c("word", "freq")
    tfidfs<-tfidfs[order(tfidfs$freq, decreasing=TRUE),]
    iranks <- max_words - c(1:nrow(tfidfs))
    tfidfs$freq<-iranks
    plotWordCloud(tfidfs[1:nrow(tfidfs),], min.freq=0, scale=c(3, 1), max.words=max_words)
    dev.off()
}

plot_word_cloud(uniacc)