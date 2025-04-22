##### Formatting training, testing, train/test.idx, clinical data
rm(list = ls())
setwd('...') # setup directory
load('CESC_LUAD_combined2472.rdata')
library(glmnet)
library(cluster)
library(cmmp)
library(nlme)
library(xtable)
library(gplots)
library(corrplot)

#corrplot(cor(snp.data[train.idx,])) #
remvindx=unique(which(abs(cor(snp.data[train.idx,]))>.7,arr.ind = T))
excludeindx=NULL
for(i in 1:dim(remvindx)[1]){
  if(remvindx[i,1]==remvindx[i,2]){
    excludeindx=c(excludeindx,i)
  }
  remvindx[i,]=sort(remvindx[i,])
}
remvindx=remvindx[-excludeindx,]
remvindx=unique(remvindx)
remvindx=unique(remvindx[,1])
#corrplot(cor(snp.data[train.idx,-remvindx]))

snp.data_filt<-snp.data[,-remvindx]
dim(snp.data_filt)

which(is.na(met.data),arr.ind = T)
met.data=met.data[,-which(is.na(met.data),arr.ind = T)[,2]]
dim(met.data)

pam1 <- function(x,k) list(cluster = pam(x,k, cluster.only=TRUE))
set.seed(1315)
gap2 <- clusGap(met.data[train.idx,],pam1,15)
gap2 #suggesting K=9
plot(gap2)
set.seed(1111)
groups=kmeans(met.data,9)
train2 <- data.frame(met.data[train.idx,],snp.data_filt[train.idx,],race=race.data[train.idx],type=type.data[train.idx],groupid=groups$cluster[train.idx])
test2 <- data.frame(met.data[test.idx,],snp.data_filt[test.idx,],race=race.data[test.idx],type=type.data[test.idx],groupid=groups$cluster[test.idx])

#plot race distribution in clusters
traintab=table(train2$race,train2$groupid)
testtab=table(test2$race,test2$groupid)
all <- traintab+testtab
print(xtable(table(train2$race,train2$groupid)))
print(xtable(table(test2$race,test2$groupid)))
print(xtable(all))
nam =rownames(table(train2$race,train2$groupid))
nam[1]='black'
nam[2]='white'
#race by side:
barplot(table(train2$race,train2$groupid), col=colors()[c(90,89)], ylim = c(0,80),
        border="white",font.axis=2, beside=T, xlab="cluster", 
        args.legend = list(x = "topright",cex = 0.6),font.lab=1,main="Training - CESC, LUAD combined")
legend('topright',legend=c('Black','White'),lty=1,col=colors()[c(90,89)],cex=.8)


library(gplots)
library(tidyverse)
library(RColorBrewer)

#load('CESC12_clin.RData')
ls()
#save(train2, test2, groups,test.idx,train.idx,CESC12_clin,file = 'data_add_clin.RData')

