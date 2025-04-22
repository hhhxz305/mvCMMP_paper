load(".../running.rdata") # load data
ls()
 memory.limit(size=225000)
library(glmnet)
library(cluster)
library(cmmp)
library(nlme)
library(xtable)
library(gplots)
library(corrplot)
setwd("e:/realdata")
source('mvCMMP_asreml.R')
source('convert.R')

dim(met.data)
colnames(met.data)
colnames(snp.data)
split.num=20
seed=5599#1516
set.seed(seed)

collector=list()
#need to save coefficients need row is outcome, col is predictor
mvcoef=list()
mvcoef$m=list()
mvcoef$beta=list()
mvcoef$rau=list()

for (z in 1:split.num){


   
	train.idx=sample(1:dim(met.data)[1],round(.7*dim(met.data)[1]))

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
dim(snp.data_filt) #631  61

which(is.na(met.data),arr.ind = T)
if(length(which(is.na(met.data),arr.ind = T)[,2])!=0){
met.data2=met.data[,-which(is.na(met.data),arr.ind = T)[,2]]} else{met.data2=met.data}
                                                                                                                                                                                                                                                                                      
pam1 <- function(x,k) list(cluster = pam(x,k, cluster.only=TRUE))
#set.seed(1315)
gap2 <- clusGap(met.data2[train.idx,],pam1,15)
gap2 #suggesting K=5, may want to make this automatic
plot(gap2)
#set.seed(1111)
nc <- maxSE(f = gap2$Tab[, "gap"], SE.f = gap2$Tab[, "SE.sim"], method = "firstSEmax", SE.factor = 1)
groups=kmeans(met.data2,nc)


#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################


train2 <- data.frame(met.data2[train.idx,],snp.data_filt[train.idx,],race=race.data[train.idx],groupid=groups$cluster[train.idx])
test2 <- data.frame(met.data2[-train.idx,],snp.data_filt[-train.idx,],race=race.data[-train.idx],groupid=groups$cluster[-train.idx])

CESC12_clin2=as.data.frame(CESC12_clin)
CESC12_clin2$Age=as.numeric(CESC12_clin2$yearstobirth)
CESC12_clin2$Gender=as.factor(CESC12_clin2$gender)
CESC12_clin2$pathologicstage=as.factor(CESC12_clin2$pathologicstage)
CESC12_clin2$pathologyTstage=as.factor(CESC12_clin2$pathologyTstage)
CESC12_clin2$pathologyNstage=as.factor(CESC12_clin2$pathologyNstage)
CESC12_clin2$pathologyMstage=as.factor(CESC12_clin2$pathologyMstage)
sum(!complete.cases(CESC12_clin2[,8:10])) #lost 49
table(CESC12_clin2$pathologyTstage,useNA = 'always')
PathoT=rep(NA,dim(train2)[1])
PathoT[CESC12_clin2$pathologyTstage=='t1'|CESC12_clin2$pathologyTstage=='t1a'|CESC12_clin2$pathologyTstage=='t1a1'|CESC12_clin2$pathologyTstage=='t1b'|CESC12_clin2$pathologyTstage=='t1b1'|CESC12_clin2$pathologyTstage=='t1b2'|CESC12_clin2$pathologyTstage=='tis']='t1'
PathoT[CESC12_clin2$pathologyTstage=='t2'|CESC12_clin2$pathologyTstage=='t2a'|CESC12_clin2$pathologyTstage=='t2a1'|CESC12_clin2$pathologyTstage=='t2a2'|CESC12_clin2$pathologyTstage=='t2b']='t2'
PathoT[CESC12_clin2$pathologyTstage=='t3'|CESC12_clin2$pathologyTstage=='t3a'|CESC12_clin2$pathologyTstage=='t3b']='t3'
PathoT[CESC12_clin2$pathologyTstage=='t4']='t4'
PathoT[CESC12_clin2$pathologyTstage=='tx']='tx'

PathoM=rep(NA,dim(train2)[1])
PathoM[CESC12_clin2$pathologyMstage=='m0']='m0'
PathoM[CESC12_clin2$pathologyMstage=='m1'|CESC12_clin2$pathologyMstage=='m1a'|CESC12_clin2$pathologyMstage=='m1b']='m0'
PathoM[CESC12_clin2$pathologyMstage=='mx']='mx'

PathoN=CESC12_clin2$pathologyNstage
PathoN[PathoN=='n3']=NA #exlucding single person category
PathoN=droplevels(PathoN)

CESC12_clin3=data.frame(Age=CESC12_clin2$Age,Gender=CESC12_clin2$Gender,PathoT=PathoT,PathoN=PathoN,PathoM=PathoM)

comptrain=complete.cases(train2)
comptest=complete.cases(test2)


###################################CESC stage variable##############################
dim(CESC12_clin2)
dim(all)
table(types)
table(CESC12_clin2$pathologicstage)
table(types,CESC12_clin2$pathologyTstage)
table(types,CESC12_clin2$pathologyNstage)
table(types,CESC12_clin2$pathologyMstage)
table(CESC12_clin2$pathologyTstage,CESC12_clin2$pathologyMstage)
newstage=rep(NA,dim(train2)[1])
str(CESC12_clin2$pathologicstage[types=='LUAD'])
newstage[types=='LUAD']=as.character(CESC12_clin2$pathologicstage[types=='LUAD'])
newstage[types=='CESC'&(CESC12_clin2$pathologyMstage=='mx'|
                          CESC12_clin2$pathologyNstage=='nx'|
                          CESC12_clin2$pathologyTstage=='tx'|
                          is.na(CESC12_clin2$pathologyMstage)|
                          is.na(CESC12_clin2$pathologyTstage))]='sx'
newstage[types=='CESC'&is.na(newstage)&(CESC12_clin2$pathologyMstage=='m1'|
                                          CESC12_clin2$pathologyMstage=='m1a'|
                                          CESC12_clin2$pathologyMstage=='m1b'|
                                          CESC12_clin2$pathologyTstage=='t4')]='s4'
newstage[types=='CESC'&is.na(newstage)&(CESC12_clin2$pathologyTstage=='t3'|
                                          CESC12_clin2$pathologyTstage=='t3a'|
                                          CESC12_clin2$pathologyTstage=='t3b'|
                                          CESC12_clin2$pathologyNstage=='n1')]='s3'
newstage[types=='CESC'&is.na(newstage)&(CESC12_clin2$pathologyTstage=='t2'|
                                          CESC12_clin2$pathologyTstage=='t2a'|
                                          CESC12_clin2$pathologyTstage=='t2a1'|
                                          CESC12_clin2$pathologyTstage=='t2a2'|
                                          CESC12_clin2$pathologyTstage=='t2b')]='s2'
newstage[types=='CESC'&is.na(newstage)&(CESC12_clin2$pathologyTstage=='t1'|
                                          CESC12_clin2$pathologyTstage=='t1a'|
                                          CESC12_clin2$pathologyTstage=='t1a1'|
                                          CESC12_clin2$pathologyTstage=='t1b'|
                                          CESC12_clin2$pathologyTstage=='t1b1'|
                                          CESC12_clin2$pathologyTstage=='t1b2')]='s1'
New.Stage=newstage
New.Stage[newstage=='s1'|newstage=='stage i'|newstage=='stage ia'|newstage=='stage ib']='Stage.1'
New.Stage[newstage=='s2'|newstage=='stage ii'|newstage=='stage iia'|newstage=='stage iib']='Stage.2'
New.Stage[newstage=='s3'|newstage=='stage iiia'|newstage=='stage iiib']='Stage.3'
New.Stage[newstage=='s4'|newstage=='stage iv']='Stage.4'
New.Stage[newstage=='sx'|is.na(newstage)]='Stage.NA'


CESC12_clin3=data.frame(Age=CESC12_clin2$Age,Gender=CESC12_clin2$Gender,Stage=New.Stage)

train.clic=CESC12_clin3[train.idx[comptrain],]
train2=data.frame(train2,train.clic)
train2=train2[complete.cases(train2),]
test2.clic=CESC12_clin3[-train.idx,][types[-train.idx]=='CESC',]
test2=data.frame(test2[-train.idx,][types[-train.idx]=='CESC',],test2.clic)
test3=test2#[test2$type=='CESC',]
test2=test2[complete.cases(test2),]
test3=test3[complete.cases(test3),]
testrace=test3$race

races=unique(train2$race)
races #races[1]=white,races[2]=black

ty_train2=unlist(lapply(train2,class))
ty_test3=unlist(lapply(test3,class))

for(i in which(ty_train2=='character')){
train2[,i]=as.factor(train2[,i])
}

for(i in which(ty_test3=='character')){
test3[,i]=as.factor(test3[,i])
}


#note matrix is overall, white, black
cmmpcol=matrix(NA,nrow=2402,ncol=3)
mvcmmpcol=matrix(NA,nrow=2402,ncol=3)
lmcol=matrix(NA,nrow=2402,ncol=3)
rfcol=matrix(NA,nrow=2402,ncol=3)

cmmpcol_m=matrix(NA,nrow=2402,ncol=3)
mvcmmpcol_m=matrix(NA,nrow=2402,ncol=3)
lmcol_m=matrix(NA,nrow=2402,ncol=3)
rfcol_m=matrix(NA,nrow=2402,ncol=3)

cmmpcol_rau=matrix(NA,nrow=2402,ncol=3)
mvcmmpcol_rau=matrix(NA,nrow=2402,ncol=3)
lmcol_rau=matrix(NA,nrow=2402,ncol=3)
rfcol_rau=matrix(NA,nrow=2402,ncol=3)

colnames(cmmpcol)=colnames(mvcmmpcol)=colnames(lmcol)=colnames(rfcol)=
colnames(cmmpcol_m)=colnames(mvcmmpcol_m)=colnames(lmcol_m)=colnames(rfcol_m)=
colnames(cmmpcol_rau)=colnames(mvcmmpcol_rau)=colnames(lmcol_rau)=colnames(rfcol_rau)=c('all','white','black')


#blup_cmmp=as.data.frame(matrix(NA,nrow=36,ncol=5))
#names(blup_cmmp)=unique(train2$groupid)
#memcol=list()
#rancol=list()

convert_m=function(p){
p2=log(p/(1-p),2)
p2[p2==Inf|p2==-Inf]=NA
p2
}

convert_rau=function(p){
 t = 2 * asin(sqrt(p))
(46.47324337*t) - 23
}

# part mv_cmmp



module.hold = list()
module.hold_3 = list()
module.hold_5 = list()


module.hold_m = list()
module.hold_3_m = list()
module.hold_5_m = list()


module.hold_r = list()
module.hold_3_r = list()
module.hold_5_r = list()

#transforming data
train2r=train2
train2r[,1:2402]=apply(train2r[,1:2402],2,convert_rau)
train2r=train2r[complete.cases(train2r),]

test3r=test3
test3r[,1:2402]=apply(test3r[,1:2402],2,convert_rau)
test3r=test3r[complete.cases(test3r),]

train2m=train2
train2m[,1:2402]=apply(train2m[,1:2402],2,convert_m)
train2m=train2m[complete.cases(train2m),]

test3m=test3
test3m[,1:2402]=apply(test3m[,1:2402],2,convert_m)
testrace_m=testrace[complete.cases(test3m)]
test3m=test3m[complete.cases(test3m),]


#beta

tt = train2[,1:2402]
#names(tt) <- substring(names(tt), 4)
mems=cutree(hclust(1-abs(as.dist(cor((tt))))),k=75)
avgcor=rep(NA,length(table(mems)))
table(mems)
#clus=53



	mvcoef$beta[[z]]=list()

for (jj in 1:dim(table(mems))) {#
  clus = jj
  
  avgcor[jj]=mean(abs(cor(train2[,which(mems==clus),drop=F])))
  trainc=data.frame((train2[,which(mems==clus),drop=F]),train2[,2403:dim(train2)[2]])
  indxout=which(mems==clus)
  

  
  #if more that 8 outcomes need to split it
  if(length(indxout)>8){
    
    mixindx=sample(indxout,length(indxout))
    
    outs2=list()
    for(i in 1:(ceiling(length(indxout)/6))){
      if((i*6)<=length(indxout)){
        outs2[[i]]=mixindx[((i-1)*6+1):(i*6)]}else{
          outs2[[i]]=mixindx[((i-1)*6+1):(length(indxout))]
        }
    }
    if(length(outs2[[length(outs2)]])<4){
      outs2[[length(outs2)-1]]=c(outs2[[length(outs2)-1]],outs2[[length(outs2)]])
      outs2[[length(outs2)]]=NULL
    }
    
    outs=which(mems==clus)
    holdermvcmmp=holdercmmp=data.frame(matrix(NA,nrow=nrow(test3),ncol=length(which(mems==clus))))
   names(holdercmmp)=names(holdermvcmmp)=names(outs)
    
    
    for(i in 1:length(outs2)){
      
      f.in2 = as.formula(paste0('cbind(',paste0(names(mems[outs2[[i]]]),collapse = ','),')~', paste(names(train2)[(2403:dim(train2)[2])[-which(names(train2[,2403:dim(train2)[2]])=='groupid')]],collapse = '+')))
      w<-NULL
      try(w<-mvCMMP(f.in2,~groupid,trainc,y.new=test3[,outs2[[i]]],x.new=test3[,2403:(dim(test3)[2]),drop=F],tolparinv=1e-5,match.train =F,clusts=1,comps=2))
      if(!is.null(w)){
        holdermvcmmp[,names(outs2[[i]])]=w$thetahat
		mvcoef$beta[[z]][[length(mvcoef$beta[[z]])+1]]=w$mod[[1]]}

	
    }
    
    for(i in 1:length(which(mems==clus))){
      f.in=as.formula(paste0(names(mems)[outs[i]],'~', paste(names(train2)[(2403:dim(train2)[2])[-which(names(train2[,2403:dim(train2)[2]])=='groupid')]],collapse = '+')))
      h=cmmp(f.in = f.in,r.in = ~1|groupid
             ,train=trainc,y.new=test3[,outs[i]],x.new=test3[,2403:dim(test3)[2]],match.train = F,x.fut = NULL,interval = T)
      holdercmmp[,i]=h$mixed.pred
    }
    
    
        lasthold=list(cmmp=colMeans((test3[,outs]-holdercmmp)^2),
                  mvcmmp=colMeans((test3[,outs]-holdermvcmmp)^2))
    module.hold[[jj]]=lasthold
    
    lasthold_3=list(cmmp=colMeans((test3[testrace==3,outs]-holdercmmp[testrace==3,])^2),
                  mvcmmp=colMeans((test3[testrace==3,outs]-holdermvcmmp[testrace==3,])^2))
    module.hold_3[[jj]]=lasthold_3
    
    lasthold_5=list(cmmp=colMeans((test3[testrace==5,outs]-holdercmmp[testrace==5,])^2),
                  mvcmmp=colMeans((test3[testrace==5,outs]-holdermvcmmp[testrace==5,])^2))
    module.hold_5[[jj]]=lasthold_5
    
  }else if (length(indxout)<=8&length(indxout)>1){
    
    outs=which(mems==clus)
    holdermvcmmp=holdercmmp=data.frame(matrix(NA,nrow=nrow(test3),ncol=length(which(mems==clus))))
   names(holdercmmp)=names(holdermvcmmp)=names(outs)
    
    outs2=list()
    outs2=indxout
    f.in2 = as.formula(paste0('cbind(',paste0(names(mems[mems==clus]),collapse = ','),')~', paste(names(train2)[(2403:dim(train2)[2])[-which(names(train2[,2403:dim(train2)[2]])=='groupid')]],collapse = '+')))
    w<-NULL
    try(w<-mvCMMP(f.in2,~groupid,trainc,y.new=test3[,which(mems==clus)],x.new=test3[,2403:dim(test3)[2],drop=F],tolparinv=1e-5,match.train =F,clusts=1,comps=2))
    if(!is.null(w)){
      holdermvcmmp[,names(mems[mems==clus])]=w$thetahat
	mvcoef$beta[[z]][[length(mvcoef$beta[[z]])+1]]=w$mod[[1]]}
    
	

	



	#holdercmmp=matrix(NA,nrow=nrow(test3),ncol=length(which(mems==clus)))
    #outs=which(mems==clus)
    for(i in 1:length(which(mems==clus))){
      f.in=as.formula(paste0(names(mems)[outs[i]],'~', paste(names(train2)[2403:(dim(train2)[2]-1)],collapse = '+')))
      h=cmmp(f.in = f.in,r.in = ~1|groupid
             ,train=trainc,y.new=test3[,outs[i]],x.new=test3[,2403:dim(test3)[2]],match.train = F,x.fut = NULL,interval = T)
      holdercmmp[,i]=h$mixed.pred
    }
    
    lasthold=list(cmmp=colMeans((test3[,outs]-holdercmmp)^2),
                  mvcmmp=colMeans((test3[,outs]-holdermvcmmp)^2))
    module.hold[[jj]]=lasthold
    
    lasthold_3=list(cmmp=colMeans((test3[testrace==3,outs]-holdercmmp[testrace==3,])^2),
                  mvcmmp=colMeans((test3[testrace==3,outs]-holdermvcmmp[testrace==3,])^2))
    module.hold_3[[jj]]=lasthold_3
    
    lasthold_5=list(cmmp=colMeans((test3[testrace==5,outs]-holdercmmp[testrace==5,])^2),
                  mvcmmp=colMeans((test3[testrace==5,outs]-holdermvcmmp[testrace==5,])^2))
    module.hold_5[[jj]]=lasthold_5
    
    
  }
  
}


#to combine
cmmphold=NULL
mvcmmphold=NULL


cmmphold_3=NULL
mvcmmphold_3=NULL


cmmphold_5=NULL
mvcmmphold_5=NULL

for(i in 1:length(module.hold)){
  cmmphold=c(cmmphold,module.hold[[i]]$cmmp)
  mvcmmphold=c(mvcmmphold,module.hold[[i]]$mvcmmp)
  
  
  cmmphold_3=c(cmmphold_3,module.hold_3[[i]]$cmmp)
  mvcmmphold_3=c(mvcmmphold_3,module.hold_3[[i]]$mvcmmp)
  
  
  cmmphold_5=c(cmmphold_5,module.hold_5[[i]]$cmmp)
  mvcmmphold_5=c(mvcmmphold_5,module.hold_5[[i]]$mvcmmp)
  
}
errhold=rbind(cmmphold,mvcmmphold)

errhold_3=rbind(cmmphold_3,mvcmmphold_3)

errhold_5=rbind(cmmphold_5,mvcmmphold_5)

library(randomForestSRC)
f.rf=as.formula(paste0('cbind(',paste0(names(train2)[1:2402],collapse = ','),')~', paste(names(train2)[(2403:dim(train2)[2])[-which(names(train2[,2403:dim(train2)[2]])=='groupid')]],collapse = '+')))
rf=rfsrc(f.rf,data=train2)
lmod=lm(f.rf,data=train2)

predrf=predict(rf,newdata=test3)
predrf_3=predict(rf,newdata=test3[testrace==3,])
predrf_5=predict(rf,newdata=test3[testrace==5,])
prelm=predict(lmod,newdata=test3)
rferrcol=NULL
rferrcol_3=NULL
rferrcol_5=NULL
for(i in 1:length(predrf$regrOutput)){
  rferrcol=c(rferrcol,predrf$regrOutput[[i]]$err.rate[length(predrf$regrOutput[[i]]$err.rate)])
  rferrcol_3=c(rferrcol_3,predrf$regrOutput[[i]]$err.rate[length(predrf$regrOutput[[i]]$err.rate)])
  rferrcol_5=c(rferrcol_5,predrf$regrOutput[[i]]$err.rate[length(predrf$regrOutput[[i]]$err.rate)])
}



lmerrcol=colMeans((prelm-test3[,1:2402])^2)
lmerrcol3=colMeans((prelm-test3[,1:2402])[testrace==3,]^2)
lmerrcol5=colMeans((prelm-test3[,1:2402])[testrace==5,]^2)


errhold2=data.frame(name=colnames(errhold),cmmp=errhold[1,],mvcmmp=errhold[2,])
errhold2_3=data.frame(name=colnames(errhold_3),cmmp=errhold_3[1,],mvcmmp=errhold_3[2,])
errhold2_5=data.frame(name=colnames(errhold_5),cmmp=errhold_5[1,],mvcmmp=errhold_5[2,])
rferrcol2=data.frame(name=names(predrf$regrOutput),rf=rferrcol)
rferrcol2_3=data.frame(name=names(predrf$regrOutput),rf=rferrcol_3)
rferrcol2_5=data.frame(name=names(predrf$regrOutput),rf=rferrcol_5)
lmerrcol2=data.frame(name=names(train2)[1:2402],mvlm=lmerrcol)
lmerrcol2_3=data.frame(name=names(train2)[1:2402],mvlm=lmerrcol3)
lmerrcol2_5=data.frame(name=names(train2)[1:2402],mvlm=lmerrcol5)

errhold3=merge(merge(errhold2,rferrcol2,'name',all=T),lmerrcol2,'name',all=T)
names(errhold3)[5]='lm'



errhold3_3=merge(merge(errhold2_3,rferrcol2_3,'name',all=T),lmerrcol2_3,'name',all=T)
names(errhold3_3)[5]='lm'

errhold3_5=merge(merge(errhold2_5,rferrcol2_5,'name',all=T),lmerrcol2_5,'name',all=T)
names(errhold3_5)[5]='lm'


##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################

#m

tt = train2m[,1:2402]
#names(tt) <- substring(names(tt), 4)
mems=cutree(hclust(1-abs(as.dist(cor((tt))))),k=75)
avgcor=rep(NA,length(table(mems)))
table(mems)
#clus=53



	mvcoef$m[[z]]=list()
for (jj in 1:dim(table(mems))) {
  clus = jj
  
  avgcor[jj]=mean(abs(cor(train2m[,which(mems==clus),drop=F])))
  trainc=data.frame((train2m[,which(mems==clus),drop=F]),train2m[,2403:dim(train2m)[2]])
  indxout=which(mems==clus)
  
  #if more that 8 outcomes need to split it
  if(length(indxout)>8){
    
    mixindx=sample(indxout,length(indxout))
    
    outs2=list()
    for(i in 1:(ceiling(length(indxout)/6))){
      if((i*6)<=length(indxout)){
        outs2[[i]]=mixindx[((i-1)*6+1):(i*6)]}else{
          outs2[[i]]=mixindx[((i-1)*6+1):(length(indxout))]
        }
    }
    if(length(outs2[[length(outs2)]])<4){
      outs2[[length(outs2)-1]]=c(outs2[[length(outs2)-1]],outs2[[length(outs2)]])
      outs2[[length(outs2)]]=NULL
    }
    
    outs=which(mems==clus)
    holdermvcmmp=holdercmmp=data.frame(matrix(NA,nrow=nrow(test3m),ncol=length(which(mems==clus))))
   names(holdercmmp)=names(holdermvcmmp)=names(outs)
    
    
    for(i in 1:length(outs2)){
      
      f.in2 = as.formula(paste0('cbind(',paste0(names(mems[outs2[[i]]]),collapse = ','),')~', paste(names(train2m)[(2403:dim(train2m)[2])[-which(names(train2m[,2403:dim(train2m)[2]])=='groupid')]],collapse = '+')))
      w<-NULL
      try(w<-mvCMMP(f.in2,~groupid,trainc,y.new=test3m[,outs2[[i]]],x.new=test3m[,2403:(dim(test3m)[2]),drop=F],tolparinv=1e-5,match.train =F,clusts=1,comps=2))
      if(!is.null(w)){
        holdermvcmmp[,names(outs2[[i]])]=w$thetahat
mvcoef$m[[z]][[length(mvcoef$m[[z]])+1]]=w$mod[[1]]}

	
    }
    
    for(i in 1:length(which(mems==clus))){
      f.in=as.formula(paste0(names(mems)[outs[i]],'~', paste(names(train2m)[(2403:dim(train2m)[2])[-which(names(train2m[,2403:dim(train2m)[2]])=='groupid')]],collapse = '+')))
      h=cmmp(f.in = f.in,r.in = ~1|groupid
             ,train=trainc,y.new=test3m[,outs[i]],x.new=test3m[,2403:dim(test3m)[2]],match.train = F,x.fut = NULL,interval = T)
      holdercmmp[,i]=h$mixed.pred
    }
    
        lasthold=list(cmmp=colMeans((test3m[,outs]-holdercmmp)^2),
                  mvcmmp=colMeans((test3m[,outs]-holdermvcmmp)^2))
    module.hold_m[[jj]]=lasthold
    
    lasthold_3=list(cmmp=colMeans((test3m[testrace_m==3,outs]-holdercmmp[testrace_m==3,])^2),
                  mvcmmp=colMeans((test3m[testrace_m==3,outs]-holdermvcmmp[testrace_m==3,])^2))
    module.hold_3_m[[jj]]=lasthold_3
    
    lasthold_5=list(cmmp=colMeans((test3m[testrace_m==5,outs]-holdercmmp[testrace_m==5,])^2),
                  mvcmmp=colMeans((test3m[testrace_m==5,outs]-holdermvcmmp[testrace_m==5,])^2))
    module.hold_5_m[[jj]]=lasthold_5
    
    
  }else if (length(indxout)<=8&length(indxout)>1){
    
    outs=which(mems==clus)
    holdermvcmmp=holdercmmp=data.frame(matrix(NA,nrow=nrow(test3m),ncol=length(which(mems==clus))))
   names(holdercmmp)=names(holdermvcmmp)=names(outs)
    
    outs2=list()
    outs2=indxout
    f.in2 = as.formula(paste0('cbind(',paste0(names(mems[mems==clus]),collapse = ','),')~', paste(names(train2m)[(2403:dim(train2m)[2])[-which(names(train2m[,2403:dim(train2m)[2]])=='groupid')]],collapse = '+')))
    w<-NULL
    try(w<-mvCMMP(f.in2,~groupid,trainc,y.new=test3m[,which(mems==clus)],x.new=test3m[,2403:dim(test3m)[2],drop=F],tolparinv=1e-5,match.train =F,clusts=1,comps=2))
    if(!is.null(w)){
      holdermvcmmp[,names(mems[mems==clus])]=w$thetahat
mvcoef$m[[z]][[length(mvcoef$m[[z]])+1]]=w$mod[[1]]}
    
	

	



	#holdercmmp=matrix(NA,nrow=nrow(test3m),ncol=length(which(mems==clus)))
    #outs=which(mems==clus)
    for(i in 1:length(which(mems==clus))){
      f.in=as.formula(paste0(names(mems)[outs[i]],'~', paste(names(train2m)[2403:(dim(train2m)[2]-1)],collapse = '+')))
      h=cmmp(f.in = f.in,r.in = ~1|groupid
             ,train=trainc,y.new=test3m[,outs[i]],x.new=test3m[,2403:dim(test3m)[2]],match.train = F,x.fut = NULL,interval = T)
      holdercmmp[,i]=h$mixed.pred
    }
    
    lasthold=list(cmmp=colMeans((test3m[,outs]-holdercmmp)^2),
                  mvcmmp=colMeans((test3m[,outs]-holdermvcmmp)^2))
    module.hold_m[[jj]]=lasthold
    
    lasthold_3=list(cmmp=colMeans((test3m[testrace_m==3,outs]-holdercmmp[testrace_m==3,])^2),
                  mvcmmp=colMeans((test3m[testrace_m==3,outs]-holdermvcmmp[testrace_m==3,])^2))
    module.hold_3_m[[jj]]=lasthold_3
    
    lasthold_5=list(cmmp=colMeans((test3m[testrace_m==5,outs]-holdercmmp[testrace_m==5,])^2),
                  mvcmmp=colMeans((test3m[testrace_m==5,outs]-holdermvcmmp[testrace_m==5,])^2))
    module.hold_5_m[[jj]]=lasthold_5
    
    
  }
  
}


#to combine
cmmphold=NULL
mvcmmphold=NULL


cmmphold_3=NULL
mvcmmphold_3=NULL


cmmphold_5=NULL
mvcmmphold_5=NULL

for(i in 1:length(module.hold_m)){
  cmmphold=c(cmmphold,module.hold_m[[i]]$cmmp)
  mvcmmphold=c(mvcmmphold,module.hold_m[[i]]$mvcmmp)
  
  
  cmmphold_3=c(cmmphold_3,module.hold_3_m[[i]]$cmmp)
  mvcmmphold_3=c(mvcmmphold_3,module.hold_3_m[[i]]$mvcmmp)
  
  
  cmmphold_5=c(cmmphold_5,module.hold_5_m[[i]]$cmmp)
  mvcmmphold_5=c(mvcmmphold_5,module.hold_5_m[[i]]$mvcmmp)
  
}
errhold=rbind(cmmphold,mvcmmphold)

errhold_3=rbind(cmmphold_3,mvcmmphold_3)

errhold_5=rbind(cmmphold_5,mvcmmphold_5)

library(randomForestSRC)
f.rf=as.formula(paste0('cbind(',paste0(names(train2m)[1:2402],collapse = ','),')~', paste(names(train2m)[(2403:dim(train2m)[2])[-which(names(train2m[,2403:dim(train2m)[2]])=='groupid')]],collapse = '+')))
rf=rfsrc(f.rf,data=train2m)
lmod=lm(f.rf,data=train2m)

predrf=predict(rf,newdata=test3m)
predrf_3=predict(rf,newdata=test3m[testrace_m==3,])
predrf_5=predict(rf,newdata=test3m[testrace_m==5,])
prelm=predict(lmod,newdata=test3m)
rferrcol=NULL
rferrcol_3=NULL
rferrcol_5=NULL
for(i in 1:length(predrf$regrOutput)){
  rferrcol=c(rferrcol,predrf$regrOutput[[i]]$err.rate[length(predrf$regrOutput[[i]]$err.rate)])
  rferrcol_3=c(rferrcol_3,predrf$regrOutput[[i]]$err.rate[length(predrf$regrOutput[[i]]$err.rate)])
  rferrcol_5=c(rferrcol_5,predrf$regrOutput[[i]]$err.rate[length(predrf$regrOutput[[i]]$err.rate)])
}



lmerrcol=colMeans((prelm-test3m[,1:2402])^2)
lmerrcol3=colMeans((prelm-test3m[,1:2402])[testrace_m==3,]^2)
lmerrcol5=colMeans((prelm-test3m[,1:2402])[testrace_m==5,]^2)


errhold2=data.frame(name=colnames(errhold),cmmp=errhold[1,],mvcmmp=errhold[2,])
errhold2_3=data.frame(name=colnames(errhold_3),cmmp=errhold_3[1,],mvcmmp=errhold_3[2,])
errhold2_5=data.frame(name=colnames(errhold_5),cmmp=errhold_5[1,],mvcmmp=errhold_5[2,])
rferrcol2=data.frame(name=names(predrf$regrOutput),rf=rferrcol)
rferrcol2_3=data.frame(name=names(predrf$regrOutput),rf=rferrcol_3)
rferrcol2_5=data.frame(name=names(predrf$regrOutput),rf=rferrcol_5)
lmerrcol2=data.frame(name=names(train2m)[1:2402],mvlm=lmerrcol)
lmerrcol2_3=data.frame(name=names(train2m)[1:2402],mvlm=lmerrcol3)
lmerrcol2_5=data.frame(name=names(train2m)[1:2402],mvlm=lmerrcol5)

errhold3_m=merge(merge(errhold2,rferrcol2,'name',all=T),lmerrcol2,'name',all=T)
names(errhold3_m)[5]='lm'



errhold3_3_m=merge(merge(errhold2_3,rferrcol2_3,'name',all=T),lmerrcol2_3,'name',all=T)
names(errhold3_3_m)[5]='lm'

errhold3_5_m=merge(merge(errhold2_5,rferrcol2_5,'name',all=T),lmerrcol2_5,'name',all=T)
names(errhold3_5_m)[5]='lm'

##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################

#rau

tt = train2r[,1:2402]
#names(tt) <- substring(names(tt), 4)
mems=cutree(hclust(1-abs(as.dist(cor((tt))))),k=75)
avgcor=rep(NA,length(table(mems)))
table(mems)
#clus=53



	mvcoef$rau[[z]]=list()
for (jj in 1:dim(table(mems))) {
  clus = jj
  
  avgcor[jj]=mean(abs(cor(train2r[,which(mems==clus),drop=F])))
  trainc=data.frame((train2r[,which(mems==clus),drop=F]),train2r[,2403:dim(train2m)[2]])
  indxout=which(mems==clus)
  
  #if more that 8 outcomes need to split it
  if(length(indxout)>8){
    
    mixindx=sample(indxout,length(indxout))
    
    outs2=list()
    for(i in 1:(ceiling(length(indxout)/6))){
      if((i*6)<=length(indxout)){
        outs2[[i]]=mixindx[((i-1)*6+1):(i*6)]}else{
          outs2[[i]]=mixindx[((i-1)*6+1):(length(indxout))]
        }
    }
    if(length(outs2[[length(outs2)]])<4){
      outs2[[length(outs2)-1]]=c(outs2[[length(outs2)-1]],outs2[[length(outs2)]])
      outs2[[length(outs2)]]=NULL
    }
    
    outs=which(mems==clus)
    holdermvcmmp=holdercmmp=data.frame(matrix(NA,nrow=nrow(test3r),ncol=length(which(mems==clus))))
   names(holdercmmp)=names(holdermvcmmp)=names(outs)
    
    
    for(i in 1:length(outs2)){
      
      f.in2 = as.formula(paste0('cbind(',paste0(names(mems[outs2[[i]]]),collapse = ','),')~', paste(names(train2r)[(2403:dim(train2r)[2])[-which(names(train2r[,2403:dim(train2r)[2]])=='groupid')]],collapse = '+')))
      w<-NULL
      try(w<-mvCMMP(f.in2,~groupid,trainc,y.new=test3r[,outs2[[i]]],x.new=test3r[,2403:(dim(test3r)[2]),drop=F],tolparinv=1e-5,match.train =F,clusts=1,comps=2))
      if(!is.null(w)){
        holdermvcmmp[,names(outs2[[i]])]=w$thetahat
mvcoef$rau[[z]][[length(mvcoef$rau[[z]])+1]]=w$mod[[1]]}

	
    }
    
    for(i in 1:length(which(mems==clus))){
      f.in=as.formula(paste0(names(mems)[outs[i]],'~', paste(names(train2r)[(2403:dim(train2r)[2])[-which(names(train2r[,2403:dim(train2r)[2]])=='groupid')]],collapse = '+')))
      h=cmmp(f.in = f.in,r.in = ~1|groupid
             ,train=trainc,y.new=test3r[,outs[i]],x.new=test3r[,2403:dim(test3r)[2]],match.train = F,x.fut = NULL,interval = T)
      holdercmmp[,i]=h$mixed.pred
    }
    
     lasthold=list(cmmp=colMeans((test3r[,outs]-holdercmmp)^2),
                  mvcmmp=colMeans((test3r[,outs]-holdermvcmmp)^2))
    module.hold_r[[jj]]=lasthold
    
    lasthold_3=list(cmmp=colMeans((test3r[testrace==3,outs]-holdercmmp[testrace==3,])^2),
                  mvcmmp=colMeans((test3r[testrace==3,outs]-holdermvcmmp[testrace==3,])^2))
    module.hold_3_r[[jj]]=lasthold_3
    
    lasthold_5=list(cmmp=colMeans((test3r[testrace==5,outs]-holdercmmp[testrace==5,])^2),
                  mvcmmp=colMeans((test3r[testrace==5,outs]-holdermvcmmp[testrace==5,])^2))
    module.hold_5_r[[jj]]=lasthold_5
    
    
  }else if (length(indxout)<=8&length(indxout)>1){
    
    outs=which(mems==clus)
    holdermvcmmp=holdercmmp=data.frame(matrix(NA,nrow=nrow(test3r),ncol=length(which(mems==clus))))
   names(holdercmmp)=names(holdermvcmmp)=names(outs)
    
    outs2=list()
    outs2=indxout
    f.in2 = as.formula(paste0('cbind(',paste0(names(mems[mems==clus]),collapse = ','),')~', paste(names(train2r)[(2403:dim(train2r)[2])[-which(names(train2r[,2403:dim(train2r)[2]])=='groupid')]],collapse = '+')))
    w<-NULL
    try(w<-mvCMMP(f.in2,~groupid,trainc,y.new=test3r[,which(mems==clus)],x.new=test3r[,2403:dim(test3r)[2],drop=F],tolparinv=1e-5,match.train =F,clusts=1,comps=2))
    if(!is.null(w)){
      holdermvcmmp[,names(mems[mems==clus])]=w$thetahat
mvcoef$rau[[z]][[length(mvcoef$rau[[z]])+1]]=w$mod[[1]]}
    
	

	



	#holdercmmp=matrix(NA,nrow=nrow(test3m),ncol=length(which(mems==clus)))
    #outs=which(mems==clus)
    for(i in 1:length(which(mems==clus))){
      f.in=as.formula(paste0(names(mems)[outs[i]],'~', paste(names(train2r)[2403:(dim(train2r)[2]-1)],collapse = '+')))
      h=cmmp(f.in = f.in,r.in = ~1|groupid
             ,train=trainc,y.new=test3r[,outs[i]],x.new=test3r[,2403:dim(test3r)[2]],match.train = F,x.fut = NULL,interval = T)
      holdercmmp[,i]=h$mixed.pred
    }
    
    lasthold=list(cmmp=colMeans((test3r[,outs]-holdercmmp)^2),
                  mvcmmp=colMeans((test3r[,outs]-holdermvcmmp)^2))
    module.hold_r[[jj]]=lasthold
    
    lasthold_3=list(cmmp=colMeans((test3r[testrace==3,outs]-holdercmmp[testrace==3,])^2),
                  mvcmmp=colMeans((test3r[testrace==3,outs]-holdermvcmmp[testrace==3,])^2))
    module.hold_3_r[[jj]]=lasthold_3
    
    lasthold_5=list(cmmp=colMeans((test3r[testrace==5,outs]-holdercmmp[testrace==5,])^2),
                  mvcmmp=colMeans((test3r[testrace==5,outs]-holdermvcmmp[testrace==5,])^2))
    module.hold_5_r[[jj]]=lasthold_5
    
    
  }
  
}


#to combine
cmmphold=NULL
mvcmmphold=NULL


cmmphold_3=NULL
mvcmmphold_3=NULL


cmmphold_5=NULL
mvcmmphold_5=NULL

for(i in 1:length(module.hold_r)){
  cmmphold=c(cmmphold,module.hold_r[[i]]$cmmp)
  mvcmmphold=c(mvcmmphold,module.hold_r[[i]]$mvcmmp)
  
  
  cmmphold_3=c(cmmphold_3,module.hold_3_r[[i]]$cmmp)
  mvcmmphold_3=c(mvcmmphold_3,module.hold_3_r[[i]]$mvcmmp)
  
  
  cmmphold_5=c(cmmphold_5,module.hold_5_r[[i]]$cmmp)
  mvcmmphold_5=c(mvcmmphold_5,module.hold_5_r[[i]]$mvcmmp)
  
}
errhold=rbind(cmmphold,mvcmmphold)

errhold_3=rbind(cmmphold_3,mvcmmphold_3)

errhold_5=rbind(cmmphold_5,mvcmmphold_5)

library(randomForestSRC)
f.rf=as.formula(paste0('cbind(',paste0(names(train2r)[1:2402],collapse = ','),')~', paste(names(train2r)[(2403:dim(train2r)[2])[-which(names(train2r[,2403:dim(train2r)[2]])=='groupid')]],collapse = '+')))
rf=rfsrc(f.rf,data=train2r)
lmod=lm(f.rf,data=train2r)

predrf=predict(rf,newdata=test3r)
predrf_3=predict(rf,newdata=test3r[testrace==3,])
predrf_5=predict(rf,newdata=test3r[testrace==5,])
prelm=predict(lmod,newdata=test3r)
rferrcol=NULL
rferrcol_3=NULL
rferrcol_5=NULL
for(i in 1:length(predrf$regrOutput)){
  rferrcol=c(rferrcol,predrf$regrOutput[[i]]$err.rate[length(predrf$regrOutput[[i]]$err.rate)])
  rferrcol_3=c(rferrcol_3,predrf$regrOutput[[i]]$err.rate[length(predrf$regrOutput[[i]]$err.rate)])
  rferrcol_5=c(rferrcol_5,predrf$regrOutput[[i]]$err.rate[length(predrf$regrOutput[[i]]$err.rate)])
}



lmerrcol=colMeans((prelm-test3r[,1:2402])^2)
lmerrcol3=colMeans((prelm-test3r[,1:2402])[testrace==3,]^2)
lmerrcol5=colMeans((prelm-test3r[,1:2402])[testrace==5,]^2)


errhold2=data.frame(name=colnames(errhold),cmmp=errhold[1,],mvcmmp=errhold[2,])
errhold2_3=data.frame(name=colnames(errhold_3),cmmp=errhold_3[1,],mvcmmp=errhold_3[2,])
errhold2_5=data.frame(name=colnames(errhold_5),cmmp=errhold_5[1,],mvcmmp=errhold_5[2,])
rferrcol2=data.frame(name=names(predrf$regrOutput),rf=rferrcol)
rferrcol2_3=data.frame(name=names(predrf$regrOutput),rf=rferrcol_3)
rferrcol2_5=data.frame(name=names(predrf$regrOutput),rf=rferrcol_5)
lmerrcol2=data.frame(name=names(train2r)[1:2402],mvlm=lmerrcol)
lmerrcol2_3=data.frame(name=names(train2r)[1:2402],mvlm=lmerrcol3)
lmerrcol2_5=data.frame(name=names(train2r)[1:2402],mvlm=lmerrcol5)

errhold3_rau=merge(merge(errhold2,rferrcol2,'name',all=T),lmerrcol2,'name',all=T)
names(errhold3_rau)[5]='lm'



errhold3_3_rau=merge(merge(errhold2_3,rferrcol2_3,'name',all=T),lmerrcol2_3,'name',all=T)
names(errhold3_3_rau)[5]='lm'

errhold3_5_rau=merge(merge(errhold2_5,rferrcol2_5,'name',all=T),lmerrcol2_5,'name',all=T)
names(errhold3_5_rau)[5]='lm'


print(z)

collector[[z]]=list(errors_beta=list(all=errhold3,black=errhold3_3,white=errhold3_5),
                    errors_m=list(all=errhold3_m,black=errhold3_3_m,white=errhold3_5_m),
                    errors_rau=list(all=errhold3_rau,black=errhold3_3_rau,white=errhold3_5_rau))

}

save.image('mvrun_10232021_v2.rdata')

##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################

##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
#collector structure 
#first is each run
#within each run is
	#coef col has 36 items for each outcome
	#errors contain one object for each method, each has a matrix with 3 columns, errors for all, white, black

#save.image('C:/Users/monster/Desktop/100 times/run100_CESC.rdata')

rownames(collector[[1]]$coefcol[[1]])==rownames(collector[[86]]$coefcol[[1]])

#need names of all snps used
nam=rownames(collector[[1]]$coefcol[[1]])
for(i in c(1:17,19:20)){
  for(j in 1:length(collector[[i]]$coefcol)){
    nam=union(nam,rownames(collector[[i]]$coefcol[[j]]))
  }
}


#coeficient plot
coefmats=list()

for( i in 1:length(collector[[1]]$coefcol)){
  coefmats[[i]]=matrix(NA,nrow=length(nam),ncol=100)
  rownames(coefmats[[i]])=nam
}

for(j in 1:length(collector[[1]]$coefcol)){
  for(i in c(1:17,19:20)){
    
    coefmats[[j]][,i]=matrix(collector[[i]]$coefcol[[j]])[match(rownames(coefmats[[j]]),rownames(collector[[i]]$coefcol[[j]])),1]
    
  }
  
}

se <- function(x) sqrt(var(x,na.rm=T)/sum(!is.na(x)))

#collecting mean and se 
meanholder=matrix(NA,nrow=length(nam),ncol=length(collector[[1]]$coefcol))
#seholder=matrix(NA,nrow=length(nam),nrow=length(collector[[1]]$coefcol))
rownames(meanholder)=nam
#rownames(seholder)=nam

for(i in 1:length(coefmats)){
  meanholder[,i]=rowMeans(coefmats[[i]],na.rm=T)
  #seholder[,i]=apply(coefmats[[i]],1)
}

meanholder2=meanholder[c(1,2:59,134,135,138,139,142,60:68,69:126,136,137,140,141,143,127:133),]#[c(1:62,137,63:130,135,136,138,131:134),]

vars=apply(meanholder2[-1,],1,var)
mvar=which(vars>quantile(vars,.9))

pdf('coefplot_CESC.pdf',width=25,height=12)
par(mar=c(6.5,4.1,2.5,1))
plot(1:(dim(meanholder2)[1]-1),meanholder2[-1,1],col=c(rep(2,63),rep(4,9),rep(3,70)),ylim=range(meanholder2),
     xlab='',ylab='Coefficent',xaxt='n',pch=20,cex=.7,main='Elastic Net Coefficients, mean over 20 samples, shaded are interactions - CESC only')
for(i in 2:length(coefmats)){
  points(1:(dim(meanholder2)[1]-1),meanholder2[-1,i],col=c(rep(2,63),rep(4,9),rep(3,70)),pch=20,cex=.7)
}
rect(72.5, -5, 150, 5, density = 8)
legend('topleft',legend=c('sCNA','Clinical','sCNA*Race'),lty=1,col=c(2,4,3,'purple'),bg='transparent')
axis(1,at=mvar,rownames(meanholder2)[-1][mvar],las=2,cex.axis=.46)
dev.off()

#error plot
vars=2472

errcol=list()
errcol$cmmpcol=matrix(NA,ncol=100,nrow=vars)
errcol$lmcol=matrix(NA,ncol=100,nrow=vars)
errcol$lascol=matrix(NA,ncol=100,nrow=vars)
errcol$las_cmmpcol=matrix(NA,ncol=100,nrow=vars)
errcol$rfcol=matrix(NA,ncol=100,nrow=vars)

errcol$cmmpcolw=matrix(NA,ncol=100,nrow=vars)
errcol$lmcolw=matrix(NA,ncol=100,nrow=vars)
errcol$lascolw=matrix(NA,ncol=100,nrow=vars)
errcol$las_cmmpcolw=matrix(NA,ncol=100,nrow=vars)
errcol$rfcolw=matrix(NA,ncol=100,nrow=vars)

errcol$cmmpcolb=matrix(NA,ncol=100,nrow=vars)
errcol$lmcolb=matrix(NA,ncol=100,nrow=vars)
errcol$lascolb=matrix(NA,ncol=100,nrow=vars)
errcol$las_cmmpcolb=matrix(NA,ncol=100,nrow=vars)
errcol$rfcolb=matrix(NA,ncol=100,nrow=vars)



errcol$cmmpcol_m=matrix(NA,ncol=100,nrow=vars)
errcol$lmcol_m=matrix(NA,ncol=100,nrow=vars)
errcol$lascol_m=matrix(NA,ncol=100,nrow=vars)
errcol$las_cmmpcol_m=matrix(NA,ncol=100,nrow=vars)
errcol$rfcol_m=matrix(NA,ncol=100,nrow=vars)

errcol$cmmpcolw_m=matrix(NA,ncol=100,nrow=vars)
errcol$lmcolw_m=matrix(NA,ncol=100,nrow=vars)
errcol$lascolw_m=matrix(NA,ncol=100,nrow=vars)
errcol$las_cmmpcolw_m=matrix(NA,ncol=100,nrow=vars)
errcol$rfcolw_m=matrix(NA,ncol=100,nrow=vars)

errcol$cmmpcolb_m=matrix(NA,ncol=100,nrow=vars)
errcol$lmcolb_m=matrix(NA,ncol=100,nrow=vars)
errcol$lascolb_m=matrix(NA,ncol=100,nrow=vars)
errcol$las_cmmpcolb_m=matrix(NA,ncol=100,nrow=vars)
errcol$rfcolb_m=matrix(NA,ncol=100,nrow=vars)



errcol$cmmpcol_rau=matrix(NA,ncol=100,nrow=vars)
errcol$lmcol_rau=matrix(NA,ncol=100,nrow=vars)
errcol$lascol_rau=matrix(NA,ncol=100,nrow=vars)
errcol$las_cmmpcol_rau=matrix(NA,ncol=100,nrow=vars)
errcol$rfcol_rau=matrix(NA,ncol=100,nrow=vars)

errcol$cmmpcolw_rau=matrix(NA,ncol=100,nrow=vars)
errcol$lmcolw_rau=matrix(NA,ncol=100,nrow=vars)
errcol$lascolw_rau=matrix(NA,ncol=100,nrow=vars)
errcol$las_cmmpcolw_rau=matrix(NA,ncol=100,nrow=vars)
errcol$rfcolw_rau=matrix(NA,ncol=100,nrow=vars)

errcol$cmmpcolb_rau=matrix(NA,ncol=100,nrow=vars)
errcol$lmcolb_rau=matrix(NA,ncol=100,nrow=vars)
errcol$lascolb_rau=matrix(NA,ncol=100,nrow=vars)
errcol$las_cmmpcolb_rau=matrix(NA,ncol=100,nrow=vars)
errcol$rfcolb_rau=matrix(NA,ncol=100,nrow=vars)






for(i in c(1:20)){
  errcol$cmmpcol[,i]=collector[[i]]$errors_beta$cmmpcol[,1]
  errcol$cmmpcolw[,i]=collector[[i]]$errors$cmmpcol[,2]
  errcol$cmmpcolb[,i]=collector[[i]]$errors$cmmpcol[,3]
  
  errcol$lmcol[,i]=collector[[i]]$errors_beta$lmcol[,1]
  errcol$lmcolw[,i]=collector[[i]]$errors_beta$lmcol[,2]
  errcol$lmcolb[,i]=collector[[i]]$errors_beta$lmcol[,3]
  
  errcol$lascol[,i]=collector[[i]]$errors_beta$lascol[,1]
  errcol$lascolw[,i]=collector[[i]]$errors_beta$lascol[,2]
  errcol$lascolb[,i]=collector[[i]]$errors_beta$lascol[,3]
  
  errcol$las_cmmpcol[,i]=collector[[i]]$errors_beta$las_cmmpcol[,1]
  errcol$las_cmmpcolw[,i]=collector[[i]]$errors_beta$las_cmmpcol[,2]
  errcol$las_cmmpcolb[,i]=collector[[i]]$errors_beta$las_cmmpcol[,3]
  
  errcol$rfcol[,i]=collector[[i]]$errors_beta$rfcol[,1]
  errcol$rfcolw[,i]=collector[[i]]$errors_beta$rfcol[,2]
  errcol$rfcolb[,i]=collector[[i]]$errors_beta$rfcol[,3]
  
  
  
  errcol$cmmpcol_m[,i]=collector[[i]]$errors_m$cmmpcol[,1]
  errcol$cmmpcolw_m[,i]=collector[[i]]$errors_m$cmmpcol[,2]
  errcol$cmmpcolb_m[,i]=collector[[i]]$errors_m$cmmpcol[,3]
  
  errcol$lmcol_m[,i]=collector[[i]]$errors_m$lmcol[,1]
  errcol$lmcolw_m[,i]=collector[[i]]$errors_m$lmcol[,2]
  errcol$lmcolb_m[,i]=collector[[i]]$errors_m$lmcol[,3]
  
  errcol$lascol_m[,i]=collector[[i]]$errors_m$lascol[,1]
  errcol$lascolw_m[,i]=collector[[i]]$errors_m$lascol[,2]
  errcol$lascolb_m[,i]=collector[[i]]$errors_m$lascol[,3]
  
  errcol$las_cmmpcol_m[,i]=collector[[i]]$errors_m$las_cmmpcol[,1]
  errcol$las_cmmpcolw_m[,i]=collector[[i]]$errors_m$las_cmmpcol[,2]
  errcol$las_cmmpcolb_m[,i]=collector[[i]]$errors_m$las_cmmpcol[,3]
  
  errcol$rfcol_m[,i]=collector[[i]]$errors_m$rfcol[,1]
  errcol$rfcolw_m[,i]=collector[[i]]$errors_m$rfcol[,2]
  errcol$rfcolb_m[,i]=collector[[i]]$errors_m$rfcol[,3]
  
  
  errcol$cmmpcol_rau[,i]=collector[[i]]$errors_rau$cmmpcol[,1]
  errcol$cmmpcolw_rau[,i]=collector[[i]]$errors_rau$cmmpcol[,2]
  errcol$cmmpcolb_rau[,i]=collector[[i]]$errors_rau$cmmpcol[,3]
  
  errcol$lmcol_rau[,i]=collector[[i]]$errors_rau$lmcol[,1]
  errcol$lmcolw_rau[,i]=collector[[i]]$errors_rau$lmcol[,2]
  errcol$lmcolb_rau[,i]=collector[[i]]$errors_rau$lmcol[,3]
  
  errcol$lascol_rau[,i]=collector[[i]]$errors_rau$lascol[,1]
  errcol$lascolw_rau[,i]=collector[[i]]$errors_rau$lascol[,2]
  errcol$lascolb_rau[,i]=collector[[i]]$errors_rau$lascol[,3]
  
  errcol$las_cmmpcol_rau[,i]=collector[[i]]$errors_rau$las_cmmpcol[,1]
  errcol$las_cmmpcolw_rau[,i]=collector[[i]]$errors_rau$las_cmmpcol[,2]
  errcol$las_cmmpcolb_rau[,i]=collector[[i]]$errors_rau$las_cmmpcol[,3]
  
  errcol$rfcol_rau[,i]=collector[[i]]$errors_rau$rfcol[,1]
  errcol$rfcolw_rau[,i]=collector[[i]]$errors_rau$rfcol[,2]
  errcol$rfcolb_rau[,i]=collector[[i]]$errors_rau$rfcol[,3]
  
}

error=list()
error$cmmp=rowMeans(errcol$cmmpcol,na.rm=T)
error$cmmpw=rowMeans(errcol$cmmpcolw,na.rm=T)
error$cmmpb=rowMeans(errcol$cmmpcolb,na.rm=T)

error$lm=rowMeans(errcol$lmcol,na.rm=T)
error$lmw=rowMeans(errcol$lmcolw,na.rm=T)
error$lmb=rowMeans(errcol$lmcolb,na.rm=T)

error$las=rowMeans(errcol$lascol,na.rm=T)
error$lasw=rowMeans(errcol$lascolw,na.rm=T)
error$lasb=rowMeans(errcol$lascolb,na.rm=T)

error$las_cmmp=rowMeans(errcol$las_cmmpcol,na.rm=T)
error$las_cmmpw=rowMeans(errcol$las_cmmpcolw,na.rm=T)
error$las_cmmpb=rowMeans(errcol$las_cmmpcolb,na.rm=T)


error$rf=rowMeans(errcol$rfcol,na.rm=T)
error$rfw=rowMeans(errcol$rfcolw,na.rm=T)
error$rfb=rowMeans(errcol$rfcolb,na.rm=T)


error$cmmp_m=rowMeans(errcol$cmmpcol_m,na.rm=T)
error$cmmpw_m=rowMeans(errcol$cmmpcolw_m,na.rm=T)
error$cmmpb_m=rowMeans(errcol$cmmpcolb_m,na.rm=T)

error$lm_m=rowMeans(errcol$lmcol_m,na.rm=T)
error$lmw_m=rowMeans(errcol$lmcolw_m,na.rm=T)
error$lmb_m=rowMeans(errcol$lmcolb_m,na.rm=T)

error$las_m=rowMeans(errcol$lascol_m,na.rm=T)
error$lasw_m=rowMeans(errcol$lascolw_m,na.rm=T)
error$lasb_m=rowMeans(errcol$lascolb_m,na.rm=T)

error$las_cmmp_m=rowMeans(errcol$las_cmmpcol_m,na.rm=T)
error$las_cmmpw_m=rowMeans(errcol$las_cmmpcolw_m,na.rm=T)
error$las_cmmpb_m=rowMeans(errcol$las_cmmpcolb_m,na.rm=T)

error$rf_m=rowMeans(errcol$rfcol_m,na.rm=T)
error$rfw_m=rowMeans(errcol$rfcolw_m,na.rm=T)
error$rfb_m=rowMeans(errcol$rfcolb_m,na.rm=T)


error$cmmp_rau=rowMeans(errcol$cmmpcol_rau,na.rm=T)
error$cmmpw_rau=rowMeans(errcol$cmmpcolw_rau,na.rm=T)
error$cmmpb_rau=rowMeans(errcol$cmmpcolb_rau,na.rm=T)

error$lm_rau=rowMeans(errcol$lmcol_rau,na.rm=T)
error$lmw_rau=rowMeans(errcol$lmcolw_rau,na.rm=T)
error$lmb_rau=rowMeans(errcol$lmcolb_rau,na.rm=T)

error$las_rau=rowMeans(errcol$lascol_rau,na.rm=T)
error$lasw_rau=rowMeans(errcol$lascolw_rau,na.rm=T)
error$lasb_rau=rowMeans(errcol$lascolb_rau,na.rm=T)

error$las_cmmp_rau=rowMeans(errcol$las_cmmpcol_rau,na.rm=T)
error$las_cmmpw_rau=rowMeans(errcol$las_cmmpcolw_rau,na.rm=T)
error$las_cmmpb_rau=rowMeans(errcol$las_cmmpcolb_rau,na.rm=T)

error$rf_rau=rowMeans(errcol$rfcol_rau,na.rm=T)
error$rfw_rau=rowMeans(errcol$rfcolw_rau,na.rm=T)
error$rfb_rau=rowMeans(errcol$rfcolb_rau,na.rm=T)


errors=list()
errors$cmmp=apply(errcol$cmmpcol,1,se)
errors$cmmpw=apply(errcol$cmmpcolw,1,se)
errors$cmmpb=apply(errcol$cmmpcolb,1,se)

errors$lm=apply(errcol$lmcol,1,se)
errors$lmw=apply(errcol$lmcolw,1,se)
errors$lmb=apply(errcol$lmcolb,1,se)

errors$las=apply(errcol$lascol,1,se)
errors$lasw=apply(errcol$lascolw,1,se)
errors$lasb=apply(errcol$lascolb,1,se)

errors$las_cmmp=apply(errcol$las_cmmpcol,1,se)
errors$las_cmmpw=apply(errcol$las_cmmpcolw,1,se)
errors$las_cmmpb=apply(errcol$las_cmmpcolb,1,se)

errors$rf=apply(errcol$rfcol,1,se)
errors$rfw=apply(errcol$rfcolw,1,se)
errors$rfb=apply(errcol$rfcolb,1,se)


errors$cmmp_m=apply(errcol$cmmpcol_m,1,se)
errors$cmmpw_m=apply(errcol$cmmpcolw_m,1,se)
errors$cmmpb_m=apply(errcol$cmmpcolb_m,1,se)

errors$lm_m=apply(errcol$lmcol_m,1,se)
errors$lmw_m=apply(errcol$lmcolw_m,1,se)
errors$lmb_m=apply(errcol$lmcolb_m,1,se)

errors$las_m=apply(errcol$lascol_m,1,se)
errors$lasw_m=apply(errcol$lascolw_m,1,se)
errors$lasb_m=apply(errcol$lascolb_m,1,se)

errors$las_cmmp_m=apply(errcol$las_cmmpcol_m,1,se)
errors$las_cmmpw_m=apply(errcol$las_cmmpcolw_m,1,se)
errors$las_cmmpb_m=apply(errcol$las_cmmpcolb_m,1,se)

errors$rf_m=apply(errcol$rfcol_m,1,se)
errors$rfw_m=apply(errcol$rfcolw_m,1,se)
errors$rfb_m=apply(errcol$rfcolb_m,1,se)


errors$cmmp_rau=apply(errcol$cmmpcol_rau,1,se)
errors$cmmpw_rau=apply(errcol$cmmpcolw_rau,1,se)
errors$cmmpb_rau=apply(errcol$cmmpcolb_rau,1,se)

errors$lm_rau=apply(errcol$lmcol_rau,1,se)
errors$lmw_rau=apply(errcol$lmcolw_rau,1,se)
errors$lmb_rau=apply(errcol$lmcolb_rau,1,se)

errors$las_rau=apply(errcol$lascol_rau,1,se)
errors$lasw_rau=apply(errcol$lascolw_rau,1,se)
errors$lasb_rau=apply(errcol$lascolb_rau,1,se)

errors$las_cmmp_rau=apply(errcol$las_cmmpcol_rau,1,se)
errors$las_cmmpw_rau=apply(errcol$las_cmmpcolw_rau,1,se)
errors$las_cmmpb_rau=apply(errcol$las_cmmpcolb_rau,1,se)

errors$rf_rau=apply(errcol$rfcol_rau,1,se)
errors$rfw_rau=apply(errcol$rfcolw_rau,1,se)
errors$rfb_rau=apply(errcol$rfcolb_rau,1,se)

pdf('cesc_beta.pdf',height = 8,width = 20)
library(vioplot)
par(mar=c(4.5, 4.1, 6.5, 1.0))
labs=apply(expand.grid(c('CMMP','LM','ENET','CMMP+ENET','RF'),c('All','White','Black')),1,paste,collapse='\n')
labs=labs[c(1,6,11,2,7,12,3,8,13,4,9,14,5,10,15)]
vioplot(error$cmmp,error$cmmpw,error$cmmpb,error$lm,error$lmw,error$lmb,error$las,error$lasw,error$lasb,
        error$las_cmmp,error$las_cmmpw,error$las_cmmpb,error$rf,error$rfw,error$rfb,names=labs,ylab='MSPE',main='MSPE averaged over 20 samples, Beta values - CESC only',cex.axis=.8,col=rep(c('coral3','lightblue','orchid'),5))

#overlay points
holder=c(error$cmmp,error$cmmpw,error$cmmpb,error$lm,error$lmw,error$lmb,error$las,error$lasw,error$lasb,
         error$las_cmmp,error$las_cmmpw,error$las_cmmpb,error$rf,error$rfw,error$rfb)
holders=c(errors$cmmp,errors$cmmpw,errors$cmmpb,errors$lm,errors$lmw,errors$lmb,errors$las,errors$lasw,errors$lasb,
          errors$las_cmmp,errors$las_cmmpw,errors$las_cmmpb,errors$rf,errors$rfw,errors$rfb)

holders_scaled=holders/max(holders,na.rm=T)*.25

for(i in 1:(15*vars)){
  arrows(floor(i/vars)+1-holders_scaled[i],holder[i],floor(i/vars)+1+holders_scaled[i],holder[i],code=1,angle=90,length=.0,col='grey50',lwd=.5)
}

pos=NULL
sca=NULL
for(i in 1:15){
  pos=c(pos,((1:15)[i]+c(-.25,-.125,.125,.25)))
  sca=c(sca,0+c(-max(holders,na.rm=T),-.5*max(holders,na.rm=T),.5*max(holders,na.rm=T),max(holders,na.rm=T)))
}
mtext('SE',at=c(0.1))
axis(3,at=pos,round(sca,3),cex.axis=.6,las=2)
dev.off()


pdf('cesc_m.pdf',height = 8,width = 20)
par(mar=c(4.5, 4.1, 6.5, 1.0))
vioplot(error$cmmp_m,error$cmmpw_m,error$cmmpb_m,error$lm_m,error$lmw_m,error$lmb_m,error$las_m,error$lasw_m,error$lasb_m,
        error$las_cmmp_m,error$las_cmmpw_m,error$las_cmmpb_m,error$rf_m,error$rfw_m,error$rfb_m,names=labs,ylab='MSPE',main='MSPE averaged over 20 samples, M values - CESC only',cex.axis=.8,col=rep(c('coral3','lightblue','orchid'),5))

#overlay points
holder=c(error$cmmp_m,error$cmmpw_m,error$cmmpb_m,error$lm_m,error$lmw_m,error$lmb_m,error$las_m,error$lasw_m,error$lasb_m,
         error$las_cmmp_m,error$las_cmmpw_m,error$las_cmmpb_m,error$rf_m,error$rfw_m,error$rfb_m)
holders=c(errors$cmmp_m,errors$cmmpw_m,errors$cmmpb_m,errors$lm_m,errors$lmw_m,errors$lmb_m,errors$las_m,errors$lasw_m,errors$lasb_m,
          errors$las_cmmp_m,errors$las_cmmpw_m,errors$las_cmmpb_m,errors$rf_m,errors$rfw_m,errors$rfb_m)

holders_scaled=holders/max(holders,na.rm=T)*.25

for(i in 1:(15*vars)){
  arrows(floor(i/vars)+1-holders_scaled[i],holder[i],floor(i/vars)+1+holders_scaled[i],holder[i],code=1,angle=90,length=.0,col='grey50',lwd=.5)
}

pos=NULL
sca=NULL
for(i in 1:15){
  pos=c(pos,((1:15)[i]+c(-.25,-.125,.125,.25)))
  sca=c(sca,0+c(-max(holders,na.rm=T),-.5*max(holders,na.rm=T),.5*max(holders,na.rm=T),max(holders,na.rm=T)))
}
mtext('SE',at=c(0.1))
axis(3,at=pos,round(sca,3),cex.axis=.6,las=2)
dev.off()

pdf('cesc_rau.pdf',height = 8,width = 20)
par(mar=c(4.5, 4.1, 6.5, 1.0))
vioplot(error$cmmp_rau,error$cmmpw_rau,error$cmmpb_rau,error$lm_rau,error$lmw_rau,error$lmb_rau,error$las_rau,error$lasw_rau,error$lasb_rau,
        error$las_cmmp_rau,error$las_cmmpw_rau,error$las_cmmpb_rau,error$rf_rau,error$rfw_rau,error$rfb_rau,names=labs,ylab='MSPE',main='MSPE averaged over 20 samples, RAU values - CESC only',cex.axis=.8,col=rep(c('coral3','lightblue','orchid'),5))

#overlay points
holder=c(error$cmmp_rau,error$cmmpw_rau,error$cmmpb_rau,error$lm_rau,error$lmw_rau,error$lmb_rau,error$las_rau,error$lasw_rau,error$lasb_rau,
         error$las_cmmp_rau,error$las_cmmpw_rau,error$las_cmmpb_rau,error$rf_rau,error$rfw_rau,error$rfb_rau)
holders=c(errors$cmmp_rau,errors$cmmpw_rau,errors$cmmpb_rau,errors$lm_rau,errors$lmw_rau,errors$lmb_rau,errors$las_rau,errors$lasw_rau,errors$lasb_rau,
          errors$las_cmmp_rau,errors$las_cmmpw_rau,errors$las_cmmpb_rau,errors$rf_rau,errors$rfw_rau,errors$rfb_rau)

holders_scaled=holders/max(holders,na.rm=T)*.25

for(i in 1:(15*vars)){
  arrows(floor(i/vars)+1-holders_scaled[i],holder[i],floor(i/vars)+1+holders_scaled[i],holder[i],code=1,angle=90,length=.0,col='grey50',lwd=.5)
}

pos=NULL
sca=NULL
for(i in 1:15){
  pos=c(pos,((1:15)[i]+c(-.25,-.125,.125,.25)))
  sca=c(sca,0+c(-max(holders,na.rm=T),-.5*max(holders,na.rm=T),.5*max(holders,na.rm=T),max(holders,na.rm=T)))
}
mtext('SE',at=c(0.1))
axis(3,at=pos,round(sca,1),cex.axis=.6,las=2)

dev.off()

#points(c(rep(1,36),rep(2,36),rep(3,36),rep(4,36),rep(5,36),rep(6,36),rep(7,36),rep(8,36),rep(9,36),rep(10,36),rep(11,36),rep(12,36))[i],c(error$cmmp,error$cmmpw,error$cmmpb,error$lm,error$lmw,error$lmb,error$las,error$lasw,error$lasb,
#error$las_cmmp,error$las_cmmpw,error$las_cmmpb)[i],col='black',pch=20,cex=.8)

#for(i in 1:(12*36)){
#points(c(rep(1,36),rep(2,36),rep(3,36),rep(4,36),rep(5,36),rep(6,36),rep(7,36),rep(8,36),rep(9,36),rep(10,36),rep(11,36),rep(12,36))[i],c(error$cmmp,error$cmmpw,error$cmmpb,error$lm,error$lmw,error$lmb,error$las,error$lasw,error$lasb,
#error$las_cmmp,error$las_cmmpw,error$las_cmmpb)[i],col='black',pch=20,cex=(holders/median(holders))[i])
#}







#cut(sds[,1],breaks=seq(min(sds),max(sds),(max(sds)-min(sds))/4))



sds=data.frame(error$cmmp,error$cmmpw,error$cmmpb,error$lm,error$lmw,error$lmb,error$las,error$lasw,error$lasb,error$las_cmmp,error$las_cmmpw,error$las_cmmpb)

cols=NULL
for(i in 1:dim(sds)[2]){
cols=c(cols,cut(sds[,i],breaks=seq(min(sds),max(sds),(max(sds)-min(sds))/4)))
}

points(c(rep(1,36),rep(2,36),rep(3,36),rep(4,36),rep(5,36),rep(6,36),rep(7,36),rep(8,36),rep(9,36),rep(10,36),rep(11,36),rep(12,36))[i],c(error$cmmp,error$cmmpw,error$cmmpb,error$lm,error$lmw,error$lmb,error$las,error$lasw,error$lasb,
error$las_cmmp,error$las_cmmpw,error$las_cmmpb)[i],pch=20,cex=.8,col=cols)





error$las_cmmp,error$las_cmmpw,error$las_cmmpb)
names(sds)=labs
library(xtable)
xtable(sds)



i=collector

cmmperr=cmmperrw=cmmperrb=mvcmmperr=mvcmmperrw=mvcmmperrb=lmerr=lmerrw=lmerrb=rferr=rferrw=rferrb=list(matrix(NA,ncol=20,nrow=2402),matrix(NA,ncol=20,nrow=2402),matrix(NA,ncol=20,nrow=2402))

names(cmmperr)=names(cmmperrw)=names(cmmperrb)=names(mvcmmperr)=names(mvcmmperrw)=names(mvcmmperrb)=names(lmerr)=
names(lmerrw)=names(lmerrb)=names(rferr)=names(rferrw)=names(rferrb)=c('rau','m','beta')


for(j in 1:20){
cmmperr$beta[,j]=i[[j]]$errors_beta$all$cmmp
cmmperrw$beta[,j]=i[[j]]$errors_beta$white$cmmp
cmmperrb$beta[,j]=i[[j]]$errors_beta$black$cmmp

mvcmmperr$beta[,j]=i[[j]]$errors_beta$all$mvcmmp
mvcmmperrw$beta[,j]=i[[j]]$errors_beta$white$mvcmmp
mvcmmperrb$beta[,j]=i[[j]]$errors_beta$black$mvcmmp

rferr$beta[,j]=i[[j]]$errors_beta$all$rf
rferrw$beta[,j]=i[[j]]$errors_beta$white$rf
rferrb$beta[,j]=i[[j]]$errors_beta$black$rf

lmerr$beta[,j]=i[[j]]$errors_beta$all$lm
lmerrw$beta[,j]=i[[j]]$errors_beta$white$lm
lmerrb$beta[,j]=i[[j]]$errors_beta$black$lm



cmmperr$rau[,j]=i[[j]]$errors_rau$all$cmmp
cmmperrw$rau[,j]=i[[j]]$errors_rau$white$cmmp
cmmperrb$rau[,j]=i[[j]]$errors_rau$black$cmmp

mvcmmperr$rau[,j]=i[[j]]$errors_rau$all$mvcmmp
mvcmmperrw$rau[,j]=i[[j]]$errors_rau$white$mvcmmp
mvcmmperrb$rau[,j]=i[[j]]$errors_rau$black$mvcmmp

rferr$rau[,j]=i[[j]]$errors_rau$all$rf
rferrw$rau[,j]=i[[j]]$errors_rau$white$rf
rferrb$rau[,j]=i[[j]]$errors_rau$black$rf

lmerr$rau[,j]=i[[j]]$errors_rau$all$lm
lmerrw$rau[,j]=i[[j]]$errors_rau$white$lm
lmerrb$rau[,j]=i[[j]]$errors_rau$black$lm




cmmperr$m[,j]=i[[j]]$errors_m$all$cmmp
cmmperrw$m[,j]=i[[j]]$errors_m$white$cmmp
cmmperrb$m[,j]=i[[j]]$errors_m$black$cmmp

mvcmmperr$m[,j]=i[[j]]$errors_m$all$mvcmmp
mvcmmperrw$m[,j]=i[[j]]$errors_m$white$mvcmmp
mvcmmperrb$m[,j]=i[[j]]$errors_m$black$mvcmmp

rferr$m[,j]=i[[j]]$errors_m$all$rf
rferrw$m[,j]=i[[j]]$errors_m$white$rf
rferrb$m[,j]=i[[j]]$errors_m$black$rf

lmerr$m[,j]=i[[j]]$errors_m$all$lm
lmerrw$m[,j]=i[[j]]$errors_m$white$lm
lmerrb$m[,j]=i[[j]]$errors_m$black$lm
}

##################new plots###############################
##########################################################
##########################################################
##########################################################
##########################################################



library(dichromat)
col=c(colorRampPalette(c('tan','tan2'))(3),colorRampPalette(c('lightblue','darkblue'))(3),col1=colorRampPalette(c('purple','purple4'))(3),
col1=colorRampPalette(c('lightgreen','darkgreen'))(3))
nams=apply(expand.grid( c('all','white','black'),c('CMMP','mvCMMP','mvRF','LM'))[2:1], 1, paste, collapse="\n")
se=function(x){ sqrt(var(x,na.rm=T)/sum(!is.na(x)))}
rowSEs=function(i,na.rm=T){apply(i,1,se)}


#M value
pdf('mval_mv.pdf',width=15,height=4)
par(mar=c(4.5, 4.1, 6.5, 1.0))
vioplot(rowMeans(cmmperr$m,na.rm=T),rowMeans(cmmperrw$m,na.rm=T),rowMeans(cmmperrb$m,na.rm=T),
rowMeans(mvcmmperr$m,na.rm=T),rowMeans(mvcmmperrw$m,na.rm=T),rowMeans(mvcmmperrb$m,na.rm=T),
rowMeans(rferr$m,na.rm=T),rowMeans(rferrw$m,na.rm=T),rowMeans(rferrb$m,na.rm=T),
rowMeans(lmerr$m,na.rm=T),rowMeans(lmerrw$m,na.rm=T),rowMeans(lmerrb$m,na.rm=T),
names=nams,main='M value',col=col,ylab='MSPE',mgp=c(4,1.5,0))

holders=c(rowSEs(cmmperr$m,na.rm=T),rowSEs(cmmperrw$m,na.rm=T),rowSEs(cmmperrb$m,na.rm=T),
rowSEs(mvcmmperr$m,na.rm=T),rowSEs(mvcmmperrw$m,na.rm=T),rowSEs(mvcmmperrb$m,na.rm=T),
rowSEs(rferr$m,na.rm=T),rowSEs(rferrw$m,na.rm=T),rowSEs(rferrb$m,na.rm=T),
rowSEs(lmerr$m,na.rm=T),rowSEs(lmerrw$m,na.rm=T),rowSEs(lmerrb$m,na.rm=T))
holder=c(rowMeans(cmmperr$m,na.rm=T),rowMeans(cmmperrw$m,na.rm=T),rowMeans(cmmperrb$m,na.rm=T),
rowMeans(mvcmmperr$m,na.rm=T),rowMeans(mvcmmperrw$m,na.rm=T),rowMeans(mvcmmperrb$m,na.rm=T),
rowMeans(rferr$m,na.rm=T),rowMeans(rferrw$m,na.rm=T),rowMeans(rferrb$m,na.rm=T),
rowMeans(lmerr$m,na.rm=T),rowMeans(lmerrw$m,na.rm=T),rowMeans(lmerrb$m,na.rm=T))

holders_scaled=holders/max(holders,na.rm=T)*.25

for(i in 1:(15*2402)){
  arrows(floor(i/vars)+1-holders_scaled[i],holder[i],floor(i/vars)+1+holders_scaled[i],holder[i],code=1,angle=90,length=.0,col='grey50',lwd=.5)
}

pos=NULL
sca=NULL
for(i in 1:15){
  pos=c(pos,((1:15)[i]+c(-.25,-.125,.125,.25)))
  sca=c(sca,0+c(-max(holders,na.rm=T),-.5*max(holders,na.rm=T),.5*max(holders,na.rm=T),max(holders,na.rm=T)))
}
mtext('SE',at=c(0.1))
axis(3,at=pos,round(sca,1),cex.axis=.6,las=2)
dev.off()



#RAU value

pdf('rauval_mv.pdf',width=15,height=4)
par(mar=c(4.5, 4.1, 6.5, 1.0))
vioplot(rowMeans(cmmperr$rau,na.rm=T),rowMeans(cmmperrw$rau,na.rm=T),rowMeans(cmmperrb$rau,na.rm=T),
rowMeans(mvcmmperr$rau,na.rm=T),rowMeans(mvcmmperrw$rau,na.rm=T),rowMeans(mvcmmperrb$rau,na.rm=T),
rowMeans(rferr$rau,na.rm=T),rowMeans(rferrw$rau,na.rm=T),rowMeans(rferrb$rau,na.rm=T),
rowMeans(lmerr$rau,na.rm=T),rowMeans(lmerrw$rau,na.rm=T),rowMeans(lmerrb$rau,na.rm=T),
names=nams,main='RAU value',col=col,ylab='MSPE',mgp=c(4,1.5,0))

holders=c(rowSEs(cmmperr$rau,na.rm=T),rowSEs(cmmperrw$rau,na.rm=T),rowSEs(cmmperrb$rau,na.rm=T),
rowSEs(mvcmmperr$rau,na.rm=T),rowSEs(mvcmmperrw$rau,na.rm=T),rowSEs(mvcmmperrb$rau,na.rm=T),
rowSEs(rferr$rau,na.rm=T),rowSEs(rferrw$rau,na.rm=T),rowSEs(rferrb$rau,na.rm=T),
rowSEs(lmerr$rau,na.rm=T),rowSEs(lmerrw$rau,na.rm=T),rowSEs(lmerrb$rau,na.rm=T))
holder=c(rowMeans(cmmperr$rau,na.rm=T),rowMeans(cmmperrw$rau,na.rm=T),rowMeans(cmmperrb$rau,na.rm=T),
rowMeans(mvcmmperr$rau,na.rm=T),rowMeans(mvcmmperrw$rau,na.rm=T),rowMeans(mvcmmperrb$rau,na.rm=T),
rowMeans(rferr$rau,na.rm=T),rowMeans(rferrw$rau,na.rm=T),rowMeans(rferrb$rau,na.rm=T),
rowMeans(lmerr$rau,na.rm=T),rowMeans(lmerrw$rau,na.rm=T),rowMeans(lmerrb$rau,na.rm=T))

holders_scaled=holders/max(holders,na.rm=T)*.25

for(i in 1:(15*2402)){
  arrows(floor(i/vars)+1-holders_scaled[i],holder[i],floor(i/vars)+1+holders_scaled[i],holder[i],code=1,angle=90,length=.0,col='grey50',lwd=.5)
}

pos=NULL
sca=NULL
for(i in 1:15){
  pos=c(pos,((1:15)[i]+c(-.25,-.125,.125,.25)))
  sca=c(sca,0+c(-max(holders,na.rm=T),-.5*max(holders,na.rm=T),.5*max(holders,na.rm=T),max(holders,na.rm=T)))
}
mtext('SE',at=c(0.1))
axis(3,at=pos,round(sca,1),cex.axis=.6,las=2)
dev.off()



#beta value
pdf('betaval_mv.pdf',width=15,height=4)
par(mar=c(4.5, 4.1, 6.5, 1.0))
vioplot(rowMeans(cmmperr$beta,na.rm=T),rowMeans(cmmperrw$beta,na.rm=T),rowMeans(cmmperrb$beta,na.rm=T),
rowMeans(mvcmmperr$beta,na.rm=T),rowMeans(mvcmmperrw$beta,na.rm=T),rowMeans(mvcmmperrb$beta,na.rm=T),
rowMeans(rferr$beta,na.rm=T),rowMeans(rferrw$beta,na.rm=T),rowMeans(rferrb$beta,na.rm=T),
rowMeans(lmerr$beta,na.rm=T),rowMeans(lmerrw$beta,na.rm=T),rowMeans(lmerrb$beta,na.rm=T),
names=nams,main='Beta value',col=col,ylab='MSPE',mgp=c(4,1.5,0))

holders=c(rowSEs(cmmperr$beta,na.rm=T),rowSEs(cmmperrw$beta,na.rm=T),rowSEs(cmmperrb$beta,na.rm=T),
rowSEs(mvcmmperr$beta,na.rm=T),rowSEs(mvcmmperrw$beta,na.rm=T),rowSEs(mvcmmperrb$beta,na.rm=T),
rowSEs(rferr$beta,na.rm=T),rowSEs(rferrw$beta,na.rm=T),rowSEs(rferrb$beta,na.rm=T),
rowSEs(lmerr$beta,na.rm=T),rowSEs(lmerrw$beta,na.rm=T),rowSEs(lmerrb$beta,na.rm=T))
holder=c(rowMeans(cmmperr$beta,na.rm=T),rowMeans(cmmperrw$beta,na.rm=T),rowMeans(cmmperrb$beta,na.rm=T),
rowMeans(mvcmmperr$beta,na.rm=T),rowMeans(mvcmmperrw$beta,na.rm=T),rowMeans(mvcmmperrb$beta,na.rm=T),
rowMeans(rferr$beta,na.rm=T),rowMeans(rferrw$beta,na.rm=T),rowMeans(rferrb$beta,na.rm=T),
rowMeans(lmerr$beta,na.rm=T),rowMeans(lmerrw$beta,na.rm=T),rowMeans(lmerrb$beta,na.rm=T))

holders_scaled=holders/max(holders,na.rm=T)*.25

for(i in 1:(15*2402)){
  arrows(floor(i/vars)+1-holders_scaled[i],holder[i],floor(i/vars)+1+holders_scaled[i],holder[i],code=1,angle=90,length=.0,col='grey50',lwd=.5)
}

pos=NULL
sca=NULL
for(i in 1:15){
  pos=c(pos,((1:15)[i]+c(-.25,-.125,.125,.25)))
  sca=c(sca,0+c(-max(holders,na.rm=T),-.5*max(holders,na.rm=T),.5*max(holders,na.rm=T),max(holders,na.rm=T)))
}
mtext('SE',at=c(0.1))
axis(3,at=pos,round(sca,4),cex.axis=.6,las=2)

dev.off()


####################################

#diffplot

#M value
pdf('mval_mv_diff.pdf',width=15,height=4)
par(mar=c(4.5, 4.1, 6.5, 1.0))
vioplot(rowMeans(cmmperr$m,na.rm=T)-rowMeans(lmerr$m,na.rm=T),rowMeans(cmmperrw$m,na.rm=T)-rowMeans(lmerrw$m,na.rm=T),rowMeans(cmmperrb$m,na.rm=T)-rowMeans(lmerrb$m,na.rm=T),
        rowMeans(mvcmmperr$m,na.rm=T)-rowMeans(lmerr$m,na.rm=T),rowMeans(mvcmmperrw$m,na.rm=T)-rowMeans(lmerrw$m,na.rm=T),rowMeans(mvcmmperrb$m,na.rm=T)-rowMeans(lmerrb$m,na.rm=T),
        rowMeans(rferr$m,na.rm=T)-rowMeans(lmerr$m,na.rm=T),rowMeans(rferrw$m,na.rm=T)-rowMeans(lmerrw$m,na.rm=T),rowMeans(rferrb$m,na.rm=T)-rowMeans(lmerrb$m,na.rm=T),
        
        names=nams[-(10:12)],main='M value',col=col,ylab='MSPE Difference',mgp=c(4,1.5,0))

holders=c(rowSEs(cmmperr$m-lmerr$m,na.rm=T),rowSEs(cmmperrw$m-lmerrw$m,na.rm=T),rowSEs(cmmperrb$m-lmerrb$m,na.rm=T),
          rowSEs(mvcmmperr$m-lmerr$m,na.rm=T),rowSEs(mvcmmperrw$m-lmerrw$m,na.rm=T),rowSEs(mvcmmperrb$m-lmerrb$m,na.rm=T),
          rowSEs(rferr$m-lmerr$m,na.rm=T),rowSEs(rferrw$m-lmerrw$m,na.rm=T),rowSEs(rferrb$m-lmerrb$m,na.rm=T))
holder=c(rowMeans(cmmperr$m,na.rm=T)-rowMeans(lmerr$m,na.rm=T),rowMeans(cmmperrw$m,na.rm=T)-rowMeans(lmerrw$m,na.rm=T),rowMeans(cmmperrb$m,na.rm=T)-rowMeans(lmerrb$m,na.rm=T),
         rowMeans(mvcmmperr$m,na.rm=T)-rowMeans(lmerr$m,na.rm=T),rowMeans(mvcmmperrw$m,na.rm=T)-rowMeans(lmerrw$m,na.rm=T),rowMeans(mvcmmperrb$m,na.rm=T)-rowMeans(lmerrb$m,na.rm=T),
         rowMeans(rferr$m,na.rm=T)-rowMeans(lmerr$m,na.rm=T),rowMeans(rferrw$m,na.rm=T)-rowMeans(lmerrw$m,na.rm=T),rowMeans(rferrb$m,na.rm=T)-rowMeans(lmerrb$m,na.rm=T)
         )

holders_scaled=holders/max(holders,na.rm=T)*.25

for(i in 1:(9*2402)){
  arrows(floor(i/vars)+1-holders_scaled[i],holder[i],floor(i/vars)+1+holders_scaled[i],holder[i],code=1,angle=90,length=.0,col='grey50',lwd=.5)
}

pos=NULL
sca=NULL
for(i in 1:9){
  pos=c(pos,((1:9)[i]+c(-.25,-.125,.125,.25)))
  sca=c(sca,0+c(-max(holders,na.rm=T),-.5*max(holders,na.rm=T),.5*max(holders,na.rm=T),max(holders,na.rm=T)))
}
mtext('SE',at=c(0.1))
axis(3,at=pos,round(sca,1),cex.axis=.6,las=2)
abline(h=0)
dev.off()



#RAU value

pdf('rauval_mv_diff.pdf',width=15,height=4)
par(mar=c(4.5, 4.1, 6.5, 1.0))
vioplot(rowMeans(cmmperr$rau,na.rm=T)-rowMeans(lmerr$rau,na.rm=T),rowMeans(cmmperrw$rau,na.rm=T)-rowMeans(lmerrw$rau,na.rm=T),rowMeans(cmmperrb$rau,na.rm=T)-rowMeans(lmerrb$rau,na.rm=T),
        rowMeans(mvcmmperr$rau,na.rm=T)-rowMeans(lmerr$rau,na.rm=T),rowMeans(mvcmmperrw$rau,na.rm=T)-rowMeans(lmerrw$rau,na.rm=T),rowMeans(mvcmmperrb$rau,na.rm=T)-rowMeans(lmerrb$rau,na.rm=T),
        rowMeans(rferr$rau,na.rm=T)-rowMeans(lmerr$rau,na.rm=T),rowMeans(rferrw$rau,na.rm=T)-rowMeans(lmerrw$rau,na.rm=T),rowMeans(rferrb$rau,na.rm=T)-rowMeans(lmerrb$rau,na.rm=T),
        
        names=nams[-(10:12)],main='RAU value',col=col,ylab='MSPE Difference',mgp=c(4,1.5,0))

holders=c(rowSEs(cmmperr$rau-lmerr$rau,na.rm=T),rowSEs(cmmperrw$rau-lmerrw$rau,na.rm=T),rowSEs(cmmperrb$rau-lmerrb$rau,na.rm=T),
          rowSEs(mvcmmperr$rau-lmerr$rau,na.rm=T),rowSEs(mvcmmperrw$rau-lmerrw$rau,na.rm=T),rowSEs(mvcmmperrb$rau-lmerrb$rau,na.rm=T),
          rowSEs(rferr$rau-lmerr$rau,na.rm=T),rowSEs(rferrw$rau-lmerrw$rau,na.rm=T),rowSEs(rferrb$rau-lmerrb$rau,na.rm=T))
holder=c(rowMeans(cmmperr$rau,na.rm=T)-rowMeans(lmerr$rau,na.rm=T),rowMeans(cmmperrw$rau,na.rm=T)-rowMeans(lmerrw$rau,na.rm=T),rowMeans(cmmperrb$rau,na.rm=T)-rowMeans(lmerrb$rau,na.rm=T),
         rowMeans(mvcmmperr$rau,na.rm=T)-rowMeans(lmerr$rau,na.rm=T),rowMeans(mvcmmperrw$rau,na.rm=T)-rowMeans(lmerrw$rau,na.rm=T),rowMeans(mvcmmperrb$rau,na.rm=T)-rowMeans(lmerrb$rau,na.rm=T),
         rowMeans(rferr$rau,na.rm=T)-rowMeans(lmerr$rau,na.rm=T),rowMeans(rferrw$rau,na.rm=T)-rowMeans(lmerrw$rau,na.rm=T),rowMeans(rferrb$rau,na.rm=T)-rowMeans(lmerrb$rau,na.rm=T)
)

holders_scaled=holders/max(holders,na.rm=T)*.25

for(i in 1:(9*2402)){
  arrows(floor(i/vars)+1-holders_scaled[i],holder[i],floor(i/vars)+1+holders_scaled[i],holder[i],code=1,angle=90,length=.0,col='grey50',lwd=.5)
}

pos=NULL
sca=NULL
for(i in 1:9){
  pos=c(pos,((1:9)[i]+c(-.25,-.125,.125,.25)))
  sca=c(sca,0+c(-max(holders,na.rm=T),-.5*max(holders,na.rm=T),.5*max(holders,na.rm=T),max(holders,na.rm=T)))
}
mtext('SE',at=c(0.1))
axis(3,at=pos,round(sca,1),cex.axis=.6,las=2)
abline(h=0)
dev.off()



#beta value
pdf('betaval_mv_diff.pdf',width=15,height=4)
par(mar=c(4.5, 4.1, 6.5, 1.0))
vioplot(rowMeans(cmmperr$beta,na.rm=T)-rowMeans(lmerr$beta,na.rm=T),rowMeans(cmmperrw$beta,na.rm=T)-rowMeans(lmerrw$beta,na.rm=T),rowMeans(cmmperrb$beta,na.rm=T)-rowMeans(lmerrb$beta,na.rm=T),
        rowMeans(mvcmmperr$beta,na.rm=T)-rowMeans(lmerr$beta,na.rm=T),rowMeans(mvcmmperrw$beta,na.rm=T)-rowMeans(lmerrw$beta,na.rm=T),rowMeans(mvcmmperrb$beta,na.rm=T)-rowMeans(lmerrb$beta,na.rm=T),
        rowMeans(rferr$beta,na.rm=T)-rowMeans(lmerr$beta,na.rm=T),rowMeans(rferrw$beta,na.rm=T)-rowMeans(lmerrw$beta,na.rm=T),rowMeans(rferrb$beta,na.rm=T)-rowMeans(lmerrb$beta,na.rm=T),
        
        names=nams[-(10:12)],main='Beta value',col=col,ylab='MSPE Difference',mgp=c(4,1.5,0))

holders=c(rowSEs(cmmperr$beta-lmerr$beta,na.rm=T),rowSEs(cmmperrw$beta-lmerrw$beta,na.rm=T),rowSEs(cmmperrb$beta-lmerrb$beta,na.rm=T),
          rowSEs(mvcmmperr$beta-lmerr$beta,na.rm=T),rowSEs(mvcmmperrw$beta-lmerrw$beta,na.rm=T),rowSEs(mvcmmperrb$beta-lmerrb$beta,na.rm=T),
          rowSEs(rferr$beta-lmerr$beta,na.rm=T),rowSEs(rferrw$beta-lmerrw$beta,na.rm=T),rowSEs(rferrb$beta-lmerrb$beta,na.rm=T))
holder=c(rowMeans(cmmperr$beta,na.rm=T)-rowMeans(lmerr$beta,na.rm=T),rowMeans(cmmperrw$beta,na.rm=T)-rowMeans(lmerrw$beta,na.rm=T),rowMeans(cmmperrb$beta,na.rm=T)-rowMeans(lmerrb$beta,na.rm=T),
         rowMeans(mvcmmperr$beta,na.rm=T)-rowMeans(lmerr$beta,na.rm=T),rowMeans(mvcmmperrw$beta,na.rm=T)-rowMeans(lmerrw$beta,na.rm=T),rowMeans(mvcmmperrb$beta,na.rm=T)-rowMeans(lmerrb$beta,na.rm=T),
         rowMeans(rferr$beta,na.rm=T)-rowMeans(lmerr$beta,na.rm=T),rowMeans(rferrw$beta,na.rm=T)-rowMeans(lmerrw$beta,na.rm=T),rowMeans(rferrb$beta,na.rm=T)-rowMeans(lmerrb$beta,na.rm=T)
)

holders_scaled=holders/max(holders,na.rm=T)*.25

for(i in 1:(9*2402)){
  arrows(floor(i/vars)+1-holders_scaled[i],holder[i],floor(i/vars)+1+holders_scaled[i],holder[i],code=1,angle=90,length=.0,col='grey50',lwd=.5)
}

pos=NULL
sca=NULL
for(i in 1:9){
  pos=c(pos,((1:9)[i]+c(-.25,-.125,.125,.25)))
  sca=c(sca,0+c(-max(holders,na.rm=T),-.5*max(holders,na.rm=T),.5*max(holders,na.rm=T),max(holders,na.rm=T)))
}
mtext('SE',at=c(0.1))
axis(3,at=pos,round(sca,3),cex.axis=.6,las=2)
abline(h=0)
dev.off()


