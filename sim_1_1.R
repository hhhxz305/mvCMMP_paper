library(MASS)
library(cmmp)
library(gimme)
library(splitstackshape)

########################### General Setting ####################################
simsetting='_1_1'
dir='...' # set directory
k=4
n_i = 5
m = 100
N = n_i*m

set.seed(1518)
#XX = cbind(rnorm(N),rnorm(N))
beta = matrix(c(5,1,-1,3,2,7,5,6),nrow=k)  
beta_L = as.vector(beta)
n.train = .7
########################### LOOPING ############################################ 
n.loop = 100
result.hold_cmmp.y1 = data.frame(matrix(NA,nrow=n.loop,ncol = 1))
result.hold_cmmp.y2 = data.frame(matrix(NA,nrow=n.loop,ncol = 1))
result.hold_cmmp.y3 = data.frame(matrix(NA,nrow=n.loop,ncol = 1))
result.hold_cmmp.y4 = data.frame(matrix(NA,nrow=n.loop,ncol = 1))
result.hold_mvcmmp.y1 = data.frame(matrix(NA,nrow=n.loop,ncol = 1))
result.hold_mvcmmp.y2 = data.frame(matrix(NA,nrow=n.loop,ncol = 1))
result.hold_mvcmmp.y3 = data.frame(matrix(NA,nrow=n.loop,ncol = 1))
result.hold_mvcmmp.y4 = data.frame(matrix(NA,nrow=n.loop,ncol = 1))
result.hold_LR.y1 = data.frame(matrix(NA,nrow=n.loop,ncol = 1))
result.hold_LR.y2 = data.frame(matrix(NA,nrow=n.loop,ncol = 1))
result.hold_LR.y3 = data.frame(matrix(NA,nrow=n.loop,ncol = 1))
result.hold_LR.y4 = data.frame(matrix(NA,nrow=n.loop,ncol = 1))
result.hold_YY = data.frame(matrix(NA,nrow=n.loop,ncol = 1))
result.hold_bb = data.frame(matrix(NA,nrow=n.loop,ncol = 1))
result.hold_ee = data.frame(matrix(NA,nrow=n.loop,ncol = 1))

random_dist_mvcmmp=list()
random_dist_cmmp=list()  
random_dist_mvcmmp_r=list()
random_dist_cmmp_r=list()  
random_dist_mvcmmp_r2=list()
random_dist_cmmp_r2=list()  

mvcmmp_rtime=rep(NA,n.loop)

set.seed(10011)
for (z in 1:n.loop) {
  step=1
  
  for (k in 1:1) {
    
    
    #####################################################
    g1=1
    g3 =  1
    g2=1
    g4 = 1
    
    rho_g12 =.1
    rho_g13 =.1
    rho_g14 =.1
    rho_g23 =.1
    rho_g24 =.1
    rho_g34 =.1
    
    G = matrix(c(g1,sqrt(g1)*sqrt(g2)*rho_g12,sqrt(g1)*sqrt(g3)*rho_g13,sqrt(g1)*sqrt(g4)*rho_g14,
                 sqrt(g1)*sqrt(g2)*rho_g12,g2,sqrt(g2)*sqrt(g3)*rho_g23,sqrt(g2)*sqrt(g4)*rho_g24,
                 sqrt(g1)*sqrt(g3)*rho_g13,sqrt(g2)*sqrt(g3)*rho_g23,g3,sqrt(g3)*sqrt(g4)*rho_g34,
                 sqrt(g1)*sqrt(g4)*rho_g14,sqrt(g2)*sqrt(g4)*rho_g24,sqrt(g3)*sqrt(g4)*rho_g34,g4),byrow=T,nrow=4)
    B = mvrnorm(m,mu=c(0,0,0,0),Sigma = G) #generate m group of intercepts
    #cor(B)
    bb = expandRows(B, count=n_i,count.is.col = F)
    b_L = as.vector(bb)
    
    #######################################################
    r1=1
    r2 =  2
    r3=4
    r4 = 5
    
    rho_r12 =.2
    rho_r13 =.4
    rho_r14 =.4
    rho_r23 =.7
    rho_r24 =.5
    rho_r34 =.8
    
    R = matrix(c(r1,sqrt(r1)*sqrt(r2)*rho_r12,sqrt(r1)*sqrt(r3)*rho_r13,sqrt(r1)*sqrt(r4)*rho_r14,
                 sqrt(r1)*sqrt(r2)*rho_r12,r2,sqrt(r2)*sqrt(r3)*rho_r23,sqrt(r2)*sqrt(r4)*rho_r24,
                 sqrt(r1)*sqrt(r3)*rho_r13,sqrt(r2)*sqrt(r3)*rho_r23,r3,sqrt(r3)*sqrt(r4)*rho_r34,
                 sqrt(r1)*sqrt(r4)*rho_r14,sqrt(r2)*sqrt(r4)*rho_r24,sqrt(r3)*sqrt(r4)*rho_r34,r4),byrow=T,nrow=4)
    #ee = mvrnorm(N,mu=c(0,0),Sigma = R)
    #ee = matrix(NA,nrow=N,ncol=2)
    ee = data.frame()
    for (i in 1:m){
      b_e = mvrnorm(n_i,mu=c(0,0,0,0),Sigma = R)
      ee = rbind(ee,b_e) 
    }
    e_L = as.vector(as.matrix(ee))
    
    #######################################################
    XX = cbind(rnorm(N),rnorm(N))
    Y_L = kronecker(diag(4),XX)%*%beta_L + b_L + e_L
    YY = cbind(y1=Y_L[1:N],y2=Y_L[(N+1):(2*N)],y2=Y_L[(2*N+1):(3*N)],y2=Y_L[(3*N+1):(4*N)])
    
    group = rep(paste('group',1:m,sep = ' '),each=n_i)
    dat.total = data.frame(YY,XX,group=group)
    names(dat.total)[1:4]=c('y1','y2','y3','y4')
    ########Sampling for training and test dataset###################
    train.idx = sample(N,N*n.train)
    dat = dat.total[train.idx,]
    dat2 = dat.total[-train.idx,]
    
    ###################generating mixed effect in the testing set###############
    #dat3=dat2
    #dat3[,1:2] = dat2[,1:2] - err[-train.idx,]
    
    #extracting unique random effects, real ones are bb
    ran_eff=unique(bb)
    
    deviance_cmmp=matrix(NA,nrow=dim(dat2)[1],ncol=4)
    deviance_mvcmmp=matrix(NA,nrow=dim(dat2)[1],ncol=4)
    deviance_cmmp_r=matrix(NA,nrow=dim(dat2)[1],ncol=4)
    deviance_mvcmmp_r=matrix(NA,nrow=dim(dat2)[1],ncol=4)
    deviance_cmmp_r2=matrix(NA,nrow=dim(dat2)[1],ncol=4)
    deviance_mvcmmp_r2=matrix(NA,nrow=dim(dat2)[1],ncol=4)
    ###########################################mvCMMP############################################
    w=NULL
    f.in = as.formula(paste0('cbind(',paste0(names(dat[,1:4]),collapse = ','),')~', paste(names(dat)[5:6],collapse = '+')))
    #w<-mvCMMP(f.in,~ group,dat,y.new=dat2[,1:2],x.new=dat2[,3:4],tolparinv=.5,match.train =T,clusts=1,comps=2)
    
    try(tm<-system.time(w<-mvCMMP(f.in2=f.in,r.in=~group,train=dat,y.new=dat2[,1:4],x.new=dat2[,5:6],tolparinv=1e-6,match.train =T,clusts=1,comps=2)))
    if(is.null(w)){
      try(tm<-system.time(w<-mvCMMP(f.in,~group,dat,y.new=dat2[,1:4],x.new=dat2[,5:6],tolparinv=0.01,match.train =T,clusts=1,comps=2)))
    }
    if(is.null(w)){
      try(tm<-system.time(w<-mvCMMP(f.in,~group,dat,y.new=dat2[,1:4],x.new=dat2[,5:6],tolparinv=0.5,match.train =T,clusts=1,comps=2)))
    }
    #w$thetahat
    if(!is.null(w)){
      ran_eff_holder=matrix(NA,nrow=dim(dat2)[1],ncol=4)
      for(e in 1:4){
        deviance_mvcmmp[,e]=bb[-train.idx,e]-ran_eff[as.numeric(w$group[,e]),e]
        ran_eff_holder[,e]=ran_eff[as.numeric(w$group[,e]),e]
      }
      deviance_mvcmmp_r=w$thetahat-dat2[,1:4]
      deviance_mvcmmp_r2=(w$thetahat-ran_eff_holder)-(dat2[,1:4]-bb[-train.idx,]) #trying to get at error by subtracting off random effects both sim and predicted
    }
    if(!is.null(w)){
      w.MSPE = colMeans((w$thetahat-dat2[,1:4])^2)
      #collecting run time
      mvcmmp_rtime[z]=tm[3]
    }else{
      w.MSPE=rep(NA,4)
    }
    ################################################################################################
    
    
    ############################################CMMP################################################
    r.in <- as.formula(~1|group) 
    f.in1 = as.formula(paste('y1 ~', paste(names(dat)[5:6],collapse = '+')))
    
    c1 = cmmp(f.in1, r.in, train = dat, x.new = dat2[,5:6], y.new = dat2[,1], x.fut = NULL, match.train = T, a1 = NULL, n.new = NULL, interval = TRUE)
    
    group_extract=function(group_pred){
      hold=strsplit(group_pred,' ')
      hold2=rep(NA,length(hold))
      for(i in 1:length(hold)){
        hold2[i]=hold[[i]][2]
      }
      return(as.numeric(hold2))
    }
    
    deviance_cmmp[,1]=bb[-train.idx,1]-ran_eff[group_extract(c1$group.est),1]
    
    f.in2 = as.formula(paste('y2 ~', paste(names(dat)[5:6],collapse = '+')))
    
    c2 = cmmp(f.in2, r.in, train = dat, x.new = dat2[,5:6], y.new = dat2[,2], x.fut = NULL, match.train = T, a1 = NULL, n.new = NULL, interval = TRUE)
    deviance_cmmp[,2]=bb[-train.idx,2]-ran_eff[group_extract(c2$group.est),2]
    
    f.in3 = as.formula(paste('y3 ~', paste(names(dat)[5:6],collapse = '+')))
    
    c3 = cmmp(f.in3, r.in, train = dat, x.new = dat2[,5:6], y.new = dat2[,3], x.fut = NULL, match.train = T, a1 = NULL, n.new = NULL, interval = TRUE)
    deviance_cmmp[,3]=bb[-train.idx,3]-ran_eff[group_extract(c3$group.est),3]
    
    f.in4 = as.formula(paste('y4 ~', paste(names(dat)[5:6],collapse = '+')))
    
    c4 = cmmp(f.in4, r.in, train = dat, x.new = dat2[,5:6], y.new = dat2[,4], x.fut = NULL, match.train = T, a1 = NULL, n.new = NULL, interval = TRUE)
    deviance_cmmp[,4]=bb[-train.idx,4]-ran_eff[group_extract(c4$group.est),4]
    c.MSPE = c(mean((c1$mixed.pred-dat2[,1])^2),mean((c2$mixed.pred-dat2[,2])^2),mean((c3$mixed.pred-dat2[,3])^2),mean((c4$mixed.pred-dat2[,4])^2))
    
    deviance_cmmp_r[,1]=c1$mixed.pred-dat2[,1]
    deviance_cmmp_r[,2]=c2$mixed.pred-dat2[,2]
    deviance_cmmp_r[,3]=c3$mixed.pred-dat2[,3]
    deviance_cmmp_r[,4]=c4$mixed.pred-dat2[,4]
    
    deviance_cmmp_r2[,1]=(c1$mixed.pred-ran_eff[group_extract(c1$group.est),1])-(dat2[,1]-bb[-train.idx,1])
    deviance_cmmp_r2[,2]=(c2$mixed.pred-ran_eff[group_extract(c2$group.est),2])-(dat2[,2]-bb[-train.idx,2])
    deviance_cmmp_r2[,3]=(c3$mixed.pred-ran_eff[group_extract(c3$group.est),3])-(dat2[,3]-bb[-train.idx,3])
    deviance_cmmp_r2[,4]=(c4$mixed.pred-ran_eff[group_extract(c4$group.est),4])-(dat2[,4]-bb[-train.idx,4])
    #################################################################################################
    random_dist_cmmp[[z]]=deviance_cmmp
    random_dist_mvcmmp[[z]]=deviance_mvcmmp
    random_dist_cmmp_r[[z]]=deviance_cmmp_r
    random_dist_mvcmmp_r[[z]]=deviance_mvcmmp_r
    random_dist_cmmp_r2[[z]]=deviance_cmmp_r2
    random_dist_mvcmmp_r2[[z]]=deviance_mvcmmp_r2
    
    #############################################LR###################################################
    lr1 = lm(f.in1,data=dat)
    lr.pred.1 = predict(lr1,dat2[,5:6])
    
    
    lr2 = lm(f.in2,data=dat)
    lr.pred.2 = predict(lr2,dat2)
    
    lr3 = lm(f.in3,data=dat)
    lr.pred.3 = predict(lr3,dat2)
    
    lr4 = lm(f.in4,data=dat)
    lr.pred.4 = predict(lr4,dat2)
    
    lr.MSPE = c(mean((lr.pred.1-dat2[,1])^2),mean((lr.pred.2-dat2[,2])^2),mean((lr.pred.3-dat2[,3])^2),mean((lr.pred.4-dat2[,4])^2))
    
    ##################################################################################################
    
    result.hold_cmmp.y1[z,step]=c.MSPE[1]
    result.hold_cmmp.y2[z,step]=c.MSPE[2]
    result.hold_cmmp.y3[z,step]=c.MSPE[3]
    result.hold_cmmp.y4[z,step]=c.MSPE[4]
    result.hold_mvcmmp.y1[z,step] = w.MSPE[1]
    result.hold_mvcmmp.y2[z,step] = w.MSPE[2]
    result.hold_mvcmmp.y3[z,step] = w.MSPE[3]
    result.hold_mvcmmp.y4[z,step] = w.MSPE[4]
    result.hold_LR.y1[z,step] = lr.MSPE[1]
    result.hold_LR.y2[z,step] = lr.MSPE[2]
    result.hold_LR.y3[z,step] = lr.MSPE[3]
    result.hold_LR.y4[z,step] = lr.MSPE[4]
    
    result.hold_YY[z,step] = cor(YY)[2,1]
    result.hold_bb[z,step] = cor(bb)[2,1]
    result.hold_ee[z,step] = cor(ee)[2,1]
    
    step=step+1 
  }
  
  print(z)
}


dat.R1 = as.data.frame(t(rbind(colMeans(result.hold_LR.y1,na.rm = T),colMeans(result.hold_cmmp.y1,na.rm = T),colMeans(result.hold_mvcmmp.y1,na.rm = T),
                               colMeans(result.hold_LR.y2,na.rm = T),colMeans(result.hold_cmmp.y2,na.rm = T),colMeans(result.hold_mvcmmp.y2,na.rm = T),
                               colMeans(result.hold_LR.y3,na.rm = T),colMeans(result.hold_cmmp.y3,na.rm = T),colMeans(result.hold_mvcmmp.y3,na.rm = T),
                               colMeans(result.hold_LR.y4,na.rm = T),colMeans(result.hold_cmmp.y4,na.rm = T),colMeans(result.hold_mvcmmp.y4,na.rm = T))))
colSE=function(i,na.rm=T){return(apply(i,2,function(j){return(sd(na.omit(j))/length(na.omit(j))^.5)}))}
dat.R2 = as.data.frame(t(rbind(colSE(result.hold_LR.y1,na.rm = T),colSE(result.hold_cmmp.y1,na.rm = T),colSE(result.hold_mvcmmp.y1,na.rm = T),
                               colSE(result.hold_LR.y2,na.rm = T),colSE(result.hold_cmmp.y2,na.rm = T),colSE(result.hold_mvcmmp.y2,na.rm = T),
                               colSE(result.hold_LR.y3,na.rm = T),colSE(result.hold_cmmp.y3,na.rm = T),colSE(result.hold_mvcmmp.y3,na.rm = T),
                               colSE(result.hold_LR.y4,na.rm = T),colSE(result.hold_cmmp.y4,na.rm = T),colSE(result.hold_mvcmmp.y4,na.rm = T))))

#row.names(dat.R1) = apply(RHO_hold2,1,function(i){paste0('G1=',i[1],'|G2=',i[2],'|GG=',i[3],'|R1=',i[4],'|R2=',i[5],'|RR=',i[6])})
colnames(dat.R2)=colnames(dat.R1) = c('LR.Y1','CMMP.Y1','mvCMMP.Y1','LR.Y2','CMMP.Y2','mvCMMP.Y2','LR.Y3','CMMP.Y3','mvCMMP.Y3','LR.Y4','CMMP.Y4','mvCMMP.Y4')

#will take colmeans for each case and then mean over iterations
rowMeans(data.frame(lapply(random_dist_cmmp,function(i){colMeans(abs(i))})),na.rm=T)
rowMeans(data.frame(lapply(random_dist_mvcmmp,function(i){colMeans(abs(i))})),na.rm=T)

library(xtable)
print.xtable(xtable(dat.R1),file=paste0(dir,'dat.R1',simsetting,'.txt'))
print.xtable(xtable(dat.R2),file=paste0(dir,'dat.R2',simsetting,'.txt'))
round(dat.R2,2)

disdf=data.frame(rbind(rowMeans(data.frame(lapply(random_dist_cmmp,function(i){colMeans(abs(i))})),na.rm=T),rowMeans(data.frame(lapply(random_dist_mvcmmp,function(i){colMeans(abs(i))})),na.rm=T)))
rownames(disdf)=c('CMMP','mvCMMP')
print.xtable(xtable(disdf),file=paste0(dir,'disdf',simsetting,'.txt'))

precor_cmmp=lapply(random_dist_cmmp,cor)
precor_mvcmmp=lapply(random_dist_mvcmmp,cor,use='pairwise.complete.obs')


cor_cmmp=matrix(0,nrow=4,ncol=4)
cor_mvcmmp=matrix(0,nrow=4,ncol=4)
mvints=0

for(i in 1:length(precor_mvcmmp)){
  cor_cmmp=cor_cmmp+precor_cmmp[[i]]
  if(!is.na(precor_mvcmmp[[i]][1,1])){
    cor_mvcmmp=cor_mvcmmp+precor_mvcmmp[[i]]
    mvints=mvints+1
  }
}
#take average correlation matrix
print.xtable(xtable(cor_cmmp/n.loop),file=paste0(dir,'cor_cmmp',simsetting,'.txt'))
print.xtable(xtable(cor_mvcmmp/mvints),file=paste0(dir,'cor_mvcmmp',simsetting,'.txt'))

#for r matrix
precor_cmmp_r=lapply(random_dist_cmmp_r,cov)
precor_mvcmmp_r=lapply(random_dist_mvcmmp_r,cov,use='pairwise.complete.obs')


cor_cmmp_r=matrix(0,nrow=4,ncol=4)
cor_mvcmmp_r=matrix(0,nrow=4,ncol=4)
mvints_r=0

for(i in 1:length(precor_mvcmmp_r)){
  cor_cmmp_r=cor_cmmp_r+precor_cmmp_r[[i]]
  if(!is.na(precor_mvcmmp_r[[i]][1,1])){
    cor_mvcmmp_r=cor_mvcmmp_r+precor_mvcmmp_r[[i]]
    mvints_r=mvints_r+1
  }
}
#take average correlation matrix
print.xtable(xtable(cor_cmmp_r/n.loop),file=paste0(dir,'cor_cmmp_r',simsetting,'.txt'))
print.xtable(xtable(cor_mvcmmp_r/mvints_r),file=paste0(dir,'cor_mvcmmp_r',simsetting,'.txt'))


#for r matrix
precor_cmmp_r2=lapply(random_dist_cmmp_r2,cov)
precor_mvcmmp_r2=lapply(random_dist_mvcmmp_r2,cov,use='pairwise.complete.obs')


cor_cmmp_r2=matrix(0,nrow=4,ncol=4)
cor_mvcmmp_r2=matrix(0,nrow=4,ncol=4)
mvints_r2=0

for(i in 1:length(precor_mvcmmp_r2)){
  cor_cmmp_r2=cor_cmmp_r2+precor_cmmp_r2[[i]]
  if(!is.na(precor_mvcmmp_r2[[i]][1,1])){
    cor_mvcmmp_r2=cor_mvcmmp_r2+precor_mvcmmp_r2[[i]]
    mvints_r2=mvints_r2+1
  }
}
#take average correlation matrix
print.xtable(xtable(cor_cmmp_r2/n.loop),file=paste0(dir,'avg_cor_cmmp',simsetting,'.txt'))
print.xtable(xtable(cor_mvcmmp_r2/mvints_r2),file=paste0(dir,'avg_cor_mvcmmp',simsetting,'.txt'))

tab=as.table(summary(mvcmmp_rtime))
print.xtable(xtable(rbind(tab,tab)),file=paste0(dir,'runtime_mvcmmp',simsetting,'.txt')) #summarizing run time

save.image(file=paste0(dir,'image',simsetting,'.rdata'))