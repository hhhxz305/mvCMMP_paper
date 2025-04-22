library(asreml)
mvCMMP = 
  
  function (f.in2, r.in, train, x.new, y.new, x.fut = NULL, match.train = F, 
            a1 = NULL, n.new = NULL, interval = NULL, tolparinv=NULL,comps=2,...) {
    
    
    
    train[,all.vars(r.in)]=as.factor(train[,all.vars(r.in)])
    
    #mod = mmer(f.in,r.in,data=train,tolparinv=tolparinv)    # fit mv linear mixed model using sommer package
    
    # trim input texts
    yhold=trimws(gsub("[\r\n]", "",as.character(f.in2))[2])
    yhold=trimws(strsplit(substr(yhold,7,nchar(yhold)-1),',')[[1]])
    xhold=trimws(gsub("[\r\n]", "",as.character(f.in2))[3])
    xhold=trimws(strsplit(xhold,'\\+')[[1]])
    
    #cut trees, no need here
    # varclust=cvtree(corclust(train[,yhold],),k = clusts)$cluster
    # modnum=length(unique(varclust))
    modlist=list()
    varclust=rep(1,length(yhold))
    singleouts=names(table(varclust))[table(varclust)==1]
    
    xmat=model.matrix(~-1+.,train[,xhold,drop=F])
    
    names(varclust)=yhold
    yvars = yhold
    # for(i in 1:length(unique(varclust))){
    #   yvars=names(varclust)[which(varclust==unique(varclust)[i])]
    #   if(length(yvars)==1){
    #     f.in=as.formula(paste(yvars,'~',paste(xhold,collapse='+')))
    #     modlist[[i]]=cmmp(f.in,r.in=as.formula(paste('~1|',all.vars(r.in))),train=data.frame(train[,c(yvars,all.vars(r.in))],xmat),y.new=y.new[,yvars],x.new=x.new,match.train = match.train,n.new = n.new,x.fut = x.fut)
    #   } else{
    if(comps<2){
      print('Error, 2 or more components are required')
      return(NULL)#error due to less than 2 components used function does not work
    }
    #test=kCompRand(as.matrix(train[,yvars]),rep('gaussian',length(yvars)),X=as.matrix(xmat),random=as.factor(train[,all.vars(r.in)]),k=comps,
      #                     method = SCGLR::methodSR("vpi", l = 4, s = 1/2, maxiter = 1000, epsilon =
      #                                                10^-3, bailout = 1000))
    #modlist[[1]]=asreml(fixed=f.in2,random=~(group),data=train)#,residual = ~id(group):us(group))
    train[,all.vars(r.in)]=as.factor(train[,all.vars(r.in)])
    
    f.new=as.formula(paste('cbind(',paste(yhold,collapse=','),')~',paste(c('trait',paste0('trait:',xhold)),collapse ='+')))
    asreml.options(tol=c(0,0),debug=F)
    modlist[[1]]=asreml(fixed=f.new,random=as.formula(paste0('~',all.vars(r.in),':us(trait)')),data=train,residual = ~units:us(trait))#as.formula(paste0('~',all.vars(r.in),':us(trait)')))
    
    #extract outcome names
    
    
    outs=yhold  # extract each outcome name  
    mgroups=c('popmean',sort(unique(train[,all.vars(r.in)]))) # holder for names of popmean and cluster names, note that new method sorts groups
    
    group=matrix(NA,nrow=dim(x.new)[1],ncol=length(outs)) # holder for the identification of groups
    yhat=matrix(NA,nrow=dim(x.new)[1],ncol=length(outs)) #holder for predicted mixed effect
    lf = as.formula(paste('~',paste(xhold,collapse='+'))) #keep only the predictors from f.in (remove outcome names)
    for(i in 1:dim(x.new)[1]){
      #xs=model.matrix(as.formula(paste('-1+',as.character(formula(delete.response(terms(f.in2)))))),data=x.new[i,xhold]) #for each new obs, formatting to make it align with multiplication with coefficients, also make it robust to variable order
      
      
      library(caret)
      dmy <- dummyVars(formula(delete.response(terms(f.in2))), data = x.new[i,xhold])
      xs <- data.frame(int=rep(1,length(i)),predict(dmy, newdata = x.new[i,xhold]))
      
      xhold2=colnames(xs)
      
      outholder=NULL
      for(j in 1:length(outs)){
        if(sum(varclust[yhold[j]]==varclust)>1){
          #matching coefs
          ma=rownames(modlist[[varclust[names(varclust)==outs[j]]]]$coefficients$fixed)[grep(outs[j],rownames(modlist[[varclust[names(varclust)==outs[j]]]]$coefficients$fixed))]
          sub=grep(outs[j],rownames(modlist[[varclust[names(varclust)==outs[j]]]]$coefficients$fixed))
          ma2=unlist(lapply(ma,function(i){holder=strsplit(i,':')[[1]][2]}))
          ma2[is.na(ma2)]='int'
          ma3=gsub('_','.',ma2)
          outholder=c(outholder,sum(xs*modlist[[varclust[names(varclust)==outs[j]]]]$coefficients$fixed[sub[match(xhold2,ma3)],1])) # fixed effect for each outcome
        }else{
          outholder=c(outholder,modlist[[varclust[names(varclust)==outs[j]]]]$mixed.pred[i])
        }
      }
      preds=matrix(NA,nrow=1+length(unique(train[,all.vars(r.in)])),ncol=length(outs)) #collect mixed effect for each outcome under each group along with popmean
      for(j in 1:length(outs)){
        if(sum(varclust[yhold[j]]==varclust)>1){
          preds[,j]=c(0,modlist[[varclust[names(varclust)==outs[j]]]]$coefficients$random[grep(yhold[j],rownames(modlist[[varclust[names(varclust)==outs[j]]]]$coefficients$random)),1])+outholder[j]
        }else{
          preds[,j]=rep(outholder[j],1+length(unique(train[,all.vars(r.in)])))
        }
        #populating predictions (plug in value in preds holder: X*beta+Z*U)
      }
      for(j in 1:length(outs)){
        if(is.null(n.new)||n.new==1){
          if(!match.train){
            #opt=which(min((preds[,j]-y.new[i,j])^2)==(preds[,j]-y.new[i,j])^2)
            opt=which(min((preds[,j]^2-2*preds[,j]*y.new[i,j]))==(preds[,j]^2-2*preds[,j]*y.new[i,j]))[1] #for each outcome, identify the optimal I(group id)
            group[i,j]=mgroups[opt] #assign group id based on opt
            if(sum(outs[j]==singleouts)==1){
              group[i,j]=modlist[[varclust[names(varclust)==outs[j]]]]$group.est[i] #currently not used any more
            }
            yhat[i,j]=preds[opt,j] #assign predicted value based on opt
          }else{ #force match to only clusters, popmean is no longer possibility
            mgroupsh=mgroups[-1]
            predsh=preds[-1,]
            opt=which(min((predsh[,j]^2-2*predsh[,j]*y.new[i,j]))==(predsh[,j]^2-2*predsh[,j]*y.new[i,j]))[1] #for each outcome, identify the optimal I(group id)
            group[i,j]=mgroupsh[opt] #assign group id based on opt
            if(sum(outs[j]==singleouts)==1){
              group[i,j]=modlist[[varclust[names(varclust)==outs[j]]]]$group.est[i] #not used anymore
            }
            yhat[i,j]=predsh[opt,j] #assign predicted value based on opt
          }
        }else{#runing group minimization instead of individual
          if(n.new!=length(y.new)|length(y.new)!=dim(x.new)[1]){
            print('Error only one n.new set should be submitted at a time')
            return(NULL)
          }
          if(!match.train){
            #opt=which(min((preds[,j]-y.new[i,j])^2)==(preds[,j]-y.new[i,j])^2)
            opt=which(min((preds[,j]^2-2*preds[,j]*mean(y.new[i,j])))==(preds[,j]^2-2*preds[,j]*mean(y.new[i,j])))[1] #for each outcome, identify the optimal I(group id)
            group[i,j]=mgroups[opt] #assign group id based on opt
            if(sum(outs[j]==singleouts)==1){
              group[i,j]=modlist[[varclust[names(varclust)==outs[j]]]]$group.est[i]
            }
            yhat[i,j]=preds[opt,j] #assign predicted value based on opt
          }else{ #force match to only clusters, popmean is no longer possibility
            mgroupsh=mgroups[-1]
            predsh=preds[-1,]
            opt=which(min((predsh[,j]^2-2*predsh[,j]*mean(y.new[i,j])))==(predsh[,j]^2-2*predsh[,j]*mean(y.new[i,j])))[1] #for each outcome, identify the optimal I(group id)
            group[i,j]=mgroupsh[opt] #assign group id based on opt
            if(sum(outs[j]==singleouts)==1){
              group[i,j]=modlist[[varclust[names(varclust)==outs[j]]]]$group.est[i]
            }
            yhat[i,j]=predsh[opt,j] #assign predicted value based on opt
            break()
          }
        }
        
      }
      
    }
    
    return(list(thetahat=yhat,group=group,mod=modlist,modgroups=varclust)) #return predicted mixed effect and group id 
  }


