####
# PrOCTOR.R
# Description: Run PrOCTOR for any given structure/target pair
# Last Updated: October 20, 2016
# Update: Fixed bug in plotting function. Added predictions for target/structure specific models and a plot all function.
# Functions:
#   1) PrOCTOR(SMILE,targets,type)
#         - applies PrOCTOR to inputted drug
#         - type = type of PrOCTOR model to use: full, target, structure (default=full)
#         - example: Dexamethasone
#             --> PrOCTOR(SMILE="[H][C@@]12C[C@@H](C)[C@](O)(C(=O)CO)[C@@]1(C)C[C@H](O)[C@@]1(F)[C@@]2([H])CCC2=CC(=O)C=C[C@]12C",targets=c("NR3C1","NR0B1","ANXA1","NOS2"))
#
#   2) PrOCTOR_plot(SMILE,targetstype)
#         - applies PrOCTOR and displays a graph that compares the inputted drug to the training set
#
#   2) PrOCTOR_plotall(SMILE,targets)
#         - applies PrOCTOR and displays graphs that compares the inputted drug to the training set
#           for full, target, and structural version of PrOCTOR
#

library(ChemmineR)
library(ChemmineOB)
library(rcdk)
library(Rcpi)
library(randomForest)
library(ggplot2)
require(gridExtra)
library(data.table)

##### Main Functions

PrOCTOR=function(SMILE,targets,type="full"){
  featureset=getAllFeatures(SMILE,targets)
  if(grepl("struct|Struct",type)){
    prob=1-median(sapply(data$tree$struct,function(this_model) predict(this_model,featureset[names(featureset) %in% rownames(data$tree$struct[[1]]$importance)],"prob")[1]  ))
    return(list(pred=c("Safe","Toxic")[as.numeric(prob<0.5)+1],prob=prob,score=log2(prob/(1-prob))))
  }else if(grepl("targ|Targ",type)){
    prob=1-median(sapply(data$tree$target,function(this_model) predict(this_model,featureset[names(featureset) %in% rownames(data$tree$target[[1]]$importance)],"prob")[1]  ))
    return(list(pred=c("Safe","Toxic")[as.numeric(prob<0.5)+1],prob=prob,score=log2(prob/(1-prob))))    
  }
  prob=1-median(sapply(data$tree$full,function(this_model) predict(this_model,featureset[!names(featureset) %in% c("Ro5","Ghose","Veber")],"prob")[1]  ))
  return(list(pred=c("Safe","Toxic")[as.numeric(prob<0.5)+1],prob=prob,score=log2(prob/(1-prob))))
}


PrOCTOR_plot=function(SMILE,targets,type="full"){
  titlemap=data.table(val=c("struct","Struct","target","Target","full"),map=c("Structural Model","Structural Model","Target Model","Target Model","Full Model"))
  res=PrOCTOR(SMILE,targets,type)
  if(grepl("struct|Struct",type)){
    trainingset_probs=apply(sapply(1:5,function(i) data$tree$struct[[i]]$votes[,2]),1,median)
  }else if(grepl("targ|Targ",type)){
    trainingset_probs=apply(sapply(1:5,function(i) data$tree$target[[i]]$votes[,2]),1,median)
  }else{
    trainingset_probs=apply(sapply(1:5,function(i) data$tree$full[[i]]$votes[,2]),1,median)
  }
  
  trainingset_scores=data.table(score=log2(trainingset_probs/(1-trainingset_probs)),class=data$tree$full[[1]]$y)
  g=ggplot(trainingset_scores,aes(x=(score),fill=class))+geom_density(alpha=0.5)+xlab("")+scale_fill_manual(values=c("red","skyblue"),name="Training Set Class")+xlab("PrOCTOR Score (- log2 Odds of Toxicity)")+theme_bw()+ theme(legend.position=c(0.25,0.8),title =element_text(size=12, face='bold'),axis.title = element_text(size=11,face='bold'),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ geom_segment(aes(x = res$score, y = 0.33, xend = res$score, yend = 0), arrow = arrow(length = unit(0.5, "cm")))+ggtitle(paste0(titlemap$map[which(sapply(titlemap$val,function(x) grepl(x,type)))],"\nPredicted ",res$pred))+ annotate("text", x = res$score, y = 0.35, label = paste0("PrOCTOR Score = ",round(abs(res$score),5)))+annotate("text",x=-6.75,y=-0.01,label="toxic",col="red")+annotate("text",x=4.25,y=-0.01,label="safe",col="skyblue3")+ theme(legend.background = element_rect( fill = NA)) 
  if(grepl("struct|Struct|target|Target",type)){
   g=g+theme(legend.position = "none") 
  }
  return(g)
}

PrOCTOR_plotall=function(SMILE,targets){
    g1=PrOCTOR_plot(SMILE,targets,"full")
    g2=PrOCTOR_plot(SMILE,targets,"target")
    g3=PrOCTOR_plot(SMILE,targets,"structure")
    grid.arrange(g1,g2,g3, ncol=3)
}

##### Sub-Functions

# QED Sub-Functions
ads=function(x,params){
  a=params$A
  b=params$B
  c=params$C
  d=params$D
  e=params$E
  f=params$F
  dx_max=params$DMAX
  return((a+(b/(1+exp(-1*(x-c+d/2)/e))*(1-1/(1+exp(-1*(x-c-d/2)/f)))))/dx_max)
}

ads_params=function(var){
  if(var=="MW"){
    return(list(A=2.817065973,B=392.5754953,C=290.7489764,D=2.419764353,E=49.22325677,F=65.37051707,DMAX=104.98055614))
  }
  if(var=="ALOGP"){
    return(list(A=3.172690585,B=137.8624751,C=2.534937431,D=4.581497897,E=0.822739154,F=0.576295591,DMAX=131.31866035))
  }
  if(var=="HBA"){
    return(list(A=2.948620388,B=160.4605972,C=3.615294657,D=4.435986202,E=0.290141953,F=1.300669958,DMAX=148.77630464))
  }
  if(var=="HBD"){
    return(list(A=1.618662227,B=1010.051101,C=0.985094388,D=0.000000000001,E=0.713820843,F=0.920922555,DMAX=258.16326158))
  }
  if(var=="PSA"){
    return(list(A=1.876861559,B=125.2232657,C=62.90773554,D=87.83366614,E=12.01999824,F=28.51324732,DMAX=104.56861672))
  }
  if(var=="ROTB"){
    return(list(A=0.01,B=272.4121427,C=2.558379970,D=1.565547684,E=1.271567166,F=2.758063707,DMAX=105.44204028))
  }
  if(var=="AROM"){
    return(list(A=3.217788970,B=957.7374108,C=2.274627939,D=0.000000000001,E=1.317690384,F=0.375760881,DMAX=312.33726097))
  }
  if(var=="ALERTS"){
    return(list(A=0.01,B=1199.094025,C=-0.09002883,D=0.000000000001,E=0.185904477,F=0.875193782,DMAX=417.72531400))
  } 
  return(NA)
}

calculate_QED=function(drug){
  weights=c(MW=0.66,ALOGP=0.46,HBA=0.05,HBD=0.61,PSA=0.06,ROTB=0.65,AROM=0.48,ALERTS=0.95)
  return(list(uQED=exp(sum(sapply(1:length(drug),function(i) log(ads(drug[i],ads_params(names(drug)[i])))))/length(weights)),wQED=exp(sum(weights * sapply(1:length(drug),function(i) log(ads(drug[i],ads_params(names(drug)[i])))))/sum(weights))))  
}

getStructuralFeatures=function(SMILE){
  mol=smiles2sdf(as(as.vector(SMILE),"SMIset"))
  props=propOB(mol)
  
  MW=props$MW # molecular weight
  XlogP=props$logP # octanol-water partition coefficient log P
  HBD=props$HBD #hydrogen bond donor count
  HBA=props$HBA2 #hydrogen bond acceptor count
  PSA=props$TPSA #polar surface area
  FC=sum(bonds(mol, type="bonds")[[1]]$charge) #formal charge
  RBC=extractDrugRotatableBondsCount(parse.smiles(as.vector(SMILE)))[[1]] #rotatable bonds count
  refr=props$MR # refractivity
  alogP=NA #
  nA=sum(atomcountMA(mol,addH=FALSE)) #number atoms
  AROMs=as.vector(rings(mol, type="count", arom=TRUE)[2])
  nALERTS=sum(sapply(data$unwantedALERTS,function(x) as.numeric(smartsSearchOB(mol,x,uniqueMatches = TRUE)!=0)))
  
  Ro5=as.numeric(MW<500&HBD<5&HBA<10&XlogP<5)
  Veber=as.numeric(RBC<=10&PSA<=140)
  Ghose=as.numeric(PSA<140&(findInterval(XlogP,c(-0.4,5.6))==1)&(findInterval(MW,c(160,480))==1)&(findInterval(nA,c(20,70))==1))
  QED=calculate_QED(c(MW=MW,ALOGP=XlogP,HBA=HBA,HBD=HBD,PSA=PSA,ROTB=RBC,AROM=AROMs,ALERTS=nALERTS)) 
  
  return(c(MolecularWeight=MW,XLogP=XlogP,HydrogenBondDonorCount=HBD,HydrogenBondAcceptorCount=HBA,PolarSurfaceArea=PSA,FormalCharge=FC,NumRings=AROMs,RotatableBondCount=RBC,Refractivity=refr,Ro5=Ro5,Ghose=Ghose,Veber=Veber,wQED=QED$wQED))
}
get_PC_value<-function(this.targets){
  expr=data$target_data$expr[rownames(data$target_data$expr) %in% this.targets,]
  if(!is.null(nrow(expr))){
    expr_features=apply(expr,2,median)
  }else{
    expr_features=expr
  }
  expr_features2=matrix(expr_features,nrow=1);colnames(expr_features2)=names(expr_features)
  pc_preds=predict(data$target_data$pc,expr_features2)
  pc_features=c(PC1=pc_preds[1],PC2=pc_preds[2],PC3=pc_preds[3])
  return(pc_features)
}


getTargetFeatures=function(this.targets){
  return(c(lossFreq=max(data$target_data$lof$deleterious[rownames(data$target_data$lof) %in% this.targets]/data$target_data$lof$Total[rownames(data$target_data$lof) %in% this.targets]),maxBtwn=max(data$target_data$btwn[names(data$target_data$btwn) %in% this.targets]),maxDegree=max(data$target_data$degree[names(data$target_data$degree) %in% this.targets]),get_PC_value(this.targets)))
}

getRotatableBondCount=function(SMILE,mol){
  out=tryCatch({RBC=extractDrugRotatableBondsCount(parse.smiles(as.vector(SMILE)))[[1]] },error=function(cond){return(as.numeric(smartsSearchOB(mol,"[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]")))})
  return(out)
}

getAllFeatures=function(SMILE,targets){
  mol=smiles2sdf(as(as.vector(SMILE),"SMIset"))
  props=propOB(mol)
  
  MW=props$MW # molecular weight
  XlogP=props$logP # octanol-water partition coefficient log P
  HBD=props$HBD #hydrogen bond donor count
  HBA=props$HBA2 #hydrogen bond acceptor count
  PSA=props$TPSA #polar surface area
  FC=sum(bonds(mol, type="bonds")[[1]]$charge) #formal charge
  RBC=getRotatableBondCount(SMILE,mol) #rotatable bonds count
  refr=props$MR # refractivity
  alogP=NA #
  nA=sum(atomcountMA(mol,addH=FALSE)) #number atoms
  AROMs=as.vector(rings(mol, type="count", arom=TRUE)[2])
  nALERTS=sum(sapply(data$unwantedALERTS,function(x) as.numeric(smartsSearchOB(mol,x,uniqueMatches = TRUE)!=0)))
  
  Ro5=as.numeric(MW<500&HBD<5&HBA<10&XlogP<5)
  Veber=as.numeric(RBC<=10&PSA<=140)
  Ghose=as.numeric(PSA<140&(findInterval(XlogP,c(-0.4,5.6))==1)&(findInterval(MW,c(160,480))==1)&(findInterval(nA,c(20,70))==1))
  QED=calculate_QED(c(MW=MW,ALOGP=XlogP,HBA=HBA,HBD=HBD,PSA=PSA,ROTB=RBC,AROM=AROMs,ALERTS=nALERTS)) 
  
  return(c(MolecularWeight=MW,XLogP=XlogP,HydrogenBondDonorCount=HBD,HydrogenBondAcceptorCount=HBA,PolarSurfaceArea=PSA,FormalCharge=FC,NumRings=AROMs,RotatableBondCount=RBC,Refractivity=refr,lossFreq=max(data$target_data$lof$deleterious[rownames(data$target_data$lof) %in% targets]/data$target_data$lof$Total[rownames(data$target_data$lof) %in% targets]),maxBtwn=max(data$target_data$btwn[names(data$target_data$btwn) %in% targets]),maxDegree=max(data$target_data$degree[names(data$target_data$degree) %in% targets]),Ro5=Ro5,Ghose=Ghose,Veber=Veber,wQED=QED$wQED,get_PC_value(targets)))
}

data=readRDS(paste0(dirname(sys.frame(1)$ofile),"/PrOCTOR.rds"))
