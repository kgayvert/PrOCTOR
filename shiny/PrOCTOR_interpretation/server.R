library(shiny)
library(shinyjs)
library(data.table)
library(plyr)
library(htmlwidgets)
library(ggplot2)
library(randomForest)
library(grid)
require(gridExtra)
library(Cairo)
load("initial_values.RData")

shinyServer(function(input, output) {
  
  ######################
  # Reactive DataFrame #
  ######################    
  df_data=reactiveValues(features=initial_feature_df,predictions=initial_prediction_df,structureFeatures=initial_feature_structure,targetFeatures=initial_feature_target,corFlag=FALSE,currentDrug="Median Values")
  
  ######################
  # Reactive Functions #
  ######################    
  observeEvent(input$correlates, {
    df_data$corFlag=input$correlates
  })
  
  observeEvent(input$goButton, {
    df_data$features=changeValueButton()
    df_data$predictions=updatePred() 
    df_data$structureFeatures=updateStructureFeature()
    df_data$targetFeatures=updateTargetFeature()
  })
  
  observeEvent(input$plot_click, {
    df_data$features=changeValueClick()
    df_data$predictions=updatePred() 
    df_data$structureFeatures=updateStructureFeature()
    df_data$targetFeatures=updateTargetFeature()
  })
  
  observeEvent(input$preselectButton, {
    df_data$features=userPreSelection()
    df_data$predictions=updatePred() 
    df_data$structureFeatures=updateStructureFeature()
    df_data$targetFeatures=updateTargetFeature()
  })
  
  changeValueButton <- function(){
    return(changeValue(this_var=input$var,newval=as.numeric(input$newval)))
  }

  changeValueClick <- function(){
    click_update=click_xy(input$plot_click)
    if(!is.null(click_update)){
      featureInd=round(click_update[1])
      new_q=click_update[2]
      if(new_q>1){new_q=1}
      if(new_q<0){new_q=0}      
      return(changeValue(this_var=colnames(sample_drugs)[featureInd],newval=quantile(sample_drugs[,featureInd],new_q)))      
    }else{
      return(NULL)
    }
  }
  
  userPreSelection <- function(){
    if((input$preselectButton == 0)|(input$drug==df_data$currentDrug)) {return (NULL)}
    feature_df=df_data$features
    df_data$currentDrug=input$drug
    current_drug=df_data$currentDrug
    feature_df$Feature_Value=sample_drugs[grep(current_drug,rownames(sample_drugs))[1],]
    feature_df$Feature_Quantile=sapply(1:ncol(sample_drugs),function(i){Fn=ecdf(sample_drugs[,i]);Fn(feature_df$Feature_Value[i])})
    #upate # tissues
    feature_df$Feature_Value[feature_df$Feature=="nTargets"]=sum(feature_df$Feature_Value[feature_df$Feature_Type=="Tissue-Specific Expression of Target"]>4)
    #update Lipinski Ro5
    if((feature_df$Feature_Value[feature_df$Feature=="MolecularWeight"]<500)&(feature_df$Feature_Value[feature_df$Feature=="XLogP"]<5)&(feature_df$Feature_Value[feature_df$Feature=="HydrogenBondDonorCount"]<5)&(feature_df$Feature_Value[feature_df$Feature=="HydrogenBondAcceptorCount"]<10)){
      feature_df$Feature_Value[feature_df$Feature=="Ro5"]=1
      feature_df$Feature_Quantile[feature_df$Feature=="Ro5"]=0   
    }else{
      feature_df$Feature_Value[feature_df$Feature=="Ro5"]=0
      feature_df$Feature_Quantile[feature_df$Feature=="Ro5"]=1   
    }
    #update Veber
    if((feature_df$Feature_Value[feature_df$Feature=="RotatableBondCount"]<12)&(feature_df$Feature_Value[feature_df$Feature=="PolarSurfaceArea"]<140)){
      feature_df$Feature_Value[feature_df$Feature=="Veber"]=1
      feature_df$Feature_Quantile[feature_df$Feature=="Veber"]=0    
    }else{
      feature_df$Feature_Value[feature_df$Feature=="Veber"]=0
      feature_df$Feature_Quantile[feature_df$Feature=="Veber"]=1 
    }
    #Update Ghose  - note: missing # atoms
    if((feature_df$Feature_Value[feature_df$Feature=="PolarSurfaceArea"]<140)&(findInterval(feature_df$Feature_Value[feature_df$Feature=="XLogP"], c(-0.4,5.6) )==1)&(findInterval(feature_df$Feature_Value[feature_df$Feature=="MolecularWeight"], c(160,480) ) ==1)&(findInterval(feature_df$Feature_Value[feature_df$Feature=="Refractivity"], c(40,130) ) ==1)&(findInterval(feature_df$Feature_Value[feature_df$Feature=="Refractivity"], c(40,130) ) ==1)){
      feature_df$Feature_Value[feature_df$Feature=="Ghose"]=1
      feature_df$Feature_Quantile[feature_df$Feature=="Ghose"]=0  
    }else{
      feature_df$Feature_Value[feature_df$Feature=="Ghose"]=0
      feature_df$Feature_Quantile[feature_df$Feature=="Ghose"]=1    
    }
    feature_df$Feature_Quantile[feature_df$Feature=="wQED"]=1-feature_df$Feature_Quantile[feature_df$Feature=="wQED"]
    return(feature_df)
  }
  
  #################
  # Change Value  #
  # Sub-Functions #
  #################
  changeValue <- function(this_var,newval){
    feature_df=df_data$features
    use_correlation_flag=df_data$corFlag
    
    Fn=ecdf(sample_drugs[,which(colnames(sample_drugs)==this_var)]);
    feature_df$Feature_Value[feature_df$Feature==this_var]=newval;
    feature_df$Feature_Quantile[feature_df$Feature==this_var]=Fn(newval);
    
    if(use_correlation_flag){
      f=which(colnames(sample_drugs)==this_var)
      q=Fn(newval)
      this_quantile_val=quantile(sample_drugs[,f],q,na.rm=T)
      low_prct=q-0.1
      hi_prct=q+0.1
      if(low_prct<=0){low_prct=0.01}
      if(hi_prct>=1){hi_prct=0.99}
      
      low_q=quantile(sample_drugs[,f],low_prct,na.rm=T)
      hi_q=quantile(sample_drugs[,f],hi_prct,na.rm=T)
      sub_sample_drugs=sample_drugs[sample_drugs[,f]>=low_q&sample_drugs[,f]<=hi_q,]
      
      cor_inds=setdiff(which(abs(feature_cor[f,])>0.5),f)    
      for(ind in cor_inds){
        replacement_val=median(sub_sample_drugs[,ind],na.rm=T)
        replacement_Fn=ecdf(sample_drugs[,ind])
        feature_df$Feature_Value[feature_df$Feature==colnames(sample_drugs)[ind]]=replacement_val
        feature_df$Feature_Quantile[feature_df$Feature==colnames(sample_drugs)[ind]]=replacement_Fn(replacement_val); 
      } 
    }
    #upate # tissues
    feature_df$Feature_Value[feature_df$Feature=="nTargets"]=sum(feature_df$Feature_Value[feature_df$Feature_Type=="Tissue-Specific Expression of Target"]>4)
    #update Lipinski Ro5
    if((feature_df$Feature_Value[feature_df$Feature=="MolecularWeight"]<500)&(feature_df$Feature_Value[feature_df$Feature=="XLogP"]<5)&(feature_df$Feature_Value[feature_df$Feature=="HydrogenBondDonorCount"]<5)&(feature_df$Feature_Value[feature_df$Feature=="HydrogenBondAcceptorCount"]<10)){
      feature_df$Feature_Value[feature_df$Feature=="Ro5"]=0
      feature_df$Feature_Quantile[feature_df$Feature=="Ro5"]=0      
    }else{
      feature_df$Feature_Value[feature_df$Feature=="Ro5"]=1
      feature_df$Feature_Quantile[feature_df$Feature=="Ro5"]=1
    }
    #update Veber
    if((feature_df$Feature_Value[feature_df$Feature=="RotatableBondCount"]<12)&(feature_df$Feature_Value[feature_df$Feature=="PolarSurfaceArea"]<140)){
      feature_df$Feature_Value[feature_df$Feature=="Veber"]=0
      feature_df$Feature_Quantile[feature_df$Feature=="Veber"]=0      
    }else{
      feature_df$Feature_Value[feature_df$Feature=="Veber"]=1
      feature_df$Feature_Quantile[feature_df$Feature=="Veber"]=1      
    }
    #Update Ghose - note: missing # atoms
    if((feature_df$Feature_Value[feature_df$Feature=="PolarSurfaceArea"]<140)&(findInterval(feature_df$Feature_Value[feature_df$Feature=="XLogP"], c(-0.4,5.6) )==1)&(findInterval(feature_df$Feature_Value[feature_df$Feature=="MolecularWeight"], c(160,480) ) ==1)&(findInterval(feature_df$Feature_Value[feature_df$Feature=="Refractivity"], c(40,130) ) ==1)&(findInterval(feature_df$Feature_Value[feature_df$Feature=="Refractivity"], c(40,130) ) ==1)){
      feature_df$Feature_Value[feature_df$Feature=="Ghose"]=0
      feature_df$Feature_Quantile[feature_df$Feature=="Ghose"]=0      
    }else{
      feature_df$Feature_Value[feature_df$Feature=="Ghose"]=1
      feature_df$Feature_Quantile[feature_df$Feature=="Ghose"]=1      
    }
    feature_df$Feature_Quantile[feature_df$Feature=="wQED"]=1-feature_df$Feature_Quantile[feature_df$Feature=="wQED"]
    
    feature_df
  }
  
  updateTargetFeature <- function(){
    feature_df=df_data$features
    feature_target=feature_df[!(feature_df$Feature_Type %in% c("Chemical Features","Drug-Likeness Rules")),]
  }
  
  updateStructureFeature <- function(){
    feature_df=df_data$features
    feature_structure=feature_df[feature_df$Feature_Type %in% c("Chemical Features","Drug-Likeness Rules"),]
    if(feature_df$Feature_Value[feature_df$Feature=="Ro5"]==1){feature_structure$Feature_Value[feature_structure$Feature=="Ro5"]="Fail"}else{feature_structure$Feature_Value[feature_structure$Feature=="Ro5"]="Pass"}
    if(feature_df$Feature_Value[feature_df$Feature=="Veber"]==1){feature_structure$Feature_Value[feature_structure$Feature=="Veber"]="Fail"}else{feature_structure$Feature_Value[feature_structure$Feature=="Veber"]="Pass"}
    if(feature_df$Feature_Value[feature_df$Feature=="Ghose"]==1){feature_structure$Feature_Value[feature_structure$Feature=="Ghose"]="Fail"}else{feature_structure$Feature_Value[feature_structure$Feature=="Ghose"]="Pass"}
    feature_structure=feature_structure[order(feature_structure$Feature_Type),]

    feature_structure
  }

  updatePred <- function(){
    feature_df=df_data$features
    prediction_df=df_data$predictions
    this_pred=median(sapply(rf.list,function(this_model) predict(this_model,get_PC_value(feature_df),"prob")[1]  ))
    this_pred_target=median(sapply(target.rf.list,function(this_model) predict(this_model,get_target_value(feature_df),"prob")[1]  ))
    this_pred_structure=median(sapply(structure.rf.list,function(this_model) predict(this_model,get_structure_value(feature_df),"prob")[1]  ))

    prediction_df$Prob[prediction_df$Model=="Overall"]=this_pred
    prediction_df$Prob[prediction_df$Model=="Structure"]=this_pred_structure
    prediction_df$Prob[prediction_df$Model=="Target"]=this_pred_target

    prediction_df$Log2_Odds[prediction_df$Model=="Overall"]=log2(this_pred/(1-this_pred))
    prediction_df$Log2_Odds[prediction_df$Model=="Structure"]=log2(this_pred_structure/(1-this_pred_structure))
    prediction_df$Log2_Odds[prediction_df$Model=="Target"]=log2(this_pred_target/(1-this_pred_target))

    if(this_pred>0.5){
      prediction_df$Prediction[prediction_df$Model=="Overall"]="Toxic"
    }else{
      prediction_df$Prediction[prediction_df$Model=="Overall"]="Safe"      
    }
    
    if(this_pred_structure>0.5){
      prediction_df$Prediction[prediction_df$Model=="Structure"]="Toxic"
    }else{
      prediction_df$Prediction[prediction_df$Model=="Structure"]="Safe"      
    }
    
    
    if(this_pred_target>0.5){
      prediction_df$Prediction[prediction_df$Model=="Target"]="Toxic"
    }else{
      prediction_df$Prediction[prediction_df$Model=="Target"]="Safe"      
    }
    
    prediction_df
    
  }

  ################
  # Plots OUTPUT #
  ################

  output$featurePlot <- renderPlot({
    feature_df=df_data$features
    plot_vars=setdiff(feature_df$Feature,c("nTargets"))
    ggplot(feature_df[feature_df$Feature %in% plot_vars,],aes(x=factor(Feature),y=Feature_Quantile,fill=Feature_Type))+geom_bar(stat = "identity")+scale_fill_manual(values=c("darkblue","forestgreen","darkorange1","red"))+theme_bw()+ theme(axis.text.x = element_text(size=12,angle = 60, vjust = 1, hjust=1),legend.position = "none")+xlab("")+ylab("Quantile Value")
  })

  output$predPlot <- renderPlot({
    prediction_df=df_data$predictions
    feature_df=df_data$features
    ggplot(output_scores_df2,aes(x=(score),fill=class))+geom_density(alpha=0.5,colour="black")+scale_fill_manual(values=c("forestgreen","red"))+xlab("PrOCTOR Score (log2 Odds of Toxicity)")+theme_bw()+theme(legend.position = "none")+ geom_segment(aes(x = prediction_df$Log2_Odds[prediction_df$Model=="Overall"], y = 0.3, xend = prediction_df$Log2_Odds[prediction_df$Model=="Overall"], yend = 0), arrow = arrow(length = unit(0.5, "cm")))+ggtitle(paste0("Overall Predicted ",prediction_df$Prediction[prediction_df$Model=="Overall"],"\nlog2 odds = ",round(abs(prediction_df$Log2_Odds[prediction_df$Model=="Overall"]),5)))    
  })

  output$predictionPanel<- renderPlot({
    prediction_df=df_data$predictions
    feature_df=df_data$features
    p1=ggplot(output_scores_df2,aes(x=(score),fill=class))+geom_density(alpha=0.5,colour="black")+scale_fill_manual(values=c("forestgreen","red"))+xlab("PrOCTOR Score (log2 Odds of Toxicity)")+theme_bw()+theme(legend.position = "none")+ geom_segment(aes(x = prediction_df$Log2_Odds[prediction_df$Model=="Overall"], y = 0.3, xend = prediction_df$Log2_Odds[prediction_df$Model=="Overall"], yend = 0), arrow = arrow(length = unit(0.5, "cm")))+ggtitle(paste0("Overall Predicted ",prediction_df$Prediction[prediction_df$Model=="Overall"],"\nlog2 odds = ",round(abs(prediction_df$Log2_Odds[prediction_df$Model=="Overall"]),5)))    
    p2=ggplot(output_scores_df2,aes(x=(score),fill=class))+geom_density(alpha=0.5,colour="black")+scale_fill_manual(values=c("forestgreen","red"))+xlab("PrOCTOR Score (log2 Odds of Toxicity)")+theme_bw()+theme(legend.position = c(0.8,0.875))+ geom_segment(aes(x = prediction_df$Log2_Odds[prediction_df$Model=="Target"], y = 0.3, xend = prediction_df$Log2_Odds[prediction_df$Model=="Target"], yend = 0), arrow = arrow(length = unit(0.5, "cm")))+ggtitle(paste0("Target Predicted ",prediction_df$Prediction[prediction_df$Model=="Target"],"\nlog2 odds = ",round(abs(prediction_df$Log2_Odds[prediction_df$Model=="Target"]),5)))
    p3=ggplot(output_scores_df2,aes(x=(score),fill=class))+geom_density(alpha=0.5,colour="black")+scale_fill_manual(values=c("forestgreen","red"))+xlab("PrOCTOR Score (log2 Odds of Toxicity)")+theme_bw()+theme(legend.position ="none")+ geom_segment(aes(x = prediction_df$Log2_Odds[prediction_df$Model=="Structure"], y = 0.3, xend = prediction_df$Log2_Odds[prediction_df$Model=="Structure"], yend = 0), arrow = arrow(length = unit(0.5, "cm")))+ggtitle(paste0("Structure Predicted ",prediction_df$Prediction[prediction_df$Model=="Structure"],"\nlog2 odds = ",round(abs(prediction_df$Log2_Odds[prediction_df$Model=="Structure"]),5)))    
    grid.arrange(p1,p2,p3, ncol=3)
  })
  
  ###############
  # Text OUTPUT #
  ###############

  # Feature Values 
  # Type = Table    
  output$featureValues <- renderDataTable({
    try({feature_visual=df_data$featureTable});
    feature_visual
  }, options = list(pageLength = 10))

  # Structure Feature Values 
  # Type = Table    
  output$structureFeatureValues <- renderDataTable({
    try({feature_structure=df_data$structureFeatures});
    feature_structure
  }, options = list(pageLength = 14))
  
  # Feature Values 
  # Type = Table    
  output$targetFeatureValues <- renderDataTable({
    try({feature_target=df_data$targetFeatures});
    feature_target
  }, options = list(pageLength = 10))

  # Output Predictions Values
  # Type = Table
  output$predictionValues <- renderTable({
    prediction_df=df_data$predictions
    prediction_df
  })

  # Output Predictions Values
  # Type = Text
  output$predictionInfo <- renderText({
    prediction_df=df_data$predictions
    paste(sapply(1:3,function(i) paste0(prediction_df$Model[i]," Model\n   Predicted ",prediction_df$Prediction[i],"\n   log2 odds = ",round(prediction_df$Log2_Odds[i],4))),collapse="\n\n")
  })

  output$overallPredictionInfo <- renderText({
    prediction_df=df_data$predictions
    paste("Predicted",prediction_df$Prediction[1],"\nProbability of Toxicity = ",round(prediction_df$Prob[1],4),"\nPrOCTOR Score = ",round(prediction_df$Log2_Odds[1],4))
  })

  output$structuralPredictionInfo <- renderText({
    prediction_df=df_data$predictions
    paste("Predicted",prediction_df$Prediction[2],"\nProbability of Toxicity = ",round(prediction_df$Prob[2],4),"\nPrOCTOR Score = ",round(prediction_df$Log2_Odds[2],4))
  })
  
  output$targetPredictionInfo <- renderText({
    prediction_df=df_data$predictions
    paste("Predicted",prediction_df$Prediction[3],"\nProbability of Toxicity = ",round(prediction_df$Prob[3],4),"\nPrOCTOR Score = ",round(prediction_df$Log2_Odds[3],4))
  })


})