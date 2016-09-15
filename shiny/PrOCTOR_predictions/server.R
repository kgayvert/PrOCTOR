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
library(ChemmineR)
library(ChemmineOB)
library(rcdk)
library(Rcpi)
load("initial_values.RData")

shinyServer(function(input, output) {
  
  ######################
  # Reactive DataFrame #
  ######################    
  df_data=reactiveValues(features=initial_feature_df,predictions=initial_prediction_df,structureFeatures=initial_feature_structure,targetFeatures=initial_feature_target,corFlag=FALSE,currentDrug="Median Values",SMILES=NA,target=NA)
  
  ######################
  # Reactive Functions #
  ######################    
  observeEvent(input$goButton, {
    df_data$SMILES=input$smiles
    df_data$target=input$target
    
    df_data$features=userInputButton()
    df_data$predictions=updatePred() 
    df_data$structureFeatures=updateStructureFeature()
    df_data$targetFeatures=updateTargetFeature()
  })
  
  userInputButton <- function(){
    return(userInput(this_targets=input$target,this_SMILES=input$smiles))
  }

  
  #################
  # Change Value  #
  # Sub-Functions #
  #################
  userInput <- function(this_targets,this_SMILES){
    feature_df=df_data$features

    this_features=c(getAllFeatures(this_SMILES,this_targets),nTargets=length(this_targets))
    feature_df$Feature_Value=as.vector(sapply(as.vector(unlist(feature_df$Feature)),function(x) this_features[names(this_features)==x]))
    
    feature_df
  }
  
  
  updateTargetFeature <- function(){
    feature_df=df_data$features
    feature_target=feature_df[(feature_df$Feature_Type %in% c("Tissue-Specific Expression of Target","Other Target-Based Features")),]
  }
  
  updateStructureFeature <- function(){
    feature_df=df_data$features
    feature_structure=feature_df[feature_df$Feature_Type %in% c("Chemical Features","Drug-Likeness Rules"),]
    if(feature_df$Feature_Value[feature_df$Feature=="Ro5"]==1){feature_structure$Feature_Value[feature_structure$Feature=="Ro5"]="Pass"}else{feature_structure$Feature_Value[feature_structure$Feature=="Ro5"]="Fail"}
    if(feature_df$Feature_Value[feature_df$Feature=="Veber"]==1){feature_structure$Feature_Value[feature_structure$Feature=="Veber"]="Pass"}else{feature_structure$Feature_Value[feature_structure$Feature=="Veber"]="Fail"}
    if(feature_df$Feature_Value[feature_df$Feature=="Ghose"]==1){feature_structure$Feature_Value[feature_structure$Feature=="Ghose"]="Pass"}else{feature_structure$Feature_Value[feature_structure$Feature=="Ghose"]="Fail"}
    feature_structure=feature_structure[order(feature_structure$Feature_Type),]

    feature_structure
  }

  updatePred <- function(){
    feature_df=df_data$features
    prediction_df=df_data$predictions
    this_pred=median(sapply(rf.list,function(this_model) predict(this_model,get_all_value(feature_df),"prob")[1]  ))
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
    ggplot(output_scores_df2,aes(x=-1*(score),fill=class))+geom_density(alpha=0.5,colour="black")+scale_fill_manual(values=rev(c("skyblue","red")))+xlab("PrOCTOR Score (log2 Odds of Toxicity)")+theme_bw()+theme(legend.position = "none")+ geom_segment(aes(x = prediction_df$Log2_Odds[prediction_df$Model=="Overall"], y = 0.3, xend = prediction_df$Log2_Odds[prediction_df$Model=="Overall"], yend = 0), arrow = arrow(length = unit(0.5, "cm")))+ggtitle(paste0("Overall Predicted ",prediction_df$Prediction[prediction_df$Model=="Overall"],"\nlog2 odds = ",round(abs(prediction_df$Log2_Odds[prediction_df$Model=="Overall"]),5)))    
  })

  output$predictionPanel<- renderPlot({
    prediction_df=df_data$predictions
    feature_df=df_data$features
    if(!(is.na(df_data$target))){
      p1=ggplot(output_scores_df,aes(x=-1*(score),fill=class))+geom_density(alpha=0.5,colour="black")+scale_fill_manual(values=(c("skyblue3","red")))+xlab("PrOCTOR Score (log2 Odds of Toxicity)")+theme_bw()+theme(legend.position = "none")+ geom_segment(aes(x = -1*prediction_df$Log2_Odds[prediction_df$Model=="Overall"], y = 0.3, xend = -1*prediction_df$Log2_Odds[prediction_df$Model=="Overall"], yend = 0), arrow = arrow(length = unit(0.5, "cm")))+ggtitle(paste0("Overall Predicted ",prediction_df$Prediction[prediction_df$Model=="Overall"],"\nlog2 odds = ",round(abs(prediction_df$Log2_Odds[prediction_df$Model=="Overall"]),5)))    
      p2=ggplot(output_target_scores_df,aes(x=-1*(score),fill=class))+geom_density(alpha=0.5,colour="black")+scale_fill_manual(values=(c("skyblue3","red")))+xlab("PrOCTOR Score (log2 Odds of Toxicity)")+theme_bw()+theme(legend.position = c(0.1,0.875))+ geom_segment(aes(x = -1*prediction_df$Log2_Odds[prediction_df$Model=="Target"], y = 0.3, xend = -1*prediction_df$Log2_Odds[prediction_df$Model=="Target"], yend = 0), arrow = arrow(length = unit(0.5, "cm")))+ggtitle(paste0("Target Predicted ",prediction_df$Prediction[prediction_df$Model=="Target"],"\nlog2 odds = ",round(abs(prediction_df$Log2_Odds[prediction_df$Model=="Target"]),5)))
      p3=ggplot(output_struct_scores_df,aes(x=-1*(score),fill=class))+geom_density(alpha=0.5,colour="black")+scale_fill_manual(values=(c("skyblue3","red")))+xlab("PrOCTOR Score (log2 Odds of Toxicity)")+theme_bw()+theme(legend.position ="none")+ geom_segment(aes(x = -1*prediction_df$Log2_Odds[prediction_df$Model=="Structure"], y = 0.3, xend = -1*prediction_df$Log2_Odds[prediction_df$Model=="Structure"], yend = 0), arrow = arrow(length = unit(0.5, "cm")))+ggtitle(paste0("Structure Predicted ",prediction_df$Prediction[prediction_df$Model=="Structure"],"\nlog2 odds = ",round(abs(prediction_df$Log2_Odds[prediction_df$Model=="Structure"]),5)))          
    }else{
      p1=ggplot(output_scores_df,aes(x=-1*(score),fill=class))+geom_density(alpha=0.5,colour="black")+scale_fill_manual(values=(c("skyblue3","red")))+xlab("PrOCTOR Score (log2 Odds of Toxicity)")+theme_bw()+theme(legend.position = "none")+ggtitle(paste0("Overall log2 odds Distribution"))    
      p2=ggplot(output_target_scores_df,aes(x=-1*(score),fill=class))+geom_density(alpha=0.5,colour="black")+scale_fill_manual(values=(c("skyblue3","red")))+xlab("PrOCTOR Score (log2 Odds of Toxicity)")+theme_bw()+theme(legend.position = c(0.1,0.875))+ggtitle(paste0("Target log2 odds Distribution"))
      p3=ggplot(output_struct_scores_df,aes(x=-1*(score),fill=class))+geom_density(alpha=0.5,colour="black")+scale_fill_manual(values=(c("skyblue3","red")))+xlab("PrOCTOR Score (log2 Odds of Toxicity)")+theme_bw()+theme(legend.position ="none")+ ggtitle(paste0("Structure Predicted log2 odds Distribution"))    
    }
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
    feature_structure=feature_structure[,-4,with=FALSE]

    feature_structure
  }, options = list(pageLength = 14))
  
  # Feature Values 
  # Type = Table    
  output$targetFeatureValues <- renderDataTable({
    try({feature_target=df_data$targetFeatures});
    feature_target=feature_target[,-4,with=FALSE]
  }, options = list(pageLength = 10))

  # Output Predictions Values
  # Type = Table
  output$predictionValues <- renderTable({
    prediction_df=df_data$predictions
    prediction_df
  })

  # Output Predictions Values
  # Type = Text
  output$drugInfo <- renderText({
    this_SMILES=df_data$predictions
    paste0("SMILES = ",df_data$SMILES,"\nTarget(s) = ",paste(df_data$target,collapse=","))
  })
  
  
  output$predictionInfo <- renderText({
    prediction_df=df_data$predictions
    if(!is.na(df_data$SMILES)){
      paste(sapply(1:3,function(i) paste0(prediction_df$Model[i]," Model\n   Predicted ",prediction_df$Prediction[i],"\n   log2 odds = ",-1*round(prediction_df$Log2_Odds[i],4))),collapse="\n\n")
    }else{
      "N/A - Please input drug"
    }
  })

  output$overallPredictionInfo <- renderText({
    prediction_df=df_data$predictions
    if(!is.na(prediction_df$Log2_Odds[1])){
      paste("Predicted",prediction_df$Prediction[1],"\nProbability of Toxicity = ",round(prediction_df$Prob[1],4),"\nPrOCTOR Score = ",-1*round(prediction_df$Log2_Odds[1],4))
    }else{
      "N/A - Please input drug"
    }
  })

  output$structuralPredictionInfo <- renderText({
    prediction_df=df_data$predictions
    if(!is.na(df_data$SMILES)){
      paste("Predicted",prediction_df$Prediction[2],"\nProbability of Toxicity = ",round(prediction_df$Prob[2],4),"\nPrOCTOR Score = ",-1*round(prediction_df$Log2_Odds[2],4))
    }else{
      "N/A - Please input drug structure"
    }

  })
  
  output$targetPredictionInfo <- renderText({
    prediction_df=df_data$predictions
    if(!is.na(df_data$target)){
      paste("Predicted",prediction_df$Prediction[3],"\nProbability of Toxicity = ",round(prediction_df$Prob[3],4),"\nPrOCTOR Score = ",-1*round(prediction_df$Log2_Odds[3],4))
    }else{
      "N/A - Please input drug target"
    }
    
  })


})