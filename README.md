# An Interactive Tool for Interpretation of PrOCTOR Predictions
###### PrOCTOR = Predicting Odds of Clinical Trial Outcomes using a Random Forest Classifier

- Last Updated: 1/30/2016
- Correspondence to:  Katie Gayvert, kmg257 [at] cornell [dot] edu

################
# Requirements #
################
- R - tested on version  3.2.2 (2015-08-14) -- "Fire Safety"
- R dependecies: shiny, shinyjs, data.table, plyr, htmlwidgets, ggplot2, randomForest, grid, gridExtra, Cairo

To Install:
```
install.packages(c("shiny", "shinyjs", "data.table", "plyr", "htmlwidgets", "ggplot2", "randomForest", "grid", "gridExtra", "Cairo"))
```
########################
# To Run (in R Studio) #
########################
```
library(shiny)
library(shinyjs)
load("/path/to/initial_values.RData")
runApp("/path/to/PrOCTOR")
```

###########
# Visuals #
###########
- Feature Quantile Plot - quantile values for each feature
- Predictions - where each model's current prediction compares to distributions of failed and approved drugs in training set
- Structure Feature Values - Current structural feature values 
- Target Feature Values - Current target-based feature values

###########
# Options #
###########
- Pre-load features existing drug
- Change individual features using manually entered values *
- Change individual features by clicking new value on barplot *
- * Change other correlated feature values along with selected feature (set to default)

