# An Interactive Tool for Interpretation and Testing of PrOCTOR Features and Predictions
#### PrOCTOR = Predicting Odds of Clinical Trial Outcomes using a Random Forest Classifier

- Last Updated: 9/16/2016
- Correspondence to:  Katie Gayvert, kmg257 [at] cornell [dot] edu

##################
### Requirements #
##################
- R - tested on version  **3.2.2** (2015-08-14) -- **"Fire Safety"**
- R dependecies (Feature Interpretation Tool): shiny, shinyjs, data.table, plyr, htmlwidgets, ggplot2, randomForest, grid, gridExtra, Cairo
- Additional R dependecies (Prediction Tool): ChemmineR, ChemmineOB, rcdk, Rcpi


To Install R Dependencies:
```
install.packages(c("shiny", "shinyjs", "data.table", "plyr", "htmlwidgets", "ggplot2", "randomForest", "grid", "gridExtra", "Cairo","rcdk"))
source("https://bioconductor.org/biocLite.R")
biocLite("Rcpi")
biocLite("ChemmineR")
biocLite("ChemmineOB")

```

## PrOCTOR - Prediction
##########################
### To Run (in R) #
##########################
```
source("/path/to/PrOCTOR/R/PrOCTOR.R")
PrOCTOR(SMILE,target_list)
```
Example (Dexamethasone):
```
PrOCTOR(SMILE="[H][C@@]12C[C@@H](C)[C@](O)(C(=O)CO)[C@@]1(C)C[C@H](O)[C@@]1(F)[C@@]2([H])CCC2=CC(=O)C=C[C@]12C",
        targets=c("NR3C1","NR0B1","ANXA1","NOS2"))
```

##########################
### To Run (in R Studio) #
##########################
```
library(shiny)
load("/path/to/PrOCTOR/shiny/PrOCTOR_prediction/initial_values.RData")
runApp("/path/to/PrOCTOR/shiny/PrOCTOR_prediction")
```


## PrOCTOR - Feature Interpretation
##########################
### To Run (in R Studio) #
##########################
```
library(shiny)
load("/path/to/PrOCTOR/shiny/PrOCTOR_interpretation/initial_values.RData")
runApp("/path/to/PrOCTOR/shiny/PrOCTOR_interpretation")
```

#############
### Visuals #
#############
- **Feature Quantile Plot** - quantile values for each feature
- **Predictions** - compares each model's current prediction to distributions of failed and approved drugs in training set
- **Structure Feature Values** - Current structural feature values 
- **Target Feature Values** - Current target-based feature values

#############
### Options #
#############
- Pre-load features existing drug
- Change individual features using manually entered values *
- Change individual features by clicking new value on barplot *

\* Change other correlated feature values along with selected feature (default set off)

#######################
### Tasks In-Progress #
#######################
- [x] Create GitHub repository
- [x] Make functional interface where all changes update
- [x] Allow user inputted structures and targets
- [ ] Finish commenting code
- [ ] Finish model fine-tuning
