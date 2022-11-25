# This is code to run ANOVA for memory errors durcing the WCST for tremor and control subjects;
# Code developed by David Pedrosa

## First specify packages of interest
packages = c("tidyverse", "ez", "rstatix", "ggpubr") #

## Now load or install & load those packages
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)


## In case of multiple people working on one project, this helps to create an automatic script
username = Sys.info()["login"]
if (username == "dpedr") {
wdir = "D:/skripte/lambda/results"
} else if (username == "david") {
wdir = "/media/storage/skripte/lambda/data"
}
setwd(wdir)

# Load data into workspace
df <- read.csv2(file.path(wdir,"anova_data.txt"), sep="\t")
ModelErrorsWCST <- ezANOVA(data=df, dv=error, wid=ID, within=c(condition, error_type), between=group, detailed=TRUE, type=3)

