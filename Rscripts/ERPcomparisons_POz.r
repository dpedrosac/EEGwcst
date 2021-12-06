# This is code to run analyses for the changes between groups and conditions;
# Code developed by David Pedrosa

## First specify the packages of interest
packages = c("tidyverse", "ez", "rstatix", "ggpubr") #

## Now load or install&load all
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
wdir = "/media/storage/skripte/lambda/results"
}
setwd(wdir)

# Read data and tidy up inputs
df <- read.csv("ERP_resuts_per_group_Pz_p600-800_shift.csv")
df$group[df$group==1] <- "ET patients" 
df$group[df$group==2] <- "CTRL subjects" 
df$group <- as.factor(df$group)

df$condition[df$condition==1] <- "wo alcohol" 
df$condition[df$condition==2] <- "with alcohol"
df$condition <- as.factor(df$condition)

# Run analyses to determine the prerequisites for an ANOVA
df %>%
  group_by(group, condition) %>%
  get_summary_stats(meanERP, type = "mean_sd")

bxp <- ggboxplot(
  df, x = "group", y = "meanERP",
  color = "condition", palette = "jco"
  )
bxp

# Identify outliers
df %>%
  group_by(group, condition) %>%
  identify_outliers(meanERP)
  
df %>%
  group_by(group, condition) %>%
  shapiro_test(meanERP)
shapiro:
_plot <- ggqqplot(df, "meanERP", ggtheme = theme_bw()) +
  facet_grid(condition ~ group, labeller = "label_both")


# Run repeated measures ANOVA with unbalanced groups
newModel <- ezANOVA(data=df, dv=meanERP, wid=ID, within=condition, between=group, detailed=TRUE, type=3)

# Pairwise comparisons between treatment groups
pwc <- df %>%
  group_by(group) %>%
  pairwise_t_test(meanERP ~ condition,
    p.adjust.method = "bonferroni"
    )
pwc