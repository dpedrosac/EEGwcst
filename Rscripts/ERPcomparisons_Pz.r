# This is code to run analyses for the changes between groups and conditions;
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
wdir = "/media/storage/skripte/lambda/results"
}
setwd(wdir)

# Read data and tidy up inputs
df <- read.csv("ERP_resuts_per_group_Pz_p600-800_shift.csv")
df$group[df$group==1] <- "ET patients" 
df$group[df$group==2] <- "CTRL subjects" 
df$group <- as.factor(df$group)

df$condition[df$condition==1] <- "no alcohol" 
df$condition[df$condition==2] <- "with alcohol"
df$condition <- as.factor(df$condition)

# Run analyses to determine the prerequisites for an ANOVA
df %>%
  group_by(group, condition) %>%
  get_summary_stats(meanERP, type = "mean_sd")

# bxp <- ggboxplot(
#  df, x = "group", y = "meanERP",
#  color = "condition", palette = "jco"
#  )
# bxp

# Run repeated measures ANOVA with unbalanced groups
newModel <- ezANOVA(data=df, dv=meanERP, wid=ID, within=condition, between=group, detailed=TRUE, type=3)

df_aggr <- with(df, aggregate(meanERP, list(group=group, condition=condition), mean))
df_aggr$se <- with(df , aggregate(meanERP, list(group=group, condition=condition), 
              function(x) sd(x)/sqrt(10)))[,3]

pd <- position_dodge(0.1)
ggplot(df_aggr, aes(x=condition, y=x, group=group)) + 
    scale_y_continuous(limits = c(0, 4.5))+
    geom_line(aes(linetype=group), size=.6, position = pd) + 
    geom_point(size = 4, position = pd) + 
    geom_point(size = 3, color="white", position = pd) +
    geom_errorbar(aes(ymax=x+se, ymin=x-se), width=.1, position = pd) +
    guides(linetype = guide_legend("Group")) +
    labs(title = paste("ERP amplitude (Pz) depending on,",
                       "condition (alcohol intake) and group.",
                       "Error bars represent 95% Confidence Intervals",
                       sep = "\n"),
        x = "",
        y = "meanERP") +
    theme(
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
         panel.border = element_blank(),
        legend.key  = element_rect(fill = "white"),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        text = element_text(size = 20))


## Pairwise comparisons between treatment groups
#pwc <- df %>%
#  group_by(group) %>%
#  pairwise_t_test(meanERP ~ condition,
#    p.adjust.method = "bonferroni"
#    )
#pwc
