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
df <- read.csv("ERPampl_reactionTime_Pz_p3B.csv")
df$group[df$group==1] <- "ET patients" 
df$group[df$group==2] <- "CTRL subjects" 
df$group <- as.factor(df$group)

df$condition[df$condition==1] <- "no alcohol" 
df$condition[df$condition==2] <- "with alcohol"
df$condition <- as.factor(df$condition)

dfERP <- df %>% tidyr::pivot_longer(cols = starts_with("ERP"), names_to = "repetition1", values_to = "amplitude", values_drop_na=TRUE) 
# %>% pivot_longer(cols = starts_with("rt"), names_to="repetition2", values_to="rt")

# Run analyses to determine the prerequisites for an ANOVA
dfERP %>%
  group_by(group, condition, repetition1) %>%
  get_summary_stats(rt, type = "mean_sd")

dfERP_aggr <- with(dfERP, aggregate(amplitude, list(group=group, repetition1=repetition1, condition=condition), mean))
dfERP_aggr$se1 <- with(dfERP , aggregate(amplitude, list(group=group, condition=condition), 
              function(x) sd(x)/sqrt(10)))[,3]
dfERP_aggr <- dfERP_aggr %>% rename(x1=x)

dfRT <- df %>% tidyr::pivot_longer(cols = starts_with("rt"), names_to = "repetition2", values_to = "rtimes", values_drop_na=TRUE) 			  
dfRT_aggr <- with(dfRT, aggregate(rtimes, list(group=group, repetition2=repetition2, condition=condition), mean))
dfRT_aggr$se2 <- with(dfRT , aggregate(rtimes, list(group=group, condition=condition), 
              function(x) sd(x)/sqrt(10)))[,3]
dfRT_aggr <- dfRT_aggr %>% rename(x2=x)

dfTotal_aggr <- cbind(dfERP_aggr, x2=dfRT_aggr$x2, se2=dfRT_aggr$se2)
								  
								  # bxp <- ggboxplot(
# Run repeated measures ANOVA with unbalanced groups
# newModel <- ezANOVA(data=df, dv=meanERP, wid=ID, within=condition, between=group, detailed=TRUE, type=3)

pd <- position_dodge(0.1)
ggplot(dfTotal_aggr, aes(x=repetition1, y=x1, group=group)) + 
    scale_y_continuous(limits = c(-.8, 0.8)) +
    geom_line(aes(linetype=group), size=.6, position = pd) + 
    geom_errorbar(aes(ymax=x1+se1, ymin=x1-se1), width=.1, position = pd) +
	geom_point(size = 4, position = pd) + 
    geom_point(size = 3, color="white", position = pd) +
	facet_grid(condition ~ .)

dev.new()									   
pd <- position_dodge(0.1)
ggplot(dfTotal_aggr, aes(x=repetition1, y=x2, group=group)) + 
    scale_y_continuous(limits = c(0, 2)) +
    geom_line(aes(linetype=group), size=.6, position = pd) + 
    geom_errorbar(aes(ymax=x2+se2, ymin=x2-se2), width=.1, position = pd) +
	geom_point(size = 4, position = pd) + 
    geom_point(size = 3, color="white", position = pd) +
	facet_grid(condition ~ .)

								  +
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