#require(nlme)
#require(car)
require(ez)
require(nparLD)
options(contrasts=c("contr.sum","contr.poly"))

# set constants & paths
if(Sys.info()["user"] == "urs") rootdir = "/home/urs/sync/projects/wcst_eeg"

# read data, convert to long
df = read.table(file.path(rootdir,"data","anova_mixedmodel.txt"), header = TRUE)

anovafun <- function(df, dv, npermute = 1000, debug = TRUE){

df.err = df[, c("ID", "group", paste(dv,"wo", sep = "_"), paste(dv, "alc", sep = "_"))]
df.err = reshape( df.err, direction = "long",
                , varying = c(paste(dv,"wo", sep = "_"), paste(dv, "alc", sep = "_"))
                , v.names = dv
                , idvar = c("ID", "group"), timevar = "condition", times = c("wo","alc"))


# ANALYSIS 1 - run non-parametric mixed anova-like model from nparLD package
t.np = nparLD( formula(paste(dv, "group*condition", sep = "~"))
             , subject = "ID", data = df.err, description = FALSE)


# ANALYSIS 2 - run non-parametric tests on main and interaction effects separately

# test the main effect of group using Mann-Whitney U Test
df.err.ag.group = aggregate( formula(paste(dv, "group+ID", sep = "~"))
                           , data = df.err, mean)
twg = wilcox.test(formula(paste(dv, "~group")), data = df.err.ag.group)

# test the main effect of alcohol
df[paste(dv, "delta", sep = "_")] = df[paste(dv,"alc",sep="_")] - df[paste(dv,"wo",sep ="_")]
twc = wilcox.test(df[[paste(dv, "delta", sep = "_")]])

# test the interaction between group and alcohol
twi = wilcox.test(formula(paste(dv, "delta~group", sep = "_")), data = df)


# ANALYSIS 3 - do a permutation test
set.seed(42)

res = list( group       = rep(NA, npermute)
          , condition   = rep(NA, npermute)
          , interaction = rep(NA, npermute))

# calculate true values
eff.group       = diff(aggregate(formula(paste(dv,"group", sep="~")), df.err, mean)[[dv]])
eff.condition   = mean(df[[paste(dv, "delta", sep = "_")]])
eff.interaction = diff(aggregate(formula(paste(dv, "delta~group", sep = "_")), df
                                , mean)[[paste(dv, "delta", sep = "_")]])

for(i in c(1:npermute)){
  df.err$rand = sample(df.err[[dv]])
  res$group[i]    = diff(aggregate(rand~group, df.err, mean)$rand)
  dft =  reshape( df.err, direction = "wide", v.names = "rand", idvar = "ID"
                , timevar = "condition")
  dft$delta = dft$rand.alc - dft$rand.wo
  res$condition[i] = mean(dft$delta)
  res$interaction[i] = diff(aggregate(delta~group, dft, mean)$delta)
}

# ANALYSIS 4 - do rm Anova with ezanova
    
df.err$dv = df.err[[dv]]
df.err$dvatan = atan(df.err$dv)     
df.err$dvlog  = log(df.err$dv)   
df.err$dvlog[df.err$dv == 0] = 0  
df.err$group = as.factor(df.err$group)
eza = ezANOVA(data=df.err, dv= dvlog, wid=ID, within=.(condition), between= .(group), type = 3) 

# compare true values to permutation results. assumption: all means of permuted
# values are 0 (which is true, I checked it)
sg = sum(abs(res$group)       >= abs(eff.group))       / npermute
sc = sum(abs(res$condition)   >= abs(eff.condition))   / npermute
si = sum(abs(res$interaction) >= abs(eff.interaction)) / npermute

print(summary(t.np))

cat("\n\n\nWilcoxon Test for group:\n\n")    
print(twg)
cat("\n\n\nWilcoxon Test for condition (alc/wo):\n\n")    
print(twc)
cat("\n\n\nWilcoxon Test for group by condition interaction:\n\n")
print(twi)

cat("\n\nPermutation test results, probability of effect greater than or equal to:\n")
cat(paste("group:      ", sg, "\n"))    
cat(paste("condition:  ", sc, "\n"))    
cat(paste("interaction:", si, "\n"))    

cat("\n\nRepeated measures anova on log normalized data results (ezanova):\n\n")
print(eza)

}

sink("error_anova_results.txt")
anovafun(df, "memserr")
anovafun(df, "sserr")
sink(file = NULL)
