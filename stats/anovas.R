#require(nlme)
#require(car)
#require(ez)
require(nparLD)

# set constants & paths
if(Sys.info()["user"] == "urs") rootdir = "/home/urs/sync/projects/wcst_eeg"

# read data, convert to long
df = read.table(file.path(rootdir,"data","anova_mixedmodel.txt"), header = TRUE)

df.memserr = subset(df, select = -c(subj, sserr_wo, sserr_alc))
df.memserr = reshape(df.memserr, direction = "long",
              ,varying = c("memserr_wo","memserr_alc"), v.names = "memserr"
              ,idvar = c("ID", "group"), timevar = "condition", times = c("wo","alc"))



# ANALYSIS 1 - run non-parametric mixed anova-like model from nparLD package
nparLD(memserr~group*condition, subject = "ID", data = df.memserr, description = FALSE)


# ANALYSIS 2 - run non-parametric tests on main and interaction effects separately

# test the main effect of group using Mann-Whitney U Test
df.memserr.ag.group = aggregate(memserr~group    +ID, data = df.memserr, mean)
wilcox.test(memserr~group    , data = df.memserr.ag.group)

# test the main effect of alcohol
df$memserr_delta = df$memserr_alc - df$memserr_wo
wilcox.test(df$memserr_delta)

# test the interaction between group and alcohol
wilcox.test(memserr_delta~group, data = df)

# ANALYSIS 3 - do a permutation test
set.seed(42)
npermute = 1000
res = list( group = rep(NA, npermute)
          , condition = rep(NA, npermute)
          , interaction = rep(NA, npermute))

# calculate true values
eff.group       = diff(aggregate(memserr~group, df.memserr, mean)$memserr)
eff.condition   = mean(df$memserr_delta)
eff.interaction = diff(aggregate(memserr_delta~group, df, mean)$memserr_delta)

for(i in c(1:npermute)){
  df.memserr$rand = sample(df.memserr$memserr)
  res$group[i]    = diff(aggregate(rand~group, df.memserr, mean)$rand)
  dft =  reshape( df.memserr, direction = "wide", v.names = "rand", idvar = "ID"
                , timevar = "condition")
  dft$delta = dft$rand.alc - dft$rand.wo
  res$condition[i] = mean(dft$delta)
  res$interaction[i] = diff(aggregate(delta~group, dft, mean)$delta)
}

sum(abs(res$group)       >= abs(eff.group))       / npermute
sum(abs(res$condition)   >= abs(eff.condition))   / npermute
sum(abs(res$interaction) >= abs(eff.interaction)) / npermute
