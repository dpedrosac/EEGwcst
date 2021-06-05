require(nlme)
require(car)
require(ez)

# set constants & paths
rootdir = D:\skripte\lambda
if(Sys.info()["user"] == "urs") rootdir = "/home/urs/sync/projects/wcst"

df = read.table(file.path(rootdir,"data","anova_rtimes.txt")


