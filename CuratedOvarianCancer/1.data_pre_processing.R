source("functions.R")
library(curatedOvarianData)
library(limma) ## needed for matrix manipulation
library(glmnet)



data(GSE9891_eset)
data(GSE26712_eset)

GSE98_eset = t(as.matrix(GSE9891_eset))
GSE26_eset = t(as.matrix(GSE26712_eset))


GSE98_vars = apply(GSE98_eset, 2, var)
GSE26_vars = apply(GSE26_eset, 2, var)

id_sel98 = which(GSE9891_eset@phenoData$vital_status =='deceased')
id_sel26 = which(GSE26712_eset@phenoData$vital_status  =='deceased')
## LASSO variable selection

## target
X = GSE98_eset[id_sel98,]
y = log(GSE9891_eset@phenoData$days_to_death[id_sel98])

target_cv = cv.glmnet(X,y, alpha =1) 
plot(target_cv)
target_sel = coef(target_cv, s = tail(target_cv$lambda))
sum(target_sel[,1]!=0)

target_26_var = colnames(GSE26_eset)[which(colnames(GSE26_eset) %in% rownames(target_sel)[target_sel[,1] !=0])]
write.csv(target_26_var, file = 'target_26_var.csv')

## proxy variables
X = GSE26_eset[id_sel26,]
y = log(GSE26712_eset@phenoData$days_to_death[id_sel26])

proxy_cv = cv.glmnet(X,y, alpha =1) 
plot(proxy_cv)

proxy_sel = coef(proxy_cv, s = tail(proxy_cv$lambda))
sum(proxy_sel[,1]!=0)
proxy_98_var = unique(c(rownames(proxy_sel)[proxy_sel[,1]!=0],target_26_var), target_26_var)
write.csv(proxy_98_var, file = 'proxy_98_var.csv')