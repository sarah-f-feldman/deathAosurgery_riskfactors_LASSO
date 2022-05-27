##################################################
#   LASSO analysis                               #
# TNT as a marker for death after Aorta surgery  #
##################################################


#---------------------------------
#Chargement objets et packages
#---------------------------------

library(glmnet)
library(boot)
library(MASS)

final <- read.csv2("data/aorte.csv")


#===============================
#evenements precoces
#===============================

#-----------------
#DATA MANAGEMENT
#-----------------
prec <- final
prec$CJPtotal <- prec$CJPhosp
table(prec$CJPtotal)
sub <- prec[, c("tnt", "creatpre","age","imc", "met", "sex", "asa", "lee", "hta", "coronaroP","icc","aomi","diab", "acfa", "avc", "CJPtotal")]
#Je supprime les NA
work <- na.omit(sub)

work$sex <- as.numeric(work$sex)-1
work$CJPtotal <- as.factor(work$CJPtotal)
var <- colnames(work)[colnames(work) != "CJPtotal"]

#-----------------
#Donnees manquantes
#-----------------
#nombre de patients avec au moins une donnee manquante : 15 patients soit <5%
round(prop.table(table(apply(sub,1, function(x)sum(is.na(x)))>0)),2)
#nombre de NA par colonnes
apply(sub,2, function(x)sum(is.na(x)))


#----------------
#virer les extremes
#----------------

formula.vec <- paste0("CJPtotal ~ ", paste(var, collapse="+"))
mod.lg <- glm(as.formula(formula.vec), family = "binomial", data = work)
a.obj <- stepAIC(mod.lg)
#plot(glm(formula = CJPtotal ~ tnt + creatpre + age + sex + hta + aomi + acfa + avc, family = "binomial", data = work))
plot(mod.lg)

#a supprimer pour evenements precoces
#work <- work[! rownames(work) %in% c(110, 168, 322), ]
work <- work[! rownames(work) %in% c(110, 168, 213, 10, 20, 322), ]

#CJPtotal ~ tnt + age + hta + aomi + acfa + avc

#2e tour de validation
mod.lg <- glm(as.formula(formula.vec), family = "binomial", data = work)
a.obj <- stepAIC(mod.lg)
plot( glm(as.formula(formula.vec), family = "binomial", data = work))

formula.vec <- paste0("CJPtotal ~ ", paste(names(coef(a.obj))[-1], collapse="+"))
g <- glm(as.formula(formula.vec), family = "binomial", data=work)
summary(g)
# ce n'est pas plus stable en les supprimant donc je garde :
# work <- work[! rownames(work) %in% c(50, 10, 20, 213), ]
# work[rownames(work) %in% 10, ]

#-----------------
#Lasso
#-----------------
x <- as.matrix(work[ , var])
y <- work$CJPtotal


fit <- glmnet(x, y, alpha = 1, family = "binomial") #1 is lasso (it's the default anyway)
plot(fit,label = TRUE)
set.seed(12345)
#set.seed(123)
cvfit = cv.glmnet(x, y, family = "binomial")
#cvfit = cv.glmnet(x, y, family = "binomial", type.measure = "auc")
plot(cvfit)
#recuperer le lambda
lambda1 <- cvfit$lambda.1se
#lambda1 <- cvfit$lambda.min
a <- coef(cvfit, s = lambda1)
a

# Pour recuperer le tableau des beta obtenus par lasso
vec <- rep(NA, a@Dim[1])
vec[a@i] <- a@x[-1]
names.a <- a@Dimnames[[1]]
var.df <- data.frame(variable = as.character(names.a[(a@i)+1])[-1], beta = as.numeric(vec[(a@i)+1])[-length(a@i)], stringsAsFactors = FALSE)
#noms des variables selectionnees par lasso
boot.var <- var.df[!is.na(var.df$beta), "variable"]


#------------------------------------
# intervalle de confiance des coefficients selectionnes par bootstrap 

#NB : on garde le lambda1 selectionne a l'etape d'avant

data <- work

set.seed(12345)
boot.lasso <- function(data, indices){
  data <- data[indices, ]
  data <- na.omit(data)
  x <- as.matrix(data[ , boot.var])
  y <- data$CJPtotal
  cvfit <- cv.glmnet(x, y, family = "binomial")
  a <- coef(cvfit, s = lambda1)
  res <- as.numeric(a)[-1]
  return(res)
}
.bootres <- boot(data = work, statistic = boot.lasso, R = 1000)

#intervalle de confiance
.n <- length (.bootres$t0)
.list.ci <- lapply(1:n, function(x) boot.ci(res,index=x,type="bca"))
.res <- data.frame (t (sapply (.list.ci, function (x) x[["bca"]][4:5]))) #selectionne les valeur de IC
rownames (.res) <- names (.bootres$t0)
colnames (.res) <- c ("CI_L", "CI_U")
.res$est <- as.numeric (.bootres$t0)
#.res$n <- nrow(work)
.ans <- round (.res, 4) #fait un arrondi sur chaque valeur
.res$variable <- boot.var
.res <- .res[ , c("variable", "est", "CI_L", "CI_U")]
.ans <- data.frame (variable=.res$variable, beta_CI=paste0 (.ans$est, " [", .ans$CI_L, "-", .ans$CI_U, "]")) #met en forme les valeurs
.ans 




#===============================
#evenements tardifs
#===============================

#---------------------------------
#DATA MANAGEMENT
#---------------------------------
tardif <- final[final$CJPhosp != 1, ]
table(tardif$CJPtotal)
sub <- tardif[, c("tnt", "creatpre","age","imc", "met", "sex", "asa", "lee", "hta", "coronaroP","icc","aomi","diab", "acfa", "avc", "CJPtotal")]
#Je supprime les NA
work <- na.omit(sub)
work$sex <- as.numeric(work$sex)-1
work$CJPtotal <- as.factor(work$CJPtotal)
var <- colnames(work)[colnames(work) != "CJPtotal"]
#var <- c("tnt", "creatpre","age","imc", "met", "sex", "asa", "lee", "hta", "coronaroP","icc","aomi","diab", "acfa", "avc")

#-----------------
#Donnees manquantes
#-----------------
#nombre de patients avec au moins une donnee manquante : 15 patients soit <5%
round(prop.table(table(apply(sub,1, function(x)sum(is.na(x)))>0)),2)
#nombre de NA par colonnes
apply(sub,2, function(x)sum(is.na(x)))


#----------------
#virer les extremes
#----------------
formula.vec <- paste0("CJPtotal ~ ", paste(var, collapse="+"))
mod.lg <- glm(as.formula(formula.vec), family = "binomial", data = work)
a.obj <- stepAIC(mod.lg)
#faut il utiliser toutes les variables(l.1) ou seulement celles selectionnees par stepAIC (l.2)?
plot(mod.lg)
plot(glm(formula = CJPtotal ~ tnt + creatpre + age + imc, family = "binomial", data = work))
#meme resultat
#a supprimer pour evenements tardifs
work <- work[! rownames(work) %in% c(74, 110), ] 

#2e tour de validation
mod.lg <- glm(as.formula(formula.vec), family = "binomial", data = work)
a.obj <- stepAIC(mod.lg) #je regarde si les variables selectionnees sont stable (cad toujour selectionnees malgre retrait de sujets)
plot(mod.lg)
#ok, plus d'individus a supprimer

#-----------------
#Lasso
#-----------------
x <- as.matrix(work[ , var])
y <- work$CJPtotal


fit <- glmnet(x, y, alpha = 1, family = "binomial") #1 is lasso (it's the default anyway)
plot(fit,label = TRUE)
set.seed(12345)
#set.seed(123)
cvfit = cv.glmnet(x, y, family = "binomial")
#cvfit = cv.glmnet(x, y, family = "binomial", type.measure = "auc")
plot(cvfit)
#recuperer le lambda
lambda1 <- cvfit$lambda.1se
#lambda1 <- cvfit$lambda.min
a <- coef(cvfit, s = lambda1)
a

# Pour recuperer le tableau des beta obtenus par lasso
vec <- rep(NA, a@Dim[1])
vec[a@i] <- a@x[-1]
names.a <- a@Dimnames[[1]]
var.df <- data.frame(variable = as.character(names.a[(a@i)+1])[-1], beta = as.numeric(vec[(a@i)+1])[-length(a@i)], stringsAsFactors = FALSE)
#noms des variables selectionnees par lasso
boot.var <- var.df[!is.na(var.df$beta), "variable"]


#------------------------------------
# intervalle de confiance des coefficients selectionnes par bootstrap 

#NB : on garde le lambda1 selectionne a l'etape d'avant

data <- work

set.seed(12345)
boot.lasso <- function(data, indices){
  data <- data[indices, ]
  data <- na.omit(data)
  x <- as.matrix(data[ , boot.var])
  y <- data$CJPtotal
  cvfit <- cv.glmnet(x, y, family = "binomial")
  a <- coef(cvfit, s = lambda1)
  res <- as.numeric(a)[-1]
  return(res)
}
.bootres <- boot(data = work, statistic = boot.lasso, R = 1000)

#intervalle de confiance
.n <- length (.bootres$t0)
.list.ci <- lapply(1:n, function(x) boot.ci(res,index=x,type="bca"))
.res <- data.frame (t (sapply (.list.ci, function (x) x[["bca"]][4:5]))) #selectionne les valeur de IC
rownames (.res) <- names (.bootres$t0)
colnames (.res) <- c ("CI_L", "CI_U")
.res$est <- as.numeric (.bootres$t0)
#.res$n <- nrow(work)
.ans <- round (.res, 4) #fait un arrondi sur chaque valeur
.ans$variable <- boot.var
.ans <- .ans[ , c("variable", "est", "CI_L", "CI_U")]
.ans <- data.frame (variable=.ans$variable, beta_CI=paste0 (.ans$est, " [", .ans$CI_L, "-", .ans$CI_U, "]")) #met en forme les valeurs
.ans 
write.table(.ans, file="clipboard", sep= "\t", row.names=F)
#OR
.ans <- .res
.ans[ ,c("est", "CI_L", "CI_U")] <- round(exp(.ans[ ,c("est", "CI_L", "CI_U")]), 3)
.ans$variable <- boot.var
.ans <- data.frame (variable=.ans$variable, beta_CI=paste0 (.ans$est, " [", .ans$CI_L, "-", .ans$CI_U, "]")) #met en forme les valeurs
.ans 
write.table(.ans, file="clipboard", sep= "\t", row.names=F)

