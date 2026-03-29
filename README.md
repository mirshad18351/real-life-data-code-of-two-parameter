# real-life-data-code-of-two-parameter
On Some Improved Two-Parameter Ridge Regression Estimators for the Heteroscedastic Linear Regression Models: Simulation and Application
#R-Code for Ridge Regression for non-symmetric errors (For real data)
#################Muhammad Irshad########################
#############Roll no P2-14########################
rm(list=ls())
set.seed(1990)
library(MASS)
p=5                                     #No of explanatory variables
I=diag(p)
setwd("D:/Phd work/Hetrocedasticity/Paper Code")      #To change working directory
data=read.csv("live stok.csv",header=TRUE)
data=data[1:12,]
x1=data$x1                # use $ to add variables
x2=data$x2
x3=data$x3
x4=data$x4
x5=data$x5                # use $ to add variables
#x6=data$x6
#x7=data$x7
#x8=data$x8
#x9=data$x9

n=length(x1)
x=cbind(x1,x2,x3,x4,x5)
x=scale(x,center = TRUE,scale = TRUE)   #Scaling x,,, standardizing x and y
y=data$y
y=(y-mean(y))/sd(y)

#y=scale(x,center = TRUE,scale = TRUE)
eigen(cor(x))
ev=eigen(cor(x))$values
vec=eigen(cor(x))$vectors
which.max(ev)
CN=max(ev)/min(ev)        #Computing condition number
CI=sqrt(max(ev)/min(ev))
beta=vec[,which.max(ev)]
C=t(x)%*%(x)
D=eigen(C)$vectors
Z=x%*%D
lam=t(Z)%*%Z
lam=round(lam,4)
lambda=diag(lam)
alpha=t(D)%*%beta

betahat=solve(t(x)%*%x)%*%t(x)%*%y
yhat=x%*%betahat
sigma.hat=(sum((y-yhat)^2))/(n-p)
alphahat=solve(lam)%*%t(Z)%*%y
#y[5]=y[5]+10*sigmahat
#y[15]=y[15]-5*sigmahat
#y[20]=y[20]+20*sigmahat
#y[35]=y[35]-750*sigmahat
#y[40]=y[40]+1000*sigmahat

#Estimation of ridge parameter
#OLS estimator
alphahat.ols=alphahat
alphahat.ols=c(alphahat.ols)

#Ridge Regression estimator of HK method (1970)

HK=sigma.hat/(max(alphahat.ols^2))
K1=HK
alphahatk1=solve(lam+I*K1)%*%lam%*%alphahat.ols
alphahatk1=c(alphahatk1)
alphahatk.HK = alphahatk1

###---------kibria method AM, GM, Median ------#####

KAM<-((sigma.hat/sum(alphahat.ols^2)))/p
K2=KAM
alphahatk2=solve(lam+I*K2)%*%lam%*%alphahat.ols
alphahatk2=c(alphahatk2)
alphahatk.KAM = alphahatk2

KGM<-sigma.hat/((prod(alphahat.ols^2))^1/p)
K3=KGM
alphahatk3=solve(lam+I*K3)%*%lam%*%alphahat.ols
alphahatk3=c(alphahatk3)
alphahatk.KGM = alphahatk3


#KMS 2013  Ridge estimator

W=max(ev)/sum(abs(alphahat.ols))
KMS=W*HK
K5=KMS
alphahatk5=solve(lam+I*K5)%*%lam%*%alphahat.ols
alphahatk5=c(alphahatk5)
alphahatk.KMS = alphahatk5

#########################-------  two parameter ridge estimator (TRE)-------################## 



### Toker and Kaciranlar 2013 Estimator
TK=sigma.hat/(max(alphahat.ols))^2
K6=TK
q.hat.opt.num<-sum(((alphahat.ols^2)*lambda)/(lambda+K6))
q.hat.opt.den<-sum(((sigma.hat*lambda^2)+(alphahat.ols^2*lambda))/(lambda+K6)^2)
q.hat.opt<-q.hat.opt.num/q.hat.opt.den
q=q.hat.opt
q=c(q)
K_opt=(q*sum(sigma.hat/lambda)+(q-1)*sum(lambda^2*alphahat.ols^2))/(sum(lambda*alphahat.ols^2))
K7=K_opt
q.num1<-(t(y)%*%Z)%*%solve(lam+I*K7)%*%(t(Z)%*%y)
q.den1<-(t(y)%*%Z)%*%solve(lam+I*K7)%*%lam%*%solve(lam+I*K7)%*%(t(Z)%*%y)
q1<-q.num1/q.den1
q1<-c(q1)
alphahatk6=q1*(solve(lam+I*K7)%*%lam%*%alphahat.ols)
alphahatk6=c(alphahatk6)
alphahatk.TK = alphahatk6

###------- RIRE two parameter ridge estimator (TRE)-------###
RIRE1<-log(1 + sum(ev * abs(alphahat.ols)) / sigma.hat)
K8=RIRE1
q.num2<-(t(y)%*%Z)%*%solve(lam+I*K8)%*%(t(Z)%*%y)
q.den2<-(t(y)%*%Z)%*%solve(lam+I*K8)%*%lam%*%solve(lam+I*K8)%*%(t(Z)%*%y)
q2<-q.num2/q.den2
q2<-c(q2)
alphahatk7=q2*(solve(lam+I*K8)%*%lam%*%alphahat.ols)
alphahatk7=c(alphahatk7)
alphahatk.RIRE1 = alphahatk7


RIRE3<-(sum(ev * abs(alphahat.ols)))^2/(p*sigma.hat)
K9=RIRE3
q.num3<-(t(y)%*%Z)%*%solve(lam+I*K9)%*%(t(Z)%*%y)
q.den3<-(t(y)%*%Z)%*%solve(lam+I*K9)%*%lam%*%solve(lam+I*K9)%*%(t(Z)%*%y)
q3<-q.num3/q.den3
q3<-c(q3)
alphahatk8=q3*(solve(lam+I*K9)%*%lam%*%alphahat.ols)
alphahatk8=c(alphahatk8)
alphahatk.RIRE3 = alphahatk8

#_____________ Heteroscedasticity ___________________#

beta.hat=solve(t(x)%*%x)%*%t(x)%*%y
y.hat=x%*%beta.hat

e.res=y-y.hat
tah = sqrt((p * sd(e.res)) / IQR(e.res))
sigma.hatro2=tah*sigma.hat

#Ridge Regression estimator of HK method (1970)

HKN=sigma.hatro2/(max(alphahat.ols^2))
K10=HKN
q.num8<-(t(y)%*%Z)%*%solve(lam+I*K10)%*%(t(Z)%*%y)
q.den8<-(t(y)%*%Z)%*%solve(lam+I*K10)%*%lam%*%solve(lam+I*K10)%*%(t(Z)%*%y)
q8<-q.num8/q.den8
q8<-c(q8)
alphahatk9=q8*(solve(lam+I*K10)%*%lam%*%alphahat.ols)
alphahatk9=c(alphahatk9)
alphahatk.HKN = alphahatk9



###---------kibria method AM, GM, Median ------#####

KAMN<-((sigma.hatro2/sum(alphahat.ols^2)))/p
K11=KAMN
q.num9<-(t(y)%*%Z)%*%solve(lam+I*K11)%*%(t(Z)%*%y)
q.den9<-(t(y)%*%Z)%*%solve(lam+I*K11)%*%lam%*%solve(lam+I*K11)%*%(t(Z)%*%y)
q9<-q.num9/q.den9
q9<-c(q9)
alphahatk10=q9*(solve(lam+I*K11)%*%lam%*%alphahat.ols)
alphahatk10=c(alphahatk10)
alphahatk.KAMN = alphahatk10


KGMN<-sigma.hatro2/((prod(alphahat.ols^2))^1/p)
K12=KGMN
q.num10<-(t(y)%*%Z)%*%solve(lam+I*K12)%*%(t(Z)%*%y)
q.den10<-(t(y)%*%Z)%*%solve(lam+I*K12)%*%lam%*%solve(lam+I*K12)%*%(t(Z)%*%y)
q10<-q.num10/q.den10
q10<-c(q10)
alphahatk11=q10*(solve(lam+I*K12)%*%lam%*%alphahat.ols)
alphahatk11=c(alphahatk11)
alphahatk.KGMN = alphahatk11




#KMS 2013  Ridge estimator

W=max(ev)/sum(abs(alphahat))
KMSN=W*HKN
K14=KMSN
q.num11<-(t(y)%*%Z)%*%solve(lam+I*K14)%*%(t(Z)%*%y)
q.den11<-(t(y)%*%Z)%*%solve(lam+I*K14)%*%lam%*%solve(lam+I*K14)%*%(t(Z)%*%y)
q11<-q.num10/q.den11
q11<-c(q11)
alphahatk12=q11*(solve(lam+I*K14)%*%lam%*%alphahat.ols)
alphahatk12=c(alphahatk12)
alphahatk.KMSN = alphahatk12


#########################-------  two parameter ridge estimator (TRE)-------################## 



### Toker and Kaciranlar 2013 Estimator
TKN=sigma.hatro2/(max(alphahat.ols))^2
K15=TKN
q.hat.opt.num4<-sum(((alphahat.ols^2)*lambda)/(lambda+K15))
q.hat.opt.den4<-sum(((sigma.hatro2*lambda^2)+(alphahat.ols^2*lambda))/(lambda+K15)^2)
q.hat.opt4<-q.hat.opt.num4/q.hat.opt.den4
q4=q.hat.opt4
q4=c(q4)
K_opt1=(q4*sum(sigma.hatro2/lambda)+(q4-1)*sum(lambda^2*alphahat.ols^2))/(sum(lambda*alphahat.ols^2))
K16=K_opt1
q.num5<-(t(y)%*%Z)%*%solve(lam+I*K16)%*%(t(Z)%*%y)
q.den5<-(t(y)%*%Z)%*%solve(lam+I*K16)%*%lam%*%solve(lam+I*K16)%*%(t(Z)%*%y)
q5<-q.num5/q.den5
q5<-c(q5)
alphahatk14=q5*(solve(lam+I*K16)%*%lam%*%alphahat.ols)
alphahatk14=c(alphahatk14)
alphahatk.TKN = alphahatk14

###------- RIRE two parameter ridge estimator (TRE)-------###
RIREN1<-log(1 + sum(ev * abs(alphahat.ols)) / sigma.hatro2)
K17=RIREN1
q.num6<-(t(y)%*%Z)%*%solve(lam+I*K17)%*%(t(Z)%*%y)
q.den6<-(t(y)%*%Z)%*%solve(lam+I*K17)%*%lam%*%solve(lam+I*K17)%*%(t(Z)%*%y)
q6<-q.num6/q.den6
q6<-c(q6)
alphahatk15=q6*(solve(lam+I*K17)%*%lam%*%alphahat.ols)
alphahatk15=c(alphahatk15)
alphahatk.RIREN1 = alphahatk15


RIREN3<-(sum(ev * abs(alphahat.ols)))^2/p*sigma.hatro2
K18=RIREN3
q.num7<-(t(y)%*%Z)%*%solve(lam+I*K18)%*%(t(Z)%*%y)
q.den7<-(t(y)%*%Z)%*%solve(lam+I*K18)%*%lam%*%solve(lam+I*K18)%*%(t(Z)%*%y)
q7<-q.num7/q.den7
q7<-c(q7)
alphahatk16=q7*(solve(lam+I*K18)%*%lam%*%alphahat.ols)
alphahatk16=c(alphahatk16)
alphahatk.RIREN3 = alphahatk16

comb.coeff=rbind(alphahat.ols,alphahatk.HK,alphahatk.KAM,alphahatk.KGM,alphahatk.KMS,
                 alphahatk.TK,alphahatk.RIRE1,alphahatk.RIRE3,alphahatk.HKN,alphahatk.KAMN,alphahatk.KGMN,alphahatk.KMSN,
                 alphahatk.TKN,alphahatk.RIREN1,alphahatk.RIREN3)
comb.coeff=round(comb.coeff,10)

###########################################################################################################
#MSE of all estimators
#MSE of OLS and RR

KR1<-c(0,K1,K2,K3,K5,K7,K8,K9) 
a1=rep(0,length(KR1)); b1=rep(0,length(KR1))                            #to obtain MSE of ridge estimators
for(i in 1:length(KR1)){
  
  a1[i]=sum((ev)/(ev+KR1[i])^2)
  b1[i]=sum(((KR1[i]^2)*(alphahat.ols^2))/(ev+KR1[i])^2)
  
}
MSE1=sigma.hat*a1+b1

# MSE of Hetero Ridge

KRH1<-c(K10,K11,K12,K14,K16,K17,K18) 
ah1=rep(0,length(KRH1)); bh1=rep(0,length(KRH1))                            #to obtain MSE of ridge estimators
for(i in 1:length(KRH1)){
  
  ah1[i]=sum((ev)/(ev+KRH1[i])^2)
  bh1[i]=sum(((KRH1[i]^2)*(alphahat.ols^2))/(ev+KRH1[i])^2)
  
}
MSE2=sigma.hatro2*ah1+bh1

MSE=c(MSE1,MSE2)
MSE=as.matrix(MSE)
MSE=round(MSE,4)
K.hat=c(K1,K2,K3,K5,K7,K8,K9,K10,K11,K12,K14,K16,K17,K18)
write.csv(K.hat,file = "D:\\Phd work\\Hetrocedasticity\\two parameter\\Two_parameter_adjusted\\TWO_PARAMETER-REAL-LIFE.csv")
row.names(MSE)=c("OLS", "HK", "KAM", "KGM", "KMS", "TK", "RIRE1", "RIRE3", "HKN", "KAMN", "KGMN", "KMSN", "TKN", "RIREN1", "RIREN3")
colnames(MSE)=c("MSE")
write.csv(MSE,file = "D:\\Phd work\\Hetrocedasticity\\two parameter\\Two_parameter_adjusted\\TWO_PARAMETER-REAL-LIFE.csv")

# Fit OLS model
ols.model <- lm(y ~ x1 + x2 + x3 + x4, data = data)

# Load lmtest package
library(lmtest)

# Breusch-Pagan Test for heteroscedasticity
bptest(ols.model)

write.csv(comb.coeff,file = "D:\\Phd work\\Hetrocedasticity\\two parameter\\Two_parameter_adjusted\\TWO_PARAMETER-REAL-LIFE-coeff.csv")

cor.x=round(cor(x),4)
cor.y=round(cor(x,y),4)
cor.xy=cbind(cor.x,cor.y)
write.csv(cor.xy,file = "D:\\Phd work\\Hetrocedasticity\\two parameter\\Two_parameter_adjusted\\TWO_PARAMETER-REAL-LIFE-corr.csv")
CN
CI
ev

library(corrplot)
# Plot the correlation matrix
corrplot(cor.x, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, addCoef.col = "black", 
         number.cex = 0.7, col = colorRampPalette(c("blue", "white", "red"))(200))

plot(ols.model$fitted.values, resid(ols.model),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")

