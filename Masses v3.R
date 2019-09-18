
# Juillet 2019
# Turgeon Roxanne 
# Initiation à la recherche / Migration assistee / Broutement / Stress_hydrique / Serre

dev.off()
remove(list = ls())

library(tidyverse)
library(plyr)
library(dplyr)
library(readr)
library(purrr)
library(ggplot2)
library(readxl)
library(car)
library(lsmeans)
library(multcompView)

Masses = read.table("C:/Users/Roxanne Turgeon/Dropbox/Initiation ? la recherche/Analyses R/Masses.txt",
                    header=TRUE,na.string = "",dec = ",")
Masses$Brout <- replace(Masses$Brout,Masses$Brout=="0","NoBrout")
Masses$Brout <- replace(Masses$Brout,Masses$Brout=="1","Brout")
Masses$Stress <- replace(Masses$Stress,Masses$Stress=="0","NoStress")
Masses$Stress <- replace(Masses$Stress,Masses$Stress=="1","Stress1")
Masses$Stress <- replace(Masses$Stress,Masses$Stress=="2","Stress2")

Masses$Prov <- as.factor(as.character(Masses$Prov))
Masses$Bloc <- as.factor(as.character(Masses$Bloc))
Masses$Brout <- as.factor(as.character(Masses$Brout))
Masses$Stress <- as.factor(as.character(Masses$Stress))

str(Masses)

CET<-filter(Masses,Esp=="CET")



#################TOUTES ESP?CES###################

boxplot(Mtot~Brout+Esp,data=Masses)
boxplot(Mtot~Stress+Esp,data=Masses)

#################CERISIERS########################

#BIOMASSE TOTALE CERISIERS 
boxplot(Mtot~Brout+Stress+Prov,data=CET,cex.axis=0.5, col = rep(c("red", "green"),6))


model1<-lm(Mtot ~ Brout + Stress + Brout*Stress, data=CET)
anova1<-anova(model1) # Pas d'int?ractions significatives mais effet du Stress2 
anova1

model1<-lm(Mtot ~ Brout + Stress + Brout*Stress + Prov + Brout*Prov + Stress*Prov, data=CET)
anova1<-anova(model1) # Pas d'int?ractions significatives mais effet du Stress2 et Provenance
anova1 # Tiens pas compte du mod?le en tiroirs...

#Conditions d'application ANOVA biomasse totale Cerisiers
plot(model1) #Normalit? & Homog?n?it? des variances

leveneTest(y=CET$Mtot,group=CET$Brout) #Variances homog?nes
leveneTest(y=CET$Mtot,group=CET$Stress) #P-value pr?s de 0.05, on v?rifie avec le boxplot
boxplot(model1$residuals~CET$Stress+CET$Brout)



#A?RIEN/RACINAIRE BROUTEMENT CERISIERS
plot(CET_Broute$`Masse_seche_aerienne(g)`,CET_Broute$`Masse_seche_racinaire(g)`, main="Broutement cerisiers", xlab="Masse aerienne (g)", ylab="Masse racinaire (g)", pch=16,col="green")
abline(lm(CET_Broute$`Masse_seche_racinaire(g)`~ CET_Broute$`Masse_seche_aerienne(g)`-1),col="green")
points(CET_Temoin$`Masse_seche_aerienne(g)`,CET_Temoin$`Masse_seche_racinaire(g)`, pch=16, col="grey")
abline(lm(CET_Temoin$`Masse_seche_racinaire(g)`~ CET_Temoin$`Masse_seche_aerienne(g)`-1),col="grey")
points(CET_Stress2$`Masse_seche_aerienne(g)`,CET_Stress2$`Masse_seche_racinaire(g)`, pch=16,col="dark red")
abline(lm(CET_Stress2$`Masse_seche_racinaire(g)`~ CET_Stress2$`Masse_seche_aerienne(g)`-1),col="dark red")

#Ratio CERISIERS
boxplot(Ratio_biomasse_aerienne_racinaire~Broutement, at=1, data=CET18_Temoin,col=555, main="Ratio a?rien racinaire cerisiers", ylab="Ratio a?rien racinaire",xlim=c(1,12),ylim=c(0,6))
boxplot(Ratio_biomasse_aerienne_racinaire~Broutement, at=2, data=CET18_Broute, add=TRUE,col=555,axis(1,at=1:2,labels=c("T18","B18")))
boxplot(Ratio_biomasse_aerienne_racinaire~Stress_hydrique, at=3, data=CET18_Stress2, add=TRUE,col=555,axis(1,at=3,labels=c("S18")))
boxplot(Ratio_biomasse_aerienne_racinaire~`Stress_hydrique`, at=4, data=CET18_Broute_Stress2, add=TRUE,col=555,axis(1, at=4, labels=c("BS18")))
boxplot(Ratio_biomasse_aerienne_racinaire~Broutement, at=5, data=CET50_Temoin, add=TRUE,col=655, axis(1,at=5,labels=c("T50")))
boxplot(Ratio_biomasse_aerienne_racinaire~Broutement, at=6, data=CET50_Broute, add=TRUE,col=655, axis(1,at=6,labels=c("B50")))
boxplot(Ratio_biomasse_aerienne_racinaire~Stress_hydrique, at=7, data=CET50_Stress2, add=TRUE,col=655, axis(1,at=7,labels=c("S50")))
boxplot(Ratio_biomasse_aerienne_racinaire~`Stress_hydrique`, at=8, data=CET50_Broute_Stress2, add=TRUE,col=655, axis(1, at=8, labels=c("BS50")))
boxplot(Ratio_biomasse_aerienne_racinaire~Broutement, at=9, data=CET80_Temoin, add=TRUE, col=554, axis(1,at=9,labels=c("T80")))
boxplot(Ratio_biomasse_aerienne_racinaire~Broutement, at=10, data=CET80_Broute, add=TRUE,col=554, axis(1,at=10,labels=c("B80")))
boxplot(Ratio_biomasse_aerienne_racinaire~Stress_hydrique, at=11, data=CET80_Stress2, add=TRUE,col=554, axis(1,at=11,labels=c("S80")))
boxplot(Ratio_biomasse_aerienne_racinaire~`Stress_hydrique`, at=12, data=CET80_Broute_Stress2, add=TRUE,col=554,axis(1, at=12, labels=c("BS80")))

model2<-lm(Ratio_biomasse_aerienne_racinaire~Broutement*Stress_hydrique*Provenance*Bloc,data=CET)
anova2<-anova(model2) # Pas d'int?ractions significatives mais Stress2 et Provenance signficatifs
summary(anova2)

#Conditions d'application ANOVA ratio Cerisiers
plot(model2,2) #Normalit?
res_anova2<-residuals(object=model2)
shapiro.test(res_anova2) #Rejette H0, distribution pas normale

plot(model2,1)#Homog?n?it? des variances 
leveneTest(Ratio_biomasse_aerienne_racinaire~Broutement*Stress_hydrique*Provenance*Bloc,data=CET) #Rejette H0, variances pas homog?nes

#Stress_hydrique CERISIERS LOG
plot(CET_Stress2$`log Masse racinaire`, CET_Stress2$`log Masse a?rienne`, main="Stress_hydrique cerisiers", xlab="Masse aerienne (g)", ylab="Masse racinaire (g)", pch=16,col="dark red", xlim = c(0,2), ylim= c(0,2.5))
regression_broute<-lm(CET_Stress2$`log Masse a?rienne`~ CET_Stress2$`log Masse racinaire`+0) # +0 Permet d'?liminer l'ordonn?e ? l'origine
regression_broute # ?quation : Racinaire = 1,465 a?rien
abline(regression_broute,col="dark red")
points(CET_Temoin$`log Masse racinaire`,CET_Temoin$`log Masse a?rienne`, pch=16, col="grey")
regression_temoin<-lm(CET_Temoin$`log Masse a?rienne`~ CET_Temoin$`log Masse racinaire`+0)
regression_temoin # ?quation : Racinaire = 1,266 a?rien
abline(regression_temoin,col="grey")
equation0 = function(x0){x0^1.465}
curve(equation0,add=TRUE)
equation1 = function(x1){0.6413113*x1}
curve(equation1,add=TRUE,col="dark red")

t.test(CET_Temoin$`Pente de la droite LOG`, CET_Stress2$`Pente de la droite LOG`) #Significatif! T?moin 0.7837636 Stress2 0.6413113


#################CHENES########################

#BIOMASSE TOTALE CHENES BROUTEMENT
boxplot(`Biomasse_totale`~Broutement, at=1, data=CHR_Temoin, main="Biomasse ch?mes", ylab="Biomasse",xlim=c(0,3),ylim=c(0,160))
boxplot(`Biomasse_totale`~Broutement, at=2, data=CHR_Broute, add=TRUE,axis(1, at=1:2, labels=c("T?moin", "Brout?")))

#BIOMASSE TOTALE CHENES BROUTEMENT
boxplot(`Biomasse_totale`~ `Stress_hydrique`, at=1, data=CHR_Temoin, main="Biomasse ch?nes", ylab="Biomasse",xlim=c(0,3),ylim=c(0,120))
boxplot(`Biomasse_totale`~ `Stress_hydrique`, at=2, data=CHR_Stress2, add=TRUE,axis(1, at=1:2, labels=c("T?moin", "Stress ?lev?")))

#BROUTEMENT CHENES
plot(CHR_Broute$`Masse_seche_aerienne(g)`,CHR_Broute$`Masse_seche_racinaire(g)`, main="Broutement ch?nes",  xlim=c(0,60), ylim=c(0,60), xlab="Masse aerienne (g)", ylab="Masse racinaire (g)", pch=16,col="green")
abline(lm(CHR_Broute$`Masse_seche_racinaire(g)`~ CHR_Broute$`Masse_seche_aerienne(g)`-1),col="green")
points(CHR_Temoin$`Masse_seche_aerienne(g)`,CHR_Temoin$`Masse_seche_racinaire(g)`, pch=16, col="grey")
abline(lm(CHR_Temoin$`Masse_seche_racinaire(g)`~ CHR_Temoin$`Masse_seche_aerienne(g)`-1),col="grey")

#Ratio broutement CHENES
t.test(CHR_Temoin$`Masse_seche_racinaire(g)`, CHR_Broute$`Masse_seche_racinaire(g)`) # Non significatif
t.test(CHR_Temoin$`Masse_seche_aerienne(g)`, CHR_Broute$`Masse_seche_aerienne(g)`) #Non significatif
t.test(CHR_Temoin$'Ratio_biomasse_aerienne_racinaire',CHR_Broute$'Ratio_biomasse_aerienne_racinaire') #Non significatif

#Stress_hydrique CHENES
plot(CHR_Stress2$`Masse_seche_aerienne(g)`,CHR_Stress2$`Masse_seche_racinaire(g)`,main="Stress_hydrique ch?nes", xlim=c(0,60), ylim=c(0,60),xlab="Masse aerienne (g)", ylab="Masse racinaire (g)", pch=16,col="dark red")
abline(lm(CHR_Stress2$`Masse_seche_racinaire(g)`~ CHR_Stress2$`Masse_seche_aerienne(g)`-1),col="dark red")
points(CHR_Temoin$`Masse_seche_aerienne(g)`,CHR_Temoin$`Masse_seche_racinaire(g)`, pch=16, col="grey")
abline(lm(CHR_Temoin$`Masse_seche_racinaire(g)`~ CHR_Temoin$`Masse_seche_aerienne(g)`-1),col="grey")

#Ratio Stress_hydrique CHENES
t.test(CHR_Temoin$`Masse_seche_racinaire(g)`, CHR_Stress2$`Masse_seche_racinaire(g)`) #Non significatif
t.test(CHR_Temoin$`Masse_seche_aerienne(g)`, CHR_Stress2$`Masse_seche_aerienne(g)`)  #Non significatif
t.test(CHR_Temoin$'Ratio_biomasse_aerienne_racinaire',CHR_Stress2$'Ratio_biomasse_aerienne_racinaire') #Presque significatif

t.test(CHR18_Temoin$'Ratio_biomasse_aerienne_racinaire',CHR18_Stress2$'Ratio_biomasse_aerienne_racinaire') #Non significatif
t.test(CHR50_Temoin$'Ratio_biomasse_aerienne_racinaire',CHR50_Stress2$'Ratio_biomasse_aerienne_racinaire') #Non significatif
t.test(CHR80_Temoin$'Ratio_biomasse_aerienne_racinaire',CHR80_Stress2$'Ratio_biomasse_aerienne_racinaire') #Non significatif


#################?RABLES########################

#BIOMASSE TOTALE ERABLES BROUTEMENT
boxplot(`Biomasse_totale`~Broutement, at=1, data=ERS_Temoin, main="Biomasse ?rables", ylab="Biomasse",xlim=c(0,3),ylim=c(0,160))
boxplot(`Biomasse_totale`~Broutement, at=2, data=ERS_Broute, add=TRUE,axis(1, at=1:2, labels=c("T?moin", "Brout?")))

#BIOMASSE TOTALE ERABLES BROUTEMENT
boxplot(`Biomasse_totale`~ `Stress_hydrique`, at=1, data=ERS_Temoin, main="Biomasse ?rables", ylab="Biomasse",xlim=c(0,4),ylim=c(0,160))
boxplot(`Biomasse_totale`~ `Stress_hydrique`, at=2, data=ERS_Stress1, add=TRUE)
boxplot(`Biomasse_totale`~ `Stress_hydrique`, at=3, data=ERS_Stress2, add=TRUE,axis(1, at=1:3, labels=c("T?moin","Stress mod?r?", "Stress ?lev?")))

#BROUTEMENT ERABLES
plot(ERS_Broute$`Masse_seche_aerienne(g)`,ERS_Broute$`Masse_seche_racinaire(g)`, main="Broutement ?rables", xlim=c(0,60), ylim=c(0,50), xlab="Masse aerienne (g)", ylab="Masse racinaire (g)", pch=16,col="green")
abline(lm(ERS_Broute$`Masse_seche_racinaire(g)`~ ERS_Broute$`Masse_seche_aerienne(g)`-1),col="green")
points(ERS_Temoin$`Masse_seche_aerienne(g)`,ERS_Temoin$`Masse_seche_racinaire(g)`, pch=16, col="grey")
abline(lm(ERS_Temoin$`Masse_seche_racinaire(g)`~ ERS_Temoin$`Masse_seche_aerienne(g)`-1),col="grey")

#Ratio broutement ERABLES
t.test(ERS_Temoin$`Masse_seche_racinaire(g)`, ERS_Broute$`Masse_seche_racinaire(g)`) #Non significatif
t.test(ERS_Temoin$`Masse_seche_aerienne(g)`, ERS_Broute$`Masse_seche_aerienne(g)`) #Non significatif
t.test(ERS_Temoin$'Ratio_biomasse_aerienne_racinaire',ERS_Broute$'Ratio_biomasse_aerienne_racinaire') #Non significatif

#Stress_hydrique ERABLES
plot(ERS_Stress2$`Masse_seche_aerienne(g)`,ERS_Stress2$`Masse_seche_racinaire(g)`,main="Stress_hydrique ?rables" , xlim = c(0,60), ylim=c(0,50),xlab="Masse aerienne (g)",ylab="Masse racinaire (g)", pch=16,col="dark red")
abline(lm(ERS_Stress2$`Masse_seche_racinaire(g)`~ ERS_Stress2$`Masse_seche_aerienne(g)`-1),col="dark red")
points(ERS_Temoin$`Masse_seche_aerienne(g)`,ERS_Temoin$`Masse_seche_racinaire(g)`, pch=16, col="grey")
abline(lm(ERS_Temoin$`Masse_seche_racinaire(g)`~ ERS_Temoin$`Masse_seche_aerienne(g)`-1),col="grey")

#Ratio Stress_hydrique ERABLES
t.test(ERS_Temoin$`Masse_seche_racinaire(g)`, ERS_Stress2$`Masse_seche_racinaire(g)`) #Non significatif
t.test(ERS_Temoin$`Masse_seche_aerienne(g)`, ERS_Stress2$`Masse_seche_aerienne(g)`)  #Non significatif
t.test(ERS_Temoin$'Ratio_biomasse_aerienne_racinaire',ERS_Stress2$'Ratio_biomasse_aerienne_racinaire') #Nombre insuffisant


#################THUYAS########################

#BIOMASSE TOTALE THUYAS BROUTEMENT
boxplot(`Biomasse_totale`~Broutement, at=1, data=THO_Temoin, main="Biomasse thuyas", ylab="Biomasse",xlim=c(0,3),ylim=c(0,90))
boxplot(`Biomasse_totale`~Broutement, at=2, data=THO_Broute, add=TRUE,axis(1, at=1:2, labels=c("T?moin", "Brout?")))

#BIOMASSE TOTALE THUYAS BROUTEMENT
boxplot(`Biomasse_totale`~ `Stress_hydrique`, at=1, data=THO_Temoin, main="Biomasse thuyas", ylab="Biomasse",xlim=c(0,4),ylim=c(0,90))
boxplot(`Biomasse_totale`~ `Stress_hydrique`, at=2, data=THO_Stress1, add=TRUE)
boxplot(`Biomasse_totale`~ `Stress_hydrique`, at=3, data=THO_Stress2, add=TRUE,axis(1, at=1:3, labels=c("T?moin","Stress mod?r?", "Stress ?lev?")))

#BROUTEMENT THUYA
plot(THO_Broute$`Masse_seche_aerienne(g)`,THO_Broute$`Masse_seche_racinaire(g)`, main="Broutement thuyas",  xlim=c(0,80), ylim=c(0,60), xlab="Masse aerienne (g)", ylab="Masse racinaire (g)", pch=16,col="green")
abline(lm(THO_Broute$`Masse_seche_racinaire(g)`~ THO_Broute$`Masse_seche_aerienne(g)`-1),col="green")
points(THO_Temoin$`Masse_seche_aerienne(g)`,THO_Temoin$`Masse_seche_racinaire(g)`, pch=16, col="grey")
abline(lm(THO_Temoin$`Masse_seche_racinaire(g)`~ THO_Temoin$`Masse_seche_aerienne(g)`-1),col="grey")

#Ratio broutement THUYA
t.test(THO_Temoin$`Masse_seche_racinaire(g)`, THO_Broute$`Masse_seche_racinaire(g)`) #Non significatif
t.test(THO_Temoin$`Masse_seche_aerienne(g)`, THO_Broute$`Masse_seche_aerienne(g)`) #Non significatif
t.test(THO_Temoin$'Ratio_biomasse_aerienne_racinaire',THO_Broute$'Ratio_biomasse_aerienne_racinaire') #Non significatif (presque)

#Stress_hydrique THUYA
plot(THO_Stress2$`Masse_seche_aerienne(g)`,THO_Stress2$`Masse_seche_racinaire(g)`,main="Stress_hydrique thuyas", xlim=c(0,100), ylim=c(0,50),xlab="Masse aerienne (g)", ylab="Masse racinaire (g)", pch=16,col="dark red")
abline(lm(THO_Stress2$`Masse_seche_racinaire(g)`~ THO_Stress2$`Masse_seche_aerienne(g)`-1),col="dark red")
points(THO_Temoin$`Masse_seche_aerienne(g)`,THO_Temoin$`Masse_seche_racinaire(g)`, pch=16, col="grey")
abline(lm(THO_Temoin$`Masse_seche_racinaire(g)`~ THO_Temoin$`Masse_seche_aerienne(g)`-1),col="grey")

#Ratio Stress_hydrique THUYA
t.test(THO_Temoin$`Masse_seche_racinaire(g)`, THO_Stress2$`Masse_seche_racinaire(g)`) #Non significatif
t.test(THO_Temoin$`Masse_seche_aerienne(g)`, THO_Stress2$`Masse_seche_aerienne(g)`)  #Thuyas sous Stress_hydrique ?lev? plus petite biomasse a?rienne
t.test(THO_Temoin$'Ratio_biomasse_aerienne_racinaire',THO_Stress2$'Ratio_biomasse_aerienne_racinaire') #Non significatif

t.test(CHR18_Temoin$'Ratio_biomasse_aerienne_racinaire',CHR18_Stress2$'Ratio_biomasse_aerienne_racinaire') #Non significatif
t.test(CHR50_Temoin$'Ratio_biomasse_aerienne_racinaire',CHR50_Stress2$'Ratio_biomasse_aerienne_racinaire') #Non significatif
t.test(CHR80_Temoin$'Ratio_biomasse_aerienne_racinaire',CHR80_Stress2$'Ratio_biomasse_aerienne_racinaire') #Non significatif


#################PINS########################

#BIOMASSE TOTALE PINS BROUTEMENT
boxplot(`Biomasse_totale`~Broutement, at=1, data=PIB_Temoin, main="Biomasse pins", ylab="Biomasse",xlim=c(0,3),ylim=c(0,50))
boxplot(`Biomasse_totale`~Broutement, at=2, data=PIB_Broute, add=TRUE,axis(1, at=1:2, labels=c("T?moin", "Brout?")))

#BIOMASSE TOTALE PINS BROUTEMENT
boxplot(`Biomasse_totale`~ `Stress_hydrique`, at=1, data=PIB_Temoin, main="Biomasse pins", ylab="Biomasse",xlim=c(0,4),ylim=c(0,60))
boxplot(`Biomasse_totale`~ `Stress_hydrique`, at=2, data=PIB_Stress1, add=TRUE)
boxplot(`Biomasse_totale`~ `Stress_hydrique`, at=3, data=PIB_Stress2, add=TRUE,axis(1, at=1:3, labels=c("T?moin","Stress mod?r?", "Stress ?lev?")))

#BROUTEMENT PIN
plot(PIB_Broute$`Masse_seche_aerienne(g)`,PIB_Broute$`Masse_seche_racinaire(g)`, main="Broutement pins",  xlim=c(0,30), ylim=c(0,20), xlab="Masse aerienne (g)", ylab="Masse racinaire (g)", pch=16,col="green")
abline(lm(PIB_Broute$`Masse_seche_racinaire(g)`~ PIB_Broute$`Masse_seche_aerienne(g)`-1),col="green")
points(PIB_Temoin$`Masse_seche_aerienne(g)`,PIB_Temoin$`Masse_seche_racinaire(g)`, pch=16, col="grey")
abline(lm(PIB_Temoin$`Masse_seche_racinaire(g)`~ PIB_Temoin$`Masse_seche_aerienne(g)`-1),col="grey")

#Ratio broutement PINS
t.test(PIB_Temoin$`Masse_seche_racinaire(g)`, PIB_Broute$`Masse_seche_racinaire(g)`) #Pins brout?s ont moins de biomasse racinaire
t.test(PIB_Temoin$`Masse_seche_aerienne(g)`, PIB_Broute$`Masse_seche_aerienne(g)`) #Pins brout?s ont moins de biomasse racinaire
t.test(PIB_Temoin$'Ratio_biomasse_aerienne_racinaire',PIB_Broute$'Ratio_biomasse_aerienne_racinaire') #Non significatif (presque)

#Stress_hydrique PINS
plot(PIB_Stress2$`Masse_seche_aerienne(g)`,PIB_Stress2$`Masse_seche_racinaire(g)`,main="Stress_hydrique pins", xlim=c(0,40), ylim=c(0,20),xlab="Masse aerienne (g)", ylab="Masse racinaire (g)", pch=16,col="dark red")
abline(lm(PIB_Stress2$`Masse_seche_racinaire(g)`~ PIB_Stress2$`Masse_seche_aerienne(g)`-1),col="dark red")
points(PIB_Temoin$`Masse_seche_aerienne(g)`,PIB_Temoin$`Masse_seche_racinaire(g)`, pch=16, col="grey")
abline(lm(PIB_Temoin$`Masse_seche_racinaire(g)`~ PIB_Temoin$`Masse_seche_aerienne(g)`-1),col="grey")

#Ratio Stress_hydrique PINS
t.test(PIB_Temoin$`Masse_seche_racinaire(g)`, PIB_Stress2$`Masse_seche_racinaire(g)`) #Non significatif
t.test(PIB_Temoin$`Masse_seche_aerienne(g)`, PIB_Stress2$`Masse_seche_aerienne(g)`)  #Non significatif
t.test(PIB_Temoin$'Ratio_biomasse_aerienne_racinaire',PIB_Stress2$'Ratio_biomasse_aerienne_racinaire') #Non significatif

t.test(CHR18_Temoin$'Ratio_biomasse_aerienne_racinaire',CHR18_Stress2$'Ratio_biomasse_aerienne_racinaire') #Non significatif
t.test(CHR50_Temoin$'Ratio_biomasse_aerienne_racinaire',CHR50_Stress2$'Ratio_biomasse_aerienne_racinaire') #Non significatif
t.test(CHR80_Temoin$'Ratio_biomasse_aerienne_racinaire',CHR80_Stress2$'Ratio_biomasse_aerienne_racinaire') #Non significatif

#SURVIE PINS
Mortalite = read_excel(file,7)
PIB_Survie<-filter(Mortalite,Espece=="PIB")
model_Broutement <- glm(PIB_Survie$Broutement~PIB_Survie$`Mortalite`, family=binomial(link='logit'))
anova(model_Broutement,test="Chisq") 

model_Hydrique <- glm(PIB_Survie$`Stress_hydrique`~PIB_Survie$`Mortalite`, family=binomial(link='logit'))
anova(model,test="Chisq") # Peut pas le faire avec trois niveaux de traitements...
