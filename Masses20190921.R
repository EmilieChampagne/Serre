
# Spetembre-Octobre 2019
# Turgeon Roxanne 
# Initiation à la recherche / Migration assistee / Broutement / Stress_hydrique / Serre

dev.off()
remove(list = ls())

library(plyr)
library(tidyverse) 
library(readr)
library(purrr)
library(ggplot2)
library(readxl)
library(car)
library(lsmeans)
library(multcompView)
library(nlme)
library(lmerTest)
library(agricolae)


Masses = read.table("C:./Masses.txt",header=TRUE,na.string = "",dec = ",") 
Masses2 <- read.table("./Masses2.txt", header=TRUE, na.string = "", dec = ",") 

Masses2$Brout <- replace(Masses2$Brout,Masses2$Brout=="0","NoBrout")
Masses2$Brout <- replace(Masses2$Brout,Masses2$Brout=="1","Brout")
Masses2$Stress <- replace(Masses2$Stress,Masses2$Stress=="0","NoStress")
Masses2$Stress <- replace(Masses2$Stress,Masses2$Stress=="1","Stress1")
Masses2$Stress <- replace(Masses2$Stress,Masses2$Stress=="2","Stress2")

Masses2$Prov <- as.factor(as.character(Masses2$Prov))
Masses2$Bloc <- as.factor(as.character(Masses2$Bloc))
Masses2$Brout <- as.factor(as.character(Masses2$Brout))
Masses2$Stress <- as.factor(as.character(Masses2$Stress))

Masses2$Hauteurini <- Masses2$Hauteurini

str(Masses2)

CET<-filter(Masses2,Esp=="CET")
CHR<-filter(Masses2,Esp=="CHR")
ERS<-filter(Masses2,Esp=="ERS")
THO<-filter(Masses2,Esp=="THO")
PIB<-filter(Masses2,Esp=="PIB")

#################TOUTES ESPECES###################

boxplot(Mtot~Brout+Esp,data=Masses2,cex.axis=0.5, las=2, ylab="Masse tot", xlab="")
boxplot(Mtot~Stress+Esp,data=Masses2,cex.axis=0.5, las=2, ylab="Masse tot", xlab="")
boxplot(Mtot~Bloc*Esp,data=Masses2, cex.axis=0.5, las=2, ylab="Masse tot", xlab="")

#################CERISIERS########################

###BIOMASSE TOTALE CERISIERS###

#Regarder si hauteur initiale est un bon estimateur de la masse totale
#Regarder la relation entre hauteur finale et la masse a la fin de l'experience
cor.test(CET$Hauteur, CET$Mtot)
plot(CET$Hauteur, CET$Mtot) 
#Oui, bonne correlation --> Utiliser hauteur initiale pour corriger pour la biomasse initiale des plants

#Savoir le nb de replicats 
ddply(CET, c("Stress", "Brout", "Prov", "Bloc"), summarise,N = length(Mtot))

#Pour la majorite des combinaisons, nous n'avons qu'un replicat
#Pour simplifier les analyses, on fait une moyenne pour les rares cas 
#ou nous avons 2 plants par combinaison

CET_mean <- aggregate(list(CET$Mtot,CET$Hauteurini,CET$Ratio), by= data.frame(CET$Stress, CET$Brout, CET$Prov, CET$Bloc), FUN= "mean")
colnames(CET_mean)[1:7] <- c("Stress", "Brout", "Prov", "Bloc","Mtot","Hauteurini","Ratio")
CET_mean <- mutate(CET_mean, Stress=factor(Stress, levels=c("NoStress", "Stress2")))#On indique qu'il n'y a pas de niveau "Stress1" pour le facteur Stress
CET <- mutate(CET, Stress=factor(Stress, levels=c("NoStress", "Stress2")))

ddply(CET_mean, c("Stress", "Brout", "Prov", "Bloc","Hauteurini"), summarise,
      N = length(Mtot)) #Verification n = 1 pour chaque combinaison

#Boxplot masses tot 
levels(CET$Brout) <- c("NoBrout","Brout")
levels(CET$Stress) <- c("NoStress","Stress")
color = c(rep("green",4),rep("yellow",4),rep("red",4))
boxplot(Mtot~Brout+Stress+Prov,data=CET,cex.axis=0.5, col = color,las=2, 
        main = "Masses tot cerisiers", ylab="Masse tot", xlab="") 

#Modèle avec plan en tiroir
A <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CET_mean)
#Message veut dire que pour cette espèce, l'effet aléatoire du bloc = 0
#On le garde pareil dans le modèle, mais il a peu d'effet
anova(A)
# Il y a seulement l'intéraction stress*brout*Hauteurini de significatif,
# mais avec un F de 4.95 contre 5.76 pour Hini seul --> garder Hini seul 
#(mtot ~ Stress*Brout*Prov+ Hauteurini + (1|Bloc/Stress)) 

#Si les interactions avec Hauteurini ne sont pas ou peu significatif, 
#refaire l’analyse avec +Hauteurini, vérifier si la covariable est significative. 
#Si oui, on la garde et sinon on refait l’analyse sans.

#On refait sans interaction avec Hini:
A <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CET_mean)
anova(A) #Hauteurini est significatif, on garde donc le modèle ainsi. On a un effet stress
summary(A)
boxplot(Mtot~Stress,data=CET) 

##Vérifier résidus en incluant les effets aléatoires (effet du design)
plot(A) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône) --> ok
leveneTest(resid(A) ~ Stress*Brout*Prov, data = CET_mean) #Peut pas inclure Hauteurini car variable continue
  #Ni prendre en compte l'effet du bloc, correct???

qqnorm(resid(A)) #Graphique de normalité ok
qqline(resid(A))
shapiro.test(resid(A)) #On vérifie avec un test, faire des transformations au besoin --> ok

####Ratio CERISIERS####

#Boxplot ratio 
boxplot(Ratio~Brout+Stress+Prov,data=CET,cex.axis=0.5, col = color,las=2, main = "Ratio cerisiers", 
        xlab="",ylab="Masse tot")

##Modele avec plan en tiroir
B <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CET_mean)
anova(B)   #Hauteur ini pas significatif, son interaction non plus --> On retire du modèle

##Modele avec plan en tiroir sans interaction Hini
B <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov+Hauteurini + (1|Bloc/Stress), data = CET_mean)
anova(B)   #Hauteur ini pas significatif

##Modele avec plan en tiroir sans Hini
B <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CET_mean)
anova(B)   #Provenance et Stress significatif 
summary(B) #Provenance 2018 et 2080 diffère mais pas 2018 et 2050

#Savoir si 2050 et 2080 sont significativement différent (on place 2050 comme valeur de base de reference)
CET_mean$Prov <- relevel(CET_mean$Prov, ref="2050")
Test <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CET_mean)
summary(Test) #Oui, Prov 2050 diffère de 2080

#Représentation de l'effet stress
boxplot(Ratio~Stress,las=2,cex.axis=0.7,xlab="",data=CET)

#Représentation de l'effet Prov
boxplot(Ratio~Prov,las=2,cex.axis=0.7,xlab="",data=CET) 
text(c(1,2,3),pos=3, offset=10, c("a","a","b"))

#Représentation de l'effet stress et Prov
color2 = c(rep("green",2),rep("yellow",2),rep("red",2))
boxplot(Ratio~Stress+Prov,las=2,cex.axis=0.7,col=color2,xlab="",data=CET) 
text(c(1,2,3),pos = 3, offset =18, c("a","a","b"))


##Vérifier résidus 
plot(B) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône) ok
leveneTest(resid(B) ~ Stress*Brout*Prov, data = CET_mean) #ok

qqnorm(resid(B)) #Graphique de normalité ok
qqline(resid(B))
shapiro.test(resid(B)) #On vérifie avec un test, faire des transformations au besoin --> ok

#################CHENES########################

###BIOMASSE TOTALE CHENES###

#Est-ce que la hauteur initiale est un bon estimateur de la masse totale?
#Regarder la relation entre hauteur finale et la masse a la fin de l'experience
cor.test(CHR$Hauteur, CHR$Mtot)
plot(CHR$Hauteur, CHR$Mtot) 
#Oui, bonne correlation, utiliser hauteurini pour corriger biomasse initiale des plants

#Savoir le nb de replicats 
ddply(CHR, c("Stress", "Brout", "Prov", "Bloc"), summarise,N = length(Mtot))

#Pour la majorite des combinaisons, nous n'avons qu'un replicat
#Pour simplifier les analyses, on fait une moyenne pour les rares cas 
# ou nous avons 2 plants par combinaison

CHR_mean <- aggregate(list(CHR$Mtot,CHR$Hauteurini,CHR$Ratio), by= data.frame(CHR$Stress, CHR$Brout, CHR$Prov, CHR$Bloc), FUN= "mean")
colnames(CHR_mean)[1:7] <- c("Stress", "Brout", "Prov", "Bloc","Mtot","Hauteurini","Ratio")
CHR_mean <- mutate(CHR_mean, Stress=factor(Stress, levels=c("NoStress", "Stress2")))#On indique qu'il n'y a pas de niveau "Stress1" pour le facteur Stress
CHR <- mutate(CHR, Stress=factor(Stress, levels=c("NoStress", "Stress2")))

#Verification: nous avons bien un n = 1 pour chaque combinaison.
ddply(CHR_mean, c("Stress", "Brout", "Prov", "Bloc","Hauteurini"), summarise,
      N = length(Mtot))

#Boxplot masses tot 
levels(CHR$Brout) <- c("NoBrout","Brout")
levels(CHR$Stress) <- c("NoStress","Stress")
color = c(rep("green",4),rep("yellow",4),rep("red",4))
boxplot(Mtot~Brout+Stress+Prov,data=CHR,cex.axis=0.5, col = color,las=2, main = "Masses tot chenes",
        xlab="",ylab="Masse tot")

#Modèle avec plan en tiroir
C <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(C) # Garder Hini seul 

#Sans interaction avec Hini:
CHR_mean$Prov <- relevel(CHR_mean$Prov, ref="2018")
C <- lmerTest::lmer(Mtot ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(C) #Hauteurini significatif, Intéraction Stress:Provenance, tendance de Stress
summary(C) #Stress2:Prov2050 et Stress2:Prov 2080 non significatis (comparaison avec 2018)

#Avec Prov2050 comme valeur de base de reference
CHR_mean$Prov <- relevel(CHR_mean$Prov, ref="2050")
C2 <- lmerTest::lmer(Mtot ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
summary(C2) #Stress2:Prov2080 non significatif (comparaison avec 2050)

####QUESTION
#Avec l'anova on voit que l'interaction Stress*Prov est significative alors que dans le summary,
#toutes les combinaisons (2018 avec 2080, 2018 avec 2050, 2050 avec 2080) ne le sont pas...
#Est-ce parce que l'anova est plus puissante? Si oui, comment on peut savoir quelle est la/les 
#différences significatives? Est ce qu'on peut se fier au summary pour déterminer ce qui est
#significatif ou non?


#Representation tendance de l'effet stress
boxplot(Mtot~Stress,data=CHR,cex.axis=0.7,las=2,xlab="")

#Representation interaction Prov*Stress
color2 = c(rep("green",2),rep("yellow",2),rep("red",2))
boxplot(Mtot~Stress+Prov,data=CHR,col=color2,cex.axis=0.7,las=2,xlab="")

##Vérifier résidus en incluant les effets aléatoires (effet du design)
plot(C) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône) --> ok
leveneTest(resid(C) ~ Stress*Brout*Prov, data = CHR_mean) #ok

qqnorm(resid(C)) #Graphique de normalité ok
qqline(resid(C))
shapiro.test(resid(C)) #On vérifie avec un test, faire des transformations au besoin --> ok

####Ratio CHENES####

#Boxplot ratio 
boxplot(Ratio~Brout+Stress+Prov,data=CHR,cex.axis=0.5, col = color,las=2, main = "Ratio chenes", 
        xlab="",ylab="Masse tot")

##Modele avec plan en tiroir
D <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(D) #Rien significatif

##Modele avec plan en tiroir sans interaction Hini
D <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(D) #Rien significatif

##Modele avec plan en tiroir sans Hini
D <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CHR_mean)
anova(D) #Tendance Stress

boxplot(Ratio~Stress,data=CHR,cex.axis=0.7,las=2,xlab="")

##Vérifier résidus 
plot(D) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(D) ~ Stress*Brout*Prov, data = CHR_mean) #ok

qqnorm(resid(B)) #Graphique de normalité ok
qqline(resid(B))
shapiro.test(resid(B)) #ok

#################ERABLES######################## SEULEMENT PROVENANCE 2018

###BIOMASSE TOTALE ERABLES###

#Est-ce que la hauteur initiale est un bon estimateur de la masse totale?
#Regarder la relation entre hauteur finale et la masse a la fin de l'experience
cor.test(ERS$Hauteur, ERS$Mtot)
plot(ERS$Hauteur, ERS$Mtot) 
#Oui, bonne correlation, utiliser hauteurini pour corriger biomasse initiale des plants

#Savoir le nb de replicats 
ddply(ERS, c("Stress", "Brout", "Bloc"), summarise,N = length(Mtot))
#Un replicat pour toutes les combinaisons

#Boxplot masses tot 
levels(ERS$Brout) <- c("NoBrout","Brout")
levels(ERS$Stress) <- c("NoStress","Stress1","Stress2")
color = c(rep("green",6))
boxplot(Mtot~Brout+Stress,data=ERS,cex.axis=0.5, col = color,las=2, main = "Masses tot Erables", 
        xlab="",ylab="Masse tot")

#Modèle avec plan en tiroir
E <- lmerTest::lmer( Mtot ~ Stress*Brout*Hauteurini + (1|Bloc/Stress), data = ERS)
anova(E) # Rien de significatif

#Sans interaction avec Hini:
E <- lmerTest::lmer( Mtot ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = ERS)
anova(E) #Hauteurini tendance a etre significatif, InteractionStress:Brout significatif
summary(E)

#Si on retire Hini
E <- lmerTest::lmer( Mtot ~ Stress*Brout + (1|Bloc/Stress), data = ERS)
anova(E) #Interaction Stress:Brout significatif
summary(E)

##Vérifier résidus en incluant les effets aléatoires (effet du design)
plot(E) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône) --> ok
leveneTest(resid(E) ~ Stress*Brout, data = ERS) #ok

qqnorm(resid(E)) #Graphique de normalité ok
qqline(resid(E))
shapiro.test(resid(E)) #ok

####Ratio ERABLES####

#Boxplot ratio 
boxplot(Ratio~Brout+Stress,data=ERS,cex.axis=0.5, col = color,las=2, main = "Ratio Erables", 
        xlab="",ylab="Masse tot")

##Modele avec plan en tiroir
G <- lmerTest::lmer( Ratio ~ Stress*Brout*Hauteurini + (1|Bloc/Stress), data = ERS)
anova(G) #Rien significatif, tendance de Brout et Interaction Brout:Hauteurini

##Modele avec plan en tiroir sans interaction Hini
G <- lmerTest::lmer( Ratio ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = ERS)
anova(G) #Rien significatif

##Modele avec plan en tiroir sans Hini
G <- lmerTest::lmer( Ratio ~ Stress*Brout + (1|Bloc/Stress), data = ERS)
anova(G) #Rien significatif

##Vérifier résidus 
plot(G) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(G) ~ Stress*Brout, data = ERS) #ok

qqnorm(resid(G)) #Graphique de normalité ok
qqline(resid(G))
shapiro.test(resid(G)) #ok


#################THUYAS########################

###BIOMASSE TOTALE THUYAS###

#Est-ce que la hauteur initiale est un bon estimateur de la masse totale?
#Regarder la relation entre hauteur finale et la masse a la fin de l'experience
cor.test(THO$Hauteur, THO$Mtot)
plot(THO$Hauteur, THO$Mtot) 
#Oui, bonne correlation, utiliser hauteurini pour corriger biomasse initiale des plants

#Savoir le nb de replicats 
ddply(THO, c("Stress", "Brout", "Prov", "Bloc"), summarise,N = length(Mtot))
#Un replicat pour toutes les combinaisons

#Boxplot masses tot 
levels(THO$Brout) <- c("NoBrout","Brout")
levels(THO$Stress) <- c("NoStress","Stress1","Stress2")
color = c(rep("green",6),rep("yellow",6),rep("red",6))
boxplot(Mtot~Brout+Stress+Prov,data=THO,cex.axis=0.5, col = color,las=2, main = "Masses tot Thuyas",
        xlab="", ylab="Masse tot")

#Modèle avec plan en tiroir
H <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = THO)
anova(H)
# Intéraction stress*brout*Hauteurini et stress*brout de significatifs,
# mais avec des F de 3.24 et 3.44 contre 21.26 pour Hini seul --> garder Hini seul 

#On refait sans interaction avec Hini:
H <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = THO)
anova(H) #Hauteurini, Stress et Brout significatifs
summary(H) #Je n'arrive pas à voir quel niveau de stress est significativement différent duquel? 

boxplot(Mtot~Brout+Stress,data=THO,cex.axis=0.5, col = color,las=2, main = "Masses tot Thuyas",
        xlab="", ylab="Masse tot")

##Vérifier résidus en incluant les effets aléatoires (effet du design)
plot(H) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône) --> ok
leveneTest(resid(H) ~ Stress*Brout*Prov, data = THO) #ok

qqnorm(resid(H)) #Graphique de normalité ok
qqline(resid(H))
shapiro.test(resid(H)) #ok

####Ratio THUYAS####

#Boxplot ratio 
boxplot(Ratio~Brout+Stress+Prov,data=THO,cex.axis=0.5, col = color,las=2, main = "Ratio thuyas", 
        xlab="",ylab="Masse tot")

##Modele avec plan en tiroir
I <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = THO)
anova(I) #Interaction Stress*Brout*Prov significatif

##Modele avec plan en tiroir sans interaction Hini
I <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = THO)
anova(I) #Tendance Stress, Interaction Stress*Brout*Prov significatif

##Modele avec plan en tiroir sans Hini
I <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = THO)
anova(I) #Tendance Stress, Interaction Stress*Brout*Prov significatif
summary(I)

##Vérifier résidus 
plot(I) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(I) ~ Stress*Brout*Prov, data = THO) #ok

qqnorm(resid(I)) #Graphique de normalité ok
qqline(resid(I))
shapiro.test(resid(I)) #ok


#################PINS########################

###BIOMASSE TOTALE PINS###

#Est-ce que la hauteur initiale est un bon estimateur de la masse totale?
#Regarder la relation entre hauteur finale et la masse a la fin de l'experience
cor.test(PIB$Hauteur, PIB$Mtot)
plot(PIB$Hauteur, PIB$Mtot) 
#Oui, bonne correlation, utiliser hauteurini pour corriger biomasse initiale des plants

#Savoir le nb de replicats 
ddply(PIB, c("Stress", "Brout", "Prov", "Bloc"), summarise,N = length(Mtot))
#Un replicat pour toutes les combinaisons

#Boxplot masses tot 
levels(PIB$Brout) <- c("NoBrout","Brout")
levels(PIB$Stress) <- c("NoStress","Stress1","Stress2")
color = c(rep("green",6),rep("yellow",6),rep("red",6))
boxplot(Mtot~Brout+Stress+Prov,data=PIB,cex.axis=0.5, col = color,las=2, main = "Masses tot Pins", 
        xlab="", ylab="Masse tot")

#Modèle avec plan en tiroir
J <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = PIB)
anova(J) #Rien significatif

#Sans interaction avec Hini:
J <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = PIB)
anova(J) #Hauteurini, Brout significatifs
summary(J)

boxplot(Mtot~Brout,data=PIB,cex.axis=0.5, col = color,las=2, main = "Masses tot Pins", 
        xlab="", ylab="Masse tot")

##Vérifier résidus en incluant les effets aléatoires (effet du design)
plot(J) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône) --> ok
leveneTest(resid(J) ~ Stress*Brout*Prov, data = PIB) #ok

qqnorm(resid(J)) #Graphique de normalité ok
qqline(resid(J))
shapiro.test(resid(J)) #ok

####Ratio PINS####

#Boxplot ratio 
boxplot(Ratio~Brout+Stress+Prov,data=PIB,cex.axis=0.5, col = color,las=2, main = "Ratio pins", 
        xlab="",ylab="Masse tot")

##Modele avec plan en tiroir
K <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = PIB)
anova(K) #Rien significatif mais tendance Interaction Brout*Hini et Stress*Hini et tendance Brout et stress
summary(K)

##Modele avec plan en tiroir sans interaction Hini
K <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = PIB)
anova(K) #Rien significatif

##Modele avec plan en tiroir sans Hini
K <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = PIB)
anova(K) #Rien significatif

##Vérifier résidus 
plot(K) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(K) ~ Stress*Brout*Prov, data = PIB) #ok

qqnorm(resid(K)) #Graphique de normalité --> ok
qqline(resid(K))
shapiro.test(resid(K)) #ok

