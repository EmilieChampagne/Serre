
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

Masses2$Hauteurini <- Masses2$Hauteurini/10

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
levels(CET_mean$Brout) <- c("NoBrout","Brout")
levels(CET_mean$Stress) <- c("NoStress","Stress")
color = c(rep("green",4),rep("yellow",4),rep("red",4))
boxplot(Mtot~Brout+Stress+Prov,data=CET_mean,cex.axis=0.5, col = color,las=2, 
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



###QUESTIONS
#Lorsqu'on voit que Prov n'est pas significatif par exemple, est-ce qu'on peut fusionner les provenances
#et prendre en compte seulement effet stress dans un graphique par exemple?
#Même chose pour le broutement, on doit exclure les données des plants broutés ou les fusionner?
boxplot(Mtot~Stress,data=CET)

##EC: Oui, exactement. Pour un graphique, tu peux fusionner les facteurs non-significatifs pour représenter
#uniquement les effets significatifs.
#En fait, même si provenance était significatif, tu pourrais présenter l'effet du stress sans tenir compte
#des provenances.

##Vérifier résidus en incluant les effets aléatoires (effet du design)
plot(A) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône) --> ok
leveneTest(resid(A) ~ Stress*Brout*Prov, data = CET_mean) 



###QUESTIONS:
#On ne peut pas inclure Hauteurini dans le teste de Levene car c'est une variable continue
#Ni prendre en compte l'effet du bloc, 
#Est-ce correct???

##EC: Hum, bonne question. En fait, le test de levene est juste une façon de formaliser ce qu'on voit
#sur le graphique. Dans ce contexte, comme on fait 2 tests, je ne crois pas que ce soit un problème.

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

##EC: attention ici...as-tu testé le modèle suivant "Ratio ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress)"?

##Modele avec plan en tiroir sans Hini
B <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CET_mean)
anova(B)   #Provenance et Stress significatif 
summary(B) #Provenance 2018 et 2080 diffère (stress valeur de p de 0.1 ??)



##QUESTIONS
#On voit que le ratio de 2018 est différent de 2080 mais pas de 2050. Est-ce qu'il y a moyen de voir
#si 2050 est différent de 2080?? Car avec le summary, ca compare toujours avec la valeur de base donc
#provenance 2018...

  ##EC: C'est ce dont je n'étais pas sûre...suite à mes recherches:
CET_mean$Prov <- relevel(CET_mean$Prov, ref="2050")
Test <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CET_mean)
summary(Test) ##Et voilà!

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
levels(CHR_mean$Brout) <- c("NoBrout","Brout")
levels(CHR_mean$Stress) <- c("NoStress","Stress")
color = c(rep("green",4),rep("yellow",4),rep("red",4))
boxplot(Mtot~Brout+Stress+Prov,data=CHR_mean,cex.axis=0.5, col = color,las=2, main = "Masses tot chenes",
        xlab="",ylab="Masse tot")

#Modèle avec plan en tiroir
C <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(C) # Garder Hini seul 

#Sans interaction avec Hini:
C <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(C) #Hauteurini significatif, Intéraction Stress:Provenance, tendance de Stress
summary(C)



###QUESTIONS:
#Lorsqu'on regarde le summary, on ne voit pas que Stress est significativement différent de NoStress
#alors que c'est supposé l'être selon l'anova. Dans le summary, ca compare Stress à NoStress dans les 
#conditons de "base" donc pour NoBrout et Prov2018. Si on regarde le graph, pour 2018 il n'y a pas
#vraiment de différence entre Stress et No Stress pour 2018 alors qu'il y en a beaucoup pour 2050 et 2080.
#D'ailleurs, c'est ce que signifie l'intéraction Stress*Prov significative dans l'anova
#Mais comment sait-t-on s'il y a une différence significative entre Stress et No Stress pour 2050 et 2080?
#Puisque le summary compare toujours avec la variable de "base" donc 2018...
#Aussi les deux p-values des lignes Stress:Prov2050 et Stress:Prov2080 sont non-significatives 
#pourtant dans l'anova on voit que l'interaction Stress*Prov est significative....pourquoi?
#Et on ne peut pas comparer 2050 avec 2080?? Toujours seulement 2050 avec 2018 et 2080 avec 2018...

##EC: OK, plusieurs questions ici. On fera mieux d'en parler demain matin. Tu peux quand même utiliser
#le code ci-haut pour changer le niveau de référence. Ça va répondre à une partie de tes questions.
#Sinon, il faut comprendre que l'ANOVA est un test statistique plus puissant que les contrastes a posteriori.
#Si les différences sont relativement faible, on peut avoir un effet significatif dans l'ANOVA qui ne se réalise
#pas dans le test a posteriori.


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
##EC: Tu as aussi testé + Hini?

##Modele avec plan en tiroir sans Hini
D <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CHR_mean)
anova(D)   #Tendance Stress
summary(D)

boxplot(Ratio~Stress,data=CHR,cex.axis=0.7,las=2,xlab="")

##Vérifier résidus 
plot(D) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(D) ~ Stress*Brout*Prov, data = CHR_mean) #ok

##EC: Tu as testé ta normalité sur B pas sur D...
qqnorm(resid(B)) #Graphique de normalité ok
qqline(resid(B))
shapiro.test(resid(B)) #ok
qqnorm(resid(D)) #
qqline(resid(D))

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
anova(K) #Rien significatif mais tendance Interaction Brout*Hini et Stress Hini et tendance Brout et stress
summary(K)

##Modele avec plan en tiroir sans Hini
K <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = PIB)
anova(K) #Rien significatif

##Vérifier résidus 
plot(K) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(K) ~ Stress*Brout*Prov, data = PIB) #ok

qqnorm(resid(K)) #Graphique de normalité --> ok
qqline(resid(K))
shapiro.test(resid(K)) #ok

