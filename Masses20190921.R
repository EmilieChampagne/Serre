
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
library(swirl)
library(Rmisc)
library(ggplot2)
library(cowplot)


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

CET_mean <- aggregate(list(CET$Mtot,CET$Hauteurini,CET$Ratio,CET$Maerien,CET$Mracine), by= data.frame(CET$Stress, CET$Brout, CET$Prov, CET$Bloc), FUN= "mean")
colnames(CET_mean)[1:9] <- c("Stress", "Brout", "Prov", "Bloc","Mtot","Hauteurini","Ratio","Maerien","Mracine")
CET_mean <- mutate(CET_mean, Stress=factor(Stress, levels=c("NoStress", "Stress2")))#On indique qu'il n'y a pas de niveau "Stress1" pour le facteur Stress
CET <- mutate(CET, Stress=factor(Stress, levels=c("NoStress", "Stress2")))

ddply(CET_mean, c("Stress", "Brout", "Prov", "Bloc","Hauteurini"), summarise,
      N = length(Mtot)) #Verification n = 1 pour chaque combinaison

#Boxplot masses tot 
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

#Représentation effet Stress
meansCET_Mtot <- summarySE(CET_mean, measurevar = "Mtot", groupvars = c("Stress")) 
ggplot(meansCET_Mtot, aes(x = Stress, y = Mtot)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Mtot-se, ymax=Mtot+se), width=.1)+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse totale (g)")
#Mtot Stress2 plus faible

##Vérifier résidus en incluant les effets aléatoires (effet du design)
plot(A) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône) --> ok
leveneTest(resid(A) ~ Stress*Brout*Prov, data = CET_mean) #ok

qqnorm(resid(A)) #Graphique de normalité ok
qqline(resid(A))
shapiro.test(resid(A)) #On vérifie avec un test, faire des transformations au besoin --> ok

####Ratio CERISIERS####

#Boxplot ratio 
boxplot(Ratio~Brout+Stress+Prov,data=CET,cex.axis=0.5, col = color,las=2, main = "Ratio cerisiers", 
        xlab="",ylab="Masse tot")

##Modele avec plan en tiroir
B <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CET_mean)
anova(B) #Interaction Hauteur ini pas significatif 

##Modele avec plan en tiroir sans interaction Hini
B <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov+Hauteurini + (1|Bloc/Stress), data = CET_mean)
anova(B)   #Hauteur ini pas significatif

##Modele avec plan en tiroir sans Hini
B <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CET_mean)
anova(B) #Provenance et Stress significatif 
summary(B)

#Représentation de Stress
meansCET_Ratio <- summarySE(CET_mean, measurevar = "Ratio", groupvars = c("Stress")) 
ggplot(meansCET_Ratio, aes(x = Stress, y = Ratio)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Ratio-se, ymax=Ratio+se), width=.1)+
  xlab("Traitement de stress hydrique") +  
  ylab("Ratio aérien/racinaire")
#Ratio de stress2 plus grand que NoStress

#Représentation Stress et Prov
meansCET_Ratio <- summarySE(CET_mean, measurevar = "Ratio", groupvars = c("Stress", "Prov")) 
ggplot(meansCET_Ratio, aes(x = Stress, y = Ratio, shape=factor(Prov))) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=Ratio-se, ymax=Ratio+se), width=.1, position=position_dodge(width=0.2))+
  scale_x_discrete(labels=c("Sans stress", "Stress")) + xlab("Traitement de stress hydrique") +  
  ylab("Ratio aérien/racinaire")

#Savoir quelle différence de Prov est significative (valeur de reference NoStress 2018)
CET_mean$Prov <- relevel(CET_mean$Prov, ref="2018")
CET_mean$Stress <- relevel(CET_mean$Stress, ref="NoStress")
B1 <- lmerTest::lmer( Ratio ~ Stress*Prov + (1|Bloc/Stress), data = CET_mean)
summary(B1) #Prov 2080 differe de 2018 mais pas 2050-2018

#Savoir si 2050 et 2080 sont significativement différent (on place 2050 comme valeur de base de reference)
CET_mean$Prov <- relevel(CET_mean$Prov, ref="2050")
B1 <- lmerTest::lmer( Ratio ~ Stress*Prov + (1|Bloc/Stress), data = CET_mean)
summary(B1) #Oui, Prov 2050 diffère de 2080

##Vérifier résidus 
plot(B) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône) ok
leveneTest(resid(B) ~ Stress*Brout*Prov, data = CET_mean) #ok

qqnorm(resid(B)) #Graphique de normalité ok
qqline(resid(B))
shapiro.test(resid(B)) #On vérifie avec un test, faire des transformations au besoin --> ok

####MASSE RACINAIRE####

#Boxplot masse racine
color = c(rep("green",4),rep("yellow",4),rep("red",4))
boxplot(Mracine~Brout+Stress+Prov,data=CET,cex.axis=0.5, col = color,las=2, 
        main = "Masses racinaire cerisiers", ylab="Masse racinaire", xlab="") 

##Modele avec plan en tiroir
CET_mean$Brout <- relevel(CET_mean$Brout, ref="NoBrout")
CET_mean$Stress <- relevel(CET_mean$Stress, ref="NoStress")
CET_mean$Prov <- relevel(CET_mean$Prov, ref="2018")
AA <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CET_mean)
anova(AA) #Interaction avec Hauteur ini significatif --> On garde le modèle
#Effet Stress, Brout*Prov, Stress*Brout*Prov, tendance Prov
summary(AA)

#Effet Stress
CET_mean$Prov <- relevel(CET_mean$Prov, ref="2018") #Pas de difference
CET_mean$Prov <- relevel(CET_mean$Prov, ref="2050") #Stress2 plus faible que No stress 2050
CET_mean$Prov <- relevel(CET_mean$Prov, ref="2080") #Pas de difference
AA <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CET_mean)
summary(AA)

#Savoir quelle difference pour Prov
CET_mean$Prov <- relevel(CET_mean$Prov, ref="2018")
CET_mean$Prov <- relevel(CET_mean$Prov, ref="2050")
AA <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CET_mean)
summary(AA) 
#Prov 2050 tendance a etre plus élevé que 2018 et 2080 pas different de 2018
#2080 plus faible que 2050

#Effet Brout?
CET_mean$Prov <- relevel(CET_mean$Prov, ref="2018") #Brout plus eleve que noBrout 2018
CET_mean$Prov <- relevel(CET_mean$Prov, ref="2050") #Brout plus faible que noBrout 2050
CET_mean$Prov <- relevel(CET_mean$Prov, ref="2080") #Pas de difference
AA <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CET_mean)
summary(AA) 

#Effet Brout*Prov
#Différence entre Brout et NoBrout de 2050 est plus petite que 2018
#Différence entre Brout et NoBrout de 2080 est égale avec 2018
#Différence entre Brout et NoBrout de 2080 est plus grande que 2050

#Effet Stress*Brout*Prov
#Stress2*Brout 2050 est plus élevé que 2018
#Stress2*Brout 2080 pas différent de 2018
#Stress2*Brout 2080 est plus faible que 2050
 #2018 Différence entre Stress2 et No Stress est plus faible pour Brout que NoBrout
 #2050 Différence entre Stress2 et No Stress est plus élevée pour Brout que NoBrout
 #2080 Différence entre Stress2 et No Stress n'es pas différente entre Brout et NoBrout

#Representation effet stress et prov
barplot(Mracine~Stress+Prov,data=CET_mean,cex.axis=0.5, las=2, 
        main = "Masses racinaire cerisiers", ylab="Masse racinaire", xlab="") 
#Diminution de la masse racinaire avec le stress
#Masse racinaire plus faible prov2080 que 2018


  #Ce n'est pas aisé de présenter une interaction triple. Il y a plusieurs groupement différents
  #qui peuvent se faire, selon ce que tu as envie de présenter. Je t'en propose 1, mais ça peut
  #être différent:

meansCET_Mracine <- summarySE(CET_mean, measurevar = "Mracine", groupvars = c("Stress", "Brout", "Prov")) 
    #On obtient les moyennes et SE des données brutes
#J'ai copié-collé & adapté un graphique que j'avais déjà fait avec ggplot. Mais on pourrait aussi le
                  #faire un graphique de base R. Ce sont les résultats que tu vas présenter,
                  #je te suggère donc d'explorer et de décider comment tu veux coder ton graphique
                  #(si tu es serrée dans le temps, tu peux aussi prendre l'option facile et reprendre
                  #ce que j'ai fait.)


#Une option est de faire un graphique par provenance
ggplot(meansCET_Mracine[which(meansCET_Mracine$Prov == "2018"),], aes(x = Stress, y = Mracine, shape=factor(Brout))) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Mracine-se, ymax=Mracine+se), width=.1)+
  scale_x_discrete(labels=c("Sans stress", "Stress")) + xlab("Traitement de stress hydrique") +  
  ylab("Masse racinaire (g)")
ggplot(meansCET_Mracine[which(meansCET_Mracine$Prov == "2050"),], aes(x = Stress, y = Mracine, shape=factor(Brout))) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Mracine-se, ymax=Mracine+se), width=.1)+
  scale_x_discrete(labels=c("Sans stress", "Stress")) + xlab("Traitement de stress hydrique") +  
  ylab("Masse racinaire (g)")
ggplot(meansCET_Mracine[which(meansCET_Mracine$Prov == "2080"),], aes(x = Stress, y = Mracine, shape=factor(Brout))) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Mracine-se, ymax=Mracine+se), width=.1)+
  scale_x_discrete(labels=c("Sans stress", "Stress")) + xlab("Traitement de stress hydrique") +  
  ylab("Masse racinaire (g)")

#Une autre option est de faire un graphique pour les broutés et les non broutés. On dirait que c'est
#le traitement qui a le moins d'effet (j'aime bien cette option)
#J'ai mis les graphiques sur la même échelle pour qu'on voit bien les différences
ggplot(meansCET_Mracine[which(meansCET_Mracine$Brout == "NoBrout"),], aes(x = Stress, y = Mracine, shape=factor(Prov))) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=Mracine-se, ymax=Mracine+se), width=.1, position=position_dodge(width=0.2))+
  scale_x_discrete(labels=c("Sans stress", "Stress")) + xlab("Traitement de stress hydrique") +  
  ylab("Masse racinaire (g)")+
  ylim(c(5,35))+
  annotate("text", x = 2, y = 25,label = "Non broutés", size = 5)
ggplot(meansCET_Mracine[which(meansCET_Mracine$Brout == "Brout"),], aes(x = Stress, y = Mracine, shape=factor(Prov))) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=Mracine-se, ymax=Mracine+se), width=.1, position=position_dodge(width=0.2))+
  scale_x_discrete(labels=c("Sans stress", "Stress")) + xlab("Traitement de stress hydrique") +  
  ylab("Masse racinaire (g)")+
  ylim(c(5,35))+
  annotate("text", x = 2, y = 25,label = "Broutés", size = 5)
    #Diviser ainsi ton interaction te permettra de mieux la comprendre
    #Tu as peut-être remarqué, j'ai utilisé "position_dodge" pour qu'on voit mieux les points (ils se superposaient)

    ##EC: Il y a par contre un problème ici...il y a une interaction avec la hauteur initiale.
    #Il va falloir que je prenne conseil de Marie-Claude sur la présentation sauf qu'elle est absente.
    #J'ai une idée de solution, mais entre temps, cette représentation est bonne. N'oublie juste pas
    #de rapporter l'effet de ta covariable quand il y en a une.


##Vérifier résidus 
plot(AA) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône) ok
leveneTest(resid(AA) ~ Stress*Brout*Prov, data = CET_mean) #ok

qqnorm(resid(AA)) #Graphique de normalité ok
qqline(resid(AA))
shapiro.test(resid(AA)) #ok

##MASSE AERIENNE##

#Boxplot masse racine
color = c(rep("green",4),rep("yellow",4),rep("red",4))
boxplot(Maerien~Brout+Stress+Prov,data=CET,cex.axis=0.5, col = color,las=2, 
        main = "Masses aeriennes cerisiers", ylab="Masse aerienne", xlab="") 

##Modele avec plan en tiroir
BB <- lmerTest::lmer( Maerien ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CET_mean)
anova(BB) 
# Interaction stress*brout*Hauteurini de significatif,
# mais avec un F de 5.56 contre 6.86 pour Hini seul --> garder Hini seul 

##Modele avec plan en tiroir sans interaction Hauteurini
CET_mean$Prov <- relevel(CET_mean$Prov, ref="2018")
BB <- lmerTest::lmer( Maerien ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CET_mean)
anova(BB)#Effet Hauteurini et stress
summary(BB) #Masse aerienne plus petite chez Stress2

##Vérifier résidus 
plot(BB) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône) ok
leveneTest(resid(BB) ~ Stress*Brout*Prov, data = CET_mean) #ok

qqnorm(resid(BB)) #Graphique de normalité ok
qqline(resid(BB))
shapiro.test(resid(BB)) # ok


##Conclusion Ratio, aerien et racinaire##
CETStress2<-filter(CET,Stress=="Stress2") 
mean(CETStress2$Maerien) #47
mean(CETStress2$Mracine) #21
mean(CETStress2$Ratio) #3.19

CETNoStress<-filter(CET,Stress=="NoStress")
mean(CETNoStress$Maerien) #34
mean(CETNoStress$Mracine) #11
mean(CETNoStress$Ratio) # 2.41

#Les plants stressés présentent une croissance plus faible aerien et racinaire, mais nettement
#plus faible pour le racinaire donnant un ratio plus élevé que chez les plants non stressés

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

CHR_mean <- aggregate(list(CHR$Mtot,CHR$Hauteurini,CHR$Ratio,CHR$Maerien,CHR$Mracine), by= data.frame(CHR$Stress, CHR$Brout, CHR$Prov, CHR$Bloc), FUN= "mean")
colnames(CHR_mean)[1:9] <- c("Stress", "Brout", "Prov", "Bloc","Mtot","Hauteurini","Ratio","Maerien","Mracine")
CHR_mean <- mutate(CHR_mean, Stress=factor(Stress, levels=c("NoStress", "Stress2")))#On indique qu'il n'y a pas de niveau "Stress1" pour le facteur Stress
CHR <- mutate(CHR, Stress=factor(Stress, levels=c("NoStress", "Stress2")))

#Verification: nous avons bien un n = 1 pour chaque combinaison.
ddply(CHR_mean, c("Stress", "Brout", "Prov", "Bloc","Hauteurini"), summarise,
      N = length(Mtot))

#Boxplot masses tot 
color = c(rep("green",4),rep("yellow",4),rep("red",4))
boxplot(Mtot~Brout+Stress+Prov,data=CHR,cex.axis=0.5, col = color,las=2, main = "Masses tot chenes",
        xlab="",ylab="Masse tot")

#Modèle avec plan en tiroir
C <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(C) # Garder Hini seul 

#Sans interaction avec Hini:
CHR_mean$Prov <- relevel(CHR_mean$Prov, ref="2018")
CHR_mean$Brout <- relevel(CHR_mean$Brout, ref="NoBrout")
C <- lmerTest::lmer(Mtot ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(C) #Hauteurini significatif, Intéraction Stress:Provenance, tendance de Stress
summary(C)#Difference entre Stress2-NoStress plus grande chez Prov2050 que 2018

#Avec Prov2050 comme valeur de base de reference
CHR_mean$Prov <- relevel(CHR_mean$Prov, ref="2050")
C2 <- lmerTest::lmer(Mtot ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
summary(C2) #Difference entre Stress2-NoStress plus grande chez Prov2050 que 2080

#Representation tendance de l'effet stress
boxplot(Mtot~Stress,data=CHR,cex.axis=0.7,col="green",las=2,xlab="")

#Representation interaction Prov*Stress
color2 = c(rep("green",2),rep("yellow",2),rep("red",2))
boxplot(Mtot~Stress+Prov,data=CHR,col=color2,cex.axis=0.7,las=2,xlab="") #outlier 

#Enlever le outlier
CHR1 <- CHR[-29, ]
CHR_mean1 <- aggregate(list(CHR1$Mtot,CHR1$Hauteurini,CHR1$Ratio), by= data.frame(CHR1$Stress, CHR1$Brout, CHR1$Prov, CHR1$Bloc), FUN= "mean")
colnames(CHR_mean1)[1:7] <- c("Stress", "Brout", "Prov", "Bloc","Mtot","Hauteurini","Ratio")
CHR_mean1 <- mutate(CHR_mean1, Stress=factor(Stress, levels=c("NoStress", "Stress2")))#On indique qu'il n'y a pas de niveau "Stress1" pour le facteur Stress
CHR1 <- mutate(CHR1, Stress=factor(Stress, levels=c("NoStress", "Stress2")))

#Modèle avec plan en tiroir sans le outlier
C <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CHR_mean1)
anova(C) # Garder Hini seul 

#Sans interaction avec Hini (sans outlier)
CHR_mean1$Prov <- relevel(CHR_mean1$Prov, ref="2018")
C <- lmerTest::lmer(Mtot ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean1)
anova(C) #Hauteurini significatif, Intéraction Stress:Provenance, tendance de Stress
summary(C) #Pas d'interaction Stress2*Prov singnificatif (en comparaison avec 2018)

#Avec Prov2050 comme valeur de base de reference (sans outlier)
CHR_mean1$Prov <- relevel(CHR_mean1$Prov, ref="2050")
C2 <- lmerTest::lmer(Mtot ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean1)
summary(C2) #Stress2:Prov2080 significatif (comparaison avec 2050), Stress significatif

#Representation tendance de l'effet stress (sans outlier)
boxplot(Mtot~Stress,data=CHR1,col="green",cex.axis=0.7,las=2,xlab="")

#Representation interaction Prov*Stress (sans outlier)
color2 = c(rep("green",2),rep("yellow",2),rep("red",2))
boxplot(Mtot~Stress+Prov,data=CHR1,col=color2,cex.axis=0.7,las=2,xlab="") 
#Conclusion: Difference entre Stress2 et NoStress est plus grande chez 2050 que 2080
#Surtout Prov2050 qui contribue à la tendance de l'effet stress...

##Vérifier résidus en incluant les effets aléatoires (effet du design)
plot(C) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône) --> ok
leveneTest(resid(C) ~ Stress*Brout*Prov, data = CHR_mean1) #ok

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

##MASSE RACINAIRE##

#Boxplot masse racine
color = c(rep("green",4),rep("yellow",4),rep("red",4))
boxplot(Mracine~Brout+Stress+Prov,data=CHR,cex.axis=0.5, col = color,las=2, 
        main = "Masses racinaire chenes", ylab="Masse racinaire", xlab="") 

##Modele avec plan en tiroir
CC <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(CC) #Hauteur ini significatif --> On garde Hiniseul

##Modele avec plan en tiroir sans interaction Hini
CC <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(CC) #Effet Stress et Stress*Prov

#Effet stress
CC <- lmerTest::lmer( Mracine ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
summary(CC) #Mracine de Stress2 est plus faible que NoStress

#Savoir quelle difference pour Stress*Prov
CHR_mean$Prov <- relevel(CHR_mean$Prov, ref="2018")
CC <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
summary(CC) #Stress2*Prov2050 significatif, pas 2080 (comparaison avec 2018)

#2050 comme reference
CHR_mean$Prov <- relevel(CHR_mean$Prov, ref="2050")
CC <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
summary(CC) #Stress2*Prov2080 significatif(comparaison avec 2050)

#Representation effet stress et prov
boxplot(Mracine~Stress+Prov,data=CHR,cex.axis=0.5, las=2, 
        main = "Masses racinaire chenes", ylab="Masse racinaire", xlab="") 

#Difference entre Stress2 et No Stress est plus importante chez Prov 2050 (meme pattern que pour Mtot)

##Vérifier résidus 
plot(CC) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(CC) ~ Stress*Brout*Prov, data = CHR_mean) #ok

qqnorm(resid(CC)) #Graphique de normalité ok
qqline(resid(CC))
shapiro.test(resid(CC)) #Essayer une transformation log


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


##MASSE AERIENNE##

#Boxplot masse racine
boxplot(Maerien~Brout+Stress+Prov,data=CHR,cex.axis=0.5, col = color,las=2, 
        main = "Masses aeriennes chenes", ylab="Masse aerienne", xlab="") 

##Modele avec plan en tiroir
DD <- lmerTest::lmer( Maerien ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(DD) #pas interaction avec Hini

##Modele avec plan en tiroir sans interaction Hauteurini
CHR_mean$Prov <- relevel(CHR_mean$Prov, ref="2018")
DD <- lmerTest::lmer(Maerien ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(DD)#Tendance Prov et Tendance Stress:Prov
summary(DD)#Prov 2050 significativement plus elevé et Stress2*Prov2050 significatif

#Reference Prov2080
CHR_mean$Prov <- relevel(CHR_mean$Prov, ref="2080")
DD <- lmerTest::lmer(Maerien ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(DD)#Tendance Prov et Tendance Stress:Prov
summary(DD)#Prov 2050 significativement plus élevé et tendance Stress2*Prov2050

#Prov 2050 a tendance a avoir Mracinaire plus elevé et une plus grande difference entre Stress2 et
#NoStress que les autres provenances

#Representation tendance Stress*Prov
boxplot(Maerien~Stress+Prov,data=CHR,cex.axis=0.5, col = color2,las=2, 
        main = "Masses aeriennes chenes", ylab="Masse aerienne", xlab="") 

##Conclusion Ratio, aerien et racinaire##
CHRStress2<-filter(CHR,Stress=="Stress2") 
mean(CHRStress2$Maerien) #20
mean(CHRStress2$Mracine) #18
mean(CHRStress2$Ratio) #1.21

CHRNoStress<-filter(CHR,Stress=="NoStress")
mean(CHRNoStress$Maerien) #23
mean(CHRNoStress$Mracine) #24
mean(CHRNoStress$Ratio) # 1.00

#Croissance aussi rapide partie aerienne pour Stress2 que NoStress (pas significatif)
#Croissance ralentie partie racinaire (surtout Prov2050) = Ratio plus élevé pour Stress2

##Vérifier résidus 
plot(DD) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(DD) ~ Stress*Brout*Prov, data = CHR_mean) #ok

qqnorm(resid(DD)) #Graphique de normalité ok
qqline(resid(DD))
shapiro.test(resid(DD)) #ok

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
color = c(rep("green",6))
boxplot(Mtot~Brout+Stress,data=ERS,cex.axis=0.5, col = color,las=2, main = "Masses tot Erables", 
        xlab="",ylab="Masse tot")

#Modèle avec plan en tiroir
E <- lmerTest::lmer( Mtot ~ Stress*Brout*Hauteurini + (1|Bloc/Stress), data = ERS)
anova(E) # Rien de significatif

#Sans interaction avec Hini:
E <- lmerTest::lmer( Mtot ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = ERS)
anova(E) #Hauteurini tendance a etre significatif, Interaction Stress*Brout significatif
summary(E)

#Si on retire Hini
ERS$Brout <- relevel(ERS$Brout, ref="NoBrout")
ERS$Stress <- relevel(ERS$Stress, ref="NoStress")
E <- lmerTest::lmer( Mtot ~ Stress*Brout + (1|Bloc/Stress), data = ERS)
anova(E) #Interaction Stress*Brout significatif
summary(E) #Brout*Stress1 significatif et tendance Brout*Stress2 (comparaison avec NoStress)

#Avec Stress1 comme valeur de base de reference
ERS$Stress <- relevel(ERS$Stress, ref="Stress1")
E2 <- lmerTest::lmer(Mtot ~ Stress*Brout + (1|Bloc/Stress), data = ERS)
summary(E2) #Tendance Brout*Stress2 (comparaison avec stress1)

##Conclusion: Différence entre Stress1 et NoStress est plus grande chez Brout que NoBrout
#Mais la tendance est inversée --> Mtot diminue avec Stress pour NonBrout mais augmente avec Brout

##Vérifier résidus en incluant les effets aléatoires (effet du design)
plot(E) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône) --> ok
leveneTest(resid(E) ~ Stress*Brout, data = ERS) #ok

qqnorm(resid(E)) #Graphique de normalité ok
qqline(resid(E))
shapiro.test(resid(E)) #ok

####Ratio ERABLES####

#Boxplot ratio 
boxplot(Ratio~Brout+Stress,data=ERS,cex.axis=0.5,las=2, main = "Ratio Erables", 
        xlab="",ylab="Masse tot")

##Modele avec plan en tiroir
G <- lmerTest::lmer( Ratio ~ Stress*Brout*Hauteurini + (1|Bloc/Stress), data = ERS)
anova(G) #Rien significatif, tendance de Brout et Interaction Brout:Hauteurini
##Si on garde ce modele en raison de la tendance de l'interaction avec Hauteurini,
#on voit une tendance de ratio plus petit pour les plants broutés 

##Modele avec plan en tiroir sans interaction Hini
G <- lmerTest::lmer( Ratio ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = ERS)
anova(G) #Rien significatif

##Modele avec plan en tiroir sans Hini
G <- lmerTest::lmer( Ratio ~ Stress*Brout + (1|Bloc/Stress), data = ERS)
anova(G) #Rien significatif

#Visualiser tendance de l'effet Brout 
boxplot(Ratio~Brout,data=ERS,cex.axis=0.5, col = color,las=2, main = "Ratio Erables", 
        xlab="",ylab="Ratio") #Outlier, essayer de l'enlever

#Enlever le outlier
ERS1 <- ERS[-13, ]

#Modèle avec plan en tiroir sans le outlier
G <- lmerTest::lmer( Mtot ~ Stress*Brout*Hauteurini + (1|Bloc/Stress), data = ERS1)
anova(G) #Rien significatif 

##Modele avec plan en tiroir sans interaction Hini
G <- lmerTest::lmer( Ratio ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = ERS1)
anova(G) #Rien significatif

##Modele avec plan en tiroir sans Hini
G <- lmerTest::lmer( Ratio ~ Stress*Brout + (1|Bloc/Stress), data = ERS1)
anova(G) #Rien significatif

##Vérifier résidus 
plot(G) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(G) ~ Stress*Brout, data = ERS) #ok

qqnorm(resid(G)) #Graphique de normalité ok
qqline(resid(G))
shapiro.test(resid(G)) #ok


##MASSE RACINAIRE##

#Boxplot masse racine
boxplot(Mracine~Brout+Stress,data=ERS,cex.axis=0.5, col = color,las=2, 
        main = "Masses racinaire erables", ylab="Masse racinaire", xlab="") 

##Modele avec plan en tiroir
EE <- lmerTest::lmer( Mracine ~ Stress*Brout*Hauteurini + (1|Bloc/Stress), data = ERS)
anova(EE) #Rien significatif

##Modele avec plan en tiroir sans interaction Hini
EE <- lmerTest::lmer( Mracine ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = ERS)
anova(EE) #Hini non significatif

##Modele avec plan en tiroir sans Hini
ERS$Brout <- relevel(ERS$Brout, ref="NoBrout")
ERS$Stress <- relevel(ERS$Stress, ref="NoStress")
EE <- lmerTest::lmer( Mracine ~ Stress*Brout + (1|Bloc/Stress), data = ERS)
anova(EE) #Stress*Brout significatif
summary(EE) #Stress1*Brout significatif (comparaison avec NoStress)

##Stress1 en reference
ERS$Brout <- relevel(ERS$Brout, ref="NoBrout")
ERS$Stress <- relevel(ERS$Stress, ref="Stress1")
EE <- lmerTest::lmer( Mracine ~ Stress*Brout + (1|Bloc/Stress), data = ERS)
summary(EE) #Stress2*Brout significatif (Comparaison avec Stress1)

##Vérifier résidus 
plot(EE) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(EE) ~ Stress*Brout*Prov, data = ERS) #ok

qqnorm(resid(EE)) #Graphique de normalité ok
qqline(resid(EE))
shapiro.test(resid(EE)) #ok

  
##MASSE AERIENNE##
  
#Boxplot masse racine
boxplot(Maerien~Brout+Stress,data=ERS,cex.axis=0.5, col = color,las=2, 
          main = "Masses aeriennes erables", ylab="Masse aerienne", xlab="") 

##Modele avec plan en tiroir
FF <- lmerTest::lmer( Maerien ~ Stress*Brout*Hauteurini + (1|Bloc/Stress), data = ERS)
anova(FF) #pas interaction avec Hini

##Modele sans interaction Hini
FF <- lmerTest::lmer( Maerien ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = ERS)
anova(FF) #Hauteurini significatif, effet Stress*Brout

##Reference NoStress
ERS$Brout <- relevel(ERS$Brout, ref="NoBrout")
ERS$Stress <- relevel(ERS$Stress, ref="NoStress")
ERS$Prov <- relevel(ERS$Prov, ref="2018")
FF <- lmerTest::lmer(Maerien ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = ERS)
summary(FF) #Stress1*Brout significatif et Stress2*Brout (comparaison avec NoStress)

#Reference Stress1
ERS$Stress <- relevel(ERS$Stress, ref="Stress1")
FF <- lmerTest::lmer(Maerien ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = ERS)
summary(FF) #Stress2*Brout pas significativement different de Stress1

##Conclusion Ratio, aerien et racinaire##
ERSNoBrout<-filter(ERS,Brout=="NoBrout") 
mean(ERSNoBrout$Maerien) #26
mean(ERSNoBrout$Mracine) #17
mean(ERSNoBrout$Ratio) #1.69

ERSBrout<-filter(ERS,Brout=="Brout")
mean(ERSBrout$Maerien) #29
mean(ERSBrout$Mracine) #21
mean(ERSBrout$Ratio) # 1.5

# Pas d'effet significatif de Brout sur Maerien ni Mracinaire
# Tendance que le ratio est plus faible pour les plants broutés non-expliquée...

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
color = c(rep("green",6),rep("yellow",6),rep("red",6))
boxplot(Mtot~Brout+Stress+Prov,data=THO,cex.axis=0.5, col = color,las=2, main = "Masses tot Thuyas",
        xlab="", ylab="Masse tot")

#Modèle avec plan en tiroir
H <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = THO)
anova(H)
# Intéraction stress*brout*Hauteurini et stress*brout de significatifs,
# mais avec des F de 3.24 et 3.44 contre 21.26 pour Hini seul --> garder Hini seul 

#On refait sans interaction avec Hini
H <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = THO)
anova(H) #Hauteurini, Stress et Brout significatifs

#Regarder quel niveau de stress differe
THO$Stress <- relevel(THO$Stress, ref="NoStress")
THO$Brout <- relevel(THO$Brout, ref="NoBrout")
H <- lmerTest::lmer( Mtot ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = THO)
summary(H) #Tendance de Stress1 et Stress2 significatif (comparaison avec NoStress)

#Stress1 comme reference
THO$Stress <- relevel(THO$Stress, ref="Stress1")
H <- lmerTest::lmer( Mtot ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = THO)
summary(H) #Stress2 significatif (comparaison avec Stress1)

#Representation effet Stress
boxplot(Mtot~Stress,data=THO,cex.axis=0.5, col = color,las=2, main = "Masses tot Thuyas",
        xlab="", ylab="Masse tot")

#Representation effet Brout
boxplot(Mtot~Brout,data=THO,cex.axis=0.5, col = color,las=2, main = "Masses tot Thuyas",
        xlab="", ylab="Masse tot")

#Representation effet stress et brout
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

#Quelle difference tendance stress
THO$Stress <- relevel(THO$Stress, ref="NoStress")
I1 <- lmerTest::lmer( Ratio ~ Stress + (1|Bloc/Stress), data = THO)
summary(I1) #Stress1 different et Tendance Stress2
THO$Stress <- relevel(THO$Stress, ref="Stress1")
I1 <- lmerTest::lmer( Ratio ~ Stress + (1|Bloc/Stress), data = THO)
summary(I1) #Stress2 pas different de Stress1

#Representation effet Stress
boxplot(Ratio~Stress,data=THO,cex.axis=0.5, col = color,las=2, main = "Masses tot Thuyas",
        xlab="", ylab="Masse tot")
#Tendance que la ratio augmente avec stress

##Analyse interaction Stress*Brout*Prov

#Avec Prov2018 comme valeur de base de reference
THO$Stress <- relevel(THO$Stress, ref="NoStress")
THO$Prov <- relevel(THO$Prov, ref="2018")
THO$Brout <- relevel(THO$Brout, ref="NoBrout")
I2 <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = THO)
summary(I2) #Stress2*Brout*Prov2080 significatif (comparaison avec 2018)
#Stress2*Brout significatif pour 2018

#Avec Prov2080 comme valeur de base de reference
THO$Prov <- relevel(THO$Prov, ref="2080")
I2 <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = THO)
summary(I2) #Rien significatif

#Representation Interaction Stress*Brout Prov 2018
THO2018<-filter(Masses2,Esp=="THO",Prov=="2018")
boxplot(Ratio~Brout+Stress,data=THO2018,cex.axis=0.5, col = "green",las=2, main = "Ratio thuyas 2018", 
        xlab="",ylab="Masse tot")

#Representation Interaction Stress*Brout Prov 2080
THO2080<-filter(Masses2,Esp=="THO",Prov=="2080")
boxplot(Ratio~Brout+Stress,data=THO2080,cex.axis=0.5, col = "red",las=2, main = "Ratio thuyas 2080", 
        xlab="",ylab="Masse tot")

##Interpretation de l'interaction Stress*Brout*Prov
#Effet Stress2*Brout significativement différent entre 2018 et 2080 =
#Différence entre Stress2 et NoStress est plus grande pour NoBrout que Brout pour Prov2018
#Alors que cette différence n'est pas significative pour Prov2080

##Vérifier résidus 
plot(I) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(I) ~ Stress*Brout*Prov, data = THO) #ok

qqnorm(resid(I)) #Graphique de normalité ok
qqline(resid(I))
shapiro.test(resid(I)) #ok

##MASSE RACINAIRE##

#Boxplot masse racine
boxplot(Mracine~Brout+Stress+Prov,data=THO,cex.axis=0.5, col = color,las=2, 
        main = "Masses racinaire thuyas", ylab="Masse racinaire", xlab="") 

##Modele avec plan en tiroir
II <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = THO)
anova(II) #Hauteur ini significatif --> On garde Hiniseul

##Modele avec plan en tiroir sans interaction Hini
II <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = THO)
anova(II) #Effet Stress et Brout

#Effet Stress et Brout
THO$Stress <- relevel(THO$Stress, ref="NoStress")
II <- lmerTest::lmer( Mracine ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = THO)
summary(II) #NoBrout plus élevé que Brout

#Stress1 comme reference
THO$Stress <- relevel(THO$Stress, ref="Stress1")
II <- lmerTest::lmer( Mracine ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = THO)
summary(II) #Aucune combinaison de stress n'est significative, on ne peut pas tirer de conclusion

#Representation effet stress et brout
boxplot(Mracine~Stress+Brout,data=THO,cex.axis=0.5, las=2, 
        main = "Masses racinaire thuyas", ylab="Masse racinaire", xlab="") 

##Vérifier résidus 
plot(II) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(II) ~ Stress*Brout*Prov, data = THO) #ok

qqnorm(resid(II)) #Graphique de normalité ok
qqline(resid(II))
shapiro.test(resid(II)) #Non Transformation à faire

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TRANSFORMATION A FAIRE
  
##MASSE AERIENNE##

#Boxplot masse aerien
boxplot(Maerien~Brout+Stress+Prov,data=THO,cex.axis=0.5, col = color,las=2, 
          main = "Masses aeriennes thuyas", ylab="Masse aerienne", xlab="") 

##Modele avec plan en tiroir
JJ <- lmerTest::lmer( Maerien ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = THO)
anova(JJ) #On garde Hiniseul

##Modele avec plan en tiroir sans interaction Hini
THO$Stress <- relevel(THO$Stress, ref="NoStress")
JJ <- lmerTest::lmer( Maerien ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = THO)
anova(JJ) #Effet Stress et Brout

#Effet Brout et Stress
JJ <- lmerTest::lmer( Maerien ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = THO)
summary(JJ) #NoBrout plus élevé que Brout
#Stress2 tendance a être plus faible que NoStress

#Stress1 comme reference
THO$Stress <- relevel(THO$Stress, ref="Stress1")
JJ <- lmerTest::lmer( Maerien ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = THO)
summary(JJ) #Stress2 pas significatif

#Representation effet stress et brout
boxplot(Maerien~Stress+Brout,data=THO,cex.axis=0.5, las=2, 
        main = "Masses aeriennes thuyas", ylab="Masse aerienne", xlab="") 

##Vérifier résidus 
plot(JJ) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(JJ) ~ Stress*Brout*Prov, data = THO) #ok

qqnorm(resid(JJ)) #Graphique de normalité ok
qqline(resid(JJ))
shapiro.test(resid(JJ)) #ok

##Conclusion Ratio, aerien et racinaire##
THOStress2<-filter(THO,Stress=="Stress2",Brout=="NoBrout") 
mean(THOStress2$Maerien) #32
mean(THOStress2$Mracine) #8
mean(THOStress2$Ratio) #4.1

THOStress1<-filter(THO,Stress=="Stress1",Brout=="NoBrout")  
mean(THOStress1$Maerien) #43
mean(THOStress1$Mracine) #13
mean(THOStress1$Ratio) #3.8

THONoStress<-filter(THO,Stress=="NoStress",Brout=="NoBrout") 
mean(THONoStress$Maerien) #50
mean(THONoStress$Mracine) #17
mean(THONoStress$Ratio) #3.2

THOBrout<-filter(THO,Stress=="NoStress",Brout=="Brout") 
mean(THOBrout$Maerien) #31
mean(THOBrout$Mracine) #10
mean(THOBrout$Ratio) #3.5

#MTot diminuée par le broutement et le stress
#Ratio ne change pas ni pour le stress (tendance à augmenter), ni pourle broutement
#Alors que Maerien et racinaire est diminuée par le stress et le broutement
# = Maerien est davantage ralentie avec stress que racinaire
# = Maerien diminué par le broutement ET Mracinaire pour atteindre un meme ratio que non-broute

##Vérifier résidus 
plot(JJ) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(JJ) ~ Stress*Brout*Prov, data = THO) #ok

qqnorm(resid(JJ)) #Graphique de normalité ok
qqline(resid(JJ))
shapiro.test(resid(JJ)) #ok


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

#Representation effet Brout
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


##MASSE RACINAIRE##

#Boxplot masse racine
boxplot(Mracine~Brout+Stress+Prov,data=PIB,cex.axis=0.5, col = color,las=2, 
        main = "Masses racinaire pins", ylab="Masse racinaire", xlab="") 

##Modele avec plan en tiroir
KK <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = PIB)
anova(KK) #Rien significatif

##Modele avec plan en tiroir sans interaction Hini
KK <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = PIB)
anova(KK) #Hauteurini significatif, Effet Brout

#Effet Brout
PIB$Brout <- relevel(PIB$Brout, ref="NoBrout")
KK <- lmerTest::lmer( Mracine ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = PIB)
summary(KK) #NoBrout plus élevé que Brout

#Representation effet brout
boxplot(Mracine~Brout,data=PIB,cex.axis=0.5, las=2, 
        main = "Masses racinaire pins", ylab="Masse racinaire", xlab="") 

##Vérifier résidus 
plot(KK) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(KK) ~ Stress*Brout*Prov, data = PIB) #ok

qqnorm(resid(KK)) #Graphique de normalité ok
qqline(resid(KK))
shapiro.test(resid(KK)) #Non transformation à faire

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Transformation à faire
  
##MASSE AERIENNE##
  
#Boxplot masse aerien
boxplot(Maerien~Brout+Stress+Prov,data=PIB,cex.axis=0.5, col = color,las=2, 
          main = "Masses aeriennes pins", ylab="Masse aerienne", xlab="") 

##Modele avec plan en tiroir
LL <- lmerTest::lmer( Maerien ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = PIB)
anova(LL) #Rien significatif

##Modele avec plan en tiroir sans interaction Hini
LL <- lmerTest::lmer( Maerien ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = PIB)
anova(LL) #Hauteurini significatif, Effet Brout

#Effet Brout
PIB$Brout <- relevel(PIB$Brout, ref="NoBrout")
LL <- lmerTest::lmer( Maerien ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = PIB)
summary(LL) #NoBrout plus élevé que Brout

#Representation effet Brout
boxplot(Maerien~Brout,data=PIB,cex.axis=0.5, las=2, 
        main = "Masses aeriennes thuyas", ylab="Masse aerienne", xlab="") 

##Vérifier résidus 
plot(LL) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(LL) ~ Stress*Brout*Prov, data = PIB) #ok

qqnorm(resid(LL)) #Graphique de normalité ok
qqline(resid(LL))
shapiro.test(resid(LL)) #ok

##Conclusion Ratio, aerien et racinaire##

PIBNoStress<-filter(PIB,Stress=="NoStress",Brout=="NoBrout") 
mean(PIBNoStress$Maerien) #21
mean(PIBNoStress$Mracine) #11
mean(PIBNoStress$Ratio) #1.9

PIBBrout<-filter(PIB,Stress=="NoStress",Brout=="Brout") 
mean(PIBBrout$Maerien) #8
mean(PIBBrout$Mracine) #5
mean(PIBBrout$Ratio) #1.6

#MTot diminuée par le broutement
#Ratio ne change pas
#Alors que Maerien et racinaire sont diminuées par le broutement

##Vérifier résidus 
plot(LL) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(LL) ~ Stress*Brout*Prov, data = PIB) #ok

qqnorm(resid(LL)) #Graphique de normalité ok
qqline(resid(LL))
shapiro.test(resid(LL)) #ok

