
# Septembre-Octobre 2019
# Turgeon Roxanne 
# Initiation a la recherche / Migration assistee / Broutement / Stress_hydrique / Serre

#EC: Ce commentaire est ajouté pour faire un exemple de l'utilisation de GitHub...

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
library(Rmisc)
library(ggplot2)
library(cowplot)

Masses3 <- read.table("./Masses3.txt", header=TRUE, na.string = "", dec = ",") 

Masses3$Brout <- replace(Masses3$Brout,Masses3$Brout=="0","NoBrout")
Masses3$Brout <- replace(Masses3$Brout,Masses3$Brout=="1","Brout")
Masses3$Stress <- replace(Masses3$Stress,Masses3$Stress=="0","NoStress")
Masses3$Stress <- replace(Masses3$Stress,Masses3$Stress=="1","Stress1")
Masses3$Stress <- replace(Masses3$Stress,Masses3$Stress=="2","Stress2")

Masses3$Prov <- as.factor(as.character(Masses3$Prov))
Masses3$Bloc <- as.factor(as.character(Masses3$Bloc))
Masses3$Brout <- as.factor(as.character(Masses3$Brout))
Masses3$Stress <- as.factor(as.character(Masses3$Stress))

Masses3$Hauteurini <- Masses3$Hauteurini

str(Masses3)

CET<-filter(Masses3,Esp=="CET")
CHR<-filter(Masses3,Esp=="CHR")
ERS<-filter(Masses3,Esp=="ERS")
THO<-filter(Masses3,Esp=="THO")
PIB<-filter(Masses3,Esp=="PIB")

Masses3Temoin<-filter(Masses3,Brout=="NoBrout",Stress=="NoStress",Prov=="2018")

#################TOUTES ESPECES###################

boxplot(Mtot~Brout+Esp,data=Masses3,cex.axis=0.5, las=2, ylab="Masse tot", xlab="")
boxplot(Mtot~Stress+Esp,data=Masses3,cex.axis=0.5, las=2, ylab="Masse tot", xlab="")
boxplot(Mtot~Bloc*Esp,data=Masses3, cex.axis=0.5, las=2, ylab="Masse tot", xlab="")

##Ratio Temoin##
boxplot(Ratio~Esp, data=Masses3Temoin, ylab="Ratio", xlab="Espèces")
meansMasses3_Temoin <- summarySE(Masses3Temoin, measurevar = "Ratio", groupvars = c("Esp")) 
ggplot(meansMasses3_Temoin, aes(x = Esp, y = Ratio)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Ratio-se, ymax=Ratio+se), width=.1)+
  xlab("Espèces") +  
  ylab("Ratio aérien/racinaire")


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
meansCET_Mtot <- mutate(meansCET_Mtot, Stress=factor(Stress, levels=c("NoStress", "Stress2")))
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

#summary
ls_means(B, which = "Prov", pairwise = TRUE) #Le test pour toutes les comparaisons, sans changer le niveau de référence
#Prov 2050-2080 et 2018-2080 différents
ls_means(B, which = "Prov")#Tu pourrais utiliser les valeurs d'estimés ici pour la représentation graphique
ls_means(B, which = "Stress")#Pour aller voir les valeurs du modèle et leur intervalle de confiance

#Représentation de Stress
meansCET_Ratio <- summarySE(CET_mean, measurevar = "Ratio", groupvars = c("Stress")) 
ggplot(meansCET_Ratio, aes(x = Stress, y = Ratio)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Ratio-se, ymax=Ratio+se), width=.1)+
  xlab("Traitement de stress hydrique") +  
  ylab("Ratio aérien/racinaire")
#Ratio de stress2 plus grand que NoStress

#Représentation de Prov
meansCET_Ratio <- summarySE(CET_mean, measurevar = "Ratio", groupvars = c("Prov")) 
ggplot(meansCET_Ratio, aes(x = Prov, y = Ratio)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Ratio-se, ymax=Ratio+se), width=.1)+
  xlab("Provenance") +  
  ylab("Ratio aérien/racinaire")

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
CET_mean$Brout <- relevel(CET_mean$Brout, ref="Brout")
CET_mean$Stress <- relevel(CET_mean$Stress, ref="Stress2")
CET_mean$Prov <- relevel(CET_mean$Prov, ref="2018")
AA <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CET_mean)
anova(AA) #Interaction avec Hini significatif --> On garde le modèle
#Effet Stress, Brout*Prov, Stress*Brout*Prov, tendance Prov
#Lorsqu'il y a une intéraction significative on laisse tomber les effets simples
#Ici on regarde seulement Stress*Brout*Prov

CET_mean$Prov <- relevel(CET_mean$Prov, ref="2018")
AA <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CET_mean)
summary(AA)
CET_mean$Prov <- relevel(CET_mean$Prov, ref="2050")
AA <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CET_mean)
summary(AA)
CET_mean$Prov <- relevel(CET_mean$Prov, ref="2080") 
AA <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CET_mean)
summary(AA)

#Effet Stress*Brout
#2018 Différence entre NoStress et Stress2 est plus faible pour NoBrout que Brout
#2050 Différence entre NoStress et Stress2 est plus élevée pour NoBrout que Brout
#2080 Différence entre NoStress et Stress2 n'est pas différente entre NoBrout et Brout

#Effet Stress*Brout*Prov
#Stress2*Brout 2050 est plus élevé que 2018
#Stress2*Brout 2080 pas différent de 2018 = La différence de 2018 est peu significative (0.04)
#Stress2*Brout 2080 est plus faible que 2050

#Effet Stress*Brout plus important chez 2050

#Representation effet Stress*Brout*Prov
meansCET_Mracine <- summarySE(CET_mean, measurevar = "Mracine", groupvars = c("Stress", "Brout", "Prov")) 
#Une option est de faire un graphique pour les broutés et les non broutés. 
#On dirait que c'est le traitement qui a le moins d'effet
ggplot(meansCET_Mracine[which(meansCET_Mracine$Brout == "NoBrout"),], aes(x = Stress, y = Mracine, shape=factor(Prov))) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=Mracine-se, ymax=Mracine+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse racinaire (g)")+
  ylim(c(5,35))+
  annotate("text", x = 2, y = 25,label = "Non broutés", size = 5)
ggplot(meansCET_Mracine[which(meansCET_Mracine$Brout == "Brout"),], aes(x = Stress, y = Mracine, shape=factor(Prov))) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=Mracine-se, ymax=Mracine+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse racinaire (g)")+
  ylim(c(5,35))+
  annotate("text", x = 2, y = 25,label = "Broutés", size = 5)
   
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

#Boxplot masse aerienne
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
anova(BB) #Effet Hauteurini et stress
summary(BB) #Masse aerienne plus petite chez Stress2

#Représentation de Stress
meansCET_Maerien <- summarySE(CET_mean, measurevar = "Maerien", groupvars = c("Stress")) 
meansCET_Maerien <- mutate(meansCET_Maerien, Stress=factor(Stress, levels=c("NoStress", "Stress2")))
ggplot(meansCET_Maerien, aes(x = Stress, y = Maerien)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Maerien-se, ymax=Maerien+se), width=.1)+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse aérienne (g)")
#Maerien plus petit chez Stress2

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
#On regarde seulement Stress*Prov
summary(C)#Difference entre Stress2-NoStress plus grande chez Prov2050 que 2018

#Enlever le outlier
CHR1 <- CHR[-29, ]
CHR_mean1 <- aggregate(list(CHR1$Mtot,CHR1$Hauteurini,CHR1$Ratio,CHR1$Mracine,CHR1$Maerien), by= data.frame(CHR1$Stress, CHR1$Brout, CHR1$Prov, CHR1$Bloc), FUN= "mean")
colnames(CHR_mean1)[1:9] <- c("Stress", "Brout", "Prov", "Bloc","Mtot","Hauteurini","Ratio","Mracine","Maerien")
CHR_mean1 <- mutate(CHR_mean1, Stress=factor(Stress, levels=c("NoStress", "Stress2")))#On indique qu'il n'y a pas de niveau "Stress1" pour le facteur Stress
CHR1 <- mutate(CHR1, Stress=factor(Stress, levels=c("NoStress", "Stress2")))

#Modèle avec plan en tiroir sans le outlier
C <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CHR_mean1)
anova(C) # Garder Hini seul 

#Sans interaction avec Hini (sans outlier)
CHR_mean1$Prov <- relevel(CHR_mean1$Prov, ref="2018")
CHR_mean1$Brout <- relevel(CHR_mean1$Brout, ref="NoBrout")
CHR_mean1$Stress <- relevel(CHR_mean1$Stress, ref="Stress2")
C <- lmerTest::lmer(Mtot ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean1)
anova(C) #Hauteurini significatif, Intéraction Stress*Provenance, tendance de Stress
summary(C) #Stress*Prov2050 plus grand que 2018 (2018 pas significatif)
#Stress*Prov égal entre 2018 et 2080 (pas significatifs les deux)

#Avec Prov2050 comme valeur de base de reference (sans outlier)
CHR_mean1$Prov <- relevel(CHR_mean1$Prov, ref="2050")
C <- lmerTest::lmer(Mtot ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean1)
summary(C) #Stress:Prov2080 plus petit que 2050


  ###QUESTIONS!!!!!!
  #Si on met Brout comme référence, on obtient aucune différence significative
  #J'ai trouvé ça pour faire toutes les comparaisons pour une interaction:
  lsmeans(C, pairwise~Stress*Prov)
  # Seulement Stress2-NoStress de 2050 significatif


#Representation effet Stress*Prov
meansCHR_Mtot <- summarySE(CHR_mean1, measurevar = "Mtot", groupvars = c("Stress", "Prov")) 
meansCHR_Mtot <- mutate(meansCHR_Mtot, Stress=factor(Stress, levels=c("NoStress", "Stress2")))
#Une option est de faire un graphique pour les broutés et les non broutés. 
#On dirait que c'est le traitement qui a le moins d'effet
ggplot(meansCHR_Mtot, aes(x = Stress, y = Mtot, shape=factor(Prov))) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=Mtot-se, ymax=Mtot+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse totale (g)")

##Vérifier résidus
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
D <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CHR_mean1)
anova(D) #Rien significatif

##Modele avec plan en tiroir sans interaction Hini
D <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean1)
anova(D) #Rien significatif

##Modele avec plan en tiroir sans Hini
D <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CHR_mean1)
anova(D) #Tendance Stress

##Vérifier résidus 
plot(D) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(D) ~ Stress*Brout*Prov, data = CHR_mean) #ok

qqnorm(resid(D)) #Graphique de normalité ok
qqline(resid(D))
shapiro.test(resid(D)) #Non transformation à faire

#Ajouter variable Transformation log de Ratio
LCHR_mean1 <- mutate(CHR_mean1, LRatio=log(Ratio))

##Modele avec plan en tiroir (log)
DD <- lmerTest::lmer( LRatio ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = LCHR_mean1)
anova(DD) #Rien significatif

##Modele avec plan en tiroir sans interaction Hini (log)
DD <- lmerTest::lmer( LRatio ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = LCHR_mean1)
anova(DD) #Rien significatif

##Modele avec plan en tiroir sans Hini (log)
DD <- lmerTest::lmer( LRatio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = LCHR_mean1)
anova(DD) #Tendance Stress, tendance Prov

#Représentation tendance effet Stress
meansLCHR_Ratio <- summarySE(LCHR_mean1, measurevar = "LRatio", groupvars = c("Stress")) 
meansLCHR_Ratio <- mutate(meansLCHR_Ratio, Stress=factor(Stress, levels=c("NoStress", "Stress2")))
ggplot(meansLCHR_Ratio, aes(x = Stress, y = LRatio)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=LRatio-se, ymax=LRatio+se), width=.1)+
  xlab("Traitement de stress hydrique") +  
  ylab("Log Ratio aérien/racinaire")
#Ratio tendance à être plus élevée pour Stress2

#Tendance Prov
ls_means(DD, which = "Prov", pairwise = TRUE)
#Prov 2050 plus grand ratio que 2018 

#Représentation tendance effet Prov
meansLCHR_Ratio <- summarySE(LCHR_mean1, measurevar = "LRatio", groupvars = c("Prov"))
meansLCHR_Ratio <- mutate(meansLCHR_Ratio, Prov=factor(Prov, levels=c("2018", "2050", "2080")))
ggplot(meansLCHR_Ratio, aes(x = Prov, y = LRatio)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=LRatio-se, ymax=LRatio+se), width=.1)+
  xlab("Provenance") +  
  ylab("Log Ratio aérien/racinaire")

#Représentation tendance Stress et Prov
meansLCHR_Ratio <- summarySE(LCHR_mean1, measurevar = "LRatio", groupvars = c("Stress", "Prov")) 
meansLCHR_Ratio <- mutate(meansLCHR_Ratio, Stress=factor(Stress, levels=c("NoStress", "Stress2")))
ggplot(meansLCHR_Ratio, aes(x = Stress, y = LRatio, shape=factor(Prov))) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=LRatio-se, ymax=LRatio+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de stress hydrique") +  
  ylab("log Ratio")

#!!!!!!!!!!Pourquoi ça dit pas effet Stress*Prov vu qu'il y a juste 2050 de significatif pour stress

##Vérifier résidus 
plot(DD) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(DD) ~ Stress*Brout*Prov, data = LCHR_mean1) #ok

qqnorm(resid(DD)) #Graphique de normalité ok
qqline(resid(DD))
shapiro.test(resid(DD)) #ok

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

#Effet Stress*Prov
CHR_mean$Prov <- relevel(CHR_mean$Prov, ref="2018")
CHR_mean$Brout <- relevel(CHR_mean$Brout, ref="NoBrout")
CHR_mean$Stress <- relevel(CHR_mean$Stress, ref="Stress2")
CC <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
summary(CC) #Stress*Prov2050 plus grande que 2018
#Stress*Prov2080 égal à 2018 (Pas significatif les deux)

#2050 comme reference
CHR_mean$Prov <- relevel(CHR_mean$Prov, ref="2050")
CC <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
summary(CC) #Stress*Prov2080 (pas significatif) est plus petit que 2050

#!!!!!!Si on met Brout comme référence on obtient rien de significatif

#Representation effet Stress*Prov
meansCHR_Mracine <- summarySE(CHR_mean, measurevar = "Mracine", groupvars = c("Stress", "Prov")) 
meansCHR_Mracine <- mutate(meansCHR_Mracine, Stress=factor(Stress, levels=c("NoStress", "Stress2")))
#Une option est de faire un graphique pour les broutés et les non broutés. 
#On dirait que c'est le traitement qui a le moins d'effet
ggplot(meansCHR_Mracine, aes(x = Stress, y = Mracine, shape=factor(Prov))) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=Mracine-se, ymax=Mracine+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse racinaire (g)")

##Vérifier résidus 
plot(CC) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(CC) ~ Stress*Brout*Prov, data = CHR_mean) #ok

qqnorm(resid(CC)) #Graphique de normalité ok
qqline(resid(CC))
shapiro.test(resid(CC)) #Essayer une transformation log

#Ajouter variable Transformation log de Mracine
LCHR_mean <- mutate(CHR, LMracine=log(Mracine))

##Modele avec plan en tiroir Log (Transformation log pour Mracine)
CC1 <- lmerTest::lmer( LMracine ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = LCHR_mean)
anova(CC1) #Garder Hini seul

##Modele avec Hini seul (Transfo log)
CC1 <- lmerTest::lmer( LMracine ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = LCHR_mean)
anova(CC1) #Effet Hini, Stress*Prov

#Effet Stress*Prov (Log)
LCHR_mean$Prov <- relevel(LCHR_mean$Prov, ref="2018")
LCHR_mean$Brout <- relevel(LCHR_mean$Brout, ref="NoBrout")
LCHR_mean$Stress <- relevel(LCHR_mean$Stress, ref="Stress2")
CC <- lmerTest::lmer( LMracine ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = LCHR_mean)
summary(CC) #NoStress-Stress de Prov2050 tendance plus grande que 2018
#NoStress-Stress de Prov2080 égal à 2018 (Pas significatif les deux)

#2050 comme reference (Log)
LCHR_mean$Prov <- relevel(LCHR_mean$Prov, ref="2050")
CC <- lmerTest::lmer( LMracine ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = LCHR_mean)
summary(CC) #NoStress-Stress de Prov2080 (pas significatif) est plus petit que 2050

#Conclusion : Seul 2050 a une diminution de biomasse racinaire avec le stress hydrique

#Representation effet Stress*Prov (Transfo log)
meansCHR_LMracine <- summarySE(LCHR_mean, measurevar = "LMracine", groupvars = c("Stress", "Prov")) 
meansCHR_LMracine <- mutate(meansCHR_LMracine, Stress=factor(Stress, levels=c("NoStress", "Stress2")))
#Une option est de faire un graphique pour les broutés et les non broutés. 
#On dirait que c'est le traitement qui a le moins d'effet
ggplot(meansCHR_LMracine, aes(x = Stress, y = LMracine, shape=factor(Prov))) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=LMracine-se, ymax=LMracine+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de stress hydrique") +  
  ylab("Log Masse racinaire (g)")

##Vérifier résidus 
plot(CC1) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(CC1) ~ Stress*Brout*Prov, data = CHR_mean) #ok

qqnorm(resid(CC1)) #Graphique de normalité ok
qqline(resid(CC1))
shapiro.test(resid(CC1)) #ok


##MASSE AERIENNE##

#Boxplot masse racine
boxplot(Maerien~Brout+Stress+Prov,data=CHR,cex.axis=0.5, col = color,las=2, 
        main = "Masses aeriennes chenes", ylab="Masse aerienne", xlab="") 

##Modele avec plan en tiroir
DD <- lmerTest::lmer( Maerien ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(DD) #pas interaction avec Hini

##Modele avec plan en tiroir sans interaction Hauteurini
CHR_mean$Prov <- relevel(CHR_mean$Prov, ref="2018")
CHR_mean$Brout <- relevel(CHR_mean$Brout, ref="NoBrout")
CHR_mean$Stress <- relevel(CHR_mean$Stress, ref="Stress2")
DD <- lmerTest::lmer(Maerien ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(DD)#Tendance Prov et Tendance Stress:Prov
summary(DD)
#Stress*Prov2050 plus elevé que 2018
#Stress2*Prov2080 égal à 2018 (pas significatif les deux)

#Reference Prov2050
CHR_mean$Prov <- relevel(CHR_mean$Prov, ref="2050")
DD <- lmerTest::lmer(Maerien ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
summary(DD)#Stress*Prov2080 plus faible que 2050

#Representation Tendance Stress*Prov
meansCHR_Maerien <- summarySE(CHR_mean1, measurevar = "Maerien", groupvars = c("Stress", "Prov")) 
meansCHR_Maerien <- mutate(meansCHR_Maerien, Stress=factor(Stress, levels=c("NoStress", "Stress2")))
#Une option est de faire un graphique pour les broutés et les non broutés. 
#On dirait que c'est le traitement qui a le moins d'effet
ggplot(meansCHR_Maerien, aes(x = Stress, y = Maerien, shape=factor(Prov))) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=Maerien-se, ymax=Maerien+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse aerienne (g)")

##Conclusion 
CHRStress2018<-filter(CHR,Stress=="Stress2",Prov=="2018") 
mean(CHRStress2018$Maerien) #20,97
mean(CHRStress2018$Mracine) #21,19
mean(CHRStress2018$Ratio) #1,05
CHRNoStress2018<-filter(CHR,Stress=="NoStress",Prov=="2018") 
mean(CHRNoStress2018$Maerien) #17,71
mean(CHRNoStress2018$Mracine) #20,23
mean(CHRNoStress2018$Ratio) #0,927

CHRStress2050<-filter(CHR,Stress=="Stress2",Prov=="2050") 
mean(CHRStress2050$Maerien) #20,92
mean(CHRStress2050$Mracine) #15,32
mean(CHRStress2050$Ratio) #1,41
CHRNoStress2050<-filter(CHR,Stress=="NoStress",Prov=="2050") 
mean(CHRNoStress2050$Maerien) #30,9
mean(CHRNoStress2050$Mracine) #30,48
mean(CHRNoStress2050$Ratio) #1,05

CHRStress2080<-filter(CHR,Stress=="Stress2",Prov=="2080") 
mean(CHRStress2080$Maerien) #18,1
mean(CHRStress2080$Mracine) #17,2
mean(CHRStress2080$Ratio) #1,12
CHRNoStress2080<-filter(CHR,Stress=="NoStress",Prov=="2080") 
mean(CHRNoStress2080$Maerien) #20,07
mean(CHRNoStress2080$Mracine) #21,44
mean(CHRNoStress2080$Ratio) #1,02

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

#Effet Stress*Brout
summary(E) 
#Différence entre Stress1 et NoStress de Brout plus grande que NoBrout (sens inverse)
#Différence entre Stress2 et NoStress Brout tendance a plus grande que Brout (sens inverse)

#Avec Stress1 comme valeur de base de reference
ERS$Stress <- relevel(ERS$Stress, ref="Stress1")
E2 <- lmerTest::lmer(Mtot ~ Stress*Brout + (1|Bloc/Stress), data = ERS)
summary(E2)
#Différence entre Stress2 et Stress 1 de Brout tendance a plus petit que NoBrout (sens inverse)

#Representation Stress*Brout
meansERS_Mtot <- summarySE(ERS, measurevar = "Mtot", groupvars = c("Stress", "Brout")) 
meansERS_Mtot <- mutate(meansERS_Mtot, Stress=factor(Stress, levels=c("NoStress","Stress1","Stress2")))
ggplot(meansERS_Mtot, aes(x = Stress, y = Mtot, shape=factor(Brout))) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=Mtot-se, ymax=Mtot+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse totale (g)")

#COnclusions effets inverses entre Brout et NoBrout
#Si on regarde les effets séparément, différence significative seulement entre NoStress
#et Stress1 autant pour Brout que NoBrou, mais effets inverses

##Vérifier résidus en incluant les effets aléatoires (effet du design)
plot(E) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône) --> ok
leveneTest(resid(E) ~ Stress*Brout, data = ERS) #ok

qqnorm(resid(E)) #Graphique de normalité ok
qqline(resid(E))
shapiro.test(resid(E)) #ok

####Ratio ERABLES####

#Boxplot ratio 
boxplot(Ratio~Brout+Stress,data=ERS,cex.axis=0.5,las=2, main = "Ratio Erables", 
        xlab="",ylab="Ratio")

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
summary(EE) 
#Différence Stress1-NoStress de Brout plus grande que NoBrout (sens inverse)
#Différence Stress2-NoStress de Brout égal à NoBrout

##Stress1 en reference
ERS$Stress <- relevel(ERS$Stress, ref="Stress1")
EE <- lmerTest::lmer( Mracine ~ Stress*Brout + (1|Bloc/Stress), data = ERS)
summary(EE) 
#Différence Stress2-Stress1 de Brout est plus petite pour NoBrout (sens inverse)

#Representation Stress*Brout
meansERS_Mracine <- summarySE(ERS, measurevar = "Mracine", groupvars = c("Stress", "Brout")) 
meansERS_Mracine <- mutate(meansERS_Mracine, Stress=factor(Stress, levels=c("NoStress","Stress1","Stress2")))
ggplot(meansERS_Mracine, aes(x = Stress, y = Mracine, shape=factor(Brout))) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=Mracine-se, ymax=Mracine+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse racine (g)")

##Vérifier résidus 
plot(EE) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(EE) ~ Stress*Brout*Prov, data = ERS) #ok

qqnorm(resid(EE)) #Graphique de normalité ok
qqline(resid(EE))
shapiro.test(resid(EE)) #ok

  
##MASSE AERIENNE##
  
#Boxplot masse aerien
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
FF <- lmerTest::lmer(Maerien ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = ERS)
summary(FF) 
#Différence Stress1-NoStress plus grand pour Brout que NoBrout (sens inverse)
#Différence Stress2-NoStress plus grand pour Brout que NoBrout (sens inverse)

#Reference Stress1
ERS$Stress <- relevel(ERS$Stress, ref="Stress1")
FF <- lmerTest::lmer(Maerien ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = ERS)
summary(FF) 
#Stress2-Stress1 pas différent entre Brout et NoBrout ???????????

#Representation Stress*Brout
meansERS_Maerien <- summarySE(ERS, measurevar = "Maerien", groupvars = c("Stress", "Brout"))
meansERS_Maerien <- mutate(meansERS_Maerien, Stress=factor(Stress, levels=c("NoStress","Stress1","Stress2")))
ggplot(meansERS_Maerien, aes(x = Stress, y = Maerien, shape=factor(Brout))) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=Maerien-se, ymax=Maerien+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse aerienne (g)")


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

#Effet Stress
ls_means(H, which = "Stress", pairwise = TRUE) 
#NoStress pas different de Stress1
#NoStress plus grand que Stress2
#Stress1 plus grand que Stress2

#Effet Brout
  #Masse totale plus faible pour plants broutés

#Representation effet stress et brout
meansTHO_Mtot <- summarySE(THO, measurevar = "Mtot", groupvars = c("Stress", "Brout")) 
#Une option est de faire un graphique pour les broutés et les non broutés. 
#On dirait que c'est le traitement qui a le moins d'effet
ggplot(meansTHO_Mtot, aes(x = Stress, y = Mtot, shape=factor(Brout))) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=Mtot-se, ymax=Mtot+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse totale (g)")

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

#Representation Stress*Brout*Prov
meansTHO_ratio <- summarySE(THO, measurevar = "Ratio", groupvars = c("Stress", "Brout", "Prov")) 
meansTHO_ratio <- mutate(meansTHO_ratio, Stress=factor(Stress, levels=c("NoStress","Stress1","Stress2")))
#Faire un graphique pour les broutés et les non broutés. 
ggplot(meansTHO_ratio[which(meansTHO_ratio$Brout == "NoBrout"),], aes(x = Stress, y = Ratio, shape=factor(Prov))) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=Ratio-se, ymax=Ratio+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de stress hydrique") +  
  ylab("Ratio aérienne/racinaire")+
  ylim(c(2,5.5))+
  annotate("text", x = 1, y = 5,label = "Non broutés", size = 5)
ggplot(meansTHO_ratio[which(meansTHO_ratio$Brout == "Brout"),], aes(x = Stress, y = Ratio, shape=factor(Prov))) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=Ratio-se, ymax=Ratio+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de stress hydrique") +  
  ylab("Ratio aérienne/racinaire")+
  ylim(c(2,5.5))+
  annotate("text", x = 1, y = 5,label = "Broutés", size = 5)


##Analyse interaction Stress*Brout*Prov
#Avec Prov2018 comme valeur de base de reference
THO$Stress <- relevel(THO$Stress, ref="NoStress")
THO$Prov <- relevel(THO$Prov, ref="2018")
THO$Brout <- relevel(THO$Brout, ref="NoBrout")
I <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = THO)
summary(I) #Stress2*Brout*Prov2080 significatif (comparaison avec 2018)
#Stress2*Brout significatif pour 2018

#Avec Prov2080 comme valeur de base de reference
THO$Prov <- relevel(THO$Prov, ref="2080")
I <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = THO)
summary(I) #Rien significatif

##Interpretation de l'interaction Stress*Brout*Prov
#Effet Stress2*Brout significativement différent entre 2018 et 2080 =
#Différence entre Stress2 et NoStress est plus grande pour NoBrout que Brout pour Prov2018 (sens inverse)
#Alors que cette différence n'est pas significative pour Prov2080

#Selon graph
#Effet surtout pour 2018
#Ratio augmente avec le stress pour non broutés
#Vs ne change pas avec le stress pour broutés
#Autres Prov ne changent pas


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
#Stress2 plus petit que NoStress
#Stress1 plus petit que NoStress

#Stress1 comme reference
THO$Stress <- relevel(THO$Stress, ref="Stress1")
II <- lmerTest::lmer( Mracine ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = THO)
summary(II) #Stress2 plus petit que NoStress

#Representation effet stress et brout
meansTHO_Mracine <- summarySE(THO, measurevar = "Mracine", groupvars = c("Stress", "Brout")) 
meansTHO_Mracine <- mutate(meansTHO_Mracine, Stress=factor(Stress, levels=c("NoStress","Stress1","Stress2")))
ggplot(meansTHO_Mracine, aes(x = Stress, y = Mracine, shape=factor(Brout))) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=Mracine-se, ymax=Mracine+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse racinaire (g)")

##Vérifier résidus 
plot(II) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(II) ~ Stress*Brout*Prov, data = THO) #ok

qqnorm(resid(II)) #Graphique de normalité ok
qqline(resid(II))
shapiro.test(resid(II)) #Non Transformation à faire

#Ajouter variable Transformation log de Mracine
THO <- mutate(THO, LMracine=log(Mracine))

##Modele avec plan en tiroir Log (Transformation log pour Mracine)
II1 <- lmerTest::lmer( LMracine ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = THO)
anova(II1) #On garde le modele, effet Brout, tendance Stress*Brout (0.9) 

#Representation effet stress et brout
meansTHO_LMracine <- summarySE(THO, measurevar = "LMracine", groupvars = c("Brout")) 
meansTHO_Mracine <- mutate(meansTHO_Mracine, Stress=factor(Stress, levels=c("NoStress","Stress1","Stress2")))
ggplot(meansTHO_LMracine, aes(x = Brout, y = LMracine)) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=LMracine-se, ymax=LMracine+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de broutement") +  
  ylab("log Masse racinaire (g)")
#Plants broutés ont masse racinaire plus faible

##Vérifier résidus 
plot(II1) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(II1) ~ Stress*Brout*Prov, data = THO) #ok

qqnorm(resid(II1)) #Graphique de normalité ok
qqline(resid(II1))
shapiro.test(resid(II1)) #ok

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
#Masse aérienne plus faible pour les plants boutés

#Effet Stress
ls_means(JJ, which = "Stress", pairwise = TRUE) #Le test pour toutes les comparaisons, sans changer le niveau de référence
#Stress2 plus faible que NoStress
#Stress2 plus faible que Stress1

#Representation effet stress
meansTHO_Maerien <- summarySE(THO, measurevar = "Maerien", groupvars = c("Stress")) 
meansTHO_Maerien <- mutate(meansTHO_Maerien, Stress=factor(Stress, levels=c("NoStress","Stress1","Stress2")))
ggplot(meansTHO_Maerien, aes(x = Stress, y = Maerien)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Maerien-se, ymax=Maerien+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse aerienne (g)")

#Representation effet brout
meansTHO_Maerien <- summarySE(THO, measurevar = "Maerien", groupvars = c("Brout"))
meansTHO_Maerien <- mutate(meansTHO_Maerien, Brout=factor(Brout, levels=c("NoBrout","Brout")))
ggplot(meansTHO_Maerien, aes(x = Brout, y = Maerien)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Maerien-se, ymax=Maerien+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de broutement") +  
  ylab("Masse aerienne (g)")

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

##????????????
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

##COMPARAISON AVEC BIOMASSE RETIRÉE##
THOBROUT<-filter(THO,Brout=="Brout") 
mean(THOBROUT$Mtot) #37.66
mean(THOBROUT$Mret) #2.10 g retire
mean(THOBROUT$Maerien) #29.00 + 2.10 = 31.10 (avec Mretiree)
mean(THOBROUT$Mracine) #8.66

THONOBROUT<-filter(THO,Brout=="NoBrout") 
mean(THONOBROUT$Mtot) #54.44
mean(THONOBROUT$Maerien) #41.53
mean(THONOBROUT$Mracine) #12.91

THOavecMret <- mutate(THO, Mtotret=Mtot+Mret, Maerien=Maerien+Mret)

##Modele avec plan en tiroir pour Maerien (ajout de Mret)
JJ1 <- lmerTest::lmer( Maerien ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = THOavecMret)
anova(JJ1) #Pas d'interaction significative

##Modele avec plan en tiroir sans interaction Hini pour Maerien (ajout de Mret)
JJ1 <- lmerTest::lmer( Maerien ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = PIBavecMret)
anova(JJ1) #Brout significatif
summary(JJ1)
#Malgre l'ajout de Mret, ceux broutés ont une Maerienne significativement plus basse
#Leur croissance a donc été ralentie suite au traitement de broutement


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
anova(J) #Hauteurini, Brout significatif
summary(J)

#Representation effet brout
meansPIB_Mtot <- summarySE(PIB, measurevar = "Mtot", groupvars = c("Brout")) 
meansPIB_Mtot <- mutate(meansPIB_Mtot, Brout=factor(Brout, levels=c("NoBrout","Brout")))
ggplot(meansPIB_Mtot, aes(x = Brout, y = Mtot)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Mtot-se, ymax=Mtot+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de broutement") +  
  ylab("Masse totale (g)")
#Biomasse totale diminue avec le broutement

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

#Representation effet brout
meansPIB_Mracine <- summarySE(PIB, measurevar = "Mracine", groupvars = c("Brout")) 
meansPIB_Mracine <- mutate(meansPIB_Mracine, Brout=factor(Brout, levels=c("NoBrout","Brout")))
ggplot(meansPIB_Mracine, aes(x = Brout, y = Mracine)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Mracine-se, ymax=Mracine+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de broutement") +  
  ylab("Masse racinaire (g)")
#Biomasse racinaire diminue avec le broutement

##Vérifier résidus 
plot(KK) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(KK) ~ Stress*Brout*Prov, data = PIB) #ok

qqnorm(resid(KK)) #Graphique de normalité ok
qqline(resid(KK))
shapiro.test(resid(KK)) #Non transformation à faire

#Ajouter variable Transformation log de Mracine
PIB <- mutate(PIB, LMracine=log(Mracine))

##Modele avec plan en tiroir Log (Transformation log pour Mracine)
KK1 <- lmerTest::lmer( LMracine ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = PIB)
anova(KK1) #On garde ce modele, Interactions avec Hini significatives

##Vérifier résidus (log)
plot(KK1) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(KK1) ~ Stress*Brout*Prov, data = PIB) #ok

qqnorm(resid(KK1)) #Graphique de normalité 
qqline(resid(KK1))
shapiro.test(resid(KK1)) #Non, essayer une autre transformation

#Ajouter variable Transformation racine carrée
PIB <- mutate(PIB, rMracine=sqrt(Mracine))

##Modele avec plan en tiroir (Transformation racine carrée pour Mracine)
KK2 <- lmerTest::lmer( rMracine ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = PIB)
anova(KK2) #Rien significatif

##Modele avec plan en tiroir sans interaction Hini(Transformation racine carrée pour Mracine)
KK2 <- lmerTest::lmer( rMracine ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = PIB)
anova(KK2) #Tendance Hini, Brout significatif

##Si on enleve Hini(Transformation racine carrée pour Mracine)
KK2 <- lmerTest::lmer( rMracine ~ Stress*Brout*Prov + (1|Bloc/Stress), data = PIB)
anova(KK2) #Brout significatif

#Representation effet brout
meansPIB_rMracine <- summarySE(PIB, measurevar = "rMracine", groupvars = c("Brout")) 
meansPIB_rMracine <- mutate(meansPIB_rMracine, Brout=factor(Brout, levels=c("NoBrout","Brout")))
ggplot(meansPIB_rMracine, aes(x = Brout, y = rMracine)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=rMracine-se, ymax=rMracine+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de broutement") +  
  ylab("sqrt Masse racinaire (g)")
#Biomasse racinaire diminue avec le broutement

##Vérifier résidus (racine carée)
plot(KK2) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(KK2) ~ Stress*Brout*Prov, data = PIB) #ok

qqnorm(resid(KK2)) #Graphique de normalité 
qqline(resid(KK2))
shapiro.test(resid(KK2)) #Ok
  
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

#Representation effet brout
meansPIB_Maerien <- summarySE(PIB, measurevar = "Maerien", groupvars = c("Brout")) 
meansPIB_Maerien <- mutate(meansPIB_Maerien, Brout=factor(Brout, levels=c("NoBrout","Brout")))
ggplot(meansPIB_Maerien, aes(x = Brout, y = Maerien)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Maerien-se, ymax=Maerien+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de broutement") +  
  ylab("Masse aerienne (g)")
#Biomasse aérienne diminue avec le broutement

##Vérifier résidus 
plot(LL) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(LL) ~ Stress*Brout*Prov, data = PIB) #ok

qqnorm(resid(LL)) #Graphique de normalité ok
qqline(resid(LL))
shapiro.test(resid(LL)) #ok

##Conclusion Ratio, aerien et racinaire##
#Ratio ne change pas
#Alors que Maerien et racinaire sont diminuées par le broutement

##Vérifier résidus 
plot(LL) #Homogénéité des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(LL) ~ Stress*Brout*Prov, data = PIB) #ok

qqnorm(resid(LL)) #Graphique de normalité ok
qqline(resid(LL))
shapiro.test(resid(LL)) #ok

##COMPARAISON AVEC BIOMASSE RETIRÉE##
PIBBROUT<-filter(PIB,Brout=="Brout") 
mean(PIBBROUT$Mtot) #13.05
mean(PIBBROUT$Mret) #4.06 (Mretiree)
mean(PIBBROUT$Maerien) #8.29 + 4.06 = 12.35 (avec Mretiree)
mean(PIBBROUT$Mracine) #4.77

PIBNOBROUT<-filter(PIB,Brout=="NoBrout") 
mean(PIBNOBROUT$Mtot) #32.03
mean(PIBNOBROUT$Maerien) #20.82
mean(PIBNOBROUT$Mracine) #11.22

PIBavecMret <- mutate(PIB, Mtotret=Mtot+Mret, Maerien=Maerien+Mret)

##Modele avec plan en tiroir pour Maerien (ajout de Mret)
LL1 <- lmerTest::lmer( Maerien ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = PIBavecMret)
anova(LL1) #Rien significatif

##Modele avec plan en tiroir sans interaction Hini pour Maerien (ajout de Mret)
LL1 <- lmerTest::lmer( Maerien ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = PIBavecMret)
anova(LL1) #Hauteurini significatif, Effet Brout
summary(LL1)
#Malgre l'ajout de Mret, ceux broutés ont une Maerienne significativement plus basse
#Leur croissance a donc été ralentie suite au traitement de broutement

##SURVIE PINS##
MortPIB <- read.table("./MortPIB.txt", header=TRUE, na.string = "", dec = ",") 
MortPIB$Prov <- as.factor(as.character(MortPIB$Prov))
MortPIB$Brout <- as.factor(as.character(MortPIB$Brout))
MortPIB$Stress <- as.factor(as.character(MortPIB$Stress))
MortPIB$Bloc <- as.factor(as.character(MortPIB$Bloc))

Z <- glm(Mort ~ Stress + Brout + Prov + Bloc, data = MortPIB, family = binomial(link='logit'))
summary(Z)
anova(Z, test="Chisq")
#Mortalite plus grande chez les pins broutés
#Tendance de l'effet du bloc 

MortPIB = data.frame(Brout = "1")
predict(Z, MortPIB, type="response")
