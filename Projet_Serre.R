#2019-2020
#Turgeon Roxanne 
#Initiation a la recherche / Migration assistee / Broutement / Stress_hydrique / Serre
#----

remove(list = ls())

library(plyr)
library(tidyverse) 
library(dplyr)
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
library(cowplot)
library(MASS)
library(gridExtra)
library(gridGraphics)


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

#################TOUTES ESPECES##########

boxplot(Mtot~Brout+Esp,data=Masses3,cex.axis=0.5, las=2, ylab="Masse tot", xlab="")
boxplot(Mtot~Stress+Esp,data=Masses3,cex.axis=0.5, las=2, ylab="Masse tot", xlab="")
boxplot(Mtot~Bloc*Esp,data=Masses3, cex.axis=0.5, las=2, ylab="Masse tot", xlab="")

##Ratio Toutes espèces sans traitement (Temoin) ##
boxplot(Ratio~Esp, data=Masses3Temoin, ylab="Ratio", xlab="Espèces")
meansMasses3_Temoin <- summarySE(Masses3Temoin, measurevar = "Ratio", groupvars = c("Esp")) 
ggplot(meansMasses3_Temoin, aes(x = Esp, y = Ratio)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Ratio-se, ymax=Ratio+se), width=.1)+
  xlab("Espèces") +  
  ylab("Ratio aerien/racinaire")



#################CERISIERS########################
####BIOMASSE TOTALE CERISIERS####

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

#Modele avec plan en tiroir
A <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CET_mean)
#Message veut dire que pour cette espece, l'effet aleatoire du bloc = 0
#On le garde pareil dans le modele, mais il a peu d'effet
anova(A)
anova_A <- anova(A) # Interaction significative avec Hini, on garde ce modele
#Hini significatif, effet   Stress
summary(lsmeans(A, ~Stress))


#Representation effet Stress
meansCET_Mtot <- summarySE(CET_mean, measurevar = "Mtot", groupvars = c("Stress")) 
meansCET_Mtot <- mutate(meansCET_Mtot, Stress=factor(Stress, levels=c("NoStress", "Stress2")))
ggplot(meansCET_Mtot, aes(x = Stress, y = Mtot)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Mtot-se, ymax=Mtot+se), width=.1)+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse totale (g)")
#Mtot Stress2 plus faible

#Representation Stress*Prov avec les estimations du modele
CET_Mtot_mod <- summary(lsmeans(A, ~Stress))
ggplot(CET_Mtot_mod, aes(x = Stress, y = lsmean, col=Stress)) +
  geom_point(size=4, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  scale_color_manual(values=c("limegreen", "lightcoral"),labels=c("Non stress?", "Stress?"))+
  scale_x_discrete(labels=c("Non stress?","Stress ?lev?"))+
  scale_y_continuous(limits = c(0,75), expand = expand_scale()) +
  xlab("Traitement de stress hydrique") +  
  ylab("Masse (g)") +
  theme(text=element_text(size=18, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))


##Verifier residus en incluant les effets aleatoires (effet du design)
plot(A) #Homogeneite des variances, on ne doit pas voir de patron particulier (cône) --> ok
leveneTest(resid(A) ~ Stress*Brout*Prov, data = CET_mean) #ok

qqnorm(resid(A)) #Graphique de normalite ok
qqline(resid(A))
shapiro.test(resid(A)) #On verifie avec un test, faire des transformations au besoin --> ok

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
anova(B)
anova_B <- anova(B) #Provenance et Stress significatif 


#summary
ls_means(B, which = "Prov", pairwise = TRUE) #Le test pour toutes les comparaisons, sans changer le niveau de reference
#Prov 2050-2080 et 2018-2080 differents

#Representation de Stress
meansCET_Ratio <- summarySE(CET_mean, measurevar = "Ratio", groupvars = c("Stress")) 
ggplot(meansCET_Ratio, aes(x = Stress, y = Ratio)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Ratio-se, ymax=Ratio+se), width=.1)+
  xlab("Traitement de stress hydrique") +  
  ylab("Ratio aerien/racinaire")
#Ratio de stress2 plus grand que NoStress
summarySE(CET_mean, measurevar = "Ratio", groupvars = c("Stress")) 
#Stress avec les estimations du modele
summary(lsmeans(B, ~Stress))

#Representation de Prov
meansCET_Ratio <- summarySE(CET_mean, measurevar = "Ratio", groupvars = c("Prov")) 
meansCET_Ratio <- mutate(meansCET_Ratio, Prov=factor(Prov, levels=c("2018","2050","2080")))

ggplot(meansCET_Ratio, aes(x = Prov, y = Ratio)) +
  geom_point(size=3) +
  theme_classic() +
  geom_errorbar(aes(ymin=Ratio-se, ymax=Ratio+se), width=.1)+
  xlab("Provenance") +  
  ylab("Ratio a?rien/racinaire") +
  geom_text(label = c("a","a","b"), vjust= -5)

#Representation Prov avec les estimations du modele
CET_Ratio_mod <- summary(lsmeans(B, ~Prov))
ggplot(CET_Ratio_mod, aes(x = Prov, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  theme_classic() +
  xlab("Provenance") +  
  ylab("Ratio aerien/racinaire estimations modele") +
  geom_text(aes(x=Prov, y=upper.CL+0.05),label = c("a","a","b"))

##Verifier residus 
plot(B) #Homogeneite des variances, on ne doit pas voir de patron particulier (cone) ok
leveneTest(resid(B) ~ Stress*Brout*Prov, data = CET_mean) #ok

qqnorm(resid(B)) #Graphique de normalite ok
qqline(resid(B))
shapiro.test(resid(B)) #On verifie avec un test, faire des transformations au besoin --> ok

####MASSE RACINAIRE####

#Boxplot masse racine
color = c(rep("green",4),rep("yellow",4),rep("red",4))
boxplot(Mracine~Brout+Stress+Prov,data=CET,cex.axis=0.5, col = color,las=2, 
        main = "Masses racinaire cerisiers", ylab="Masse racinaire", xlab="") 

##Modele avec plan en tiroir
AA <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CET_mean)
anova(AA)
anova_AA <- anova(AA) #Interaction avec Hini significatif --> On garde le modele
#Effet Stress, Brout*Prov, Stress*Brout*Prov, tendance Prov
#Lorsqu'il y a une interaction significative on laisse tomber les effets simples
#Ici on regarde seulement Stress*Brout*Prov
lsmeans(AA, pairwise~Stress*Brout*Prov)

#Combinaisons diff?rentes significativement
#Stress2,Brout,2018 - NoStress,NoBrout,2050 *
#NoStress,Brout,2080 - NoStress,NoBrout,2050
#Stress2,NoBrout,2080 - NoStress,NoBrout,2050 *
#Stress2,Brout,2050 - NoStress,NoBrout,2050 *
#Stress2,NoBrout,2050 - NoStress,NoBrout,2050 *

#Representation effet Stress*Brout*Prov
meansCET_Mracine <- summarySE(CET_mean, measurevar = "Mracine", groupvars = c("Stress", "Brout", "Prov")) 
meansCET_Mracine <- mutate(meansCET_Mracine, Stress=factor(Stress, levels=c("NoStress", "Stress2")))
meansCET_Mracine <- mutate(meansCET_Mracine, Prov=factor(Prov,levels=c("2018","2050","2080")))

ggplot(meansCET_Mracine, aes(x = Stress, y = Mracine, shape=factor(Prov),color=Brout)) +
  geom_point(size=3, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=Mracine-se, ymax=Mracine+se), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  xlab("Traitement de stress hydrique") +  
  ylab("Masse racinaire (g)")+
  ylim(c(5,35)) +
  geom_text(aes(x=Stress, y=Mracine+se+ 1),
            label = c("ab","b","a","ab","ab","ab",   "b","ab","b","b","ab","b"),  
            position=position_dodge(width=0.4)) 

#Representation Stress*Brout*Prov avec les estimations du modele
CET_Mracine_mod <- summary(lsmeans(AA, ~Prov*Stress*Brout))
ggplot(CET_Mracine_mod, aes(x = Stress, y = lsmean, shape=factor(Prov),color=Brout)) +
  geom_point(size=3, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  xlab("Provenance") +  
  ylab("Masse racinaire (g)") +
  labs(shape="Provenance", colour="Broutement") +
  scale_x_discrete(labels=c("Non stress?","Stress ?lev?"))+
  geom_text(aes(x=Stress, y=upper.CL+1),
            label = c("ab","b","a","ab","ab","ab",   "b","ab","b","b","ab","b"),  
            position=position_dodge(width=0.4), show.legend=F) 

#Representation Stress*Brout pour Prov2050 avec les estimations du modele
Graph_2050<-ggplot(CET_Mracine_mod [which(meansCET_Mracine$Prov == "2050"),], aes(x = Stress, y = lsmean,color=Brout)) +
  geom_point(size=2.5, shape=17, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  xlab("Stress hydrique") +  
  ylab("Masse racinaire (g)") +
  labs(colour="Broutement") +
  ylim(0,48)+
  scale_x_discrete(labels=c("Non stress?","Stress ?lev?"))+
  geom_text(aes(x=Stress, y=upper.CL+2),
            label = c("a","ab","b","b"),  
            position=position_dodge(width=0.4), show.legend=F) 

#Representation Stress*Brout pour Prov2018 avec les estimations du modele
Graph_2018<-ggplot(CET_Mracine_mod [which(meansCET_Mracine$Prov == "2018"),], aes(x = Stress, y = lsmean,color=Brout)) +
  geom_point(size=2.5, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  xlab("") +  
  ylab("Masse racinaire (g)") +
  labs(colour="Broutement") +
  ylim(0,48)+
  scale_x_discrete(labels=c("Non stress?","Stress ?lev?"))+
  geom_text(aes(x=Stress, y=upper.CL+2),
            label = c("ab","ab","ab","b"),  
            position=position_dodge(width=0.4), show.legend=F) 

#Representation Stress*Brout pour Prov2080 avec les estimations du modele
Graph_2080<-ggplot(CET_Mracine_mod [which(meansCET_Mracine$Prov == "2080"),], aes(x = Stress, y = lsmean,color=Brout)) +
  geom_point(size=2.5,shape=15, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  xlab("") +  
  ylab("Masse racinaire (g)") +
  ylim(0,48)+
  labs(colour="Broutement") +
  scale_x_discrete(labels=c("Non stress?","Stress ?lev?"))+
  scale_color_discrete(labels=c("Brout?", "Non brout?")) +
  geom_text(aes(x=Stress, y=upper.CL+2),
            label = c("ab","b","b","ab"),  
            position=position_dodge(width=0.4), show.legend=F)

#Panneau avec les trois provenances
get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(Graph_2080)
Graph_2018 <- Graph_2018 + theme(legend.position="none")
Graph_2050 <- Graph_2050 + theme(legend.position="none")
Graph_2080 <- Graph_2080 + theme(legend.position="none")

plot_grid(Graph_2018, Graph_2050, Graph_2080, legend, labels = c('2018', '2050','2080'), 
          label_size = 13,label_x = 0.15, rel_widths = c(1, 1, 1, 0.2), ncol=4)

##Verifier residus 
plot(AA) #Homogeneite des variances, on ne doit pas voir de patron particulier (cône) ok
leveneTest(resid(AA) ~ Stress*Brout*Prov, data = CET_mean) #ok

qqnorm(resid(AA)) #Graphique de normalite ok
qqline(resid(AA))
shapiro.test(resid(AA)) #ok

####MASSE AERIENNE####

#Boxplot masse aerienne
color = c(rep("green",4),rep("yellow",4),rep("red",4))
boxplot(Maerien~Brout+Stress+Prov,data=CET,cex.axis=0.5, col = color,las=2, 
        main = "Masses aeriennes cerisiers", ylab="Masse aerienne", xlab="") 

##Modele avec plan en tiroir
BB <- lmerTest::lmer( Maerien ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CET_mean)
anova(BB)
anova_BB <- anova(BB) 
# Une interaction Stress*Brout*Hauteurini significatif, Hini significatif
#Stress*Brout significatif

lsmeans(BB, pairwise~Stress*Brout) #Rien significatif

#Representation de Stress*Brout
meansCET_Maerien <- summarySE(CET_mean, measurevar = "Maerien", groupvars = c("Stress","Brout")) 
meansCET_Maerien <- mutate(meansCET_Maerien, Stress=factor(Stress, levels=c("NoStress", "Stress2")))

ggplot(meansCET_Maerien, aes(x = Stress, y = Maerien, color=Brout)) +
  geom_point(size=3, position=position_dodge(width=-0.2)) +
  geom_errorbar(aes(ymin=Maerien-se, ymax=Maerien+se), width=.1, position=position_dodge(width=-0.2))+
  theme_classic() +
  xlab("Traitement de stress hydrique") +  
  ylab("Masse aerienne (g)")

#Representation Stress*Brout avec les estimations du modele 
CET_Maerien_mod <- summary(lsmeans(BB, ~Stress*Brout))
ggplot(CET_Maerien_mod, aes(x = Stress, y = lsmean, color=Brout)) +
  geom_point(size=3, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  labs(colour="Broutement") +
  scale_x_discrete(labels=c("Non stress?","Stress ?lev?"))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse a?rienne (g)") 

##Verifier residus 
plot(BB) #Homogeneite des variances, on ne doit pas voir de patron particulier (cone) ok
leveneTest(resid(BB) ~ Stress*Brout*Prov, data = CET_mean) #ok

qqnorm(resid(BB)) #Graphique de normalite ok
qqline(resid(BB))
shapiro.test(resid(BB)) # ok


##Conclusion Ratio, aerien et racinaire##
#Les plants stresses presentent une croissance plus faible aerien et racinaire, mais nettement
#plus faible pour le racinaire donnant un ratio plus eleve que chez les plants non stresses

##Tables anova cerisiers
capture.output(anova_A, anova_B, anova_AA, anova_BB, file="Cerisiers.doc")

#################CHENES########################
####BIOMASSE TOTALE CHENES####

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

#Modele avec plan en tiroir
C <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(C) # Garder Hini seul 
summary(C)

#Sans interaction avec Hini:
CHR_mean$Prov <- relevel(CHR_mean$Prov, ref="2018")
CHR_mean$Brout <- relevel(CHR_mean$Brout, ref="NoBrout")
C <- lmerTest::lmer(Mtot ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(C) #Hauteurini significatif, Interaction Stress:Provenance, tendance de Stress
anova_C <- anova(C)
#On regarde seulement Stress*Prov

lsmeans(C, pairwise~Stress*Prov)
#NoStress2018 differe de NoStress 2050 et Stress2-NoStress de 2050 significatif

#Representation effet Stress*Prov
meansCHR_Mtot <- summarySE(CHR_mean, measurevar = "Mtot", groupvars = c("Stress", "Prov")) 
meansCHR_Mtot <- mutate(meansCHR_Mtot, Stress=factor(Stress, levels=c("NoStress", "Stress2")))
#Une option est de faire un graphique pour les broutes et les non broutes. 
#On dirait que c'est le traitement qui a le moins d'effet
ggplot(meansCHR_Mtot, aes(x = Stress, y = Mtot, shape=factor(Prov))) +
  geom_point(size=3, position=position_dodge(width=-0.4)) +
  theme_classic() +
  geom_text(aes(x=Stress, y=Mtot+se+2),label = c("ab","a","b","ab","b","ab"), position=position_dodge(width=-0.4))+
  geom_errorbar(aes(ymin=Mtot-se, ymax=Mtot+se), width=.1, position=position_dodge(width=-0.4))+
  ylim(28,72) +
  xlab("Traitement de stress hydrique") +  
  ylab("Masse totale (g)")

#Representation Stress*Prov avec les estimations du modele
CHR_Mtot_mod <- summary(lsmeans(C, ~Prov*Stress))
ggplot(CHR_Mtot_mod, aes(x = Stress, y = lsmean, shape=factor(Prov),col=Stress)) +
  geom_point(size=4, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  scale_shape_manual(values=c(10, 19, 1),labels=c("Nord", "Interm?diaire","Sud"))+
  scale_color_manual(values=c("limegreen", "lightcoral"),labels=c("Non stress?", "Stress?"))+
  labs(shape="Provenance") +
  scale_x_discrete(labels=c("Non stress?","Stress ?lev?"))+
  scale_y_continuous(limits = c(0,75), expand = expand_scale()) +
  xlab("Traitement de stress hydrique") +  
  ylab("Masse (g)") +
  theme(text=element_text(size=18, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+ 
  geom_text(aes(x=Stress, y=upper.CL+1.9),
            label = c("ab","a","b","ab","b","ab"),  
            position=position_dodge(width=0.4),
            show.legend=F, size=4.5) 


##Verifier residus
plot(C) #Homogeneite des variances, on ne doit pas voir de patron particulier (cône) --> ok
leveneTest(resid(C) ~ Stress*Brout*Prov, data = CHR_mean1) #ok

qqnorm(resid(C)) #Graphique de normalite ok
qqline(resid(C))
shapiro.test(resid(C)) #On verifie avec un test, faire des transformations au besoin --> ok

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
anova_D <- anova(D)

##Verifier residus 
plot(D) #Homogeneite des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(D) ~ Stress*Brout*Prov, data = CHR_mean) #ok

qqnorm(resid(D)) #Graphique de normalite ok
qqline(resid(D))
shapiro.test(resid(D)) #Non transformation a faire

#Ajouter variable Transformation log de Ratio
LCHR_mean <- mutate(CHR_mean, LRatio=log(Ratio))

##Modele avec plan en tiroir (log)
DD <- lmerTest::lmer( LRatio ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = LCHR_mean)
anova(DD) #Rien significatif

##Modele avec plan en tiroir sans interaction Hini (log)
DD <- lmerTest::lmer( LRatio ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = LCHR_mean)
anova(DD) #Rien significatif

##Modele avec plan en tiroir sans Hini (log)
DD <- lmerTest::lmer( LRatio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = LCHR_mean)
anova(DD) #Effet Prov, Tendance Stress 

ls_means(DD, which = "Prov", pairwise = TRUE)
#Prov2018-2050 significatif

#Representation tendance effet Stress
meansLCHR_Ratio <- summarySE(LCHR_mean, measurevar = "LRatio", groupvars = c("Stress")) 
meansLCHR_Ratio <- mutate(meansLCHR_Ratio, Stress=factor(Stress, levels=c("NoStress", "Stress2")))
ggplot(meansLCHR_Ratio, aes(x = Stress, y = LRatio)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=LRatio-se, ymax=LRatio+se), width=.1)+
  xlab("Traitement de stress hydrique") +  
  ylab("Log Ratio aerien/racinaire")
#Ratio tendance a etre plus elevee pour Stress2

#Representation effet Prov
meansLCHR_Ratio <- summarySE(LCHR_mean, measurevar = "LRatio", groupvars = c("Prov"))
meansLCHR_Ratio <- mutate(meansLCHR_Ratio, Prov=factor(Prov, levels=c("2018", "2050", "2080")))
ggplot(meansLCHR_Ratio, aes(x = Prov, y = LRatio)) +
  geom_point(size=3) +
  theme_classic()+
  geom_errorbar(aes(ymin=LRatio-se, ymax=LRatio+se), width=.1)+
  xlab("Provenance") +  
  ylab("Log Ratio aerien/racinaire")

#Representation Prov avec les estimations du modele
CHR_LRatio_mod <- summary(lsmeans(DD, ~Prov))
ggplot(CHR_LRatio_mod, aes(x = Prov, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  theme_classic() +
  xlab("Provenance") +  
  ylab("Log du ratio a?rien/racinaire")+
  geom_text(aes(x=Prov, y=upper.CL+0.02),
            label = c("a","b","ab")) 

##Verifier residus 
plot(DD) #Homogeneite des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(DD) ~ Stress*Brout*Prov, data = LCHR_mean1) #ok

qqnorm(resid(DD)) #Graphique de normalite ok
qqline(resid(DD))
shapiro.test(resid(DD)) #ok

####MASSE RACINAIRE####

#Boxplot masse racine
color = c(rep("green",4),rep("yellow",4),rep("red",4))
boxplot(Mracine~Brout+Stress+Prov,data=CHR,cex.axis=0.5, col = color,las=2, 
        main = "Masses racinaire chenes", ylab="Masse racinaire", xlab="") 

##Modele avec plan en tiroir
CC <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(CC) #Seulement Hauteur ini significatif --> On garde Hiniseul

##Modele avec plan en tiroir sans interaction Hini
CC <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(CC) #Effet Stress et Stress*Prov

##Verifier residus 
plot(CC) #Homogeneite des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(CC) ~ Stress*Brout*Prov, data = CHR_mean) #ok

qqnorm(resid(CC)) #Graphique de normalite ok
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
anova_CC1 <- anova(CC1)

lsmeans(CC1, pairwise~Stress*Prov)
#Difference entre Stress et No Stress de 2050

#Representation effet Stress*Prov (Log)
meansCHR_LMracine <- summarySE(LCHR_mean, measurevar = "LMracine", groupvars = c("Stress", "Prov")) 
meansCHR_LMracine <- mutate(meansCHR_LMracine, Stress=factor(Stress, levels=c("NoStress", "Stress2")))

ggplot(meansCHR_LMracine, aes(x = Stress, y = LMracine, shape=factor(Prov))) +
  geom_point(size=3, position=position_dodge(width=-0.3)) +
  geom_errorbar(aes(ymin=LMracine-se, ymax=LMracine+se), width=.1, position=position_dodge(width=-0.3))+
  geom_text(aes(x=Stress, y=LMracine+se+ 0.05),label = c("ab","a","ab","ab","b","ab"), position=position_dodge(width=-0.3))+
  theme_classic() +
  xlab("Traitement de stress hydrique") +  
  ylab("Log Masse racinaire (g)")

#Representation Stress*Prov avec les estimations du modele 
CHR_LMracine_mod <- summary(lsmeans(CC1, ~Prov*Stress))
ggplot(CHR_LMracine_mod, aes(x = Stress, y = lsmean, shape=factor(Prov))) +
  geom_point(size=3, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  labs(shape="Provenance") +
  scale_x_discrete(labels=c("Non stress?","Stress ?lev?"))+
  xlab("Traitement de stress hydrique") +  
  ylab("Log Masse racinaire (g)") +
  geom_text(aes(x=Stress, y=upper.CL+0.05),
            label = c("ab","a","ab","ab","b","ab"),  
            position=position_dodge(width=0.4)) 

##Verifier residus 
plot(CC1) #Homogeneite des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(CC1) ~ Stress*Brout*Prov, data = CHR_mean) #ok

qqnorm(resid(CC1)) #Graphique de normalite ok
qqline(resid(CC1))
shapiro.test(resid(CC1)) #ok


####MASSE AERIENNE####

#Boxplot masse racine
boxplot(Maerien~Brout+Stress+Prov,data=CHR,cex.axis=0.5, col = color,las=2, 
        main = "Masses aeriennes chenes", ylab="Masse aerienne", xlab="") 

##Modele avec plan en tiroir
DD <- lmerTest::lmer( Maerien ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(DD) #pas interaction avec Hini

##Modele avec plan en tiroir sans interaction Hauteurini
DD <- lmerTest::lmer(Maerien ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(DD)#Effet Hauteurini Tendance Prov et Tendance Stress*Prov
anova_DD <- anova(DD)

#Representation Tendance Stress*Prov
meansCHR_Maerien <- summarySE(CHR_mean, measurevar = "Maerien", groupvars = c("Stress", "Prov")) 
meansCHR_Maerien <- mutate(meansCHR_Maerien, Stress=factor(Stress, levels=c("NoStress", "Stress2")))
ggplot(meansCHR_Maerien, aes(x = Stress, y = Maerien, shape=factor(Prov))) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=Maerien-se, ymax=Maerien+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse aerienne (g)")

##Verifier residus 
plot(DD) #Homogeneite des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(DD) ~ Stress*Brout*Prov, data = CHR_mean) #ok

qqnorm(resid(DD)) #Graphique de normalite ok
qqline(resid(DD))
shapiro.test(resid(DD)) #ok

##Tables anova chenes
capture.output(anova_C, anova_D, anova_CC1, anova_DD, file="Chenes.doc")

#################ERABLES(PROV2018)#############
####BIOMASSE TOTALE ERABLES####

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

#Modele avec plan en tiroir
E <- lmerTest::lmer( Mtot ~ Stress*Brout*Hauteurini + (1|Bloc/Stress), data = ERS)
anova(E) # Rien de significatif

#Sans interaction avec Hini:
E <- lmerTest::lmer( Mtot ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = ERS)
anova(E) #Hauteurini non significative

#On retire Hini
E <- lmerTest::lmer( Mtot ~ Stress*Brout + (1|Bloc/Stress), data = ERS)
anova(E) #Stress*Brout significatif
anova_E <- anova(E)
lsmeans(E, pairwise~Stress*Brout)
#Difference entre Brout et NoBrout pour Stress1

#Representation Stress*Brout
meansERS_Mtot <- summarySE(ERS, measurevar = "Mtot", groupvars = c("Stress", "Brout")) 
meansERS_Mtot <- mutate(meansERS_Mtot, Stress=factor(Stress, levels=c("NoStress","Stress1","Stress2")))

ggplot(meansERS_Mtot, aes(x = Stress, y = Mtot, color=Brout)) +
  geom_point(size=3, position=position_dodge(width=-0.2)) +
  geom_errorbar(aes(ymin=Mtot-se, ymax=Mtot+se), width=.1, position=position_dodge(width=-0.2))+
  theme_classic() +
  geom_text(aes(x=Stress, y=Mtot+se+2),label = c("ab","ab","a","b","ab","ab"), position=position_dodge(width=-0.2))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse totale (g)")

#Representation Stress*Brout avec les estimations du modele 
ERS_Mtot_mod <- summary(lsmeans(E, ~Stress*Brout))

ERS_Mtot_mod$Brout <- as.vector(ERS_Mtot_mod$Brout) #get rid of factors
ERS_Mtot_mod$Brout = factor(ERS_Mtot_mod$Brout, levels=c("NoBrout", "Brout"))

ggplot(ERS_Mtot_mod, aes(x = Stress, y = lsmean, color=Brout)) +
  geom_point(size=4, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  theme(legend.title=element_blank()) + #Enlever titres des l?gendes
  scale_x_discrete(labels=c("Non stress?","Stress mod?r?", "Stress ?lev?"))+
  scale_y_continuous(limits = c(0,85), expand = expand_scale()) + 
  scale_color_manual(values=c("limegreen","lightcoral"), labels=c("Non Brout?", "Brout?"))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse (g)") +
  theme(text=element_text(size=18, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+ 
  geom_text(aes(x=Stress, y=upper.CL+2),
            label = c("ab","ab","a","b","ab","ab"),  
            position=position_dodge(width=0.4), size=4.5, show.legend=F)

##Verifier residus 
plot(E) #Homogeneite des variances --> ok
leveneTest(resid(E) ~ Stress*Brout, data = ERS) #ok

qqnorm(resid(E)) #Graphique de normalite ok
qqline(resid(E))
shapiro.test(resid(E)) #ok

####Ratio ERABLES####

#Boxplot ratio 
boxplot(Ratio~Brout+Stress,data=ERS,cex.axis=0.5,las=2, main = "Ratio Erables", 
        xlab="",ylab="Ratio")

##Modele avec plan en tiroir
G <- lmerTest::lmer( Ratio ~ Stress*Brout*Hauteurini + (1|Bloc/Stress), data = ERS)
anova(G) #Rien significatif

##Modele avec plan en tiroir sans interaction Hini
G <- lmerTest::lmer( Ratio ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = ERS)
anova(G) #Rien significatif

##Modele avec plan en tiroir sans Hini
G <- lmerTest::lmer( Ratio ~ Stress*Brout + (1|Bloc/Stress), data = ERS)
anova(G) #Rien significatif
anova_G <- anova(G)

##Verifier residus 
plot(G) #Homogeneite des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(G) ~ Stress*Brout, data = ERS) #ok

qqnorm(resid(G)) #Graphique de normalite ok
qqline(resid(G))
shapiro.test(resid(G)) #ok


####MASSE RACINAIRE####

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
EE <- lmerTest::lmer( Mracine ~ Stress*Brout + (1|Bloc/Stress), data = ERS)
anova(EE) #Stress*Brout significatif
anova_EE <- anova(EE)
lsmeans(EE, pairwise~Stress*Brout)

#Representation Stress*Brout
meansERS_Mracine <- summarySE(ERS, measurevar = "Mracine", groupvars = c("Stress", "Brout")) 
meansERS_Mracine <- mutate(meansERS_Mracine, Stress=factor(Stress, levels=c("NoStress","Stress1","Stress2")))
ggplot(meansERS_Mracine, aes(x = Stress, y = Mracine, color=Brout)) +
  geom_point(size=3, position=position_dodge(width=-0.3)) +
  geom_errorbar(aes(ymin=Mracine-se, ymax=Mracine+se), width=.1, position=position_dodge(width=-0.3))+
  geom_text(aes(x=Stress, y=Mracine+se+1),label = c("ab","ab","b","a","ab","ab"), position=position_dodge(width=-0.3))+
  theme_classic() +
  xlab("Traitement de stress hydrique") +  
  ylab("Masse racine (g)")

#Representation Stress*Brout avec les estimations du modele 
ERS_Mracine_mod <- summary(lsmeans(EE, ~Stress*Brout))
ggplot(ERS_Mracine_mod, aes(x = Stress, y = lsmean, color=Brout)) +
  geom_point(size=3, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  labs(colour="Broutement") +
  scale_x_discrete(labels=c("Non stress?","Stress mod?r?", "Stress ?lev?"))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse racinaire (g)") +
  geom_text(aes(x=Stress, y=upper.CL+2),
            label = c("ab","ab","b","a","ab","ab"),  
            position=position_dodge(width=0.4), show.legend=F) 

##Verifier residus 
plot(EE) #Homogeneite des variances, on ne doit pas voir de patron particulier (cone)
leveneTest(resid(EE) ~ Stress*Brout*Prov, data = ERS) #ok

qqnorm(resid(EE)) #Graphique de normalite ok
qqline(resid(EE))
shapiro.test(resid(EE)) #ok

  
####MASSE AERIENNE####
  
#Boxplot masse aerien
boxplot(Maerien~Brout+Stress,data=ERS,cex.axis=0.5, col = color,las=2, 
          main = "Masses aeriennes erables", ylab="Masse aerienne", xlab="") 

##Modele avec plan en tiroir
FF <- lmerTest::lmer( Maerien ~ Stress*Brout*Hauteurini + (1|Bloc/Stress), data = ERS)
anova(FF) #pas interaction avec Hini

##Modele sans interaction Hini
FF <- lmerTest::lmer( Maerien ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = ERS)
anova(FF) #Hauteurini significatif, effet Stress*Brout
anova_FF <- anova(FF)
lsmeans(FF, pairwise~Stress*Brout)
#Rien significatif

#Representation Stress*Brout
meansERS_Maerien <- summarySE(ERS, measurevar = "Maerien", groupvars = c("Stress", "Brout"))
meansERS_Maerien <- mutate(meansERS_Maerien, Stress=factor(Stress, levels=c("NoStress","Stress1","Stress2")))
ggplot(meansERS_Maerien, aes(x = Stress, y = Maerien, shape=factor(Brout))) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  theme_classic() +
  geom_errorbar(aes(ymin=Maerien-se, ymax=Maerien+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse aerienne (g)")

#Representation Stress*Brout avec les estimations du modele 
ERS_Maerien_mod <- summary(lsmeans(FF, ~Stress*Brout))
ggplot(ERS_Maerien_mod, aes(x = Stress, y = lsmean, color=Brout)) +
  geom_point(size=3, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  labs(colour="Broutement") +
  scale_x_discrete(labels=c("Non stress?","Stress mod?r?", "Stress ?lev?"))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse aerienne estimations modele (g)")

##Verifier residus 
plot(FF) #Homogeneite des variances, on ne doit pas voir de patron particulier (cone)
leveneTest(resid(FF) ~ Stress*Brout*Prov, data = ERS) #ok

qqnorm(resid(FF)) #Graphique de normalite,pas certaine
qqline(resid(FF))
shapiro.test(resid(FF)) #ok

##Tables anova erables
capture.output(anova_E, anova_G, anova_EE, anova_FF, file="Erables.doc")


##Surcompensation - Nb de ramilles
ERS_Brout<-filter(ERS,Brout=="Brout")
ERS_NoBrout<-filter(ERS,Brout=="NoBrout")
mean(ERS_NoBrout$Ram) #2,5
mean(ERS_Brout$Ram) #1,5
t.test(ERS_Brout$Ram,ERS_NoBrout$Ram, var.equal=TRUE)
#NON, finalement c'est ceux non brout?s qui ont plus de ramille
#diff?rence de 1, c'est probablement juste la tige apicale compt?e de plus chez les noBrout

#################THUYAS########################
####BIOMASSE TOTALE THUYAS####

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

#Modele avec plan en tiroir
H <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = THO)
anova(H) # Interaction stress*brout*Hauteurini, Hini significatif
#et stress*brout de significatif

##Verifier residus en incluant les effets aleatoires (effet du design)
plot(H) #Homogeneite des variances, on ne doit pas voir de patron particulier (cone) --> ok
leveneTest(resid(H) ~ Stress*Brout*Prov, data = THO) #ok

qqnorm(resid(H)) #Graphique de normalite, pas certaine
qqline(resid(H))
shapiro.test(resid(H)) #Non, essayer transformation log

#Ajouter variable Transformation log de Mtot
THO <- mutate(THO, LMtot=log(Mtot))

#Modele avec plan en tiroir (log)
H1 <- lmerTest::lmer( LMtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = THO)
anova(H1) #Interactions avec Hini significatives, Hini significatif
anova_H1 <- anova(H1)
#Et Prov et Brout
lsmeans(H1, ~Brout)

#Representation effet Brout
meansTHO_LMtot <- summarySE(THO, measurevar = "LMtot", groupvars = "Brout") 
ggplot(meansTHO_LMtot, aes(x = Brout, y = LMtot)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=LMtot-se, ymax=LMtot+se), width=.1)+
  xlab("Broutement") +  
  ylab("Log Masse totale (g)")
ls_means(H1, which = "Brout", pairwise = TRUE)

#Effet Brout selon les estim?s du mod?le
summary(lsmeans(H1, ~Brout))

#Effet Prov
ls_means(H1, which = "Prov", pairwise = TRUE)
#Difference 2050 et 2080

#Representation effet Prov
meansTHO_LMtot <- summarySE(THO, measurevar = "LMtot", groupvars = "Prov") 
meansTHO_LMtot <- mutate(meansTHO_LMtot, Prov=factor(Prov, levels=c("2018", "2050","2080")))
ggplot(meansTHO_LMtot, aes(x = Prov, y = LMtot)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=LMtot-se, ymax=LMtot+se), width=.1)+
  theme_classic() +
  xlab("Provenance") +  
  ylab("Log Masse totale (g)") +
  geom_text(aes(x=Prov, y=LMtot+se+0.03),label = c("","",""))
#Pourtant ce n'est pas ce qu'on observe dans le graph (2018 semble beaucoup plus ?lev? que 2050 et 2080)
#Le test statistique n'est pas r?alis? sur les donn?es brutes, mais sur les donn?es du mod?le
#Et en plus, ce mod?le utilise les valeurs log. 
#Si on fait :
ls_means(H1, which = "Prov")
#on voit les moyennes de LS means pour chacun des niveaux de traitement
#Le r?sultat est bien moins surprenant
#C'est la diff?rence entre pr?senter les moyennes des mod?les vs les donn?es brutes. 
#C'est un cas assez extr?me, mais au lieu des donn?es brutes, 
#essayer plut?t de pr?senter les estim?s de la table ls means, 
#avec leur intervalle de confiance (colonne upper & lower).

#Representation effet Prov avec les estim?s du mod?le
THO_Mtot_mod <- summary(lsmeans(H1, ~Prov))
ggplot(THO_Mtot_mod, aes(x = Prov, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  theme_classic() +
  xlab("Provenance") +  
  ylab("Masse totale (g)") +
  geom_text(aes(x=Prov, y=upper.CL+0.03),label = c("ab","a","b"))

##Verifier residus en incluant les effets aleatoires (effet du design)
plot(H1) #Homogeneite des variances, on ne doit pas voir de patron particulier (cone) --> ok
leveneTest(resid(H1) ~ Stress*Brout*Prov, data = THO) #ok

qqnorm(resid(H1)) #Graphique de normalite, pas certaine
qqline(resid(H1))
shapiro.test(resid(H1)) #ok

####Ratio THUYAS####

#Boxplot ratio 
boxplot(Ratio~Brout+Stress+Prov,data=THO,cex.axis=0.5, col = color,las=2, main = "Ratio thuyas", 
        xlab="",ylab="Masse tot")

##Modele avec plan en tiroir
I <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = THO)
anova(I) #Pas interaction avec Hini

##Modele avec plan en tiroir sans interaction Hini
I <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = THO)
anova(I) #Hini pas significatif

##Modele avec plan en tiroir sans Hini
I <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = THO)
anova(I) #Tendance Stress, Interaction Stress*Brout*Prov significatif
anova_I <- anova(I)

lsmeans(I, pairwise~Stress*Brout*Prov)
#Rien de significatif

#Representation Stress*Brout*Prov
meansTHO_ratio <- summarySE(THO, measurevar = "Ratio", groupvars = c("Stress", "Brout", "Prov")) 
meansTHO_ratio <- mutate(meansTHO_ratio, Stress=factor(Stress, levels=c("NoStress","Stress1","Stress2")))
meansTHO_ratio <- mutate(meansTHO_ratio, Prov=factor(Prov,levels=c("2018","2050","2080")))

ggplot(meansTHO_ratio, aes(x = Stress, y = Ratio, shape=factor(Prov),color=Brout)) +
  geom_point(size=3, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=Ratio-se, ymax=Ratio+se), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  xlab("Traitement de stress hydrique") +  
  ylab("Ratio masse a?rienne/racinaire")

#Representation Stress*Brout*Prov avec les estimations du modele
THO_Ratio_mod <- summary(lsmeans(I, ~Prov*Stress*Brout))
ggplot(THO_Ratio_mod, aes(x = Stress, y = lsmean, shape=factor(Prov),color=Brout)) +
  geom_point(size=3, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  labs(shape="Provenance", colour="Broutement") +
  scale_x_discrete(labels=c("Non stress?","Stress mod?r?", "Stress ?lev?"))+
  scale_color_discrete(labels=c("Brout?", "Non brout?")) +
  xlab("Stress hydrique") +  
  ylab("Ratio de biomasse a?rienne/racinaire")

#Graphique pour les broutes et les non broutes (donnes brutes) 
ggplot(meansTHO_ratio[which(meansTHO_ratio$Brout == "NoBrout"),], aes(x = Stress, y = Ratio, shape=factor(Prov))) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=Ratio-se, ymax=Ratio+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de stress hydrique") +  
  ylab("Ratio aerienne/racinaire")+
  ylim(c(2,5.5))+
  annotate("text", x = 1, y = 5,label = "Non broutes", size = 5)
ggplot(meansTHO_ratio[which(meansTHO_ratio$Brout == "Brout"),], aes(x = Stress, y = Ratio, shape=factor(Prov))) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=Ratio-se, ymax=Ratio+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de stress hydrique") +  
  ylab("Ratio aerienne/racinaire")+
  ylim(c(2,5.5))+
  annotate("text", x = 1, y = 5,label = "Broutes", size = 5)

##Verifier residus 
plot(I) #Homogeneite des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(I) ~ Stress*Brout*Prov, data = THO) #ok

qqnorm(resid(I)) #Graphique de normalite ok
qqline(resid(I))
shapiro.test(resid(I)) #ok

####MASSE RACINAIRE####

#Boxplot masse racine
boxplot(Mracine~Brout+Stress+Prov,data=THO,cex.axis=0.5, col = color,las=2, 
        main = "Masses racinaire thuyas", ylab="Masse racinaire", xlab="") 

##Modele avec plan en tiroir
II <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = THO)
anova(II) #Pas interaction avec Hauteur ini 

##Modele avec plan en tiroir sans interaction Hini
II <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov+Hauteurini + (1|Bloc/Stress), data = THO)
#Erreur de convergence, on test un optimizer different
II <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = THO,
                      control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(II) #Hauteur ini significatif, Effet Stress et Brout
anova_II <-anova(II)

lsmeans(II, ~Stress)
ls_means(II, which = "Stress", pairwise = TRUE)
#NoStress differe de Stress1
#NoStress differe de Stress2

#Representation effet Brout
meansTHO_Mracine <- summarySE(THO, measurevar = "Mracine", groupvars = c("Brout")) 
ggplot(meansTHO_Mracine, aes(x = Brout, y = Mracine)) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=Mracine-se, ymax=Mracine+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de broutement") +  
  ylab("Masse racinaire (g)")
#Plants broutes ont masse racinaire plus faible
ls_means(II, which = "Brout", pairwise = TRUE)
summary(lsmeans(II, ~Brout))

#Representation effet Stress
meansTHO_Mracine <- summarySE(THO, measurevar = "Mracine", groupvars = c("Stress")) 
ggplot(meansTHO_Mracine, aes(x = Stress, y = Mracine)) +
  theme_classic() +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=Mracine-se, ymax=Mracine+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse racinaire (g)")+
  geom_text(aes(x=Stress, y=Mracine+se+0.5),label = c("a","b","b")) 
  
#Representation Stress avec les estim?s du mod?le
THO_Mracine_mod <- summary(lsmeans(II, ~Stress))
ggplot(THO_Mracine_mod, aes(x = Stress, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  theme_classic() +
  scale_x_discrete(labels=c("Non stress?","Stress mod?r?", "Stress ?lev?"))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse racinaire (g)") +
  geom_text(aes(x=Stress, y=upper.CL+0.5),label = c("a","b","b"))

##Verifier residus 
plot(II) #Homogeneite des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(II) ~ Stress*Brout*Prov, data = THO) #ok

qqnorm(resid(II)) #Graphique de normalite, pas certaine
qqline(resid(II))
shapiro.test(resid(II)) #0.06, ok

####MASSE AERIENNE####

#Boxplot masse aerien
boxplot(Maerien~Brout+Stress+Prov,data=THO,cex.axis=0.5, col = color,las=2, 
          main = "Masses aeriennes thuyas", ylab="Masse aerienne", xlab="") 

##Modele avec plan en tiroir
JJ <- lmerTest::lmer( Maerien ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = THO)
anova(JJ) #On garde Hiniseul

##Modele avec plan en tiroir sans interaction Hini
JJ <- lmerTest::lmer( Maerien ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = THO)
anova(JJ) #Effet Hini, Stress et Brout
anova_JJ <- anova(JJ)
lsmeans(JJ, ~Stress)
#Masse aerienne plus faible pour les plants boutes
ls_means(JJ, which = "Brout", pairwise = TRUE) 

ls_means(JJ, which = "Stress", pairwise = TRUE) 
#Stress2 plus faible que NoStress
#Stress2 plus faible que Stress1

#Representation effet stress
meansTHO_Maerien <- summarySE(THO, measurevar = "Maerien", groupvars = c("Stress")) 
meansTHO_Maerien <- mutate(meansTHO_Maerien, Stress=factor(Stress, levels=c("NoStress","Stress1","Stress2")))
ggplot(meansTHO_Maerien, aes(x = Stress, y = Maerien)) +
  geom_point(size=3) +
  theme_classic() +
  geom_errorbar(aes(ymin=Maerien-se, ymax=Maerien+se), width=.1, position=position_dodge(width=0.2))+
  geom_text(aes(x=Stress, y=Maerien+se+1),label = c("a","a","b")) +
  xlab("Traitement de stress hydrique") +  
  ylab("Masse aerienne (g)")

#Representation Stress avec les estim?s du mod?le
THO_Maerien_mod <- summary(lsmeans(JJ, ~Stress))
ggplot(THO_Maerien_mod, aes(x = Stress, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  theme_classic() +
  scale_x_discrete(labels=c("Non stress?","Stress mod?r?", "Stress ?lev?"))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse a?rienne (g)") +
  geom_text(aes(x=Stress, y=upper.CL+0.5),label = c("a","a","b"))

#Representation effet brout
meansTHO_Maerien <- summarySE(THO, measurevar = "Maerien", groupvars = c("Brout"))
meansTHO_Maerien <- mutate(meansTHO_Maerien, Brout=factor(Brout, levels=c("NoBrout","Brout")))
ggplot(meansTHO_Maerien, aes(x = Brout, y = Maerien)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Maerien-se, ymax=Maerien+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de broutement") +  
  ylab("Masse aerienne (g)")
summary(lsmeans(JJ, ~Brout))

##Verifier residus 
plot(JJ) #Homogeneite des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(JJ) ~ Stress*Brout*Prov, data = THO) #ok

qqnorm(resid(JJ)) #Graphique de normalite ok
qqline(resid(JJ))
shapiro.test(resid(JJ)) #ok


##COMPARAISON AVEC BIOMASSE RETIREE##

THOavecMret <- mutate(THO, Mtotret=Mtot+Mret, Maerien=Maerien+Mret)
summarySE(THO, measurevar = "Mtot", groupvars = c("Brout")) 

##Modele avec plan en tiroir pour Maerien (ajout de Mret)
JJ1 <- lmerTest::lmer( Maerien ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = THOavecMret)
anova(JJ1) #Pas d'interaction avec Hini significative

##Modele avec plan en tiroir sans interaction Hini pour Maerien (ajout de Mret)
JJ1 <- lmerTest::lmer( Maerien ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = THOavecMret)
anova(JJ1) #Hini significatif, Brout significatif
summary(JJ1)
meansTHO_Maerienmret <- summarySE(THOavecMret, measurevar = "Maerien", groupvars = c("Brout")) 
#Malgre l'ajout de Mret, ceux broutes ont une Maerienne significativement plus basse
#Leur croissance a donc ete ralentie suite au traitement de broutement
lsmeans(JJ1, ~Brout)

##TEST Graphique Biomasse aerienne et racinaire et totale
#Donn?es Mtot log retransformer pour comparer avec racine et a?rien
THO_Mtot_mod <- summary(lsmeans(H1, ~Stress))
#exp(1) est valeur de e (du log naturel)
THO_Mtot_mod <- mutate(THO_Mtot_mod, lsmean=exp(lsmean),upper.CL=exp(upper.CL),lower.CL=exp(lower.CL))

#Fusionner 3 dataframes
THO_Mtot_mod <- add_column(THO_Mtot_mod,Biomasse="Mtot", .after = "Stress")
THO_Maerien_mod <- add_column(THO_Maerien_mod,Biomasse="Maerien", .after = "Stress")
THO_Mracine_mod <- add_column(THO_Mracine_mod,Biomasse="Mracine", .after = "Stress")

THO_Biomasse_mod <- dplyr::union(THO_Mtot_mod,THO_Maerien_mod)
THO_Biomasse_mod <- dplyr::union(THO_Biomasse_mod,THO_Mracine_mod)

ggplot(THO_Biomasse_mod, aes(x = Stress, y = lsmean, color=Biomasse)) +
  geom_point(size=2, position=position_dodge(width=-0.2))+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2, position=position_dodge(width=-0.2))+
  theme_classic() +
  scale_color_manual(values=c("limegreen","darkorange", "grey40"),
                     labels=c("A?rienne", "Racinaire", "Totale"))+
  scale_x_discrete(labels=c("Non stress?","Stress mod?r?", "Stress ?lev?"))+
  scale_y_continuous(limits = c(0,85), expand = expand_scale()) +
  theme(text=element_text(size=12, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+ 
  xlab("Traitement de stress hydrique") +  
  ylab("Masse (g)")+
  geom_text(aes(x=c(1.06,2.06,3.06), y=upper.CL+2),data = THO_Maerien_mod,label = c("A","B","B"), show.legend=F) +
  geom_text(aes(x=Stress, y=upper.CL+2),data = THO_Mracine_mod, label = c("a","a","b"), show.legend=F)  

#Graphique Maerien et racinaire seulement
THO_Biomasse_mod2 <- dplyr::union(THO_Mracine_mod,THO_Maerien_mod)

ggplot(THO_Biomasse_mod2, aes(x = Stress, y = lsmean, color=Biomasse)) +
  geom_point(size=2, position=position_dodge(width=-0.2))+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2, position=position_dodge(width=-0.2))+
  theme_classic() +
  scale_color_manual(values=c("limegreen","darkorange"),
                     labels=c("A?rienne", "Racinaire"))+
  scale_x_discrete(labels=c("Non stress?","Stress mod?r?", "Stress ?lev?"))+
  scale_y_continuous(limits = c(0,50), expand = expand_scale()) +
  theme(text=element_text(size=12, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+ 
  xlab("Traitement de stress hydrique") +  
  ylab("Masse (g)")+
  geom_text(aes(x=c(1.05,2.05,3.05), y=upper.CL+1.5),data = THO_Maerien_mod,label = c("A","B","B"), show.legend=F) +
  geom_text(aes(x=c(0.95,1.95,2.05), y=upper.CL+1.5),data = THO_Mracine_mod, label = c("a","a","b"), show.legend=F)  



##Tables anova Thuyas
capture.output(anova_H1, anova_I, anova_II, anova_JJ, file="Thuyas.doc")


#################PINS########################
####BIOMASSE TOTALE PINS####

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

#Modele avec plan en tiroir
J <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = PIB)
anova(J) #Rien significatif

#Sans interaction avec Hini:
J <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = PIB)
anova(J) #Hauteurini, Brout significatif
summary(J)
anova_J <- anova(J)
lsmeans(J, ~Brout)

#Representation effet brout
meansPIB_Mtot <- summarySE(PIB, measurevar = "Mtot", groupvars = c("Brout")) 
meansPIB_Mtot <- mutate(meansPIB_Mtot, Brout=factor(Brout, levels=c("NoBrout","Brout")))
ggplot(meansPIB_Mtot, aes(x = Brout, y = Mtot)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Mtot-se, ymax=Mtot+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de broutement") +  
  ylab("Masse totale (g)")
#Biomasse totale diminue avec le broutement


#Representation avec donn?es du modele
PIB_Mtot_mod <- summary(lsmeans(J, ~Brout))

PIB_Mtot_mod$Brout <- as.vector(PIB_Mtot_mod$Brout) #get rid of factors
PIB_Mtot_mod$Brout = factor(PIB_Mtot_mod$Brout, levels=c("NoBrout", "Brout"))

ggplot(PIB_Mtot_mod, aes(x = Brout, y = lsmean)) +
  geom_point(size=4, colour=c("lightcoral","limegreen"))+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL),size=0.7, width=.1, colour=c("lightcoral","limegreen"))+
  theme_classic() +
  ylab("Masse (g)")+
  xlab("Traitement de broutement") +
  theme(text=element_text(size=18, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+ 
  scale_x_discrete(labels=c("Non brout?","Brout?"))+
  scale_y_continuous(limits = c(0,40), expand = expand_scale())  #Enlever espace entre z?ro et axe
  

##Verifier residus en incluant les effets aleatoires (effet du design)
plot(J) #Homogeneite des variances, on ne doit pas voir de patron particulier (cône) --> ok
leveneTest(resid(J) ~ Stress*Brout*Prov, data = PIB) #ok

qqnorm(resid(J)) #Graphique de normalite ok
qqline(resid(J))
shapiro.test(resid(J)) #ok

####Ratio PINS####

#Boxplot ratio 
boxplot(Ratio~Brout+Stress+Prov,data=PIB,cex.axis=0.5, col = color,las=2, main = "Ratio pins", 
        xlab="",ylab="Masse tot")

##Modele avec plan en tiroir
K <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = PIB)
anova(K) #Rien significatif 

##Modele avec plan en tiroir sans interaction Hini
K <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = PIB)
anova(K) #Rien significatif

##Modele avec plan en tiroir sans Hini
K <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = PIB)
anova(K) #Rien significatif
anova_K <- anova(K)
lsmeans(K, ~Brout)

##Verifier residus 
plot(K) #Homogeneite des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(K) ~ Stress*Brout*Prov, data = PIB) #ok

qqnorm(resid(K)) #Graphique de normalite --> ok
qqline(resid(K))
shapiro.test(resid(K)) #ok


####MASSE RACINAIRE####

#Boxplot masse racine
boxplot(Mracine~Brout+Stress+Prov,data=PIB,cex.axis=0.5, col = color,las=2, 
        main = "Masses racinaire pins", ylab="Masse racinaire", xlab="") 

##Modele avec plan en tiroir
KK <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = PIB)
anova(KK) #Rien significatif

##Modele avec plan en tiroir sans interaction Hini
KK <- lmerTest::lmer( Mracine ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = PIB)
anova(KK) #Hauteurini significatif, Effet Brout

##Verifier residus 
plot(KK) #Homogeneite des variances, on ne doit pas voir de patron particulier (cone)
leveneTest(resid(KK) ~ Stress*Brout*Prov, data = PIB) #ok

qqnorm(resid(KK)) #Graphique de normalite ne semble pas normal
qqline(resid(KK))
shapiro.test(resid(KK)) #Non transformation a faire

#Ajouter variable Transformation log de Mracine
PIB <- mutate(PIB, LMracine=log(Mracine))

##Modele avec plan en tiroir Log (Transformation log pour Mracine)
KK1 <- lmerTest::lmer( LMracine ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = PIB)
anova(KK1) #On garde ce modele, Interactions avec Hini significatives

##Verifier residus (log)
plot(KK1) #Homogeneite des variances, on ne doit pas voir de patron particulier (cone)
leveneTest(resid(KK1) ~ Stress*Brout*Prov, data = PIB) #ok

qqnorm(resid(KK1)) #Mieux mais ne semble pas encore normal 
qqline(resid(KK1))
shapiro.test(resid(KK1)) #Non, essayer une autre transformation

#Ajouter variable Transformation racine carree
PIB <- mutate(PIB, rMracine=sqrt(Mracine))

##Modele avec plan en tiroir (Transformation racine carree pour Mracine)
KK2 <- lmerTest::lmer( rMracine ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = PIB)
anova(KK2) #Rien significatif

##Modele avec plan en tiroir sans interaction Hini(Transformation racine carree pour Mracine)
KK2 <- lmerTest::lmer( rMracine ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = PIB)
anova(KK2) #Hini pas significatif

##On enleve Hini(Transformation racine carree pour Mracine)
KK2 <- lmerTest::lmer( rMracine ~ Stress*Brout*Prov + (1|Bloc/Stress), data = PIB)
anova(KK2) #Brout significatif
summary(KK2)
anova_KK2 <- anova(KK2)

#Representation effet brout
meansPIB_rMracine <- summarySE(PIB, measurevar = "rMracine", groupvars = c("Brout")) 
meansPIB_rMracine <- mutate(meansPIB_rMracine, Brout=factor(Brout, levels=c("NoBrout","Brout")))
ggplot(meansPIB_rMracine, aes(x = Brout, y = rMracine)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=rMracine-se, ymax=rMracine+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de broutement") +  
  ylab("sqrt Masse racinaire (g)")
#Biomasse racinaire diminue avec le broutement

#Representation Brout avec les estimations du modele (racine carree)
PIB_Mracine_mod <- summary(lsmeans(KK2, ~Brout))
PIB_Mracine_mod <- mutate(PIB_Mracine_mod, Brout=factor(Brout, levels=c("NoBrout","Brout")))
ggplot(PIB_Mracine_mod, aes(x = Brout, y = lsmean, color=Brout)) +
  geom_point(size=3, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  ylim(0,4)+
  scale_x_discrete(labels=c("Brout?","Non Brout?"))+
  xlab("Broutement") +  
  ylab("Masse (g)") 

##Verifier residus (racine caree)
plot(KK2) #Homogeneite des variances, on ne doit pas voir de patron particulier (cone)
leveneTest(resid(KK2) ~ Stress*Brout*Prov, data = PIB) #ok

qqnorm(resid(KK2)) #Graphique de normalite, pas certaine??
qqline(resid(KK2))
shapiro.test(resid(KK2)) #Ok
  
####MASSE AERIENNE####
  
#Boxplot masse aerien
boxplot(Maerien~Brout+Stress+Prov,data=PIB,cex.axis=0.5, col = color,las=2, 
          main = "Masses aeriennes pins", ylab="Masse aerienne", xlab="") 

##Modele avec plan en tiroir
LL <- lmerTest::lmer( Maerien ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = PIB)
anova(LL) #Rien significatif

##Modele avec plan en tiroir sans interaction Hini
LL <- lmerTest::lmer( Maerien ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = PIB)
anova(LL) #Hauteurini significatif, Effet Brout
summary(LL)
anova_LL <- anova(LL)
lsmeans(LL, ~Brout)

#Representation effet brout
meansPIB_Maerien <- summarySE(PIB, measurevar = "Maerien", groupvars = c("Brout")) 
meansPIB_Maerien <- mutate(meansPIB_Maerien, Brout=factor(Brout, levels=c("NoBrout","Brout")))
ggplot(meansPIB_Maerien, aes(x = Brout, y = Maerien)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=Maerien-se, ymax=Maerien+se), width=.1, position=position_dodge(width=0.2))+
  xlab("Traitement de broutement") +  
  ylab("Masse aerienne (g)")
#Biomasse aerienne diminue avec le broutement

#Representation Brout avec les estimations du modele 
PIB_Maerien_mod <- summary(lsmeans(LL, ~Brout))
PIB_Maerien_mod <- mutate(PIB_Maerien_mod, Brout=factor(Brout, levels=c("NoBrout","Brout")))
ggplot(PIB_Maerien_mod, aes(x = Brout, y = lsmean, color=Brout)) +
  geom_point(size=3, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  scale_x_discrete(labels=c("Brout?","Non Brout?"))+
  xlab("Broutement") +  
  ylab("Masse (g)") 


##Verifier residus 
plot(LL) #Homogeneite des variances, on ne doit pas voir de patron particulier (cône)
leveneTest(resid(LL) ~ Stress*Brout*Prov, data = PIB) #ok

qqnorm(resid(LL)) #Graphique de normalite ok
qqline(resid(LL))
shapiro.test(resid(LL)) #ok

##COMPARAISON AVEC BIOMASSE RETIREE##
PIBavecMret <- mutate(PIB, Mtotret=Mtot+Mret, Maerien=Maerien+Mret)
summarySE(PIB, measurevar = "Mtot", groupvars = c("Brout")) 

##Modele avec plan en tiroir pour Maerien (ajout de Mret)
LL1 <- lmerTest::lmer( Maerien ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = PIBavecMret)
anova(LL1) #Rien significatif

##Modele avec plan en tiroir sans interaction Hini pour Maerien (ajout de Mret)
LL1 <- lmerTest::lmer( Maerien ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = PIBavecMret)
anova(LL1) #Hauteurini significatif, Effet Brout
summary(LL1)
#Malgre l'ajout de Mret, ceux broutes ont une Maerienne significativement plus basse
#Leur croissance a donc ete ralentie suite au traitement de broutement
lsmeans(LL1, ~Brout)
ls_means(LL1, which = "Brout", pairwise = TRUE)
meansPIB_Maerienmret <- summarySE(PIBavecMret, measurevar = "Maerien", groupvars = c("Brout")) 


##TEST Graphique Biomasse racinaire et aerien
#Donn?es Mracine racine carree retransformer pour comparer avec a?rien
PIB_Mracine_mod <- summary(lsmeans(KK2, ~Brout))
PIB_Mracine_mod <- mutate(PIB_Mracine_mod, lsmean=(lsmean)^2,upper.CL=(upper.CL)^2,lower.CL=(lower.CL)^2,
                          Brout=factor(Brout, levels=c("NoBrout","Brout")))
PIB_Maerien_mod <- summary(lsmeans(LL, ~Brout))
PIB_Maerien_mod <- mutate(PIB_Maerien_mod, Brout=factor(Brout, levels=c("NoBrout","Brout")))


#Fusionner 2 dataframes
PIB_Maerien_mod <- add_column(PIB_Maerien_mod,Biomasse="A?rienne")
PIB_Mracine_mod <- add_column(PIB_Mracine_mod,Biomasse="Racinaire")
PIB_Biomasse_mod <- dplyr::union(PIB_Mracine_mod,PIB_Maerien_mod)


ggplot(PIB_Biomasse_mod, aes(x = Brout, y = lsmean, color=Biomasse)) +
  geom_point(size=2.6, position=position_dodge(width=-0.2))+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL),size=0.7, width=.1, position=position_dodge(width=-0.2))+
  theme_classic() +
  scale_color_manual(values=c("limegreen", "darkorange"))+
  ylab("Masse (g)")+
  xlab("Traitement de broutement") +
  theme(text=element_text(size=12, family="sans"), #Police Arial ("serif" pour Time New Roman)
                          panel.border = element_rect(colour = "black", fill=NA, size=0.5))+ 
  scale_x_discrete(labels=c("Non brout?","Brout?"))+
  scale_y_continuous(limits = c(0,24.5), expand = expand_scale()) + #Enlever espace entre z?ro et axe
  geom_text(aes(x=c(1.95,0.95), y=upper.CL+0.7),data = PIB_Mracine_mod,label = c("b","a"), show.legend=F) +
  geom_text(aes(x=c(2.05,1.05), y=upper.CL+0.7),size=3.2,data = PIB_Maerien_mod, label = c("B","A"), show.legend=F) 



##Tables anova Pins
capture.output(anova_J, anova_K, anova_KK2, anova_LL, file="Pins.doc")
