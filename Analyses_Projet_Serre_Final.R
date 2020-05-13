
#2019-2020
#Projet Serre / Broutement / Stress_hydrique / Migration assistee

remove(list = ls())

#####Importation des donnees####
library(dplyr)
library(ggplot2)
library(plyr) #pour fonction ddply 
library(car) #pour anovas type III
library(lsmeans) #fonction lsmeans
library(MASS) #boxcox

Serre <- read.table("./Masses3.txt", header=TRUE, na.string = "", dec = ",") 

Serre$Brout <- replace(Serre$Brout,Serre$Brout=="0","NoBrout")
Serre$Brout <- replace(Serre$Brout,Serre$Brout=="1","Brout")
Serre$Stress <- replace(Serre$Stress,Serre$Stress=="0","NoStress")
Serre$Stress <- replace(Serre$Stress,Serre$Stress=="1","Stress1")
Serre$Stress <- replace(Serre$Stress,Serre$Stress=="2","Stress2")

Serre$Prov <- as.factor(as.character(Serre$Prov)) #Definir comme facteur
Serre$Bloc <- as.factor(as.character(Serre$Bloc))
Serre$Brout <- as.factor(as.character(Serre$Brout))
Serre$Stress <- as.factor(as.character(Serre$Stress))

str(Serre)
Serre$Brout <- factor(Serre$Brout, levels=c("NoBrout","Brout")) #Reorganiser niveaux facteur Brout


CET<-filter(Serre,Esp=="CET")
CET$Stress <- factor(CET$Stress) #Enlever le niveau Stress 1 absent de CET
CHR<-filter(Serre,Esp=="CHR")
CHR$Stress <- factor(CHR$Stress) #Enlever le niveau Stress 1 absent de CHR
ERS<-filter(Serre,Esp=="ERS")
THO<-filter(Serre,Esp=="THO")
PIB<-filter(Serre,Esp=="PIB")

SerreTemoin<-filter(Serre,Brout=="NoBrout",Stress=="NoStress",Prov=="2018")



#####Exploration des donnees####
#Outliers

boxplot(CET$Mtot, xlab="Masse totale Cerisier")

boxplot(CHR$Mtot, xlab="Masse totale Chene")
dotchart(CHR$Mtot,xlab="Masse totale Chene", ylab = "Ordre des donnees") #Un outlier? Semble ok dans les donnees

boxplot(ERS$Mtot, xlab="Masse totale Erable")

boxplot(PIB$Mtot, xlab="Masse totale Pin")

boxplot(THO$Mtot, xlab="Masse totale Thuya")
dotchart(THO$Mtot,xlab="Masse totale Thuya", ylab = "Ordre des donnees") #Pas d'outlier important

#Normalite Mtot
hist(CET$Mtot)
hist(CHR$Mtot)
hist(ERS$Mtot)
hist(PIB$Mtot)
hist(THO$Mtot)

#################CERISIERS########################
####BIOMASSE TOTALE CERISIERS####

#Regarder si hauteur initiale est un bon estimateur de la masse totale
#Regarder la relation entre hauteur finale et la masse a la fin de l'experience
cor.test(CET$Hauteur, CET$Mtot) #0,63
plot(CET$Hauteur, CET$Mtot) 
#Oui, bonne correlation, on utilise la hauteur initiale pour corriger pour la biomasse initiale des plants

#Savoir le nb de replicats 
ddply(CET, c("Stress", "Brout", "Prov", "Bloc"), summarise, N = length(Mtot))

#Pour la majorite des combinaisons, nous n'avons qu'un replicat
#Pour simplifier les analyses, on fait une moyenne pour les rares cas 
#ou nous avons 2 plants par combinaison

CET_mean <- aggregate(list(CET$Mtot,CET$Hauteurini,CET$Ratio,CET$Maerien,CET$Mracine), by= data.frame(CET$Stress, CET$Brout, CET$Prov, CET$Bloc), FUN= "mean")
colnames(CET_mean)[1:9] <- c("Stress", "Brout", "Prov", "Bloc","Mtot","Hauteurini","Ratio","Maerien","Mracine")

ddply(CET_mean, c("Stress", "Brout", "Prov", "Bloc","Hauteurini"), summarise,
      N = length(Mtot)) #Verification n = 1 pour chaque combinaison

#Boxplot masses tot 
color1 = c(rep("green",2),rep("red",2))
boxplot(Mtot~Brout+Stress+Prov,data=CET,cex.axis=0.5, col = color1,las=2, 
        main = "Masses tot cerisiers", ylab="Masse tot", xlab="") 

#Modele avec plan en tiroir
mod_CET_mtot <- lmerTest::lmer(Mtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CET_mean)
#Message veut dire que l'effet aleatoire du bloc est tr?s pr?s de 0 ou nul
#On le garde requand meme dans le modele, mais il a peu d'effet
summary(mod_CET_mtot)
Anova(mod_CET_mtot, type=3) #Interaction significative avec Hini, on garde ce modele
#Effet Stress*Brout*Prov

#Test a posteriori
summary(lsmeans(mod_CET_mtot, pairwise~Stress*Brout*Prov))
#Rien n'est significatif dans le test ? posteriori

#Representation Stress*Brout*Prov avec les estimes du modele
CET_Mtot_mod <- summary(lsmeans(mod_CET_mtot, ~Stress*Brout*Prov))
ggplot(CET_Mtot_mod, aes(x = Stress, y = lsmean, color=Brout, shape=Prov)) +
  geom_point(size=4, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.5))+
  theme_classic() +
  scale_shape_manual(values=c(10, 19, 1))+
  scale_x_discrete(labels=c("Non stress?","Stress ?lev?"))+
  scale_y_continuous(limits = c(0,125), expand = expand_scale()) +
  xlab("Traitement de stress hydrique") +  
  ylab("Masse (g)") +
  theme(text=element_text(size=12, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+ 
  geom_text(aes(x=Stress, y=upper.CL+2),
            label = c("","","","","","","","","","","",""),  
            position=position_dodge(width=0.4), size=4.5, show.legend=F)


##Verifier assomptions de modeles
plot(mod_CET_mtot) #Homogeneite des variances ok
leveneTest(resid(mod_CET_mtot) ~ Stress*Brout*Prov, data = CET_mean) #ok

qqPlot(resid(mod_CET_mtot)) #Graphique de normalite ok
shapiro.test(resid(mod_CET_mtot)) #On verifie avec un test ok


#################CHENES########################
####BIOMASSE TOTALE CHENES####

#Est-ce que la hauteur initiale est un bon estimateur de la masse totale?
#Regarder la relation entre hauteur finale et la masse a la fin de l'experience
cor.test(CHR$Hauteur, CHR$Mtot) #0,67
plot(CHR$Hauteur, CHR$Mtot) 
#Oui, bonne correlation, utiliser hauteurini pour corriger biomasse initiale des plants

#Savoir le nb de replicats 
ddply(CHR, c("Stress", "Brout", "Prov", "Bloc"), summarise,N = length(Mtot))

#Pour la majorite des combinaisons, nous n'avons qu'un replicat
#Pour simplifier les analyses, on fait une moyenne pour les rares cas ou nous avons 2 plants par combinaison
CHR_mean <- aggregate(list(CHR$Mtot,CHR$Hauteurini,CHR$Ratio,CHR$Maerien,CHR$Mracine), by= data.frame(CHR$Stress, CHR$Brout, CHR$Prov, CHR$Bloc), FUN= "mean")
colnames(CHR_mean)[1:9] <- c("Stress", "Brout", "Prov", "Bloc","Mtot","Hauteurini","Ratio","Maerien","Mracine")

#Verification: nous avons bien un n = 1 pour chaque combinaison.
ddply(CHR_mean, c("Stress", "Brout", "Prov", "Bloc","Hauteurini"), summarise, N = length(Mtot))

#Boxplot masses tot 
boxplot(Mtot~Brout+Stress+Prov,data=CHR,cex.axis=0.5, col = color1,las=2, main = "Masses tot chenes",
        xlab="",ylab="Masse tot")

#Modele avec plan en tiroir
mod_CHR_mtot <- lmerTest::lmer(Mtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CHR_mean)
Anova(mod_CHR_mtot, type=3) # Garder Hini seul 

#Modele sans interaction avec Hini:
mod_CHR_mtot <- lmerTest::lmer(Mtot ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
summary(mod_CHR_mtot)
Anova(mod_CHR_mtot, type=3) #Hauteurini significatif, on garde ce modele
#Interaction Stress*Provenance significatif, Effet de Provenance
#On regarde seulement Stress*Prov

lsmeans(mod_CHR_mtot, pairwise~Stress*Prov)
#NoStress2018 differe de NoStress 2050 et Stress2-NoStress de 2050 significatif

#Representation Stress*Prov avec les estimes du modele
CHR_Mtot_mod <- summary(lsmeans(mod_CHR_mtot, ~Prov*Stress))
ggplot(CHR_Mtot_mod, aes(x = Stress, y = lsmean, shape=Prov)) +
  geom_point(size=4, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  scale_shape_manual(values=c(10, 19, 1))+
  labs(shape="Provenance") +
  scale_x_discrete(labels=c("Non stress?","Stress ?lev?"))+
  scale_y_continuous(limits = c(0,75), expand = expand_scale()) +
  xlab("Traitement de stress hydrique") +  
  ylab("Masse (g)") +
  theme(text=element_text(size=12, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+ 
  geom_text(aes(x=Stress, y=upper.CL+1.9),
            label = c("ab","a","b","ab","b","ab"),  
            position=position_dodge(width=0.4),
            show.legend=F, size=4.5) 


##Verifier Assomptions de modele
plot(mod_CHR_mtot) #Homogeneite des variances ok
leveneTest(resid(mod_CHR_mtot) ~ Stress*Brout*Prov, data = CHR_mean) #ok

qqPlot(resid(mod_CHR_mtot)) #Graphique de normalite ok
shapiro.test(resid(mod_CHR_mtot)) #On verifie avec un test, ok

#################ERABLES(PROV2018)#############
####BIOMASSE TOTALE ERABLES####

#Est-ce que la hauteur initiale est un bon estimateur de la masse totale?
#Regarder la relation entre hauteur finale et la masse a la fin de l'experience
cor.test(ERS$Hauteur, ERS$Mtot) #0.80
plot(ERS$Hauteur, ERS$Mtot) 
#Oui, bonne correlation, utiliser hauteurini pour corriger biomasse initiale des plants

#Savoir le nb de replicats 
ddply(ERS, c("Stress", "Brout", "Bloc"), summarise,N = length(Mtot))
#Un replicat pour toutes les combinaisons

#Boxplot masses tot 
color2 = c(rep("green",2),rep("yellow",2),rep("red",2))
boxplot(Mtot~Brout+Stress,data=ERS,cex.axis=0.5, col = color2,las=2, main = "Masses tot Erables", 
        xlab="",ylab="Masse tot")

#Modele avec plan en tiroir
mod_ERS_mtot <- lmerTest::lmer( Mtot ~ Stress*Brout*Hauteurini + (1|Bloc/Stress), data = ERS)
Anova(mod_ERS_mtot, type=3) #Rien de significatif

#On retire interaction avec Hini
mod_ERS_mtot <- lmerTest::lmer(Mtot ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = ERS)
Anova(mod_ERS_mtot, type=3) #Hauteurini significative, on garde ce modele
#Stress*Brout significatif

#Test a posteriori
summary(lsmeans(mod_ERS_mtot, pairwise~Stress*Brout))
#Rien de significatif

#Representation Stress*Brout avec les estimes du modele 
ERS_Mtot_mod <- summary(lsmeans(mod_ERS_mtot, ~Stress*Brout))

ggplot(ERS_Mtot_mod, aes(x = Stress, y = lsmean, color=Brout)) +
  geom_point(size=4, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  theme(legend.title=element_blank()) + #Enlever titres des legendes
  scale_x_discrete(labels=c("Non stress?","Stress mod?r?", "Stress ?lev?"))+
  scale_y_continuous(limits = c(0,85), expand = expand_scale()) + 
  scale_color_manual(values=c("limegreen","lightcoral"), labels=c("Non Brout?", "Brout?"))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse (g)") +
  theme(text=element_text(size=12, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+ 
  geom_text(aes(x=Stress, y=upper.CL+2),
            label = c("","","","","",""),  
            position=position_dodge(width=0.4), size=4.5, show.legend=F)

##Verifier residus 
plot(mod_ERS_mtot) #Homogeneite des variances ok
leveneTest(resid(mod_ERS_mtot) ~ Stress*Brout, data = ERS) #ok

qqPlot(resid(mod_ERS_mtot)) #Graphique de normalite NOOOOOOOON
shapiro.test(resid(mod_ERS_mtot)) #ok





#################THUYAS########################
####BIOMASSE TOTALE THUYAS####

#Est-ce que la hauteur initiale est un bon estimateur de la masse totale?
#Regarder la relation entre hauteur finale et la masse a la fin de l'experience
cor.test(THO$Hauteur, THO$Mtot) #0.77
plot(THO$Hauteur, THO$Mtot) 
#Oui, bonne correlation, utiliser hauteurini pour corriger biomasse initiale des plants

#Savoir le nb de replicats 
ddply(THO, c("Stress", "Brout", "Prov", "Bloc"), summarise, N = length(Mtot))
#Un replicat pour toutes les combinaisons

#Boxplot masses tot 
boxplot(Mtot~Brout+Stress+Prov,data=THO,cex.axis=0.5, col = color2,las=2, main = "Masses tot Thuyas",
        xlab="", ylab="Masse tot")

#Modele avec plan en tiroir
mod_THO_mtot <- lmerTest::lmer(Mtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = THO)
summary(mod_THO_mtot)
Anova(mod_THO_mtot, type=3) # Interaction avec Hauteurini, on garde ce modele
#Effet Stress

##Verifier Assomptions de modele
plot(mod_THO_mtot) #Homogeneite des variances ok
leveneTest(resid(mod_THO_mtot) ~ Stress*Brout*Prov, data = THO) #ok

qqPlot(resid(mod_THO_mtot), main="Non transform?") #Graphique de normalite, pas certaine
shapiro.test(resid(mod_THO_mtot)) #Non, essayer transformation log

#Ajouter variable Transformation log de Mtot
THO <- mutate(THO, LMtot=log(Mtot))

#Modele avec plan en tiroir (log)
mod_THO_mtotlog <- lmerTest::lmer(LMtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = THO)
Anova(mod_THO_mtotlog, type=3) #Rien significatf

#On enl?ve interaction avec Hini
mod_THO_mtotlog <- lmerTest::lmer(LMtot ~ Stress*Brout*Prov+Hauteurini + (1|Bloc/Stress), data = THO)
summary(mod_THO_mtotlog)
Anova(mod_THO_mtotlog, type=3) #Hini significatif, on garde ce modele
#Pas d'effet significatif


##Verifier Assomptions
plot(mod_THO_mtotlog) #Homogeneite des variances ok
leveneTest(resid(mod_THO_mtotlog) ~ Stress*Brout*Prov, data = THO) #ok

qqPlot(resid(mod_THO_mtotlog), main="Log") #Graphique de normalite, non
shapiro.test(resid(mod_THO_mtotlog)) #ok


#Ajouter variable Transformation sqrt de Mtot
THO <- mutate(THO, sqrtMtot=sqrt(Mtot))

#Modele avec plan en tiroir (sqrt)
mod_THO_mtot_sqrt <- lmerTest::lmer( sqrtMtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = THO)
Anova(mod_THO_mtot_sqrt, type=3) #Rien significatf

#On enl?ve interaction avec Hini
mod_THO_mtot_sqrt <- lmerTest::lmer( sqrtMtot ~ Stress*Brout*Prov+Hauteurini + (1|Bloc/Stress), data = THO)
summary(mod_THO_mtot_sqrt)
Anova(mod_THO_mtot_sqrt, type=3) #Hini significatif, on garde ce modele
#Effet de Brout et Stress


##Verifier Assomptions
plot(mod_THO_mtot_sqrt) #Homogeneite des variances ok
leveneTest(resid(mod_THO_mtot_sqrt) ~ Stress*Brout*Prov, data = THO) #ok

qqPlot(resid(mod_THO_mtot_sqrt), main="Racine carr?e") #Graphique normalite ok
shapiro.test(resid(mod_THO_mtot_sqrt)) #ok


     
#################PINS##########################
####BIOMASSE TOTALE PINS####

#Est-ce que la hauteur initiale est un bon estimateur de la masse totale?
#Regarder la relation entre hauteur finale et la masse a la fin de l'experience
cor.test(PIB$Hauteur, PIB$Mtot) #0.82
plot(PIB$Hauteur, PIB$Mtot) 
#Oui, bonne correlation, utiliser hauteurini pour corriger biomasse initiale des plants

#Savoir le nb de replicats 
ddply(PIB, c("Stress", "Brout", "Prov", "Bloc"), summarise,N = length(Mtot))
#Un replicat pour toutes les combinaisons

#Boxplot masses tot 
boxplot(Mtot~Brout+Stress+Prov,data=PIB,cex.axis=0.5, col = color2,las=2, main = "Masses tot Pins", 
        xlab="", ylab="Masse tot")

#Modele avec plan en tiroir
mod_PIB_mtot <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = PIB)
Anova(mod_PIB_mtot, type=3) #Interactions avec Hini significatifs
#Effet Stress*Brout
summary(lsmeans(mod_PIB_mtot, pairwise~Stress*Brout))
#NoBrout-NoStress differe de Brout-Stress2
#NoBrout-Stress1 differe de Brout-Stress2
#NoBrout-Stress2 differe de Brout-Stress2

#Representation avec estimes du modele
PIB_Mtot_mod <- summary(lsmeans(mod_PIB_mtot, ~Stress*Brout))

ggplot(PIB_Mtot_mod, aes(x = Stress, y = lsmean, color=Brout)) +
  geom_point(size=4, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  theme(legend.title=element_blank()) + #Enlever titres des legendes
  scale_x_discrete(labels=c("Non stress?","Stress mod?r?", "Stress ?lev?"))+
  scale_y_continuous(limits = c(0,45), expand = expand_scale()) + 
  scale_color_manual(values=c("limegreen","lightcoral"), labels=c("Non Brout?", "Brout?"))+
  xlab("Traitement de stress hydrique") +  
  ylab("Masse (g)") +
  theme(text=element_text(size=12, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+ 
  geom_text(aes(x=Stress, y=upper.CL+2),
            label = c("ab","a","ab","a","b","a"),  
            position=position_dodge(width=0.4), size=4.5, show.legend=F)


##Verifier residus en incluant les effets aleatoires (effet du design)
plot(mod_PIB_mtot) #Homogeneite des variances ok
leveneTest(resid(mod_PIB_mtot) ~ Stress*Brout*Prov, data = PIB) #ok

qqPlot(resid(mod_PIB_mtot)) #Graphique de normalite, non
shapiro.test(resid(mod_PIB_mtot)) #non



#Ajouter variable Transformation log de Mtot
PIB <- mutate(PIB, LMtot=log(Mtot))

#Modele avec plan en tiroir (log)
mod_PIB_mtotlog <- lmerTest::lmer( LMtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = PIB)
Anova(mod_PIB_mtotlog, type=3) #Interaction significative avec Hini, on garde ce modele

##Verifier Assomptions
qqPlot(resid(mod_PIB_mtotlog)) #Graphique de normalite, non
shapiro.test(resid(mod_PIB_mtotlog)) #non

#Ajouter variable Transformation sqrt de Mtot
PIB <- mutate(PIB, sqrtMtot=sqrt(Mtot))

#Modele avec plan en tiroir (sqrt)
mod_PIB_mtot_sqrt <- lmerTest::lmer( sqrtMtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = PIB)
Anova(mod_PIB_mtot_sqrt, type=3) #Interaction significative avec Hini, on garde ce modele

##Verifier Assomptions
qqPlot(resid(mod_PIB_mtot_sqrt))#Graphique de normalite, non
shapiro.test(resid(mod_PIB_mtot_sqrt)) #non

#BOXCOX pour trouver quelle puissance serait la meilleure pour rendre les donnees normales
bc<-boxcox(PIB$Mtot ~ 1)
bc$x[which.max(bc$y)] #lambda max ? 0,46 --> semblable a faire une transfo racine carree


#Modele Robuste
library(robustlmm)
mod_PIB_mtot_robust<-rlmer(Mtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = PIB)
summary(mod_PIB_mtot_robust)
Anova(mod_PIB_mtot_robust, type=3)




