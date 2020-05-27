
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
#Masse totale - Outliers

boxplot(CET$Mtot, xlab="Masse totale Cerisier")

boxplot(CHR$Mtot, xlab="Masse totale Chene")
dotchart(CHR$Mtot,xlab="Masse totale Chene", ylab = "Ordre des donnees") #Un outlier? Semble ok dans les donnees

boxplot(ERS$Mtot, xlab="Masse totale Erable")

boxplot(THO$Mtot, xlab="Masse totale Thuya")
dotchart(THO$Mtot,xlab="Masse totale Thuya", ylab = "Ordre des donnees") #Pas d'outlier important

boxplot(PIB$Mtot, xlab="Masse totale Pin")
dotchart(PIB$Mtot,xlab="Masse totale Pin", ylab = "Ordre des donnees") #Pas d'outlier important

#Ratio - Outliers

boxplot(CET$Ratio, xlab="Ratio Cerisier")
dotchart(CET$Ratio, xlab="Ratio Cerisier", ylab = "Ordre des donnees") 

boxplot(CHR$Ratio, xlab="Ratio Chene")
dotchart(CHR$Ratio, xlab="Ratio Chene", ylab = "Ordre des donnees") #Quelques donnes plus extremes

boxplot(ERS$Ratio, xlab="Ratio Erable")

boxplot(THO$Ratio, xlab="Ratio Thuya")
dotchart(THO$Ratio, xlab="Ratio Thuya", ylab = "Ordre des donnees") 

boxplot(PIB$Ratio, xlab="Ratio Pin")


#Normalite Masse totale
hist(CET$Mtot)
hist(CHR$Mtot)
hist(ERS$Mtot)
hist(PIB$Mtot)
hist(THO$Mtot)


#Normalite Ratio
hist(CET$Ratio)
hist(CHR$Ratio)
hist(ERS$Ratio)
hist(PIB$Ratio)
hist(THO$Ratio)

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
#Message veut dire que l'effet aleatoire du bloc est tres pres de 0 ou nul
#On le garde requand meme dans le modele, mais il a peu d'effet
summary(mod_CET_mtot)
anova(mod_CET_mtot) #Fait un anova de type III de Satterthwaite
#Interaction significative avec Hini, on garde ce modele
#Effet Stress

#Estimés et intervalle de confiance
lsmeans(mod_CET_mtot,~Stress)
#Rien n'est significatif

#Representation Stress avec les estimes du modele
CET_Mtot_mod <- summary(lsmeans(mod_CET_mtot, ~Stress))
ggplot(CET_Mtot_mod, aes(x = Stress, y = lsmean)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  theme_classic() +
  scale_x_discrete(labels=c("No stress","High stress"))+
  scale_y_continuous(limits = c(0,90), expand = expand_scale()) +
  xlab("Water stress treatment") +  
  ylab("Mass (g)") +
  theme(text=element_text(size=12, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  geom_text(aes(x=Stress, y=upper.CL+2),
            label = c("",""), size=4.5, show.legend=F)

##Verifier assomptions de modeles
plot(mod_CET_mtot) #Homogeneite des variances ok
qqPlot(resid(mod_CET_mtot)) #Graphique de normalite ok


####RATIO CERISIERS####

#Modele avec plan en tiroir
mod_CET_ratio <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CET_mean)
summary(mod_CET_ratio)
anova(mod_CET_ratio) #Interaction non significative avec Hini

#Modele sans interaction Hini
mod_CET_ratio <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov+Hauteurini + (1|Bloc/Stress), data = CET_mean)
summary(mod_CET_ratio)
anova(mod_CET_ratio) #Hini non significative 

#Modele sans Hini
mod_CET_ratio <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CET_mean)
summary(mod_CET_ratio)
anova(mod_CET_ratio) #Effet Stress et Prov

#Estimés et intervalle de confiance
lsmeans(mod_CET_ratio,~Stress) #Ratio plus eleve pour stress2
lsmeans(mod_CET_ratio,~Prov) #Ratio plus eleve chez 2080 que 2018 et 2050 


#Representation Stress avec les estimes du modele
CET_Ratio_mod <- summary(lsmeans(mod_CET_ratio, ~Stress))
ggplot(CET_Ratio_mod, aes(x = Stress, y = lsmean)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  theme_classic() +
  scale_x_discrete(labels=c("No stress","High stress"))+
  scale_y_continuous(limits = c(0,4), expand = expand_scale()) +
  xlab("Water stress treatment") +  
  ylab("Shoot:root ratio") +
  theme(text=element_text(size=12, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  geom_text(aes(x=Stress, y=upper.CL+0.2),
          label = c("a","b"), size=4.5, show.legend=F)

#Representation Prov 
CET_Ratio_mod <- summary(lsmeans(mod_CET_ratio, ~Prov))
ggplot(CET_Ratio_mod, aes(x = Prov, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  theme_classic() +
  xlab("Provenance") +  
  ylab("Shoot:root ratio") +
  theme(text=element_text(size=12, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  geom_text(aes(x=Prov, y=upper.CL+0.05),label = c("a","a","b"))

##Verifier assomptions de modeles
plot(mod_CET_ratio) #Homogeneite des variances ok
qqPlot(resid(mod_CET_ratio)) #Graphique de normalite ok


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
anova(mod_CHR_mtot) # Garder Hini seul 

#Modele sans interaction avec Hini:
mod_CHR_mtot <- lmerTest::lmer(Mtot ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
summary(mod_CHR_mtot)
anova(mod_CHR_mtot) #Hauteurini significatif, on garde ce modele
#Interaction Stress*Provenance significatif

#Estimes et intervalles de confiance
lsmeans(mod_CHR_mtot,~Stress*Prov)
#Stress2-NoStress de 2050 different

#Representation Stress*Prov avec les estimes du modele
CHR_Mtot_mod <- summary(lsmeans(mod_CHR_mtot, ~Prov*Stress))
ggplot(CHR_Mtot_mod, aes(x = Stress, y = lsmean, shape=Prov)) +
  geom_point(size=4, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  scale_shape_manual(values=c(10, 19, 1))+
  labs(shape="Provenance") +
  scale_x_discrete(labels=c("No stress","High stress"))+
  scale_y_continuous(limits = c(0,75), expand = expand_scale()) +
  xlab("Water stress treatment") +  
  ylab("Mass (g)") +
  theme(text=element_text(size=12, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+ 
  geom_text(aes(x=Stress, y=upper.CL+1.9),
            label = c("ab","a","ab","ab","b","ab"),  
            position=position_dodge(width=0.4),
            show.legend=F, size=4.5) 


##Verifier Assomptions de modele
plot(mod_CHR_mtot) #Homogeneite des variances ok
qqPlot(resid(mod_CHR_mtot)) #Graphique de normalite ok

####RATIO CHENES####

##Modele avec plan en tiroir
mod_CHR_ratio <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(mod_CHR_ratio) #Rien significatif

##Modele avec plan en tiroir sans interaction Hini
mod_CHR_ratio <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(mod_CHR_ratio) #Rien significatif

##Modele avec plan en tiroir sans Hini
mod_CHR_ratio <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CHR_mean)
anova(mod_CHR_ratio) #Rien significatif, tendance Stress

##Verifier assomptions de modeles
plot(mod_CHR_ratio) #Homogeneite des variances ok
qqPlot(resid(mod_CHR_ratio)) #Graphique de normalite ok?


#Essayer d'enlever trois valeurs plus extreme (30,12,11) 
CHR_mean2 <- CHR_mean[-30,]
CHR_mean2 <- CHR_mean2[-12,]
CHR_mean2 <- CHR_mean2[-11,]
dotchart(CHR_mean$Ratio,xlab="Ratio CHR", ylab = "Ordre des donnees") #Avant
dotchart(CHR_mean2$Ratio,xlab="Ratio CHR", ylab = "Ordre des donnees") #Apres

##Modele avec plan en tiroir sans outlier
mod_CHR_ratio2 <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CHR_mean2)
anova(mod_CHR_ratio2) #Rien significatif

##Modele avec plan en tiroir sans interaction Hini sans outlier
mod_CHR_ratio2 <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean2)
anova(mod_CHR_ratio2) #Rien significatif

##Modele avec plan en tiroir sans Hini sans outlier
mod_CHR_ratio2 <- lmerTest::lmer( Ratio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CHR_mean2)
anova(mod_CHR_ratio2) #Rien significatif, tendance Stress

##Verifier assomptions de modeles
plot(mod_CHR_ratio2) #Homogeneite des variances ok
qqPlot(resid(mod_CHR_ratio2)) #Graphique de normalite pas mieux!?


#Ajouter variable Transformation log de Ratio
LCHR_mean <- mutate(CHR_mean, LRatio=log(Ratio))

##Modele avec plan en tiroir (log)
mod_CHR_Lratio <- lmerTest::lmer( LRatio ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = LCHR_mean)
anova(mod_CHR_Lratio) #Rien significatif

##Modele avec plan en tiroir sans interaction Hini (log)
mod_CHR_Lratio <- lmerTest::lmer( LRatio ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = LCHR_mean)
anova(mod_CHR_Lratio) #Rien significatif

##Modele avec plan en tiroir sans Hini (log)
mod_CHR_Lratio <- lmerTest::lmer( LRatio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = LCHR_mean)
anova(mod_CHR_Lratio) #Effet Prov, Tendance Stress 

lsmeans(mod_CHR_Lratio,~Prov)
#Rien

#Test a posteriori aurait trouve que 2018-2050 significatif.......???
ls_means(mod_CHR_Lratio,  which = "Prov", pairwise = TRUE)

#Representation Prov avec les estimations du modele (log)
CHR_LRatio_mod <- summary(lsmeans(mod_CHR_Lratio, ~Prov))
ggplot(CHR_LRatio_mod, aes(x = Prov, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  theme_classic() +
  xlab("Provenance") +  
  ylab("Log shoot:root ratio")+
  theme(text=element_text(size=12, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))

#Enlever le log
CHR_Ratio_mod <- mutate(CHR_LRatio_mod, lsmean=10^(lsmean),lower.CL=10^(lower.CL), upper.CL=10^(upper.CL))
ggplot(CHR_Ratio_mod, aes(x = Prov, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  theme_classic() +
  xlab("Provenance") +  
  ylab("Log shoot:root ratio")

#Representation estime modele tendance effet Stress
CHR_LRatio_mod <- summary(lsmeans(mod_CHR_Lratio, ~Stress))
ggplot(CHR_LRatio_mod, aes(x = Stress, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  theme_classic() +
  xlab("Stress") +  
  ylab("Log shoot:root ratio")
#Ratio tendance a etre plus elevee pour Stress2


##Verifier assomptions de modeles
plot(mod_CHR_Lratio) #Homogeneite des variances ok
qqPlot(resid(mod_CHR_Lratio)) #Graphique de normalite ok


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
anova(mod_ERS_mtot) #Rien de significatif

#On retire interaction avec Hini
mod_ERS_mtot <- lmerTest::lmer(Mtot ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = ERS)
summary(mod_ERS_mtot)
anova(mod_ERS_mtot) #Hauteurini non significative, on retire du modele

#On retire Hini
mod_ERS_mtot <- lmerTest::lmer(Mtot ~ Stress*Brout + (1|Bloc/Stress), data = ERS)
summary(mod_ERS_mtot)
anova(mod_ERS_mtot) 
#Stress*Brout significatif

#Estimés et intervalles de confiance
lsmeans(mod_ERS_mtot, ~Stress*Brout)
#Différence entre NoBrout et Brout de Stress1

#Representation Stress*Brout avec les estimes du modele 
ERS_Mtot_mod <- summary(lsmeans(mod_ERS_mtot, ~Stress*Brout))
ggplot(ERS_Mtot_mod, aes(x = Stress, y = lsmean, shape=Brout)) +
  geom_point(size=3, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  theme(legend.title=element_blank()) + #Enlever titres des legendes
  scale_x_discrete(labels=c("No stress","Moderate stress", "High stress"))+
  scale_y_continuous(limits = c(0,85), expand = expand_scale()) + 
  scale_shape_manual(values=c(19,17), labels=c("Unbrowsed", "Browsed"))+
  xlab("Water stress treatment") +  
  ylab("Mass (g)") +
  theme(text=element_text(size=12, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+ 
  geom_text(aes(x=Stress, y=upper.CL+2),
            label = c("ab","a","ab","ab","b","ab"),  
            position=position_dodge(width=0.4), size=4.5, show.legend=F)

##Verifier residus 
plot(mod_ERS_mtot) #Homogeneite des variances ok
qqPlot(resid(mod_ERS_mtot)) #Graphique de normalite ok

#Reponse compensatoire pour les plants de Stress1
ERS_stress1_nobrout <- filter(ERS, Brout=="NoBrout", Stress=="Stress1")
moyERS_stress1_nobrout <- mean(ERS_stress1_nobrout$Mtot) #30.62

ERS_stress1_brout <- filter(ERS, Brout=="Brout", Stress=="Stress1")
moyERS_stress1_brout <- mean(ERS_stress1_brout$Mtot) #63.6
moyERS_stress1_Mret <- mean(ERS_stress1_brout$Mret) #0.28

((moyERS_stress1_brout + moyERS_stress1_Mret) - moyERS_stress1_nobrout)*(100/moyERS_stress1_nobrout) # -46.59
# 108.62

####RATIO ERABLES####

#Modele avec plan en tiroir
mod_ERS_ratio <- lmerTest::lmer(Ratio ~ Stress*Brout*Hauteurini + (1|Bloc/Stress), data = ERS)
anova(mod_ERS_ratio) #Rien de significatif

#On retire interaction avec Hini
mod_ERS_ratio <- lmerTest::lmer(Ratio ~ Stress*Brout + Hauteurini + (1|Bloc/Stress), data = ERS)
summary(mod_ERS_ratio)
anova(mod_ERS_ratio) #Hauteurini non significative, on retire du modele

#On retire Hini
mod_ERS_ratio <- lmerTest::lmer(Ratio ~ Stress*Brout + (1|Bloc/Stress), data = ERS)
summary(mod_ERS_ratio)
anova(mod_ERS_ratio) #Rien

##Verifier residus 
plot(mod_ERS_ratio) #Homogeneite des variances ok
qqPlot(resid(mod_ERS_ratio)) #Graphique de normalite ok

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
anova(mod_THO_mtot) # Interaction avec Hauteurini, on garde ce modele
#Effet Stress*Brout

#Estimés et intervalles de confiance
lsmeans(mod_THO_mtot, ~Stress*Brout)
#NoStressNoBrout differe de Stress2NoBrout 
#NoStressNoBrout differe de Stress2Brout 

#Representation Stress et Brout avec estimes du modele (sans outlier)
THO_Mtot_mod <- summary(lsmeans(mod_THO_mtot, ~Stress*Brout))
ggplot(THO_Mtot_mod, aes(x = Stress, y = lsmean, shape=Brout)) +
  geom_point(size=3, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  theme(legend.title=element_blank()) + #Enlever titres des legendes
  scale_x_discrete(labels=c("No stress","Moderate stress", "High stress"))+
  scale_y_continuous(limits = c(0,85), expand = expand_scale()) + 
  scale_shape_manual(values=c(19,17), labels=c("Unbrowsed", "Browsed"))+
  xlab("Water stress treatment") +  
  ylab("Mass (g)") +
  theme(text=element_text(size=12, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+ 
  geom_text(aes(x=Stress, y=upper.CL+2),
            label = c("a","ab","b","ab","ab","b"),  
            position=position_dodge(width=0.4), size=4.5, show.legend=F)

##Verifier Assomptions de modele
plot(mod_THO_mtot) #Homogeneite des variances ok
qqPlot(resid(mod_THO_mtot), main="Non transforme") #Graphique de normalite, ok

#Reponse compensatoire de BroutStress2 par rapport a NoBroutNoStress
THO_nobrout_nostress <- filter(THO, Brout=="NoBrout", Stress=="NoStress")
moyTHO_nobrout_nostress <- mean(THO_nobrout_nostress$Mtot) #66.98

THO_brout_stress2 <- filter(THO, Brout=="Brout", Stress=="Stress2")
moyTHO_brout_stress2 <- mean(THO_brout_stress2$Mtot) #34.15
moyTHO_stress2_Mret <- mean(THO_brout_stress2$Mret) #2.1

((moyTHO_brout_stress2 + moyTHO_stress2_Mret) - moyTHO_nobrout_nostress)*(100/moyTHO_nobrout_nostress) # -46.59
# -45.87

####RATIO THUYAS####

#Modele avec plan en tiroir
mod_THO_ratio <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = THO)
summary(mod_THO_ratio)
anova(mod_THO_ratio) # Pas interaction avec Hauteurini

#Modele avec plan en tiroir sans interaction Hini
mod_THO_ratio <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = THO)
summary(mod_THO_ratio)
anova(mod_THO_ratio) # Hauteurini pas significatif

#Modele sans Hini
mod_THO_ratio <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = THO)
summary(mod_THO_ratio)
anova(mod_THO_ratio) #Stress*Brout*Prov

#Estimés et intervalles de confiance
lsmeans(mod_THO_ratio, ~Stress*Brout*Prov) #Rien mais 2018 tend a avoir effet du broutement sous stress2

#Representation Stress*Brout*Prov avec estimes du modele 
THO_Ratio_mod <- summary(lsmeans(mod_THO_ratio, ~Stress*Brout*Prov))
ggplot(THO_Ratio_mod, aes(x = Stress, y = lsmean, col=Brout, shape=Prov)) +
  geom_point(size=3, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  theme_classic() +
  theme(legend.title=element_blank()) + #Enlever titres des legendes
  scale_x_discrete(labels=c("No stress","Moderate stress", "High stress"))+
  scale_y_continuous(limits = c(0,6.5), expand = expand_scale()) + 
  scale_shape_manual(values=c(10, 19, 1))+
  scale_color_manual(values=c("turquoise2","lightcoral"),
                     labels=c("Unbrowsed", "Browsed"))+
  xlab("Water stress treatment") +  
  ylab("Shoot:root ratio") +
  theme(text=element_text(size=12, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  geom_text(aes(x=Stress, y=upper.CL+2),
            label = c("a","ab","b","ab","ab","b"),  
            position=position_dodge(width=0.4), size=4.5, show.legend=F)

##Verifier Assomptions de modele
plot(mod_THO_ratio) #Homogeneite des variances ok
qqPlot(resid(mod_THO_ratio)) #Graphique de normalite, ok


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
summary(mod_PIB_mtot)
anova(mod_PIB_mtot) #Rien significatif

#Modele avec plan en tiroir sans interaction Hini
mod_PIB_mtot <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = PIB)
summary(mod_PIB_mtot)
anova(mod_PIB_mtot) #Hini significatif on garde ce modele
#Brout

#Representation effet Brout avec estimes du modele
PIB_Mtot_mod <- summary(lsmeans(mod_PIB_mtot, ~ Brout))
ggplot(PIB_Mtot_mod, aes(x = Brout, y = lsmean)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  theme_classic() +
  theme(legend.title=element_blank()) + #Enlever titres des legendes
  scale_x_discrete(labels=c("Unbrowsed","Browsed"))+
  scale_y_continuous(limits = c(0,45), expand = expand_scale()) + 
  xlab("Browsing treatment") +  
  ylab("Mass (g)") +
  theme(text=element_text(size=12, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+ 
  geom_text(aes(x=Brout, y=upper.CL+2),
            label = c("a","b"),  
            position=position_dodge(width=0.4), size=4.5, show.legend=F)

##Verifier Assomptions
plot(mod_PIB_mtot) #Homogeneite des variances ok
qqPlot(resid(mod_PIB_mtot)) #Graphique de normalite ok

#Essayer d'enlever la valeur plux extreme (31) 
PIB2 <- PIB[-31,]
dotchart(PIB$Mtot,xlab="Masse totale Pin", ylab = "Ordre des donnees") #Avant
dotchart(PIB2$Mtot,xlab="Masse totale Pin", ylab = "Ordre des donnees") #Apres

#Modele avec plan en tiroir
mod_PIB2_mtot <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = PIB2)
anova(mod_PIB2_mtot) #Rien

#Modele avec plan en tiroir sans interaction hini
mod_PIB2_mtot <- lmerTest::lmer( Mtot ~ Stress*Brout*Prov+Hauteurini + (1|Bloc/Stress), data = PIB2)
anova(mod_PIB2_mtot) #Hini significatif on garde ce modele
#Effet Brout

##Verifier Assomptions
plot(mod_PIB2_mtot) #Homogeneite des variances ok
qqPlot(resid(mod_PIB2_mtot)) #Graphique de normalite ok

#Reponse compensatoire
PIB_brout <- filter(PIB, Brout=="Brout")
moyPIB_brout <- mean(PIB_brout$Mtot)
moyPIB_Mret <- mean(PIB_brout$Mret)
PIB_nobrout <- filter(PIB, Brout=="NoBrout")
moyPIB_nobrout <- mean(PIB_nobrout$Mtot)

((moyPIB_brout + moyPIB_Mret) - moyPIB_nobrout)*(100/moyPIB_nobrout) # -46.59

#Donne le meme resultat 
####RATIO PINS####

#Modele avec plan en tiroir
mod_PIB_ratio <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = PIB)
summary(mod_PIB_ratio)
anova(mod_PIB_ratio) #Rien significatif

#Modele avec plan en tiroir sans interaction Hini
mod_PIB_ratio <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov+Hauteurini + (1|Bloc/Stress), data = PIB)
summary(mod_PIB_ratio)
anova(mod_PIB_ratio) #Rien significatif

#Modele avec plan en tiroir sans Hini
mod_PIB_ratio <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov (1|Bloc/Stress), data = PIB)
summary(mod_PIB_ratio)
anova(mod_PIB_ratio) #Rien significatif

##Verifier Assomptions
plot(mod_PIB_ratio) #Homogeneite des variances ok
qqPlot(resid(mod_PIB_ratio)) #Graphique de normalite ok

