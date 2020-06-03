
#2019-2020
#Projet Serre / Broutement / Stress_hydrique / Migration assistee

remove(list = ls())

#####IMPORTATION DES DONNEES####

#Packages
library(dplyr)
library(ggplot2)
library(plyr) #pour fonction ddply 
library(car) #pour anovas type III
library(lsmeans) #fonction lsmeans
library(MASS) #boxcox
library(plyr)
library(tidyverse)
library(corrplot)
library(Rmisc)#pour summarySE
library(cowplot)
  theme_set(theme_cowplot())
library(lme4)
library(corrplot)
library(Hmisc) #violin plot
  
#Manipulation du jeu de donnees

#Jeu de donnees des masses
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

Serre$Brout <- factor(Serre$Brout, levels=c("NoBrout","Brout")) #Reorganiser niveaux facteur Brout

#Jeux de donnees de l'analyse chimique
chemPIB <- read.table("./PIB.txt", header = TRUE, dec = ",", stringsAsFactors = FALSE)
chemERS <- read.table("./ERS.txt", header = TRUE, dec = ",", stringsAsFactors = FALSE)
chemCET <- read.table("./CET.txt", header = TRUE, dec = ",", stringsAsFactors = FALSE)
chemTHO <- read.table("./THO.txt", header = TRUE, dec = ",", stringsAsFactors = FALSE)
chemCHR <- read.table("./CHR.txt", header = TRUE, dec = ",", stringsAsFactors = FALSE)

chem <- rbind(chemPIB, chemERS, chemCET, chemTHO, chemCHR)
colnames(chem)[2] <- "ID"
chem <- subset(chem, select = -ID_lab)

#Jeu de donnees fusionnes
data <- inner_join(Serre, chem, by = "ID")

#Jeu de donnees par espece
CET<-filter(data,Esp=="CET")
CET$Stress <- factor(CET$Stress) #Enlever le niveau Stress 1 absent de CET
CHR<-filter(data,Esp=="CHR")
CHR$Stress <- factor(CHR$Stress) #Enlever le niveau Stress 1 absent de CHR
ERS<-filter(data,Esp=="ERS")
PIB<-filter(data,Esp=="PIB")
THO<-filter(data,Esp=="THO")
cTHO <- na.omit(THO) #Enlever 2 lignes avec des NA pour l'analyse chimique seulement

#Savoir le nb de replicats pour chaque espece
ddply(CET, c("Stress", "Brout", "Prov", "Bloc"), summarise, N = length(Mtot))
ddply(CHR, c("Stress", "Brout", "Prov", "Bloc"), summarise,N = length(Mtot))
ddply(ERS, c("Stress", "Brout", "Bloc"), summarise,N = length(Mtot))
ddply(THO, c("Stress", "Brout", "Prov", "Bloc"), summarise, N = length(Mtot))
ddply(PIB, c("Stress", "Brout", "Prov", "Bloc"), summarise,N = length(Mtot))


#Pour la majorite des combinaisons, nous n'avons qu'un replicat
#Pour simplifier les analyses, on fait une moyenne pour les rares cas 
#ou nous avons 2 plants par combinaison (CET et CHR seulement)

CET_mean <- aggregate(list(CET$Mtot,CET$Hauteurini,CET$Ratio,CET$Maerien,CET$Mracine, CET$N, CET$Phen, CET$Flav), by= data.frame(CET$Stress, CET$Brout, CET$Prov, CET$Bloc), FUN= "mean")
colnames(CET_mean)[1:12] <- c("Stress", "Brout", "Prov", "Bloc","Mtot","Hauteurini","Ratio","Maerien","Mracine", "N", "Phen", "Flav")

CHR_mean <- aggregate(list(CHR$Mtot,CHR$Hauteurini,CHR$Ratio,CHR$Maerien,CHR$Mracine, CHR$N, CHR$Phen, CHR$Flav), by= data.frame(CHR$Stress, CHR$Brout, CHR$Prov, CHR$Bloc), FUN= "mean")
colnames(CHR_mean)[1:12] <- c("Stress", "Brout", "Prov", "Bloc","Mtot","Hauteurini","Ratio","Maerien","Mracine", "N", "Phen", "Flav")




#####MASSE - Exploration des donnees####

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

#Regarder si hauteur initiale est un bon estimateur de la masse totale
#Regarder la relation entre hauteur finale et la masse a la fin de l'experience
cor.test(CET$Hauteur, CET$Mtot) #0,63
plot(CET$Hauteur, CET$Mtot) 
cor.test(CHR$Hauteur, CHR$Mtot) #0,67
plot(CHR$Hauteur, CHR$Mtot) 
cor.test(ERS$Hauteur, ERS$Mtot) #0.80
plot(ERS$Hauteur, ERS$Mtot) 
cor.test(THO$Hauteur, THO$Mtot) #0.77
plot(THO$Hauteur, THO$Mtot) 
cor.test(PIB$Hauteur, PIB$Mtot) #0.82
plot(PIB$Hauteur, PIB$Mtot) 
#Oui, bonne correlation, on utilise la hauteur initiale pour corriger pour la biomasse initiale des plants


#Boxplot Masses totales
color1 = c(rep("green",2),rep("red",2))
color2 = c(rep("green",2),rep("yellow",2),rep("red",2))

boxplot(Mtot~Brout+Stress+Prov,data=CET,cex.axis=0.5, col = color1,las=2, 
        main = "Masses tot cerisiers", ylab="Masse tot", xlab="") 

boxplot(Mtot~Brout+Stress+Prov,data=CHR,cex.axis=0.5, col = color1,las=2, 
        main = "Masses tot chenes", xlab="",ylab="Masse tot")

boxplot(Mtot~Brout+Stress,data=ERS,cex.axis=0.5, col = color2,las=2, 
        main = "Masses tot Erables", xlab="",ylab="Masse tot")

boxplot(Mtot~Brout+Stress+Prov,data=THO,cex.axis=0.5, col = color2,las=2, 
        main = "Masses tot Thuyas", xlab="", ylab="Masse tot")

boxplot(Mtot~Brout+Stress+Prov,data=PIB,cex.axis=0.5, col = color2,las=2, 
        main = "Masses tot Pins", xlab="", ylab="Masse tot")


#####MASSE CERISIERS#################
    ####MASSE TOTALE CERISIERS####

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
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  geom_text(aes(x=Stress, y=upper.CL+2.5),
            label = c("a","b"), size=4.5, show.legend=F)

#Donnees brutes
summarySE(CET, measurevar = "Mtot", groupvars = c("Stress")) 

#Violin plot Stress avec les donnees brutes (moyenne et intervalle de confiance)
ggplot(CET, aes(x = Stress, y = Mtot)) +
  geom_violin(fill="grey") +
  scale_x_discrete(labels=c("No stress","High stress"))+
  scale_y_continuous(limits = c(0,110), breaks= seq(0, 100, by = 20), expand = expand_scale()) +
  xlab("Water stress treatment") +  
  ylab("Mass (g)") +
  stat_summary(fun = mean, geom = "point") + 
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange") + 
  theme(text=element_text(size=12, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))

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

#Estimés et intervalle de confiance de Stress
lsmeans(mod_CET_ratio,~Stress) #Ratio plus eleve pour stress2

#Estimés et intervalle de confiance et test aposteriori de Provenance
lsmeans(mod_CET_ratio,pairwise~Prov) #Ratio plus eleve chez 2080 que 2018 et 2050 


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

#Representation Prov estimes du modele
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

#Donnees brutes
summarySE(CET, measurevar = "Ratio", groupvars = c("Stress")) 
summarySE(CET, measurevar = "Ratio", groupvars = c("Prov")) 

#Violin plot Stress avec les donnees brutes (moyenne et intervalle de confiance)
ggplot(CET, aes(x = Stress, y = Ratio)) +
  geom_violin(fill="grey") +
  scale_x_discrete(labels=c("No stress","High stress"))+
  scale_y_continuous(limits = c(0,6), breaks= seq(0, 6, by = 1), expand = expand_scale()) +
  xlab("Water stress treatment") +  
  ylab("Mass (g)") +
  stat_summary(fun = mean, geom = "point") + 
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange") + 
  theme(text=element_text(size=12, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))

#Violin plot Provenance avec les donnees brutes (moyenne et intervalle de confiance)
ggplot(CET, aes(x = Prov, y = Ratio)) +
  geom_violin(fill="grey") +
  scale_x_discrete(labels=c("2018","2050","2080"))+
  #scale_y_continuous(limits = c(0,6), breaks= seq(0, 6, by = 1), expand = expand_scale()) +
  xlab("Provenance") +  
  ylab("Shoot :Root Ratio") +
  stat_summary(fun = mean, geom = "point") + 
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange") + 
  theme(text=element_text(size=12, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))

##Verifier assomptions de modeles
plot(mod_CET_ratio) #Homogeneite des variances ok
qqPlot(resid(mod_CET_ratio)) #Graphique de normalite ok


#####MASSE CHENES###################
    ####MASSE TOTALE CHENES####

#Modele avec plan en tiroir
mod_CHR_mtot <- lmerTest::lmer(Mtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = CHR_mean)
anova(mod_CHR_mtot) # Garder Hini seul 

#Modele sans interaction avec Hini:
mod_CHR_mtot <- lmerTest::lmer(Mtot ~ Stress*Brout*Prov + Hauteurini + (1|Bloc/Stress), data = CHR_mean)
summary(mod_CHR_mtot)
anova(mod_CHR_mtot) #Hauteurini significatif, on garde ce modele
#Interaction Stress*Provenance significatif

#Estimes et intervalles de confiance et test a posteriori
lsmeans(mod_CHR_mtot, pairwise~Stress*Prov)
#NoStress de 2018 et 2050 different
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
            label = c("b","a","ab","ab","b","ab"),  
            position=position_dodge(width=0.4),
            show.legend=F, size=4.5) 

#Violin plot Stress*Prov avec les donnees brutes (moyenne et intervalle de confiance)
ggplot(CHR, aes(x = Stress, y = Mtot, shape=Prov)) +
  geom_violin(fill="grey", position=position_dodge(width=0.8)) +
  stat_summary(fun = mean, geom = "point", position=position_dodge(width=0.8)) + 
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange", position=position_dodge(width=0.8)) + 
  theme_classic() +
  scale_shape_manual(values=c(19, 15, 17))+
  labs(shape="Provenance") +
  scale_x_discrete(labels=c("No stress","High stress"))+
  scale_y_continuous(limits = c(0,120), expand = expand_scale()) +
  xlab("Water stress treatment") +  
  ylab("Mass (g)") +
  theme(text=element_text(size=12, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) 
  geom_text(aes(x=Stress, y=),
            label = c("b","a","ab","ab","b","ab"),  
            position=position_dodge(width=0.4),
            show.legend=F, size=4.5)

#Violin plot Stress*Prov avec les donnees brutes (jitter)
ggplot(CHR, aes(x = Stress, y = Mtot, shape=Prov)) +
  geom_violin(fill="grey", position=position_dodge(width=0.8), bw=8) + #trim=FALSE 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, position=position_dodge(width=0.8))+
    theme_classic() +
  scale_shape_manual(values=c(19, 15, 17))+
  labs(shape="Provenance") +
  scale_x_discrete(labels=c("No stress","High stress"))+
  scale_y_continuous(limits = c(0,120), expand = expand_scale()) +
  xlab("Water stress treatment") +  
  ylab("Mass (g)") +
  theme(text=element_text(size=12, family="sans"), #Police Arial ("serif" pour Time New Roman)
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)) 


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
qqPlot(resid(mod_CHR_ratio)) #Graphique de normalite ???


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
qqPlot(resid(mod_CHR_ratio2)) #Graphique de normalite pas beaucoup mieux?

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

#Estimes, intervalle de confiance et test a posteriori pour Provenance
lsmeans(mod_CHR_Lratio,pairwise~Prov)
#Difference entre 2018 et 2050

#Representation Prov avec les estimations du modele (log)
CHR_LRatio_mod <- summary(lsmeans(mod_CHR_Lratio, ~Prov))
ggplot(CHR_LRatio_mod, aes(x = Prov, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  theme_classic() +
  xlab("Provenance") +  
  ylab("Log shoot:root ratio")+
  theme(text=element_text(size=12, family="sans"), #Police Arial ("serif" pour Time New Roman)
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  geom_text(aes(x=Prov, y=upper.CL+0.05),label = c("a","b","ab"))

#Enlever le log
CHR_Ratio_mod <- mutate(CHR_LRatio_mod, lsmean=10^(lsmean),lower.CL=10^(lower.CL), upper.CL=10^(upper.CL))
ggplot(CHR_Ratio_mod, aes(x = Prov, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  theme_classic() +
  xlab("Provenance") +  
  ylab("Shoot:root ratio")+
  geom_text(aes(x=Prov, y=upper.CL+0.05),label = c("a","b","ab"))

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


#####MASSE ERABLES(PROV2018)########
    ####MASSE TOTALE ERABLES####

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

#Estimés et intervalles de confiance et test a posteriori
lsmeans(mod_ERS_mtot, pairwise~Stress*Brout)
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

#####MASSE THUYAS###################
    ####MASSE TOTALE THUYAS####

#Modele avec plan en tiroir
mod_THO_mtot <- lmerTest::lmer(Mtot ~ Stress*Brout*Prov*Hauteurini + (1|Bloc/Stress), data = THO)
summary(mod_THO_mtot)
anova(mod_THO_mtot) # Interaction avec Hauteurini, on garde ce modele
#Effet Stress*Brout

#Estimés et intervalles de confiance et test a posteriori
lsmeans(mod_THO_mtot, pairwise~Stress*Brout)
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

#Estimés et intervalles de confiance et test a posteriori
lsmeans(mod_THO_ratio, pairwise~Stress*Brout*Prov) 
#Rien mais 2018 tend a avoir effet du broutement sous stress2

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


#####MASSE PINS#####################
    ####BIOMASSE TOTALE PINS####

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
#Donc au final, donne le meme resultat 

#Reponse compensatoire
PIB_brout <- filter(PIB, Brout=="Brout")
moyPIB_brout <- mean(PIB_brout$Mtot)
moyPIB_Mret <- mean(PIB_brout$Mret)
PIB_nobrout <- filter(PIB, Brout=="NoBrout")
moyPIB_nobrout <- mean(PIB_nobrout$Mtot)

((moyPIB_brout + moyPIB_Mret) - moyPIB_nobrout)*(100/moyPIB_nobrout) # -46.59


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
mod_PIB_ratio <- lmerTest::lmer(Ratio ~ Stress*Brout*Prov + (1|Bloc/Stress), data = PIB)
summary(mod_PIB_ratio)
anova(mod_PIB_ratio) #Rien significatif

##Verifier Assomptions
plot(mod_PIB_ratio) #Homogeneite des variances ok
qqPlot(resid(mod_PIB_ratio)) #Graphique de normalite ok



###

######
####ANALYSES CHIMIQUES - Exploration des données####
    ####Outliers####
#CET Outliers
boxplot(CET$N, xlab="N CET") #5 donnes plus extremes un peu
dotchart(CET$N,xlab="N CET", ylab = "Ordre des donnees") 
boxplot(CET$Phen, xlab="Phen CET")
dotchart(CET$Phen, xlab="Phen CET", ylab = "Ordre des donnees") 
boxplot(CET$Flav, xlab="Flav CET")
dotchart(CET$Flav,xlab="Flav CET", ylab = "Ordre des donnees") 

#CHR Outliers
boxplot(CHR$N, xlab="N CHR") #5 donnes plus extremes un peu
dotchart(CHR$N,xlab="N CHR", ylab = "Ordre des donnees") 
boxplot(CHR$Phen, xlab="Phen CHR")
dotchart(CHR$Phen, xlab="Phen CHR", ylab = "Ordre des donnees") 
boxplot(CHR$Flav, xlab="Flav CHR") #1 donnee plus extreme un peu
dotchart(CHR$Flav,xlab="Flav CHR", ylab = "Ordre des donnees") 

#ERS Outliers
boxplot(ERS$N, xlab="N ERS") #3 donnes plus extremes un peu
dotchart(ERS$N,xlab="N ERS", ylab = "Ordre des donnees") 
boxplot(ERS$Phen, xlab="Phen ERS") #2 donnes plus extremes un peu
dotchart(ERS$Phen, xlab="Phen ERS", ylab = "Ordre des donnees") 
boxplot(ERS$Flav, xlab="Flav ERS")
dotchart(ERS$Flav,xlab="Flav ERS", ylab = "Ordre des donnees") 

#PIB Outliers
boxplot(PIB$N, xlab="N PIB") 
dotchart(PIB$N,xlab="N PIB", ylab = "Ordre des donnees") 
boxplot(PIB$Phen, xlab="Phen PIB")
dotchart(PIB$Phen, xlab="Phen PIB", ylab = "Ordre des donnees") 
boxplot(PIB$Flav, xlab="Flav PIB")
dotchart(PIB$Flav,xlab="Flav PIB", ylab = "Ordre des donnees") #3 donnes un peu plus extreme

#THO Outliers
boxplot(THO$N, xlab="N THO") #5 donnes plus extremes un peu
dotchart(THO$N,xlab="N THO", ylab = "Ordre des donnees") 
boxplot(THO$Phen, xlab="Phen THO")
dotchart(THO$Phen, xlab="Phen THO", ylab = "Ordre des donnees") 
boxplot(THO$Flav, xlab="Flav THO")
dotchart(THO$Flav,xlab="Flav THO", ylab = "Ordre des donnees") #3 donnees un peu plus extremes

    ####Correlations####
##Général
cor(data[,c(8, 11:14, 17:19)], use = "complete.obs") %>%
  corrplot(method = "shade", type = "lower", addCoef.col= "grey75", 
           tl.col = "black", number.cex = 0.8, tl.cex = 0.8, tl.srt = 45, diag = FALSE)
#Il y a peut-être quelque chose à investiguer ici...
#Un lien entre les phénoliques et la croissance?

#Hauteur initiale selon Phénolique 
ggplot(data, aes(x = Phen, y = Hauteur, color = Esp)) +
  geom_point(size=3) +
  xlab("Concentration en phénoliques (mg/g)") +  
  ylab("Taille du plant (cm)") 
#Semble surtout chez THO
ggplot(data[which(data$Esp == "THO"),], aes(x = Phen, y = Hauteur)) +
  geom_point(size=3) +
  xlab("Concentration en phénoliques (mg/g)") +  
  ylab("Taille du plant (cm)")

#Masse finale totale selon phénolique
ggplot(data, aes(x = Phen, y = Mtot, color = Esp)) +
  geom_point(size=3) +
  xlab("Concentration en phénoliques (mg/g)") +  
  ylab("Masse totale (g)") 

#Hauteur initiale selon flavonoide
ggplot(data, aes(x = Flav, y = Hauteur, color = Esp)) +
  geom_point(size=3) +
  xlab("Concentration en flavonoïdes (mg/g)") +  
  ylab("Taille du plant (cm)") 

#Masse finale totale selon flavonoide
ggplot(data, aes(x = Flav, y = Mtot, color = Esp)) +
  geom_point(size=3) +
  xlab("Concentration en flavonoïdes (mg/g)") +  
  ylab("Masse totale (g)") 

    ####Boxplot####
##CET
boxplot(N ~ Prov, data=CET)
boxplot(Phen ~ Prov, data= CET) #Differences possibles
boxplot(Flav ~ Prov, data= CET)

boxplot(N~ Brout, data=CET)
boxplot(Phen ~ Brout, data=CET)
boxplot(Flav ~ Brout, data=CET)

boxplot(N ~ Stress, data= CET) #Differences possibles
boxplot(Phen ~ Stress, data= CET)
boxplot(Flav ~ Stress, data= CET)

##CHR
boxplot(N ~ Prov, data= CHR) 
boxplot(Phen ~ Prov, data= CHR)
boxplot(Flav ~ Prov, data= CHR)

boxplot(N ~ Brout, data= CHR)
boxplot(Phen ~ Brout, data= CHR)
boxplot(Flav ~ Brout, data= CHR)

boxplot(N ~ Stress, data= CHR) 
boxplot(Phen ~ Stress, data= CHR)
boxplot(Flav ~ Stress, data= CHR)

##ERS
boxplot(N ~ Brout, data= ERS)#Variance hétérogènes?
boxplot(Phen ~ Brout, data= ERS)#Variance hétérogènes?
boxplot(Flav ~ Brout, data= ERS)#Variance hétérogènes?

boxplot(N ~ Stress, data= ERS)
boxplot(Phen ~ Stress, data= ERS)
boxplot(Flav ~ Stress, data= ERS)#Différence possible

##PIB
boxplot(N ~ Prov, data= PIB)#Différence possible
boxplot(Phen ~ Prov, data= PIB)#Variances à surveiller
boxplot(Flav ~ Prov, data= PIB)#Variances à surveiller

boxplot(N ~ Brout, data= PIB)#Différence possible
boxplot(Phen ~ Brout, data= PIB)#Différence possible
boxplot(Flav ~ Brout, data= PIB)#Différence possible

boxplot(N ~ Stress, data= PIB)
boxplot(Phen ~ Stress, data= PIB)
boxplot(Flav ~ Stress, data= PIB)

##THO
boxplot(N ~ Prov, data= cTHO)#Différence possible
boxplot(Phen ~ Prov, data= cTHO)#Différence possible
boxplot(Flav ~ Prov, data= cTHO)#Variance hétérogènes

boxplot(N ~ Brout, data= cTHO)
boxplot(Phen ~ Brout, data= cTHO)
boxplot(Flav ~ Brout, data= cTHO)

boxplot(N ~ Stress, data= cTHO
boxplot(Phen ~ Stress, data= cTHO)
boxplot(Flav ~ Stress, data= cTHO)

####ANALYSES CHIMIQUES CRISIERS#######
    ####AZOTE CERISIERS####

#Modele avec plan en tiroir Azote
CET_N <- lmerTest::lmer( N ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CET_mean)
anova(CET_N) #Effet du stress

##Graphique Azote
CET_N_mod <- summary(lsmeans(CET_N, ~Stress))
ggplot(CET_N_mod, aes(x = Stress, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  xlab("Hydric stress") +  
  ylab("Nitrogen concentration (g/kg)") + 
  ylim(0,30) 

#Conditions d'application
plot(CET_N) #Homogeneite des variances ok
qqPlot(resid(CET_N)) #Graphique de normalite ok


    ####PHENOLS CERISIERS####

##Modèle avec plan en tiroir pour les phénoliques
CET_phen <- lmerTest::lmer( Phen ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CET_mean)
anova(CET_phen) #Rien de significatif

#Conditions d'application
plot(CET_phen) #Homogeneite des variances ok
qqPlot(resid(CET_phen)) #Graphique de normalite ok

    ####FLAVONOIDES CERISIERS####

##Modèle avec plan en tiroir pour les flavonoïdes
CET_flav <- lmerTest::lmer( Flav ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CET_mean)
anova(CET_flav) #Interaction Brout*Prov
summary(lsmeans(CET_flav, pairwise~Brout*Prov))
#Brout 2018 differe de NoBrout 2080

##Graph: veut pas dire grand chose
CET_flav_mod <- summary(lsmeans(CET_flav, ~Brout*Prov))
ggplot(CET_flav_mod, aes(x = Brout, y = lsmean, col=Prov)) +
  geom_point(size=4, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, 
                position=position_dodge(width=0.4))+
  scale_color_manual(values=c("limegreen", "lightcoral", "lightblue"),
                     labels=c("2018", "2050", "2080"))+
  scale_x_discrete(labels=c("0","1"))+
  xlab("Simulated browsing") +  
  ylab("Flavonoids (mg/g)") +
  geom_text(aes(x=Brout, y=upper.CL+0.3),
            label = c("ab","b","ab","ab","a","ab"),  
            position=position_dodge(width=0.4),
            show.legend=F, size=4.5) 

#Conditions d'application
plot(CET_flav) #Homogeneite des variances ok
qqPlot(resid(CET_flav)) #Graphique de normalite ok

####ANALYSES CHIMIQUES CHENES#########
    ####AZOTE CHENES####

##Modèle avec plan en tiroir pour l'azote
CHR_N <- lmerTest::lmer(N ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CHR_mean)
anova(CHR_N) #Effet Stress*brout

#Estimes, intervalles de confiance et Test a posteriori
summary(lsmeans(CHR_N, pairwise~Stress*Brout))
#NoStress-Brout differe de Stress2-Brout 

##Graphique Stress*Brout
CHR_N_mod <- summary(lsmeans(CHR_N, ~Stress*Brout))
ggplot(CHR_N_mod, aes(x = Stress, y = lsmean, col=Brout)) +
  geom_point(size=4, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  scale_color_manual(values=c("limegreen", "lightcoral"),labels=c("Unbrowsed", "Browsed"))+
  scale_x_discrete(labels=c("No stress","High stress")) +
  xlab("Water stress treatment") +  
  ylab("Nitrogen concentration (g/kg)") +
  geom_text(aes(x=Stress, y=upper.CL+0.3),
            label = c("ab","ab","a","b"),  
            position=position_dodge(width=0.4),
            show.legend=F, size=4.5) 

#Conditions d'application
plot(CHR_N) #Homogeneite des variances --> 3 valeurs un peu extremes
qqPlot(resid(CHR_N)) #Graphique de normalite --> 3 valeurs un peu extremes

#Essayer d'enlever 3 valeurs plus extremes (20,21, 36) 
CHR2_mean <- CHR_mean [-36,]
CHR2_mean <- CHR2_mean [-21,]
CHR2_mean <- CHR2_mean[-20,]
dotchart(CHR_mean$N,xlab="N CHR", ylab = "Ordre des donnees") 
dotchart(CHR2_mean$N,xlab="N CHR", ylab = "Ordre des donnees") 

##Modèle avec plan en tiroir pour l'azote (sans outlier)
CHR2_N <- lmerTest::lmer( N ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CHR2_mean)
anova(CHR2_N) #Rien de significatif

#Conditions d'application
plot(CHR2_N) #Homogeneite des variances ok
qqPlot(resid(CHR2_N)) #Graphique de normalite ok

    ####PHENOLS CHENES####

##Modèle avec plan en tiroir pour les phénoliques
CHR_phen <- lmerTest::lmer( Phen ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CHR_mean)
anova(CHR_phen) #Rien de significatif

#Conditions d'application
plot(CHR_phen) #Homogeneite des variances ok
qqPlot(resid(CHR_phen)) #Graphique de normalite ok

    ####FLAVONOIDES CHENES####

##Modèle avec plan en tiroir pour les flavonoïdes
CHR_flav <- lmerTest::lmer( Flav ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CHR2_mean)
anova(CHR_flav) #Rien de significatif

#Conditions d'application
plot(CHR_flav) #Homogeneite des variances ok
qqPlot(resid(CHR_flav)) #Graphique de normalite ok


####ANALSYES CHIMIQUES ERABLES########
    ####AZOTE ERABLES####

##Modèle avec plan en tiroir pour l'azote
ERS_N <- lmerTest::lmer( N ~ Stress*Brout + (1|Bloc/Stress), data = ERS)
anova(ERS_N) #Rien

#Conditions d'application
plot(ERS_N) #Homogeneite des variances ok
qqPlot(resid(ERS_N)) #Graphique de normalite ok

#Essayer d'enlever 3 valeurs plus extremes (8,9, 22) 
ERS2 <- ERS [-22,]
ERS2 <- ERS2 [-9,]
ERS2 <- ERS2 [-8,]
dotchart(ERS$N,xlab="N CHR", ylab = "Ordre des donnees") 
dotchart(ERS2$N,xlab="N CHR", ylab = "Ordre des donnees") 

##Modèle avec plan en tiroir pour l'azote
ERS2_N <- lmerTest::lmer( N ~ Stress*Brout + (1|Bloc/Stress), data = ERS2)
anova(ERS2_N) #Rien

#Conditions d'application
plot(ERS2_N) #Homogeneite des variances ok
qqPlot(resid(ERS2_N)) #Graphique de normalite ok

    ####PHENOLS ERABLES####

##Modèle avec plan en tiroir pour les phénoliques
ERS_phen <- lmerTest::lmer( Phen ~ Stress*Brout + (1|Bloc/Stress), data = ERS)
anova(ERS_phen) #Rien

#Conditions d'application
plot(ERS_phen) #Homogeneite des variances ok
qqPlot(resid(ERS_phen)) #Graphique de normalite ok

    ####FLAVONOIDES ERABLES####

##Modèle avec plan en tiroir pour les flavonoïdes
ERS_flav <- lmerTest::lmer( Flav ~ Stress*Brout + (1|Bloc/Stress), data = ERS)
anova(ERS_flav) #Tendance de stress
lsmeans(ERS_flav, pairwise~Stress)

#Conditions d'application
plot(ERS_flav) #Homogeneite des variances ok
qqPlot(resid(ERS_flav)) #Graphique de normalite ok

####ANALYSES CHIMIQUES PINS###########
    ####AZOTE PINS####

##Modèle avec plan en tiroir pour l'azote
PIB_N <- lmerTest::lmer( N ~ Stress*Brout*Prov + (1|Bloc/Stress), data = PIB)
anova(PIB_N) #Effet Broutement et Provenance

#Test a posteriori Provenance
summary(lsmeans(PIB_N, pairwise~Prov))
#2018 differe de 2050 et 2080

##Graphique Broutement
PIB_N_mod <- summary(lsmeans(PIB_N, ~Brout))
ggplot(PIB_N_mod, aes(x = Brout, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin = lower.CL, ymax=upper.CL), width=.1)+
  xlab("Traitement de broutement") +  
  ylab("Concentration en azote (g/kg)") #Broutés, plus riches en azote

#Graphique Provenance
PIB_N_mod <- summary(lsmeans(PIB_N, ~Prov))
ggplot(PIB_N_mod, aes(x = Prov, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  xlab("Provenance") +  
  ylab("Concentration en azote (g/kg)") #2018, plus riche en azote

#Conditions d'application
plot(PIB_N) #Homogeneite des variances ok
qqPlot(resid(PIB_N)) #Graphique de normalite ok

    ####PHENOLS PINS####

##Modèle avec plan en tiroir pour les phénoliques
PIB_phen <- lmerTest::lmer( Phen ~ Stress*Brout*Prov + (1|Bloc/Stress),data = PIB)
anova(PIB_phen)#Brout

#Estime et intervalles de confiance
lsmeans(PIB_phen, ~Brout)

#Graphique Broutement
PIB_phen_mod <- summary(lsmeans(PIB_phen, ~Brout))
ggplot(PIB_phen_mod, aes(x = Brout, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  xlab("Traitement de broutement") +  
  ylab("Concentration en phénoliques (mg/g)") #Non broutés, plus riches en phénols

#Conditions d'application
plot(PIB_phen) #Homogeneite des variances ok
qqPlot(resid(PIB_phen)) #Graphique de normalite ok

    ####FLAVONOIDES PINS####

##Modèle avec plan en tiroir pour les flavonoïdes
PIB_flav <- lmerTest::lmer( Flav ~ Stress*Brout*Prov + (1|Bloc/Stress), data = PIB)
anova(PIB_flav) #Brout
lsmeans(PIB_flav, ~Brout)

##Graph
PIB_flav_mod <- summary(lsmeans(PIB_flav, ~Brout)) 
ggplot(PIB_flav_mod, aes(x = Brout, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  xlab("Traitement de broutement") +  
  ylab("Concentration en flavonoïdes (mg/g))") #Broutés, plus faibles en flavo
                                               #Différence biologiquement significative?

#Conditions d'application
plot(PIB_flav) #Homogeneite des variances ok
qqPlot(resid(PIB_flav)) #Graphique de normalite ok??????

#Essayer d'enlever 2 valeurs plus extremes (10,42) 
PIB2 <- PIB [-42,]
PIB2 <- PIB2 [-10,]
dotchart(PIB$Flav, xlab="flav PIB", ylab = "Ordre des donnees") 
dotchart(PIB2$Flav, xlab="flav PIB", ylab = "Ordre des donnees") 

##Modèle avec plan en tiroir pour les flavonoïdes
PIB2_flav <- lmerTest::lmer( Flav ~ Stress*Brout*Prov + (1|Bloc/Stress), data = PIB2)
anova(PIB2_flav) #Brout

#Conditions d'application
plot(PIB2_flav) #Homogeneite des variances ok
qqPlot(resid(PIB2_flav)) #Graphique de normalite ok
#Donne les memes resultats

####ANALYSES CHIMIQUES THUYAS#########
    ####AZOTE THUYAS####

##Modèle avec plan en tiroir pour l'azote
THO_N <- lmerTest::lmer( N ~ Stress*Brout*Prov + (1|Bloc/Stress), data = cTHO)
anova(THO_N) #Effet Provenance

#Estimes, intervalle de confiance et Test a posteriori
summary(lsmeans(THO_N, pairwise~Prov))
summary(lsmeans(THO_N, pairwise~Stress*Brout*Prov))

##Graphique Provenance
THO_N_mod <- summary(lsmeans(THO_N, ~Prov)) 
ggplot(THO_N_mod, aes(x = Prov, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  xlab("Provenance") +  

#Conditions d'application
plot(THO_N) #Homogeneite des variances ok
qqPlot(resid(THO_N)) #Graphique de normalite ok

    ####PHENOLS THUYAS####

##Modèle avec plan en tiroir pour les phénoliques
THO_phen <- lmerTest::lmer( Phen ~ Stress*Brout*Prov + (1|Bloc/Stress), data = cTHO)
anova(THO_phen) #Effet provenance, Effet stress

#Estimes, intervalles de confiance et Test a posteriori
lsmeans(THO_phen, pairwise~Stress)
lsmeans(THO_phen, pairwise~Prov)

#Graphique Stress
THO_phen_mod <- summary(lsmeans(THO_phen, ~Stress)) 
ggplot(THO_phen_mod, aes(x = Stress, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  xlab("Traitement de stress hydrique") +  
  ylab("Concentration en phénoliques (mg/g)") 

#Graphique Provenance
THO_phen_mod <- summary(lsmeans(THO_phen, ~Prov)) 
ggplot(THO_phen_mod, aes(x = Prov, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  xlab("Provenance") +  
  ylab("Concentration en phénoliques (mg/g)")#claire diminution

#Conditions d'application
plot(THO_phen) #Homogeneite des variances ok
qqPlot(resid(THO_phen)) #Graphique de normalite ok

    ####FLAVONOIDES THUYAS####

##Modèle avec plan en tiroir pour les flavonoïdes
THO_flav <- lmerTest::lmer( Flav ~ Stress*Brout*Prov + (1|Bloc/Stress), data = cTHO)
anova(THO_flav) #Brout*Prov

#Estimes, intervalles de confiance et Test a posteriori
lsmeans(THO_flav, pairwise~Brout*Prov)
#Rien de significatif selon le test a posteriori

##Graphique Brout*Provenance
THO_flav_mod <- summary(lsmeans(THO_flav, ~Brout*Prov))
ggplot(THO_flav_mod, aes(x = Brout, y = lsmean, col=Prov)) +
  geom_point(size=4, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, 
                position=position_dodge(width=0.4))+
  scale_color_manual(values=c("limegreen", "lightcoral", "lightblue"),
                     labels=c("2018", "2050", "2080"))+
  scale_x_discrete(labels=c("0","1"))+
  xlab("Traitement de broutement") +  
  ylab("Flavonoïdes (mg/g)") +
  theme(text=element_text(size=18, family="sans"), 
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) 

#Conditions d'application
plot(THO_flav) #Homogeneite des variances ok
qqPlot(resid(THO_flav)) #Graphique de normalite ok


####Régressions linéaires croissances et phénoliques####

##CET
summary(lmerTest::lmer(Mtot ~ Phen + (1|Bloc/Stress), data = CET_mean))#Significant
summary(lmerTest::lmer(Mtot ~ Flav + (1|Bloc/Stress), data = CET_mean))

summary(lmerTest::lmer(Maerien ~ Phen + (1|Bloc/Stress), data = CET_mean))#Significant
summary(lmerTest::lmer(Maerien ~ Flav + (1|Bloc/Stress), data = CET_mean))

summary(lmerTest::lmer(Mracine ~ Phen + (1|Bloc/Stress), data = CET_mean))
summary(lmerTest::lmer(Mracine ~ Flav + (1|Bloc/Stress), data = CET_mean))

summary(lmerTest::lmer(Hauteur ~ Phen + (1|Bloc/Stress), data = CET_mean))#Tendency
summary(lmerTest::lmer(Hauteur ~ Flav + (1|Bloc/Stress), data = CET_mean))#Significant


##CHR
summary(lmerTest::lmer(Mtot ~ Phen + (1|Bloc/Stress), data = CHR_mean))
summary(lmerTest::lmer(Mtot ~ Flav + (1|Bloc/Stress), data = CHR_mean))#Significant

summary(lmerTest::lmer(Maerien ~ Phen + (1|Bloc/Stress), data = CHR_mean))
summary(lmerTest::lmer(Maerien ~ Flav + (1|Bloc/Stress), data = CHR_mean))#Significant

summary(lmerTest::lmer(Mracine ~ Phen + (1|Bloc/Stress), data = CHR_mean))
summary(lmerTest::lmer(Mracine ~ Flav + (1|Bloc/Stress), data = CHR_mean))#Tendency

summary(lmerTest::lmer(Hauteur ~ Phen + (1|Bloc/Stress), data = CHR_mean))
summary(lmerTest::lmer(Hauteur ~ Flav + (1|Bloc/Stress), data = CHR_mean))#Significant


##ERS
summary(lmerTest::lmer(Mtot ~ Phen + (1|Bloc/Stress), data = ERS))#Significant
summary(lmerTest::lmer(Mtot ~ Flav + (1|Bloc/Stress), data = ERS))#Significant

summary(lmerTest::lmer(Maerien ~ Phen + (1|Bloc/Stress), data = ERS))#Significant
summary(lmerTest::lmer(Maerien ~ Flav + (1|Bloc/Stress), data = ERS))#Significant

summary(lmerTest::lmer(Mracine ~ Phen + (1|Bloc/Stress), data = ERS))
summary(lmerTest::lmer(Mracine ~ Flav + (1|Bloc/Stress), data = ERS))

summary(lmerTest::lmer(Hauteur ~ Phen + (1|Bloc/Stress), data = ERS))#Significant
summary(lmerTest::lmer(Hauteur ~ Flav + (1|Bloc/Stress), data = ERS))#Significant


##PIB
summary(lmerTest::lmer(Mtot ~ Phen + (1|Bloc/Stress), data = PIB))#Significant
summary(lmerTest::lmer(Mtot ~ Flav + (1|Bloc/Stress), data = PIB))#Significant

summary(lmerTest::lmer(Maerien ~ Phen + (1|Bloc/Stress), data = PIB))#Significant
summary(lmerTest::lmer(Maerien ~ Flav + (1|Bloc/Stress), data = PIB))#Significant

summary(lmerTest::lmer(Mracine ~ Phen + (1|Bloc/Stress), data = PIB))#Significant
summary(lmerTest::lmer(Mracine ~ Flav + (1|Bloc/Stress), data = PIB))#Significant

summary(lmerTest::lmer(Hauteur ~ Phen + (1|Bloc/Stress), data = PIB))#Significant
summary(lmerTest::lmer(Hauteur ~ Flav + (1|Bloc/Stress), data = PIB))#Significant


##THO
summary(lmerTest::lmer(Mtot ~ Phen + (1|Bloc/Stress), data = THO))#Significant
summary(lmerTest::lmer(Mtot ~ Flav + (1|Bloc/Stress), data = THO))

summary(lmerTest::lmer(Maerien ~ Phen + (1|Bloc/Stress), data = THO))#Significant
summary(lmerTest::lmer(Maerien ~ Flav + (1|Bloc/Stress), data = THO))

summary(lmerTest::lmer(Mracine ~ Phen + (1|Bloc/Stress), data = THO))#Significant
summary(lmerTest::lmer(Mracine ~ Flav + (1|Bloc/Stress), data = THO))

summary(lmerTest::lmer(Hauteur ~ Phen + (1|Bloc/Stress), data = THO))#Significant
summary(lmerTest::lmer(Hauteur ~ Flav + (1|Bloc/Stress), data = THO))







