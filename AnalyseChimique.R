####Analyses de composition chimique####

####1. Packages####
library(plyr)
library(tidyverse)
library(corrplot)
library(lsmeans)
library(Rmisc)
library(ggplot2)
library(cowplot)
  theme_set(theme_cowplot())
library(car)
library(lme4)

####2. Données et manipulation de données####
PIB <- read.table("./PIB.txt", header = TRUE, dec = ",", stringsAsFactors = FALSE)
ERS <- read.table("./ERS.txt", header = TRUE, dec = ",", stringsAsFactors = FALSE)
CET <- read.table("./CET.txt", header = TRUE, dec = ",", stringsAsFactors = FALSE)
THO <- read.table("./THO.txt", header = TRUE, dec = ",", stringsAsFactors = FALSE)
CHR <- read.table("./CHR.txt", header = TRUE, dec = ",", stringsAsFactors = FALSE)

chem <- rbind(PIB, ERS, CET, THO, CHR)
colnames(chem)[2] <- "ID"
  
data <- read.table("./Masses3.txt", header=TRUE, na.string = "", dec = ",", 
                      stringsAsFactors = FALSE) 

data <- inner_join(data, chem, by = "ID")

data$Esp <- as.factor(data$Esp)
data$Prov <- as.factor(as.character(data$Prov))
data$Brout <- as.factor(as.character(data$Brout))
data$Stress <- as.factor(as.character(data$Stress))

colnames(data)[17] <- "Nitro"

CET <- data[which(data$Esp == "CET"),]
CHR <- data[which(data$Esp == "CHR"),]
ERS <- data[which(data$Esp == "ERS"),]
PIB <- data[which(data$Esp == "PIB"),]
THO <- data[which(data$Esp == "THO"),]
THO <- THO[-which(is.na(THO$Nitro)),]

####3. Exploration des données####

##Général
cor(data[,c(8, 11:14, 17:19)], use = "complete.obs") %>%
  corrplot(method = "shade", type = "lower", addCoef.col= "grey75", 
           tl.col = "black", number.cex = 0.8, tl.cex = 0.8, tl.srt = 45, diag = FALSE)
  #Il y a peut-être quelque chose à investiguer ici...
  #Un lien entre les phénoliques et la croissance?
  ggplot(data, aes(x = Phen, y = Hauteur, color = Esp)) +
    geom_point(size=3) +
    xlab("Concentration en phénoliques (mg/g)") +  
    ylab("Taille du plant (cm)") 
  ggplot(data[which(data$Esp == "THO"),], aes(x = Phen, y = Hauteur)) +
    geom_point(size=3) +
    xlab("Concentration en phénoliques (mg/g)") +  
    ylab("Taille du plant (cm)") 
  ggplot(data, aes(x = Phen, y = Mtot, color = Esp)) +
    geom_point(size=3) +
    xlab("Concentration en phénoliques (mg/g)") +  
    ylab("Masse totale (g)") 
  ggplot(data, aes(x = Flav, y = Hauteur, color = Esp)) +
    geom_point(size=3) +
    xlab("Concentration en flavonoïdes (mg/g)") +  
    ylab("Taille du plant (cm)") 
  ggplot(data, aes(x = Flav, y = Mtot, color = Esp)) +
    geom_point(size=3) +
    xlab("Concentration en flavonoïdes (mg/g)") +  
    ylab("Masse totale (g)") 

##CET
boxplot(Nitro ~ Prov, data= data[which(data$Esp == "CET"),])
boxplot(Phen ~ Prov, data= data[which(data$Esp == "CET"),])
boxplot(Flav ~ Prov, data= data[which(data$Esp == "CET"),])

boxplot(Nitro ~ Brout, data= data[which(data$Esp == "CET"),])
boxplot(Phen ~ Brout, data= data[which(data$Esp == "CET"),])
boxplot(Flav ~ Brout, data= data[which(data$Esp == "CET"),])

boxplot(Nitro ~ Stress, data= data[which(data$Esp == "CET"),])#Différence possible
boxplot(Phen ~ Stress, data= data[which(data$Esp == "CET"),])
boxplot(Flav ~ Stress, data= data[which(data$Esp == "CET"),])

##CHR
boxplot(Nitro ~ Prov, data= data[which(data$Esp == "CHR"),])
boxplot(Phen ~ Prov, data= data[which(data$Esp == "CHR"),])
boxplot(Flav ~ Prov, data= data[which(data$Esp == "CHR"),])

boxplot(Nitro ~ Brout, data= data[which(data$Esp == "CHR"),])
boxplot(Phen ~ Brout, data= data[which(data$Esp == "CHR"),])
boxplot(Flav ~ Brout, data= data[which(data$Esp == "CHR"),])

boxplot(Nitro ~ Stress, data= data[which(data$Esp == "CHR"),])#Différence possible
boxplot(Phen ~ Stress, data= data[which(data$Esp == "CHR"),])
boxplot(Flav ~ Stress, data= data[which(data$Esp == "CHR"),])

##ERS
boxplot(Nitro ~ Brout, data= data[which(data$Esp == "ERS"),])#Variance hétérogènes?
boxplot(Phen ~ Brout, data= data[which(data$Esp == "ERS"),])#Variance hétérogènes?
boxplot(Flav ~ Brout, data= data[which(data$Esp == "ERS"),])#Variance hétérogènes?

boxplot(Nitro ~ Stress, data= data[which(data$Esp == "ERS"),])#Différence possible
boxplot(Phen ~ Stress, data= data[which(data$Esp == "ERS"),])
boxplot(Flav ~ Stress, data= data[which(data$Esp == "ERS"),])#Différence possible

##PIB
boxplot(Nitro ~ Prov, data= data[which(data$Esp == "PIB"),])#Différence possible
boxplot(Phen ~ Prov, data= data[which(data$Esp == "PIB"),])#Variances à surveiller
boxplot(Flav ~ Prov, data= data[which(data$Esp == "PIB"),])#Variances à surveiller

boxplot(Nitro ~ Brout, data= data[which(data$Esp == "PIB"),])#Différence possible
boxplot(Phen ~ Brout, data= data[which(data$Esp == "PIB"),])#Différence possible
boxplot(Flav ~ Brout, data= data[which(data$Esp == "PIB"),])#Différence possible

boxplot(Nitro ~ Stress, data= data[which(data$Esp == "PIB"),])
boxplot(Phen ~ Stress, data= data[which(data$Esp == "PIB"),])
boxplot(Flav ~ Stress, data= data[which(data$Esp == "PIB"),])

##THO
boxplot(Nitro ~ Prov, data= data[which(data$Esp == "THO"),])
boxplot(Phen ~ Prov, data= data[which(data$Esp == "THO"),])#Différence possible
boxplot(Flav ~ Prov, data= data[which(data$Esp == "THO"),])#Variance hétérogènes

boxplot(Nitro ~ Brout, data= data[which(data$Esp == "THO"),])
boxplot(Phen ~ Brout, data= data[which(data$Esp == "THO"),])
boxplot(Flav ~ Brout, data= data[which(data$Esp == "THO"),])

boxplot(Nitro ~ Stress, data= data[which(data$Esp == "THO"),])#Différence possible
boxplot(Phen ~ Stress, data= data[which(data$Esp == "THO"),])
boxplot(Flav ~ Stress, data= data[which(data$Esp == "THO"),])

####4. Analyses statistiques####

####CET####
ddply(CET, c("Stress", "Brout", "Prov", "Bloc"), 
      summarise, N = length(Nitro))

##Combiner pour n'avoir qu'une valeur par combinaison de traitement dans un bloc

CET_mean <- aggregate(list(CET$Nitro, CET$Phen, CET$Flav, CET$Mtot, CET$Maerien, CET$Mracine, CET$Hauteur), 
                      by= data.frame(CET$Stress, CET$Brout, CET$Prov, CET$Bloc), 
                      FUN= "mean")
colnames(CET_mean)[1:11] <- c("Stress", "Brout", "Prov", "Bloc","Nitro", "Phen", "Flav", "Mtot", "Maerien", "Mracine", "Hauteur")
CET_mean <- mutate(CET_mean, 
                   Stress=factor(Stress, levels=c("0", "2")))
  #On indique qu'il n'y a pas de niveau 1 pour le facteur Stress

ddply(CET_mean, c("Stress", "Brout", "Prov", "Bloc"), summarise,
      N = length(Nitro)) #Verification n = 1 pour chaque combinaison

##Modèle avec plan en tiroir pour l'azote
CET_N <- lmerTest::lmer( Nitro ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CET_mean)
  anova(CET_N) #Effet du stress
  summary(lsmeans(CET_N, ~Stress))
  
##Vérifier les conditions d'application
  # plot(CET_N) 
  # leveneTest(resid(CET_N) ~ Stress*Brout*Prov, data = CET_mean) #ok
  # 
  # qqnorm(resid(CET_N))
  # qqline(resid(CET_N))
  # shapiro.test(resid(CET_N)) 
  
##Graph
  temp <- summary(lsmeans(CET_N, ~Stress))
ggplot(temp, aes(x = Stress, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  xlab("Hydric stress") +  
  ylab("Nitrogen concentration (g/kg)") + 
  ylim(0,30) #Plants stressés plus riche en azote, mais pas différent selon lsmean

##Modèle avec plan en tiroir pour les phénoliques
CET_phen <- lmerTest::lmer( Phen ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CET_mean)
anova(CET_phen) #Effet non significatif du stress, de l'interaction stress*brout et de 
          #l'interaction brout*prov

##Vérifier les conditions d'application
plot(CET_phen) 
leveneTest(resid(CET_phen) ~ Stress*Brout*Prov, data = CET_mean) #ok

qqnorm(resid(CET_phen))
qqline(resid(CET_phen))
shapiro.test(resid(CET_phen)) 

##Graph: Rien d'intéressant
# temp <- summary(lsmeans(CET_phen, ~Stress))
# ggplot(temp, aes(x = Stress, y = lsmean)) +
#   geom_point(size=3) +
#   geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
#   xlab("Traitement de stress hydrique") +  
#   ylab("Composés phénoliques (mg/g)") 
# 
# temp <- summary(lsmeans(CET_phen, ~Stress*Brout))
# ggplot(temp, aes(x = Stress, y = lsmean, col=Brout)) +
#   geom_point(size=4, position=position_dodge(width=0.4)) +
#   geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
#   scale_color_manual(values=c("limegreen", "lightcoral"),labels=c("Sans broutement", "Brouté"))+
#   scale_x_discrete(labels=c("0","2"))+
#   scale_y_continuous(limits = c(0,75), expand = expand_scale()) +
#   xlab("Traitement de stress hydrique") +  
#   ylab("Composés phénoliques (mg/g)") +
#   theme(text=element_text(size=18, family="sans"), 
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5))
# 
# temp <- summary(lsmeans(CET_phen, ~Brout*Prov))
# ggplot(temp, aes(x = Brout, y = lsmean, col=Prov)) +
#   geom_point(size=4, position=position_dodge(width=0.4)) +
#   geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, 
#                 position=position_dodge(width=0.4))+
#   scale_color_manual(values=c("limegreen", "lightcoral", "lightblue"),
#                      labels=c("2018", "2050", "2080"))+
#   scale_x_discrete(labels=c("0","1"))+
#   scale_y_continuous(limits = c(0,75), expand = expand_scale()) +
#   xlab("Traitement de broutement") +  
#   ylab("Composés phénoliques (mg/g)") +
#   theme(text=element_text(size=18, family="sans"), 
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5))

##Modèle avec plan en tiroir pour les flavonoïdes
CET_flav <- lmerTest::lmer( Flav ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CET_mean)
anova(CET_flav) #Interaction Brout*Prov
summary(lsmeans(CET_flav, ~Brout*Prov))

##Vérifier les conditions d'application
plot(CET_flav) 
leveneTest(resid(CET_flav) ~ Stress*Brout*Prov, data = CET_mean) #ok

qqnorm(resid(CET_flav))
qqline(resid(CET_flav))
shapiro.test(resid(CET_flav)) 

##Graph: veut pas dire grand chose
temp <- summary(lsmeans(CET_flav, ~Brout*Prov))
ggplot(temp, aes(x = Brout, y = lsmean, col=Prov)) +
  geom_point(size=4, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, 
                position=position_dodge(width=0.4))+
  scale_color_manual(values=c("limegreen", "lightcoral", "lightblue"),
                     labels=c("2018", "2050", "2080"))+
  scale_x_discrete(labels=c("0","1"))+
  xlab("Simulated browsing") +  
  ylab("Flavonoids (mg/g)") 

####CHR####
ddply(CHR, c("Stress", "Brout", "Prov", "Bloc"), 
      summarise, N = length(Nitro))

##Combiner pour n'avoir qu'une valeur par 
#combinaison de traitement dans un bloc

CHR_mean <- aggregate(list(CHR$Nitro, CHR$Phen, CHR$Flav, CHR$Mtot, CHR$Maerien, CHR$Mracine, CHR$Hauteur), 
                      by= data.frame(CHR$Stress, CHR$Brout, CHR$Prov, CHR$Bloc), 
                      FUN= "mean")
colnames(CHR_mean)[1:11] <- c("Stress", "Brout", "Prov", "Bloc","Nitro", "Phen", "Flav", "Mtot", "Maerien", "Mracine", "Hauteur")
CHR_mean <- mutate(CHR_mean, 
                   Stress=factor(Stress, levels=c("0", "2")))
#On indique qu'il n'y a pas de niveau 1 pour le facteur Stress

ddply(CHR_mean, c("Stress", "Brout", "Prov", "Bloc"), summarise,
      N = length(Nitro)) #Verification n = 1 pour chaque combinaison

##Modèle avec plan en tiroir pour l'azote
CHR_N <- lmerTest::lmer( Nitro ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CHR_mean)
anova(CHR_N) #Effet du stress*brout
summary(lsmeans(CHR_N, ~Stress*Brout))

##Vérifier les conditions d'application
plot(CHR_N) 
leveneTest(resid(CHR_N) ~ Stress*Brout*Prov, data = CHR_mean) #ok

qqnorm(resid(CHR_N))
qqline(resid(CHR_N))
shapiro.test(resid(CHR_N)) 
#Correct d'ignorer

##Graph
temp <- summary(lsmeans(CHR_N, ~Stress*Brout))
ggplot(temp, aes(x = Stress, y = lsmean, col=Brout)) +
  geom_point(size=4, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, position=position_dodge(width=0.4))+
  scale_color_manual(values=c("limegreen", "lightcoral"),labels=c("Sans broutement", "Brouté"))+
  scale_x_discrete(labels=c("0","2")) +
  xlab("Traitement de stress hydrique") +  
  ylab("Concentration azote (g/kg)") 

##Modèle avec plan en tiroir pour les phénoliques
CHR_phen <- lmerTest::lmer( Phen ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CHR_mean)
anova(CHR_phen) #Effet non significatif de l'interaction stress*brout

##Vérifier les conditions d'application
plot(CHR_phen) 
leveneTest(resid(CHR_phen) ~ Stress*Brout*Prov, data = CHR_mean) #ok

qqnorm(resid(CHR_phen))
qqline(resid(CHR_phen))
shapiro.test(resid(CHR_phen)) 

##Modèle avec plan en tiroir pour les flavonoïdes
CHR_flav <- lmerTest::lmer( Flav ~ Stress*Brout*Prov + (1|Bloc/Stress), data = CHR_mean,
                     control = lmerControl(optimizer= "bobyqa"))
  anova(CHR_flav)#Nothing

##Vérifier les conditions d'application
plot(CHR_flav)
leveneTest(resid(CHR_flav) ~ Stress*Brout*Prov, data = CHR_mean) #ok

qqnorm(resid(CHR_flav))
qqline(resid(CHR_flav))
shapiro.test(resid(CHR_flav))


####ERS####
ddply(ERS, c("Stress", "Brout", "Bloc"), 
      summarise, N = length(Nitro))

##Modèle avec plan en tiroir pour l'azote
ERS_N <- lmerTest::lmer( Nitro ~ Stress*Brout + (1|Bloc/Stress), data = ERS)
anova(ERS_N) #Rien

##Vérifier les conditions d'application
plot(ERS_N) 
leveneTest(resid(ERS_N) ~ Stress*Brout, data = ERS) #ok...

qqnorm(resid(ERS_N))
qqline(resid(ERS_N))
shapiro.test(resid(ERS_N)) 

##Modèle avec plan en tiroir pour les phénoliques
ERS_phen <- lmerTest::lmer( Phen ~ Stress*Brout + (1|Bloc/Stress), data = ERS)
anova(ERS_phen) #Rien

##Vérifier les conditions d'application
plot(ERS_phen) 
leveneTest(resid(ERS_phen) ~ Stress*Brout, data = ERS) #ok

qqnorm(resid(ERS_phen))
qqline(resid(ERS_phen))
shapiro.test(resid(ERS_phen)) 

##Modèle avec plan en tiroir pour les flavonoïdes
ERS_flav <- lmerTest::lmer( Flav ~ Stress*Brout + (1|Bloc/Stress), 
                     control = lmerControl(optimizer= "bobyqa"), data = ERS)
anova(ERS_flav)#Tendence, stress

##Vérifier les conditions d'application
plot(ERS_flav)
leveneTest(resid(ERS_flav) ~ Stress*Brout*Prov, data = ERS) #ok

qqnorm(resid(ERS_flav))
qqline(resid(ERS_flav))
shapiro.test(resid(ERS_flav))

####PIB####
ddply(PIB, c("Stress", "Brout", "Prov", "Bloc"), 
      summarise, N = length(Nitro))

##Modèle avec plan en tiroir pour l'azote
PIB_N <- lmerTest::lmer( Nitro ~ Stress*Brout*Prov + (1|Bloc/Stress), data = PIB)
anova(PIB_N) #Effet du broutement et de la provenance
summary(lsmeans(PIB_N, ~Brout))
summary(lsmeans(PIB_N, ~Prov))

##Vérifier les conditions d'application
plot(PIB_N) 
leveneTest(resid(PIB_N) ~ Stress*Brout*Prov, data = PIB) #ok

qqnorm(resid(PIB_N))
qqline(resid(PIB_N))
shapiro.test(resid(PIB_N)) 

##Graph
temp <- summary(lsmeans(PIB_N, ~Brout))
ggplot(temp, aes(x = Brout, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin = lower.CL, ymax=upper.CL), width=.1)+
  xlab("Traitement de broutement") +  
  ylab("Concentration en azote (g/kg)") #Broutés, plus riches en azote

temp <- summary(lsmeans(PIB_N, ~Prov))
ggplot(temp, aes(x = Prov, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  xlab("Provenance") +  
  ylab("Concentration en azote (g/kg)") #2018, plus riche en azote

##Modèle avec plan en tiroir pour les phénoliques
PIB_phen <- lmerTest::lmer( Phen ~ Stress*Brout*Prov + (1|Bloc/Stress), 
                     control = lmerControl(optimizer= "bobyqa"), data = PIB)

anova(PIB_phen)#Brout
summary(lsmeans(PIB_phen, ~Brout))

##Vérifier les conditions d'application
plot(PIB_phen)
leveneTest(resid(PIB_phen) ~ Stress*Brout*Prov, data = PIB) #ok

qqnorm(resid(PIB_phen))
qqline(resid(PIB_phen))
shapiro.test(resid(PIB_phen))

temp <- summary(lsmeans(PIB_phen, ~Brout))
ggplot(temp, aes(x = Brout, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  xlab("Traitement de broutement") +  
  ylab("Concentration en phénoliques (mg/g)") #Broutés, plus riches en azote

##Modèle avec plan en tiroir pour les flavonoïdes
PIB_flav <- lmerTest::lmer( Flav ~ Stress*Brout*Prov + (1|Bloc/Stress), data = PIB)
anova(PIB_flav) #Brout
summary(lsmeans(PIB_flav, ~Brout))

##Vérifier les conditions d'application
plot(PIB_flav) 
leveneTest(resid(PIB_flav) ~ Stress*Brout*Prov, data = PIB) #ok

qqnorm(resid(PIB_flav))
qqline(resid(PIB_flav))
shapiro.test(resid(PIB_flav)) 

##Graph
temp <- summary(lsmeans(PIB_flav, ~Brout)) 
ggplot(temp, aes(x = Brout, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  xlab("Traitement de broutement") +  
  ylab("Concentration en flavonoïdes (mg/g))") #Broutés, plus faibles en flavo
  #Différence biologiquement significative?

####THO####
ddply(THO, c("Stress", "Brout", "Prov", "Bloc"), 
      summarise, N = length(Nitro))

##Modèle avec plan en tiroir pour l'azote
THO_N <- lmerTest::lmer( Nitro ~ Stress*Brout*Prov + (1|Bloc/Stress), data = THO)
anova(THO_N) #Effet de la provenance, tendance interaction triple et stress
summary(lsmeans(THO_N, ~Stress))
summary(lsmeans(THO_N, ~Prov))
summary(lsmeans(THO_N, ~Stress*Brout*Prov))

##Vérifier les conditions d'application
plot(THO_N) 
leveneTest(resid(THO_N) ~ Stress*Brout*Prov, data = THO)

qqnorm(resid(THO_N))
qqline(resid(THO_N))
shapiro.test(resid(THO_N)) 

##Graph
temp <- summary(lsmeans(THO_N, ~Prov)) 
ggplot(temp, aes(x = Prov, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  xlab("Provenance") +  
  ylab("Concentration en azote (g/kg)")

temp <- summary(lsmeans(THO_N, ~Stress*Brout*Prov))
ggplot(temp, aes(x = Stress, y = lsmean, shape=factor(Prov),color=Brout)) +
  geom_point(size=3, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax= upper.CL), width=.1, 
                position=position_dodge(width=0.4))+
  xlab("Traitement de stress hydrique") +  
  ylab("Concentration en azote (g/kg)")#Complicated

##Modèle avec plan en tiroir pour les phénoliques
THO_phen <- lmerTest::lmer( Phen ~ Stress*Brout*Prov + (1|Bloc/Stress), data = THO)
anova(THO_phen)#Effet provenance, Effet stress

##Vérifier les conditions d'application
plot(THO_phen)
leveneTest(resid(THO_phen) ~ Stress*Brout*Prov, data = THO) #ok

qqnorm(resid(THO_phen))
qqline(resid(THO_phen))
shapiro.test(resid(THO_phen))

temp <- summary(lsmeans(THO_phen, ~Stress)) 
ggplot(temp, aes(x = Stress, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  xlab("Traitement de stress hydrique") +  
  ylab("Concentration en phénoliques (mg/g)") 


temp <- summary(lsmeans(THO_phen, ~Prov)) 
ggplot(temp, aes(x = Prov, y = lsmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1)+
  xlab("Provenance") +  
  ylab("Concentration en phénoliques (mg/g)")#claire diminution

##Modèle avec plan en tiroir pour les flavonoïdes
THO_flav <- lmerTest::lmer( Flav ~ Stress*Brout*Prov + (1|Bloc/Stress), data = THO)
anova(THO_flav) #Brout*Prov
summary(lsmeans(THO_flav, ~Brout*Prov))

##Vérifier les conditions d'application
plot(THO_flav) 
leveneTest(resid(THO_flav) ~ Stress*Brout*Prov, data = THO) #ok

qqnorm(resid(THO_flav))
qqline(resid(THO_flav))
shapiro.test(resid(THO_flav)) 

##Graph
temp <- summary(lsmeans(THO_flav, ~Brout*Prov))
ggplot(temp, aes(x = Brout, y = lsmean, col=Prov)) +
  geom_point(size=4, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.1, 
                position=position_dodge(width=0.4))+
  scale_color_manual(values=c("limegreen", "lightcoral", "lightblue"),
                     labels=c("2018", "2050", "2080"))+
  scale_x_discrete(labels=c("0","1"))+
  xlab("Traitement de broutement") +  
  ylab("Flavonoïdes (mg/g)") +
  theme(text=element_text(size=18, family="sans"), 
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))#Bof

####5. Régressions linéaires croissances et phénoliques####

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
  
