
# Juillet 2019
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
#library(agricolae)##EC: Ajout pour test de Tukey
library(lme4)##Ajout par Marie-Claude
library(lmerTest)##Ajout par Marie-Claude


#Masses = read.table("C:/Users/Roxanne Turgeon/Dropbox/Initiation recherche/Analyses R/Masses.txt",header=TRUE,na.string = "",dec = ",") 
#Pour moi le point ne fonctionne pas...
##EC: C'est parce que ton fichier n'est pas dans le dossier du projet...
##Pour la suite, je vais utiliser ton nouveau fichier pour avoir les hauteurs initiales
#Masses <- read.table("./Masses.txt", header=TRUE, na.string = "", dec = ",") 
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

boxplot(Mtot~Brout+Esp,data=Masses2)
boxplot(Mtot~Stress+Esp,data=Masses2)
boxplot(Mtot~Bloc*Esp,data=Masses2)

#################CERISIERS########################

###BIOMASSE TOTALE CERISIERS###

#Est-ce que la hauteur initiale est un bon estimateur de la masse totale?
#Regarder la relation entre hauteur finale et la masse a la fin de l'experience
cor.test(CET$Hauteur, CET$Mtot)
plot(CET$Hauteur, CET$Mtot) 
#Oui, bonne correlation pour les cerisiers. 
#Utiliser hauteur initiale pourrait corriger pour la biomasse initiale des plants

#Savoir le nb de replicats 
ddply(CET, c("Stress", "Brout", "Prov", "Bloc"), summarise,N = length(Mtot))

##NOUVEAU CODE: Ajustements par Marie-Claude...il fallait faire une moyenne des hauteurs aussi!
##Ce code est aussi plus élégant que celui que j'avais fait. (Tu peux supprimer ce commentaire après lecture)
#Pour la majorite des combinaisons, nous n'avons qu'un replicat
#Pour simplifier les analyses, on fait une moyenne pour les rares cas 
#ou nous avons 2 plants par combinaison

CET_mean <- CET %>%  group_by(Esp, Prov, Bloc, Brout, Stress) %>%
  summarise(hini=mean(Hauteurini), mtot=mean(Mtot)) 

ddply(CET_mean, c("Stress", "Brout", "Prov", "Bloc"), summarise,
      N = length(mtot)) #Verification n = 1 pour chaque combinaison

#Boxplot masses tot 
levels(CET$Brout) <- c("NoBrout","Brout")
levels(CET$Stress) <- c("NoStress","Stress")
color = c(rep("green",4),rep("yellow",4),rep("red",4))
boxplot(Mtot~Brout+Stress+Prov,data=CET,cex.axis=0.5, col = color,las=2, 
        main = "Masses tot cerisiers", ylab="Masse tot") #Est-ce qu'il faudrait utiliser CET_mean ici aussi?
  ##EC: Le boxplot est surtout une façon de visualiser les données, et non un graphique
    #final. Tu peux utiliser CET_mean ou non. Belle amélioration du graphique! 

##Modele avec plan en tiroir
mod <- aov(Mtot ~ Stress*Brout*Prov +  Error(Bloc + Bloc:Stress), data=CET_mean)   
summary(mod)
model.tables(mod, type="means") #Moyennes estimées par le modèle 
 
  ##EC:Pour répondre à ta question, j'ai consulté un livre (The analysis of Biological Data, Whitlock & Schluter)
    #que j'aimerais te prêter sous peu. C'est un livre pour ceux qui débute en analyses.
    #Dans la section pertinente ('Adjusting for the effects of a covariate'), voici
    #la méthode recommandée (nommée : Analysis of covariance ou ANCOVA):
    #1-Tester un modèle où il y a une interaction avec la covariable (ici: hauteur initiale)
    #2-Si l'interaction n'est pas significative, enlever cette interaction
    #Dans tes données, les interactions ne sont pas significatives, on peut donc l'enlever
    #Le modèle choisi sera donc mod3. Il faudra refaire cette vérification pour
    #chaque espèce
##EC: J'ai mis cette section en commentaire, car nous avons une nouvelle version 
# #Ajout de la hauteur initiale
# mod2 <- aov(Mtot ~ Stress*Brout*Prov*Hauteurini +  Error(Bloc + Bloc:Stress), data=CET_mean)
# #Sans la partie Bloc:Stress car ça empêche l'affichage de l'effet du Stress seul, est-ce normal?
# 
# ##EC: Je crois que c'est parce qu'il n'y a plus de degrés de liberté à l'erreur.
#   #Vérifions si le modèle sans design aucun a une interaction significative avec la covariable
# mod2 <- aov(Mtot ~ Stress*Brout*Prov*Hauteurini, data=CET_mean)
#   summary(mod2)#Non
# 
# # mod2 <- aov(Mtot ~ Stress*Brout*Prov*Hauteurini +  Error(Bloc), data=CET_mean) 
# # summary(mod2) # Effet significatif de Prov, Stress et Hauteur ini
# ##EC: Cette option est à éviter car elle ne tient pas compte du design en tiroir,
#   #mais elle nous permet de voir qu'il n'y a pas d'interactions avec la hauteur initiale
# 
# ##EC: Modèle sélectionné:
# mod3 <- aov(Mtot ~ Stress*Brout*Prov + Hauteurini +  Error(Bloc + Bloc:Stress), data=CET_mean)
# summary(mod3) #Effet du stress, de la provenance et de la hauteur initiale
# 
# plot(CET_mean$Mtot ~ CET_mean$Hauteurini)


#Effet provenance et Stress seulement (est-ce comme ça qu'on doit faire? Car je ne peux faire ce post test
#avec toute l'équation complète)
##EC: non, on ne peux faire ce modèle, parce qu'il ne tient ni compte du design, ni compte des traitements

# boxplot(Mtot~Stress*Prov,data=CET, main = "Masses tot cerisiers", ylab="Masse tot")
# model <- aov(Mtot ~ Stress*Prov,data=CET) ##EC: N'oublie pas d'utiliser CET_mean
# summary(model)
# TukeyHSD(model) #Différence Mtot moins en moins grande avec la provenance
#Si on compare les Mtot sans stress, celle de 2080 est plus faible
#Est-ce qu'on doit corriger pour la hauteur ini? Si oui comment?
##EC: en utilisant le bon modèle et en modifiant ton code
  #J'ai trouvé la solution sur ce site: http://www.personal.psu.edu/mar36/stat_461/split_plot/split_plot.html
  #Ça ne sert à rien de faire un Tukey pour le stress...comme il n'y a que 2 niveaux!
  #Non plus pour Hauteurini, car c'est un facteur continu.
  #Il nous reste le facteur provenance, significatif et à trois niveaux
print(with(CET_mean, HSD.test(Mtot, Prov, DFerror = 39, MSerror = 2873)))
  #Pas de différence a posteriori, probablement à cause de la sévérité du test
#Ce code, plus simple, donne le même test:
lsmeans(mod3, ~ Prov) %>% cld(alpha=0.05) 
#Je ne sais pas encore lequel donne les bons estimés...il va falloir que je regarde

#Conditions d'application ANOVA biomasse totale Cerisiers
plot(mod3) #Normalite & Homogeneite des variances --> Ne marche pas!
##EC: Selon le site http://www.personal.psu.edu/mar36/stat_461/split_plot/split_plot.html
  #on vérifie les suppositions sans l'erreur. Je te laisse tester toi-même si ça fonctionne

shapiro.test(CET_mean$Mtot) #Distribution normale

leveneTest(y=CET_mean$Mtot,group=CET_mean$Brout) #Variances homogenes de Brout
leveneTest(y=CET_mean$Mtot,group=CET_mean$Stress) #Variances homogenes de Stress

#Est-ce qu'on doit vérifier sur les résidus aussi?
##EC: Toujours! Les modèles sont plus résistants à des déviations de la normale à l'hétéroscédasticité
  #À faire avec le bon modèle.
pr <- proj(mod2)                                                                  
res <- pr[["Within"]][,"Residuals"] #Extraire résidus

tapply(res,CET_mean$Brout, shapiro.test) #Normalite Brout ok
tapply(res,CET_mean$Stress, shapiro.test) #Normalite Stress ok

bartlett.test(res~CET_mean$Brout) #Variances homogenes Brout
bartlett.test(res~CET_mean$Stress) # Variance non homogene pour Stress
boxplot(res~CET_mean$Stress+CET_mean$Brout) #on peut verifier avec le boxplot --> Ok

####NOUVEAU CODE: voici le modèle proposé par Marie-Claude####
B <- lmerTest::lmer( mtot ~ Stress*Brout*Prov*hini + (1|Bloc/Stress), data = CET_mean)
#Tu vas voir un message. Il veut simplement dire que pour cette espèce, l'effet aléatoire
#du bloc = 0. On le garde pareil dans le modèle, mais il a peu d'effet
summary(B) #Dans ce summary, tu peux voir toutes les comparaisons a posteriori. 
#Je soupçonne que tu vas peut-être trouver la notation confuse.
#Appelle moi et je t'explique le tout
anova(B)
#Voici une copie d'un échange que j'ai eu avec Marie-Claude:
#Ma question: S’il n’y avait pas d’interaction significative avec hini, est-ce qu’on enlèverait l’interaction?
#Sa réponse: 3.	Dans CET, il y a seulement stress*brout*Hini, mais avec un F de 4.95 
#contre 5.76 pour Hini seul, j’aurais tendance à garder Hini seul 
#(mtot ~ Stress*Brout*Prov+ Hini + (1|Bloc/Stress)). 
#Si les interactions avec Hini ne sont pas ou peu significatif, 
#refaire l’analyse avec +Hini, vérifier si la covariable est significative. 
#Si oui, on la garde et sinon on refait l’analyse sans.

#Donc, on refait sans interaction avec Hini:
B <- lmerTest::lmer( mtot ~ Stress*Brout*Prov + hini + (1|Bloc/Stress), data = CET_mean)
anova(B) #Hini est significatif, on garde donc le modèle ainsi. On a un effet stress
summary(B)

##Pour vérifier tes résidus en incluant les effets aléatoires (soit l'effet du design)
plot(B)
##Et le graphique de normalité
qqnorm(resid(B))


####Ratio CERISIERS####

#Boxplot ratio 
boxplot(Ratio~Brout+Stress+Prov,data=CET,cex.axis=0.5, col = color,las=2, main = "Ratio cerisiers", ylab="Masse tot")

##Modele avec plan en tiroir
mod4 <- aov(Ratio ~ Stress*Brout*Prov +  Error(Bloc + Bloc:Stress), data=CET_mean)   
summary(mod4)
model.tables(mod4, type="means") #Moyennes estimées par le modèle 

#Ajout de la hauteur initiale
mod5 <- aov(Ratio ~ Stress*Brout*Prov*Hauteurini +  Error(Bloc + Bloc:Stress), data=CET_mean)
summary(mod5)
mod5 <- aov(Ratio ~ Stress*Brout*Prov*Hauteurini +  Error(Bloc), data=CET_mean) #Sans la partie Bloc:Stress car ça empêche l'affichage de l'effet du Stress seul
summary(mod5)
#OU
mod6 <- aov(Ratio ~ Stress*Brout*Prov + Hauteurini +  Error(Bloc + Bloc:Stress), data=CET_mean)
summary(mod6) 

#Effet provenance et Stress seulement
boxplot(Ratio~Stress*Prov,data=CET, main = "Ratio cerisiers", ylab="Ratio")
model2 <- aov(Ratio ~ Stress*Prov,data=CET)
summary(model2)
TukeyHSD(model2)

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
boxplot(Mtot~Brout+Stress+Prov,data=CHR,cex.axis=0.5, col = color,las=2, main = "Masses tot chenes", ylab="Masse tot")

##Modele avec plan en tiroir
mod7 <- aov(Mtot ~ Stress*Brout*Prov +  Error(Bloc + Bloc:Stress), data=CHR_mean)   
summary(mod7)
model.tables(mod7, type="means") #Moyennes estimées par le modèle 

#Ajout de la hauteur initiale
mod8 <- aov(Mtot ~ Stress*Brout*Prov*Hauteurini +  Error(Bloc + Bloc:Stress), data=CHR_mean)
summary(mod8)
#OU
mod9 <- aov(Mtot ~ Stress*Brout*Prov + Hauteurini +  Error(Bloc + Bloc:Stress), data=CHR_mean)
summary(mod9) 

#Conditions d'application ANOVA biomasse totale chenes MARCHE PAS
plot(mod7) #Normalite & Homogeneite des variances

leveneTest(y=CHR_mean$Mtot,group=CHR_mean$Brout) #Variances homogenes
leveneTest(y=CHR_mean$Mtot,group=CHR_mean$Stress) #P-value pres de 0.05, on verifie avec le boxplot
boxplot(mo7d$residuals~CHR_mean$Stress+CHR_mean$Brout)

####RATIO CHENES####

#Boxplot ratio 
boxplot(Ratio~Brout+Stress+Prov,data=CHR,cex.axis=0.5, col = color,las=2, main = "Ratio chenes", ylab="Masse tot")

##Modele avec plan en tiroir
mod10 <- aov(Ratio ~ Stress*Brout*Prov +  Error(Bloc + Bloc:Stress), data=CHR_mean)   
summary(mod10)
model.tables(mod10, type="means") #Moyennes estimées par le modèle 

#Ajout de la hauteur initiale
mod11 <- aov(Ratio ~ Stress*Brout*Prov*Hauteurini +  Error(Bloc + Bloc:Stress), data=CHR_mean)
summary(mod11)
mod11 <- aov(Ratio ~ Stress*Brout*Prov*Hauteurini +  Error(Bloc), data=CET_mean) #Sans la partie Bloc:Stress car ça empêche l'affichage de l'effet du Stress seul

#OU
mod12 <- aov(Ratio ~ Stress*Brout*Prov + Hauteurini +  Error(Bloc + Bloc:Stress), data=CHR_mean)
summary(mod12) 

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
boxplot(Mtot~Brout+Stress,data=ERS,cex.axis=0.5, col = color,las=2, main = "Masses tot Erables", ylab="Masse tot")

##Modele avec plan en tiroir
mod13 <- aov(Mtot ~ Stress*Brout +  Error(Bloc + Bloc:Stress), data=ERS)   
summary(mod13)
model.tables(mod13, type="means") #Moyennes estimées par le modèle 

#Ajout de la hauteur initiale
mod14 <- aov(Mtot ~ Stress*Brout*Hauteurini +  Error(Bloc + Bloc:Stress), data=ERS)
summary(mod14)
#OU
mod15 <- aov(Mtot ~ Stress*Brout + Hauteurini +  Error(Bloc + Bloc:Stress), data=ERS)
summary(mod15) 

#Conditions d'application ANOVA biomasse totale Erables MARCHE PAS
plot(mod13) #Normalite & Homogeneite des variances

leveneTest(y=ERS$Mtot,group=ERS$Brout) #Variances homogenes
leveneTest(y=ERS$Mtot,group=ERS$Stress) #P-value pres de 0.05, on verifie avec le boxplot
boxplot(mo13d$residuals~ERS$Stress+ERS$Brout)

####RATIO Erables####

#Boxplot ratio 
boxplot(Ratio~Brout+Stress,data=ERS,cex.axis=0.5, col = color,las=2, main = "Ratio Erables", ylab="Masse tot")

##Modele avec plan en tiroir
mod16 <- aov(Ratio ~ Stress*Brout +  Error(Bloc + Bloc:Stress), data=ERS)   
summary(mod16)
model.tables(mod10, type="means") #Moyennes estimées par le modèle 

#Ajout de la hauteur initiale
mod17 <- aov(Ratio ~ Stress*Brout*Hauteurini +  Error(Bloc + Bloc:Stress), data=ERS)
summary(mod17)
#OU
mod18 <- aov(Ratio ~ Stress*Brout + Hauteurini +  Error(Bloc + Bloc:Stress), data=ERS)
summary(mod18) 

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
boxplot(Mtot~Brout+Stress+Prov,data=THO,cex.axis=0.5, col = color,las=2, main = "Masses tot THUYAS", ylab="Masse tot")

##Modele avec plan en tiroir
mod19 <- aov(Mtot ~ Stress*Brout*Prov +  Error(Bloc + Bloc:Stress), data=THO)   
summary(mod19)
model.tables(mod19, type="means") #Moyennes estimées par le modèle 

#Ajout de la hauteur initiale
mod20 <- aov(Mtot ~ Stress*Brout*Prov*Hauteurini +  Error(Bloc + Bloc:Stress), data=THO)
summary(mod20)
#OU
mod21 <- aov(Mtot ~ Stress*Brout*Prov + Hauteurini +  Error(Bloc + Bloc:Stress), data=THO)
summary(mod21) 

#Conditions d'application ANOVA biomasse totale THUYAS MARCHE PAS
plot(mod19) #Normalite & Homogeneite des variances

leveneTest(y=THO$Mtot,group=THO$Brout) #Variances homogenes
leveneTest(y=THO$Mtot,group=THO$Stress) #P-value pres de 0.05, on verifie avec le boxplot
boxplot(mo19d$residuals~THO$Stress+THO$Brout)

####RATIO THUYAS####

#Boxplot ratio 
boxplot(Ratio~Brout+Stress+Prov,data=THO,cex.axis=0.5, col = color,las=2, main = "Ratio THUYAS", ylab="Masse tot")

##Modele avec plan en tiroir
mod22 <- aov(Ratio ~ Stress*Brout*Prov +  Error(Bloc + Bloc:Stress), data=THO)   
summary(mod22)
model.tables(mod22, type="means") #Moyennes estimées par le modèle 

#Ajout de la hauteur initiale
mod23 <- aov(Ratio ~ Stress*Brout*Prov*Hauteurini +  Error(Bloc + Bloc:Stress), data=THO)
summary(mod23)
#OU
mod24 <- aov(Ratio ~ Stress*Brout*Prov + Hauteurini +  Error(Bloc + Bloc:Stress), data=THO)
summary(mod24) 


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
boxplot(Mtot~Brout+Stress+Prov,data=PIB,cex.axis=0.5, col = color,las=2, main = "Masses tot Pins", ylab="Masse tot")

##Modele avec plan en tiroir
mod25 <- aov(Mtot ~ Stress*Brout*Prov +  Error(BlocBloc:Stress), data=PIB)   
summary(mod25)
model.tables(mod25, type="means") #Moyennes estimées par le modèle 

#Ajout de la hauteur initiale
mod26 <- aov(Mtot ~ Stress*Brout*Prov*Hauteurini +  Error(Bloc + Bloc:Stress), data=PIB)
summary(mod26)
#OU
mod27 <- aov(Mtot ~ Stress*Brout*Prov + Hauteurini +  Error(Bloc + Bloc:Stress), data=PIB)
summary(mod27) 

#Conditions d'application ANOVA biomasse totale Pins MARCHE PAS
plot(mod19) #Normalite & Homogeneite des variances

leveneTest(y=PIB$Mtot,group=PIB$Brout) #Variances homogenes
leveneTest(y=PIB$Mtot,group=PIB$Stress) #P-value pres de 0.05, on verifie avec le boxplot
boxplot(mo19d$residuals~PIB$Stress+PIB$Brout)

####RATIO Pins####

#Boxplot ratio 
boxplot(Ratio~Brout+Stress+Prov,data=PIB,cex.axis=0.5, col = color,las=2, main = "Ratio Pins", ylab="Masse tot")

##Modele avec plan en tiroir
mod28 <- aov(Ratio ~ Stress*Brout*Prov +  Error(Bloc + Bloc:Stress), data=PIB)   
summary(mod28)
model.tables(mod28, type="means") #Moyennes estimées par le modèle 

#Ajout de la hauteur initiale
mod29 <- aov(Ratio ~ Stress*Brout*Prov*Hauteurini +  Error(Bloc + Bloc:Stress), data=PIB)
summary(mod29)
#OU
mod30 <- aov(Ratio ~ Stress*Brout*Prov + Hauteurini +  Error(Bloc + Bloc:Stress), data=PIB)
summary(mod30) 
