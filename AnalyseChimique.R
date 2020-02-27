####Analyses de composition chimique####

####Packages####
library(tidyr)

####Données####
PIB <- read.table("./PIB.txt", header = TRUE, dec = ",", stringsAsFactors = FALSE)
ERS <- read.table("./ERS.txt", header = TRUE, dec = ",", stringsAsFactors = FALSE)

PIB <- separate(PIB, ID_plant, into = c("Esp", "Prov", "ID"), sep = "-")
ERS <- separate(ERS, ID_plant, into = c("Esp", "Prov", "ID"), sep = "-")

####Graph####

boxplot(N ~ Prov, data= PIB)
boxplot(ERS$N)

boxplot(Phen ~ Prov, data= PIB)
boxplot(ERS$Phen)

boxplot(Flav ~ Prov, data= PIB)
boxplot(ERS$Flav)
