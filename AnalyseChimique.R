####Analyses de composition chimique####

####Packages####
library(tidyr)

####Données####
PIB <- read.table("./PIB.txt", header = TRUE, dec = ",", stringsAsFactors = FALSE)
ERS <- read.table("./ERS.txt", header = TRUE, dec = ",", stringsAsFactors = FALSE)
CET <- read.table("./CET.txt", header = TRUE, dec = ",", stringsAsFactors = FALSE)

PIB <- separate(PIB, ID_plant, into = c("Esp", "Prov", "ID"), sep = "-")
ERS <- separate(ERS, ID_plant, into = c("Esp", "Prov", "ID"), sep = "-")
CET <- separate(CET, ID_plant, into = c("Esp", "Prov", "ID"), sep = "-")

####Graph####

boxplot(N ~ Prov, data= PIB)
boxplot(ERS$N)

boxplot(Phen ~ Prov, data= PIB)
boxplot(ERS$Phen)

boxplot(Flav ~ Prov, data= PIB)
boxplot(ERS$Flav)

boxplot(N ~ Prov, data= CET)
boxplot(Phen ~ Prov, data= CET)
boxplot(Flav ~ Prov, data= CET)
