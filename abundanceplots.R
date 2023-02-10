library(readxl)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(tidyr)


new_data<-read_excel("/home/enovik/R/Thesis/divided-level-6.xlsx")
View(new_data)

data2<-read_excel("/home/enovik/R/Thesis/V2updated-divided-level-6.xlsx")
View(data2)

data2<-V2divided_level_6

#PHYLUM
phylumdata <- select(data2, "Phylum","Genus", "Ban Max", "BLO 02", "BLO 16", "Clad", "NF F")
#View(phylumdata)

fphylumdata<-phylumdata[!grepl("g__Chloroplast", phylumdata$Genus),]
fphylumdata<-fphylumdata[!grepl("g__Mitochondria", fphylumdata$Genus),]
fphylumdata<-select(fphylumdata, "Phylum", "Ban Max", "BLO 02", "BLO 16", "Clad", "NF F")
View(fphylumdata)



#sum of total amount of sequences in each sample 
sums <- list()
for (i in 2:ncol(fphylumdata)){
  temp<-sum(fphylumdata[,i])
  sums <- c(sums, temp)
  temp<-0
}
print(sums)



for (i in 2:ncol(fphylumdata)){
  print(i)
  for (j in 1:nrow(fphylumdata)){
    fphylumdata[j,i]<-(fphylumdata[j,i]/sums[i-1])*100
  }
}
#combining the same phyla into one row
library(plyr)
fphylumdata<-ddply(fphylumdata,"Phylum",numcolwise(sum))
View(fphylumdata)


#Remove "p__"
fphylumdata$Phylum <- as.character(fphylumdata$Phylum)
fphylumdata$Phylum<- gsub("^.{0,3}", "", fphylumdata$Phylum)
View(fphylumdata)
print(sums)

#fphylumdata<-select(fphylumdata, "Phylum", "BLO 02", "BLO 16", "Clad", "NF F")

names(fphylumdata)[2] <- "BAN_MAX"
names(fphylumdata)[3] <- "BLO_02"
names(fphylumdata)[4] <- "BLO_16"
names(fphylumdata)[6] <- "NF_F"

long_phylum <- gather(fphylumdata, key = "samples", value = "rel_ab",
                     BAN_MAX, BLO_02, BLO_16, Clad, NF_F, convert=FALSE)
View(long_phylum)

phylaplot <- (ggplot(long_phylum, aes(fill=Phylum, y=rel_ab, x=samples))+
        geom_bar(position="stack", stat="identity"))+
  #ggtitle("Relative abundances of microbial phyla")+
  labs(x = "Alga Sample", y = "Relative abundance", fill = "Phyla") +
  theme(legend.title = element_text(size = 2),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.5, "cm"))+
  
  scale_colour_viridis_d(alpha=1, begin = 0, end = 1, direction = 1,
                         option = "D", aesthetics = c("colour","fill"))
phylaplot


#GENUS
genusdata<-select(data2, "Genus", "Ban Max", "BLO 02", "BLO 16", "Clad", "NF F")
#View(genusdata)

#Filter out chloroplast and mitochondria counts
fgenusdata<-genusdata[!grepl("g__Chloroplast", genusdata$Genus),]
fgenusdata<-fgenusdata[!grepl("g__Mitochondria", fgenusdata$Genus),]

View(fgenusdata)

#sum of total amount of sequences in each sample 
sums <- list()
for (i in 2:ncol(fgenusdata)){
  temp<-sum(fgenusdata[,i])
  sums <- c(sums, temp)
  temp<-0
}

for (i in 2:ncol(fgenusdata)){
  print(i)
  for (j in 1:nrow(fgenusdata)){
    fgenusdata[j,i]<-(fgenusdata[j,i]/sums[i-1])*100
  }
}

which(grepl("g__Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", fgenusdata$Genus))
print(fgenusdata[79,1])
#View(fgenusdata)

fgenusdata[79,1]="g__Allorhizobium"

fgenusdata$Genus <- as.character(fgenusdata$Genus)
fgenusdata$Genus<- gsub("^.{0,3}", "", fgenusdata$Genus)
View(fgenusdata)
print(sums)


names(fgenusdata)[2] <- "Ban_Max"
names(fgenusdata)[3] <- "BLO_02"
names(fgenusdata)[4] <- "BLO_16"
names(fgenusdata)[6] <- "NF_F"




long_genus <- gather(fgenusdata, key = "samples", value = "rel_ab",
                     Ban_Max, BLO_02, BLO_16, Clad, NF_F, convert=FALSE)
View(long_genus)


subsetlong <- long_genus[1:100,]
View(subsetlong)

print(ggplot(long_genus, aes(fill=Genus, y=rel_ab, x=samples))+
  geom_bar(position="stack", stat="identity"))+
  xlab("Alga Sample") + ylab("Relative abundance") +
  #ggtitle("Relative abundances of microbial genera")+
  theme(legend.title = element_text(size = 0.25),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.05, "cm"),
        legend.text = element_text(face = "italic"))+
    scale_colour_viridis_d(alpha=1, begin = 0, end = 1, direction = 1,
                        option = "D", aesthetics = c("colour","fill"))
  #theme(legend.position = "none")


#Rel abundance of genera from different phyla: Proteobacteria, Cyanobacteria, Bacteroidota
filterdata <- select(data2, "Phylum","Genus", "Ban Max", "BLO 02", "BLO 16", "Clad", "NF F")

filterdata<-filterdata[!grepl("g__Chloroplast", filterdata$Genus),]
filterdata<-filterdata[!grepl("g__Mitochondria", filterdata$Genus),]

View(filterdata)

proteodata<-filterdata[filterdata$Phylum=="p__Proteobacteria",]
#View(proteodata)

proteodata<-select(proteodata, "Genus", "Ban Max", "BLO 02", "BLO 16", "Clad", "NF F")
sums <- list()
for (i in 2:ncol(proteodata)){
  temp<-sum(proteodata[,i])
  sums <- c(sums, temp)
  temp<-0
}

for (i in 2:ncol(proteodata)){
  print(i)
  for (j in 1:nrow(proteodata)){
    proteodata[j,i]<-(proteodata[j,i]/sums[i-1])*100
  }
}

names(proteodata)[2] <- "Ban_Max"
names(proteodata)[3] <- "BLO_02"
names(proteodata)[4] <- "BLO_16"
names(proteodata)[6] <- "NF_F"



proteodata$Genus <- as.character(proteodata$Genus)
proteodata$Genus<- gsub("^.{0,3}", "", proteodata$Genus)
View(proteodata)

which(grepl("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", proteodata$Genus))
proteodata[13,1]="Allorhizobium"
which(grepl("FukuN57", proteodata$Genus))
proteodata[29,1]="Pseudochelatococcus"
long_proteo <- gather(proteodata, key = "samples", value = "rel_ab",
                     Ban_Max, BLO_02, BLO_16, Clad, NF_F, convert=FALSE)
View(long_proteo)

proteoplot <- (ggplot(long_proteo, aes(fill=Genus, y=rel_ab, x=samples))+
                geom_bar(position="stack", stat="identity"))+
  #ggtitle("Relative abundances of genera from phylum Proteobacteria")+
  labs(x = "Alga Sample", y = "Relative abundance", fill = "Genus") +
  theme(legend.title = element_text(size = 1),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.25, "cm"),
        legend.text = element_text(face = "italic"))+
  
  scale_colour_viridis_d(alpha=1, begin = 0, end = 1, direction = 1,
                         option = "D", aesthetics = c("colour","fill"))
proteoplot

#CYANOBACTERIA

cyanodata<-filterdata[filterdata$Phylum=="p__Cyanobacteria",]
View(cyanodata)

cyanodata<-select(cyanodata, "Genus", "Ban Max", "BLO 02", "BLO 16", "Clad", "NF F")
sums <- list()
for (i in 2:ncol(cyanodata)){
  temp<-sum(cyanodata[,i])
  sums <- c(sums, temp)
  temp<-0
}

for (i in 2:ncol(cyanodata)){
  print(i)
  for (j in 1:nrow(cyanodata)){
    cyanodata[j,i]<-(cyanodata[j,i]/sums[i-1])*100
  }
}

names(cyanodata)[2] <- "Ban_Max"
names(cyanodata)[3] <- "BLO_02"
names(cyanodata)[4] <- "BLO_16"
names(cyanodata)[6] <- "NF_F"



cyanodata$Genus <- as.character(cyanodata$Genus)
cyanodata$Genus<- gsub("^.{0,3}", "", cyanodata$Genus)
#View(proteodata)

cyanodata[is.na(cyanodata)]<-0
View(cyanodata)

#which(grepl("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", proteodata$Genus))
#proteodata[7,1]="Allorhizobium"



long_cyano <- gather(cyanodata, key = "samples", value = "rel_ab",
                      Ban_Max, BLO_02, BLO_16, Clad, NF_F, convert=FALSE)
#View(long_proteo)

cyanoplot <- (ggplot(long_cyano, aes(fill=Genus, y=rel_ab, x=samples))+
                 geom_bar(position="stack", stat="identity"))+
  ggtitle("Relative abundances of genera from phylum Cyanobacteria")+
  labs(x = "Alga Sample", y = "Relative abundance", fill = "Genus") +
  theme(legend.title = element_text(size = 1),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.25, "cm"))+
  
  scale_colour_viridis_d(alpha=1, begin = 0, end = 1, direction = 1,
                         option = "D", aesthetics = c("colour","fill"))
cyanoplot


#BACTEROIDOTA
View(filterdata)
bacterodata<-filterdata[filterdata$Phylum=="p__Bacteroidetes",]
View(bacterodata)

bacterodata<-select(bacterodata, "Genus", "Ban Max", "BLO 02", "BLO 16", "Clad", "NF F")
sums <- list()
for (i in 2:ncol(bacterodata)){
  temp<-sum(bacterodata[,i])
  sums <- c(sums, temp)
  temp<-0
}

for (i in 2:ncol(bacterodata)){
  print(i)
  for (j in 1:nrow(bacterodata)){
    bacterodata[j,i]<-(bacterodata[j,i]/sums[i-1])*100
  }
}

names(bacterodata)[2] <- "Ban_Max"
names(bacterodata)[3] <- "BLO_02"
names(bacterodata)[4] <- "BLO_16"
names(bacterodata)[6] <- "NF_F"



bacterodata$Genus <- as.character(bacterodata$Genus)
bacterodata$Genus<- gsub("^.{0,3}", "", bacterodata$Genus)
View(bacterodata)

#which(grepl("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", proteodata$Genus))
#proteodata[7,1]="Allorhizobium"

long_bactero <- gather(bacterodata, key = "samples", value = "rel_ab",
                      Ban_Max, BLO_02, BLO_16, Clad, NF_F, convert=FALSE)
View(long_bactero)

bacteroplot <- (ggplot(long_bactero, aes(fill=Genus, y=rel_ab, x=samples))+
                 geom_bar(position="stack", stat="identity"))+
  #ggtitle("Relative abundances of genera from phylum Bacteriodota")+
  labs(x = "Alga Sample", y = "Relative abundance", fill = "Genus") +
  theme(legend.title = element_text(size = 1),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.25, "cm"),
        legend.text = element_text(face = "italic"))+
  
  scale_colour_viridis_d(alpha=1, begin = 0, end = 1, direction = 1,
                         option = "D", aesthetics = c("colour","fill"))
bacteroplot

#PLANCTOMYCETES

planctodata<-filterdata[filterdata$Phylum=="p__Planctomycetota",]
View(planctodata)

planctodata<-select(planctodata, "Genus", "Ban Max", "BLO 02", "BLO 16", "Clad", "NF F")
sums <- list()
for (i in 2:ncol(planctodata)){
  temp<-sum(planctodata[,i])
  sums <- c(sums, temp)
  temp<-0
}

for (i in 2:ncol(planctodata)){
  print(i)
  for (j in 1:nrow(planctodata)){
    planctodata[j,i]<-(planctodata[j,i]/sums[i-1])*100
  }
}

View(planctodata)

planctodata[is.na(planctodata)]<-0
View(planctodata)

names(planctodata)[2] <- "Ban_Max"
names(planctodata)[3] <- "BLO_02"
names(planctodata)[4] <- "BLO_16"
names(planctodata)[6] <- "NF_F"



planctodata$Genus <- as.character(planctodata$Genus)
planctodata$Genus<- gsub("^.{0,3}", "", planctodata$Genus)


#which(grepl("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", proteodata$Genus))
#proteodata[7,1]="Allorhizobium"

long_plancto <- gather(planctodata, key = "samples", value = "rel_ab",
                      Ban_Max, BLO_02, BLO_16, Clad, NF_F, convert=FALSE)
View(long_plancto)

planctoplot <- (ggplot(long_plancto, aes(fill=Genus, y=rel_ab, x=samples))+
                 geom_bar(position="stack", stat="identity"))+
  ggtitle("Relative abundances of genera from phylum Planctomycetes")+
  labs(x = "Alga Sample", y = "Relative abundance", fill = "Genus") +
  theme(legend.title = element_text(size = 1),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.25, "cm"))+
  
  scale_colour_viridis_d(alpha=1, begin = 0, end = 1, direction = 1,
                         option = "D", aesthetics = c("colour","fill"))
planctoplot

