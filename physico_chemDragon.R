######LOAD LIBRARIES ######
require(destiny)
require(dynamicTreeCut)
require(WGCNA)
require(xlsx)
library(umap)
library(rococo)

####### DEFINE FUNCTIONS ######
require(ggplot2)
plot.std <- function(df, xname, yname,legend, title){
  
  ggplot(df, aes(x=x, y=y, color=z))+
    geom_point(size=5)+
    xlab(xname)+
    ylab(yname)+
    ggtitle(title)+
    theme(axis.title.x = element_text(face="bold", size=30, vjust=-2),
          axis.text.x  = element_text( size=25),
          axis.title.y = element_text(face="bold", size=30,vjust=2),
          axis.text.y  = element_text( size=30),
          plot.margin=unit(c(1,1,1.5,1.2),"cm"),
          legend.text=element_text(size=15),#size of legend
          legend.title=element_text(size=15),
          plot.title = element_text(lineheight=.8, face="bold",size=30)) +
    scale_color_discrete(name=legend)#it depends which variable the variable in the legend is mapped to
  
  
}


#############################################


#####LOAD DATA ##########

setwd("~/Desktop/SpatialTranscriptomics_USB/newData/current/REALIGN/3D/")
#OR data 
or.info<-read.table("../supp_tables/allORs_data.csv", sep=",", header = T, row.names = 1)
#row.names(or.info)<-as.character(or.info[,1])


#fill the maxTopic info for the non-detected genes
predicted<-row.names(or.info)[is.na(or.info$maxTopic) &
                                !is.na(or.info$DPTindex)]#find genes that were predicted

plot(or.info[predicted, "DPTindex"],
     or.info[predicted, "indexMiyamichiReal"])#doesn't look that great...

detected<-row.names(or.info)[!is.na(or.info$maxTopic)]

#correaltion dpt index and myamichi "real"
plot(or.info[detected, "DPTindex"],
     or.info[detected, "indexMiyamichiReal"])
points(or.info[predicted, "DPTindex"],
     or.info[predicted, "indexMiyamichiReal"],col="red")#THIS MIGHT BE A BETTER WAY TO VISUALIZE THE PREDICTED GENES

#correaltion dpt index and myamichi predicted in chem senses
plot(or.info[detected, "DPTindex"],
     or.info[detected, "indexMiyamichiChemS"])
points(or.info[predicted, "DPTindex"],
       or.info[predicted, "indexMiyamichiChemS"],col="red")#THIS MIGHT BE A BETTER WAY TO VISUALIZE THE PREDICTED GENES



#look at the relationship between the DPTindex (and myamichi index) and the topics
pdf(file = "../figures_pdf2/F6TopicsVsMiyamichi.pdf", width=7, height = 7)
boxplot(or.info[detected, "indexMiyamichiReal"]~or.info[detected, "maxTopic"], 
        xlab="Topic", ylab="Myamichi index")
dev.off()

pdf(file = "../figures_pdf2/F6TopicsVsDoBs.pdf", width=7, height = 7)
boxplot(or.info[detected, "DPTindex"]~or.info[detected, "maxTopic"], 
        xlab="Topic", ylab="DPT index")
dev.off()

#convert index into topic: consider the median DPT for each topic, take the closest 
median.topic<-vector(length=5)
median.topic[1]<-median(or.info[row.names(or.info)[!is.na(or.info$maxTopic)&
                                                     or.info$maxTopic==1], "DPTindex"])
median.topic[2]<-median(or.info[row.names(or.info)[!is.na(or.info$maxTopic)&
                                                     or.info$maxTopic==2], "DPTindex"])
median.topic[3]<-median(or.info[row.names(or.info)[!is.na(or.info$maxTopic)&
                                                     or.info$maxTopic==3], "DPTindex"])
median.topic[4]<-median(or.info[row.names(or.info)[!is.na(or.info$maxTopic)&
                                                     or.info$maxTopic==4], "DPTindex"])
median.topic[5]<-median(or.info[row.names(or.info)[!is.na(or.info$maxTopic)&
                                                     or.info$maxTopic==5], "DPTindex"])

predicted.topic<-vector(length=length(predicted))
names(predicted.topic)<-predicted
predicted.topic<-unlist(sapply(predicted, function(x){
  ind<-or.info[x, "DPTindex"]
  diff<-abs(median.topic-ind)
  res<-which(diff==min(diff))
  return(res)
  } ))
table(predicted.topic)

#predicted topic vs class
table(predicted.topic, or.info[names(predicted.topic), "Class"])

#fill in the maxTopic with the predicted topic for the predicted ORs; and attach a column with predicted/non-predicted OR
or.info<-cbind(or.info, 
               Predicted=row.names(or.info)%in%predicted)
or.info[predicted, "maxTopic"]<-predicted.topic

#NEW
#Physico-chem properties of odourants 
physico.chem<-read.csv(file="../supp_tables/ORLigandsDragon.txt", sep="\t")[,-1]
physico.chem<-physico.chem[-which(is.na(physico.chem$MW)),]
physico.chem<-physico.chem[!duplicated(physico.chem$NAME),]
row.names(physico.chem)<-physico.chem$NAME

#Physico-chem properties of odourants 
#physico.chem<-read.csv(file="../supp_tables/physChem4model.csv")
#row.names(physico.chem)<-physico.chem$X
physico.chem<-physico.chem[,-1]
physico.chem<-physico.chem[!duplicated(physico.chem),]
#physico.chem<-apply(physico.chem, 2, as.numeric)
temp<-apply(physico.chem, 2, sd)

id<-which(temp!=0)# remove those properties that are constant for all the odourants in the dataset
physico.chem<-physico.chem[,id]

#DOB of odourants 
all.dob<-read.csv(file="../supp_tables/odourantsDoBs4model.csv")
row.names(all.dob)<-all.dob$X
all.dob<-all.dob[,-c(1,2)]

#Table Odourants vs ORs 
#NEW
or.odourant<-read.csv(file="../supp_tables/mouse_OR_ligand_ pairs_Feb2021.csv")
id<-which(or.odourant$Odorant.CAS.Number!="N/A")
or.odourant<-or.odourant[id,]

#tmp<-read.csv("../supp_tables/cas2smilesV2.csv")
#cas2smiles<-tmp$Odorant.CAS.Number
#names(cas2smiles)<-tmp$SMILES

#cas4physC<-cas2smiles[rownames(physico.chem)]
#rownames(physico.chem)<-cas4physC

#physico.chem<-physico.chem[-which(rownames(physico.chem)%in%c("76748-69-1", "6745-75-1", "1222-05-5", "77-54-3", "81-14-1")),] #remove outliers (PC1 > 150, grey cluster 2nd and 3rd attempt or PC1 > 30 from the beginning)
#physico.chem<-physico.chem[-which(rownames(physico.chem)%in%c("60763-41-9", "6745-75-1", "1222-05-5", "77-54-3", "81-14-1")),] #remove outliers (PC1 > 150, grey cluster 2nd and 3rd attempt or PC1 > 30 from the beginning)
physico.chem<-physico.chem[!duplicated(physico.chem),]

temp<-apply(physico.chem, 2, sd)

id<-which(temp!=0)# remove those properties that are constant for all the odourants in the dataset
physico.chem<-physico.chem[,id]

#or.odourant<-read.csv(file="../supp_tables/ORs_odourants.csv",row.names=1)
#id<-which(or.odourant$data.Detected.Odorant.CAS.Number!="N/A")
#or.odourant<-or.odourant[id,]

sort(table(or.odourant[,2]), decrease=T)
sort(table(or.odourant[,3]), decrease=T)

#Solubility of odourants 
#Downloaded from http://www.vcclab.org/lab/alogps/
solub<-read.table(file="../supp_tables/solubility.tsv", sep=" ", comment.char = "")
row.names(solub)<-solub[,1]


#####PCA OF ODOURANTS #######

pca<-prcomp(physico.chem, scale. = T)

summary(pca)

df<-data.frame(x=pca$x[,1], y=pca$x[,2], z=as.factor(all.dob[row.names(physico.chem),"max_zone"]))

plot.std(df, "PC1", "PC2", legend="", title="")
plot(df$x, df$y, xlab="PC1", ylab="PC2", pch=19)

#pdf(file = "../figures_pdf2/F6PCAwOutliersV3.pdf", width=7, height = 7)
plot(df$x, df$y, xlab="PC1", ylab="PC2", pch=19)
#abline(v = 50, col="red")
#dev.off()
#plot.std(df[-which(df$x>150),], "PC1", "PC2", legend="", title="")

#plot.std(df[-which(df$x>25),], "PC1", "PC2", legend="", title="")
#####DIFFUSION MAP OF ODOURANTS#######
require(destiny)

data.scaled<-apply(physico.chem, 2, function(x) (x-mean(x))/sd(x))

dm<-DiffusionMap(data=data.scaled)

plot(dm$DC1, dm$DC2, col=as.factor(all.dob[row.names(physico.chem),"max_zone"]))
dpt<-DPT(dm)

df<-data.frame(x=dm$DC1, y=dm$DC2, z=as.factor(all.dob[row.names(physico.chem),"max_zone"]))

plot.std(df, "DC1", "DC2", legend="", title="")




######CHECK THE CLASSES OF ORS DETECTING ODOURANTS WITH HIGH PC2 #######
id<-names(sort(pca$x[,2], decreasing = T))[1:10]
id1<-which(or.odourant$Odorant.CAS.Number%in%id)
names.or<-unique(or.odourant[id1,"Olfr.Nomenclature"])
or.info[names.or,"Class"]



######HIGHLIGHT ODOURANTS DETECTED BY CLASS I ORS #######


temp<-or.odourant
temp<-cbind(temp, class=or.info[temp$Olfr.Nomenclature, "Class"])

id<-which(temp$class==1)#find odourants detected by at least 1 class-1-OR
check.od<-unique(temp[id,"Odorant.CAS.Number"])

table(all.dob[check.od,"max_zone"])#what is the maxTopic of these odourants?

#highlight "class I odourants" in the PCA
df<-data.frame(x=pca$x[,1], y=pca$x[,2], 
               z=row.names(physico.chem)%in%check.od)

plot.std(df, "PC1", "PC2", legend="", title="")

boxplot(solub[check.od,3], 
        solub[setdiff(row.names(physico.chem), check.od),3], 
        ylab="Solubility", names = c("Class I", "All"))#plot solubilities







#CHECK DISTRIBUTION OF ODOURANTS DETECTED BY ORS IN EACH TOPIC IN THE PCA SPACE #######

#select the ORs that are included in the OR/odourant table
or.included<-intersect(unique(or.odourant[,2]), row.names(or.info)[!is.na(or.info$maxTopic)])

table(or.info[or.included,"maxTopic"])

odourant.topic<-list()

topics<-1:5
for(t in topics){
  id<-or.included[or.info[or.included,"maxTopic"]==t]
  odourant.topic[[t]]<-unique(or.odourant[which(or.odourant$Olfr.Nomenclature%in%id),"Odorant.CAS.Number"])
  
}


for(topic in 1:5){## plot the odourant in the PCA or diff map colored according to whether they are detected in topic 1-5 or not

  label=vector(length=nrow(physico.chem))
  
  names(label)<-row.names(physico.chem)
  label[intersect(names(label),odourant.topic[[topic]])]<-"TRUE"
  
  df<-data.frame(x=pca$x[,1], y=pca$x[,2], z=label)
  
  
  p1<-plot.std(df, "PC1", "PC2", legend="", title=paste0("Topic ", topic))
  print(p1)
  
}






####CLUSTERING OF ODOURANTS #########

data.scaled<-apply(physico.chem, 2, function(x) (x-mean(x))/sd(x))
#standardize data
dist.mat<-dist(data.scaled)
odourant.clust<-hclust(dist.mat, method="ward.D2")#"average")
cut2<-cutreeDynamic(odourant.clust,distM=as.matrix(dist.mat),
                    minClusterSize=10, method="hybrid",deepSplit = 1)
clust.colour<-labels2colors(cut2)
names(clust.colour)<-row.names(physico.chem)
df<-data.frame(x=pca$x[,1], y=pca$x[,2], z=clust.colour)#plot clusters
#plot.std(df, "PC1", "PC2", legend="", title="Clusters")
table(clust.colour)

OdsUmap<-umap(data.scaled) 

#pdf(file = "../figures_pdf2/F6A.2.pdf", width=7, height = 7) #min clust size = 5
#pdf(file = "../figures_pdf2/F6A.2V3.pdf", width=7, height = 7) #min clust size = 5
pdf(file = "../figures_pdf2/F6A.2V4.pdf", width=7, height = 7) #min clust size = 10, method clustering = ward
#pdf(file = "../figures_pdf2/F6A.3.pdf", width=7, height = 7) #min clust size = 10
#pdf(file = "../figures_pdf2/F6A.4.pdf", width=7, height = 7) #min clust size = 3
plot.std(df, "PC1", "PC2", legend="", title="Clusters")

plot(df$x, df$y, main="Odourants clusters", xlab="PC1", ylab="PC2", col=clust.colour, pch=19)
plot(OdsUmap$layout[,1], OdsUmap$layout[,2], col=clust.colour, pch=19, xlab="UMAP1", ylab="UMAP2")
plot(dm$DC1, dm$DC2, col=clust.colour, pch=19, xlab="DC1", ylab="UDC2")

dev.off()

#exclude grey cluster
#physico.chem<-physico.chem[-which(clust.colour=="grey"),]
#clust.colour=clust.colour[-which(clust.colour=="grey")]
#####CLUSTERS OF ODOURANTS VS TOPICS ######

#create a table with cluster of odourants vs topic of odourant
clust.vs.max.topic<-table(clust.colour[row.names(all.dob)], all.dob$max_zone)
chisq.test(clust.vs.max.topic)

#create a table with cluster of odourants vs whether odourant is detected in topic X or not ##########

#Blue cluster ####
od<-names(clust.colour)[clust.colour=="blue"]

or<-unique(or.odourant[or.odourant$Odorant.CAS.Number%in%od,"Olfr.Nomenclature"])
top<-or.info[or, "maxTopic"]
top<-unique(top[!is.na(top)])
top

#Brown cluster ####
od<-names(clust.colour)[clust.colour=="brown"]

or<-unique(or.odourant[or.odourant$Odorant.CAS.Number%in%od,"Olfr.Nomenclature"])
top<-or.info[or, "maxTopic"]
top<-unique(top[!is.na(top)])
top

#green cluster #####
od<-names(clust.colour)[clust.colour=="grey"]

or<-unique(or.odourant[or.odourant$Odorant.CAS.Number%in%od,"Olfr.Nomenclature"])
top<-or.info[or, "maxTopic"]
top<-unique(top[!is.na(top)])
top

#turquoise cluster #####
od<-names(clust.colour)[clust.colour=="turquoise"]

or<-unique(or.odourant[or.odourant$Odorant.CAS.Number%in%od,"Olfr.Nomenclature"])
top<-or.info[or, "maxTopic"]
top<-unique(top[!is.na(top)])
top

#yellow cluster ######
od<-names(clust.colour)[clust.colour=="yellow"]

or<-unique(or.odourant[or.odourant$Odorant.CAS.Number%in%od,"Olfr.Nomenclature"])
top<-or.info[or, "maxTopic"]
top<-unique(top[!is.na(top)])
top

#green cluster ######
od<-names(clust.colour)[clust.colour=="green"]

or<-unique(or.odourant[or.odourant$Odorant.CAS.Number%in%od,"Olfr.Nomenclature"])
top<-or.info[or, "maxTopic"]
top<-unique(top[!is.na(top)])
top

#red cluster ######
od<-names(clust.colour)[clust.colour=="red"]

or<-unique(or.odourant[or.odourant$Odorant.CAS.Number%in%od,"Olfr.Nomenclature"])
top<-or.info[or, "maxTopic"]
top<-unique(top[!is.na(top)])
top

#pink cluster ######
od<-names(clust.colour)[clust.colour=="pink"]

or<-unique(or.odourant[or.odourant$Odorant.CAS.Number%in%od,"Olfr.Nomenclature"])
top<-or.info[or, "maxTopic"]
top<-unique(top[!is.na(top)])
top

#black cluster ######
od<-names(clust.colour)[clust.colour=="black"]

or<-unique(or.odourant[or.odourant$Odorant.CAS.Number%in%od,"Olfr.Nomenclature"])
top<-or.info[or, "maxTopic"]
top<-unique(top[!is.na(top)])
top

#magenta cluster ######
od<-names(clust.colour)[clust.colour=="magenta"]

or<-unique(or.odourant[or.odourant$Odorant.CAS.Number%in%od,"Olfr.Nomenclature"])
top<-or.info[or, "maxTopic"]
top<-unique(top[!is.na(top)])
top

####

airMucusPCtbl<-read.csv("../supp_tables/Unique_CAS_numbers_v4.csv")
airMucusPC<-airMucusPCtbl$Air.mucus.partition.coefficient..log..B..
#airMucusPC<-as.numeric(gsub("−", "-", airMucusPC))
names(airMucusPC)<-airMucusPCtbl$Odorant.CAS.Number

ods_2ORsPlus<-table(or.odourant$Odorant.CAS.Number)
ods_2ORsPlus<-names(ods_2ORsPlus)[which(ods_2ORsPlus>1)]

set.seed(10)
odsDPT<-data.frame(or.info[or.odourant$Olfr.Nomenclature,]$DPTindex_3T, 
                   or.info[or.odourant$Olfr.Nomenclature,]$DPTindex_4T, 
                   or.info[or.odourant$Olfr.Nomenclature,]$DPTindex_6T, 
                   or.info[or.odourant$Olfr.Nomenclature,]$DPTindex, sample(or.info[or.odourant$Olfr.Nomenclature,]$DPTindex), or.odourant$Olfr.Nomenclature)
colnames(odsDPT)<-c("DPTindex_3T", "DPTindex_4T", "DPTindex_6T", "DPTindex", "DPTindex_rand", "Olfr.Nomenclature")
odsDPT$Odorant.CAS.Number<-or.odourant$Odorant.CAS.Number
odsDPT_all<-odsDPT

class1ORs<-rownames(or.info)[which(or.info$Class==1)]
RFpredicted<-rownames(or.info)[which(!is.na(or.info$index_RFpredicted))]
odsDPTnoC1<-odsDPT[-which(odsDPT$Olfr.Nomenclature %in% class1ORs),]
odsDPT_pred<-odsDPT[-which(odsDPT$Olfr.Nomenclature %in% RFpredicted),]
odsDPT_NotPred<-odsDPT[which(odsDPT$Olfr.Nomenclature %in% RFpredicted),]

length(which(odsDPT$Odorant.CAS.Number %in% ods_2ORsPlus))
length(which(odsDPT_pred$Odorant.CAS.Number %in% ods_2ORsPlus))
length(which(odsDPT_NotPred$Odorant.CAS.Number %in% ods_2ORsPlus))
length(which(odsDPTnoC1$Odorant.CAS.Number %in% ods_2ORsPlus))

c1ods<-or.odourant[which(or.odourant$Olfr.Nomenclature %in% class1ORs),]
c2ods<-unique(or.odourant$Odorant.CAS.Number[-which(or.odourant$Odorant.CAS.Number %in% c1ods$Odorant.CAS.Number)])

set.seed(10)
odsDPTrand<-tapply(odsDPT$DPTindex_rand, odsDPT$Odorant.CAS.Number, function(x){mean(x, na.rm=T)})
odsDPTrand<-odsDPTrand[-which(is.na(odsDPTrand))]


odsDPT<-tapply(odsDPT$DPTindex, odsDPT$Odorant.CAS.Number, function(x){mean(x, na.rm=T)})
odsDPT<-odsDPT[-which(is.na(odsDPT))]

odsDPT_3T<-tapply(odsDPT_all$DPTindex_3, odsDPT_all$Odorant.CAS.Number, function(x){mean(x, na.rm=T)})
odsDPT_3T<-odsDPT_3T[-which(is.na(odsDPT_3T))]
odsDPT_4T<-tapply(odsDPT_all$DPTindex_4T, odsDPT_all$Odorant.CAS.Number, function(x){mean(x, na.rm=T)})
odsDPT_4T<-odsDPT_4T[-which(is.na(odsDPT_4T))]
odsDPT_6T<-tapply(odsDPT_all$DPTindex_6T, odsDPT_all$Odorant.CAS.Number, function(x){mean(x, na.rm=T)})
odsDPT_6T<-odsDPT_6T[-which(is.na(odsDPT_6T))]

odsDPTnoC1<-tapply(odsDPTnoC1$DPTindex, odsDPTnoC1$Odorant.CAS.Number, function(x){mean(x, na.rm=T)})
odsDPTnoC1<-odsDPTnoC1[-which(is.na(odsDPTnoC1))]

odsDPT_pred<-tapply(odsDPT_pred$DPTindex, odsDPT_pred$Odorant.CAS.Number, function(x){mean(x, na.rm=T)})
odsDPT_pred<-odsDPT_pred[-which(is.na(odsDPT_pred))]

odsDPT_NotPred<-tapply(odsDPT_NotPred$DPTindex, odsDPT_NotPred$Odorant.CAS.Number, function(x){mean(x, na.rm=T)})
#odsDPT_NotPred<-odsDPT_NotPred[-which(is.na(odsDPT_NotPred))]

plotXYcorrelation<-function(x, y, main, xlab, ylab, col, ylim=c(min(y, na.rm=T), max(y, na.rm=T)), xlim=c(min(x, na.rm=T), max(x, na.rm=T))){
  #x=values in x axis, y=values in y axis, xlab=x axis name, ylab=y axis name, col=vector of colors per sample
  plot(x, y, main=main, ylab=ylab, xlab=xlab, col=alpha(col, 0.3), cex.main=1, cex=0.8, pch=19, ylim=ylim, xlim=xlim)
  #plot(x, y, main=main, ylab=ylab, xlab=xlab, col=col, cex.main=1, cex=0.8, pch=19, ylim=ylim, xlim=xlim)
  abline(lm(y~x))
  mtext(c(paste("R=", round(cor.test(x, y, method="spearman")$estimate, digits=2), "\t", "\t", "p=", cor.test(x, y, method="spearman")$p.value)), side=3, cex=0.6)
}


plotXYcorrelation2<-function(x, y, main, xlab, ylab, col, ylim=c(min(y, na.rm=T), max(y, na.rm=T)), xlim=c(min(x, na.rm=T), max(x, na.rm=T))){
  #x=values in x axis, y=values in y axis, xlab=x axis name, ylab=y axis name, col=vector of colors per sample
  plot(x, y, main=main, ylab=ylab, xlab=xlab, col=alpha(col, 0.3), cex.main=1, cex=0.8, pch=19, ylim=ylim, xlim=xlim)
  abline(lm(y~x))
  mtext(c(paste("R=", round(cor.test(x, y, method="pearson")$estimate, digits=2), "\t", "\t", "p=", cor.test(x, y, method="pearson")$p.value)), side=3, cex=0.6)
}


odsMaxTopic<-data.frame(or.odourant$Olfr.Nomenclature, or.info[or.odourant$Olfr.Nomenclature,]$maxTopic)
odsMaxTopic$Odorant.CAS.Number<-or.odourant$Odorant.CAS.Number
colnames(odsMaxTopic)<-c("Olfr.Nomenclature", "maxTopic", "Odorant.CAS.Number")
odsMaxTopic<-tapply(odsMaxTopic$maxTopic, odsMaxTopic$Odorant.CAS.Number, function(x){mean(x, na.rm=T)})
odsMaxTopic<-odsMaxTopic[-which(is.na(odsMaxTopic))]
odsMaxTopic<-round(odsMaxTopic)



#pdf(file = "../figures_pdf2/F6_airMucusPC.pdf", width=7, height = 7)
pdf(file = "../figures_pdf2/F7_airMucusPCnoC1_supp.pdf", width=7, height = 7)

plotXYcorrelation(x=odsDPTnoC1, y=airMucusPC[names(odsDPTnoC1)], xlab = "avgDPT", 
                  ylab="airMucusCoeff", col=1, main="no Class1 ORs")

plotXYcorrelation(x=odsDPTnoC1[ods_2ORsPlus], y=airMucusPC[ods_2ORsPlus], xlab = "avgDPT", 
                  ylab="airMucusCoeff", col=1, main="Ods detected by 2 or more non C1 ORs")

plotXYcorrelation(x=odsDPT[c2ods], y=airMucusPC[c2ods], xlab = "avgDPT", 
                  ylab="airMucusCoeff", col=1, main="Ods detected by non C1 ORs")

plotXYcorrelation(x=odsDPT[c2ods[which(c2ods %in% ods_2ORsPlus)]], y=airMucusPC[c2ods[which(c2ods %in% ods_2ORsPlus)]], xlab = "avgDPT", 
                  ylab="airMucusCoeff", col=1, main="Ods detected by non C1 ORs and more than 2 ORs")

dev.off()

pdf(file = "../figures_pdf2/FS5_airMucusPCnoC1_rev.pdf", width=4, height = 4)
topicCol_noC1<-odsMaxTopic[ods_2ORsPlus]+1
topicCol_noC1[which(topicCol_noC1==5)]<-7
plotXYcorrelation(x=odsDPTnoC1[ods_2ORsPlus], y=airMucusPC[ods_2ORsPlus], xlab = "avgDPT", 
                  ylab="airMucusCoeff", col=topicCol_noC1, main="No C1 ORs")
legend("bottomright", c("T1", "T2", "T3", "T4", "T5"), pch=19, col=alpha(c(2, 3, 4, 7, 6), 0.3))
dev.off()

pdf(file = "../figures_pdf2/F6_airMucusPC_supp_rev.pdf", width=7, height = 7)

#plotXYcorrelation(x=odsDPT_3T, y=airMucusPC[names(odsDPT_3T)], xlab = "avgDPT", 
#                  ylab="airMucusCoeff", col=1, main="all data 3 topics")

#plotXYcorrelation(x=odsDPT_4T, y=airMucusPC[names(odsDPT_4T)], xlab = "avgDPT", 
#                  ylab="airMucusCoeff", col=1, main="all data 4 topics")

#plotXYcorrelation(x=odsDPT_6T, y=airMucusPC[names(odsDPT_6T)], xlab = "avgDPT", 
#                  ylab="airMucusCoeff", col=1, main="all data 6 topics")

#plotXYcorrelation(x=odsDPT, y=airMucusPC[names(odsDPT)], xlab = "avgDPT", 
#                  ylab="airMucusCoeff", col=1, main="all data")

#plotXYcorrelation(x=odsDPTrand, y=airMucusPC[names(odsDPTrand)], xlab = "avgDPT", 
#                  ylab="airMucusCoeff", col=1, main="all data randomized")

#plotXYcorrelation2(x=odsDPT, y=airMucusPC[names(odsDPT)], xlab = "avgDPT", 
#                  ylab="airMucusCoeff", col=1, main="all data")

topicCol<-odsMaxTopic[names(odsDPT)]+1
topicCol[which(topicCol==5)]<-7
#text(airMucusPC[names(odsDPT)]~odsDPT, labels=names(odsDPT),cex=0.3, font=2)

#plotXYcorrelation(x=odsDPT, y=airMucusPC[names(odsDPT)], xlab = "avgDPT", 
#                  ylab="airMucusCoeff", col=alpha(topicCol, 0.3), main="all data")
#legend("bottomleft", c("T1", "T2", "T3", "T4", "T5"), pch=19, col=alpha(c(2, 3, 4, 7, 6), 0.3))

### all data
#plotXYcorrelation2(x=odsDPT, y=airMucusPC[names(odsDPT)], xlab = "avgDPT", 
#                   ylab="airMucusCoeff", col=odsMaxTopic[names(odsDPT)], main="all data")

#### Just experimentally got coefficients

exp_coeffs<-airMucusPCtbl$Odorant.CAS.Number[which(airMucusPCtbl$logKow..0.experimental.database..1.KowWin.estimated.==1)]
airMucusPC_exp<-airMucusPC[exp_coeffs]
  
#plotXYcorrelation(x=odsDPT[exp_coeffs], y=airMucusPC[exp_coeffs], xlab = "avgDPT", 
#                  ylab="airMucusCoeff", col=1, main="experimentally obtained coeffs.")

#### Just DE ORs
or.info_DE<-or.info[which(is.na(or.info$index_RFpredicted)),]
odsDPT_DE<-data.frame(rownames(or.info_DE)[or.odourant$Olfr.Nomenclature], 
                        or.info_DE[or.odourant$Olfr.Nomenclature,]$DPTindex)
colnames(odsDPT_DE)[1:2]<-c("Olfr.Nomenclature", "DPTindex")

odsDPT_DE$Odorant.CAS.Number<-or.odourant$Odorant.CAS.Number
odsDPT_DE<-tapply(odsDPT_DE$DPTindex, odsDPT_DE$Odorant.CAS.Number, function(x){mean(x, na.rm=T)})
odsDPT_DE<-odsDPT_DE[-which(is.na(odsDPT_DE))]

#plotXYcorrelation(x=odsDPT_DE, y=airMucusPC[names(odsDPT_DE)], xlab = "avgDPT", 
#                  ylab="airMucusCoeff", col=1, main="Differentially expressed ORs")

#odsDPT_DE2<-odsDPT_DE[which(names(odsDPT_DE) %in% ods_2ORsPlus)]
odsDPT_DE2<-odsDPT_DE[ods_2ORsPlus]
topicCol_s<-odsMaxTopic[names(odsDPT_DE2)]+1
topicCol_s[which(topicCol_s==5)]<-7
plotXYcorrelation(x=odsDPT_DE2, y=airMucusPC[names(odsDPT_DE2)], xlab = "avgDPT", 
                  ylab="airMucusCoeff", col=alpha(topicCol_s, 0.3), main="Not predicted ORs (Ods 2+ORs)")
legend("bottomleft", c("T1", "T2", "T3", "T4", "T5"), pch=19, col=alpha(c(2, 3, 4, 7, 6), 0.3))

#### Just predicted ORs
or.info_noDE<-or.info[which(!is.na(or.info$index_RFpredicted)),]
odsDPT_noDE<-data.frame(rownames(or.info_noDE)[or.odourant$Olfr.Nomenclature], 
                      or.info_noDE[or.odourant$Olfr.Nomenclature,]$DPTindex)
colnames(odsDPT_noDE)[1:2]<-c("Olfr.Nomenclature", "DPTindex")

odsDPT_noDE$Odorant.CAS.Number<-or.odourant$Odorant.CAS.Number
odsDPT_noDE<-tapply(odsDPT_noDE$DPTindex, odsDPT_noDE$Odorant.CAS.Number, function(x){mean(x, na.rm=T)})
odsDPT_noDE<-odsDPT_noDE[-which(is.na(odsDPT_noDE))]

odsDPT_noDE2<-odsDPT_noDE[which(names(odsDPT_noDE) %in% ods_2ORsPlus)]
topicCol_s2<-odsMaxTopic[names(odsDPT_noDE2)]+1
topicCol_s2[which(topicCol_s2==5)]<-7
plotXYcorrelation(x=odsDPT_noDE2, y=airMucusPC[names(odsDPT_noDE2)], xlab = "avgDPT", 
                  ylab="airMucusCoeff", col=alpha(topicCol_s2, 0.3), main="Predicted ORs (Ods 2+ORs)")
legend("bottomright", c("T1", "T2", "T3", "T4", "T5"), pch=19, col=alpha(c(2, 3, 4, 7, 6), 0.3))



dev.off()

### Just odourants detected by more than one OR

pdf(file = "../figures_pdf2/F6_airMucusPC.pdf", width=4, height = 4)

plotXYcorrelation(x=odsDPT[ods_2ORsPlus], y=airMucusPC[ods_2ORsPlus], xlab = "avgDPT", 
                  ylab="airMucusCoeff", col=1, main="Odourants detected by 2 or more ORs")

#plotXYcorrelation2(x=odsDPT[ods_2ORsPlus], y=airMucusPC[ods_2ORsPlus], xlab = "avgDPT", 
#                   ylab="airMucusCoeff", col=1, main="Odourants detected by 2 or more ORs")

topicCol<-odsMaxTopic[ods_2ORsPlus]+1
topicCol[which(topicCol==5)]<-7
plotXYcorrelation(x=odsDPT[ods_2ORsPlus], y=airMucusPC[ods_2ORsPlus], xlab = "avgDPT", 
                  ylab="airMucusCoeff", col=alpha(topicCol, 0.3), main="Odourants detected by 2 or more ORs")
legend("bottomright", c("T1", "T2", "T3", "T4"), pch=19, col=alpha(c(2, 3, 4, 7), 0.3))
text(airMucusPC[ods_2ORsPlus]~odsDPT[ods_2ORsPlus], labels=ods_2ORsPlus,cex=0.3, font=2)

plotXYcorrelation(x=odsDPT[ods_2ORsPlus], y=airMucusPC[ods_2ORsPlus], xlab = "average DPT", 
                  ylab="air/mucus partition Coeff. (log10)", col=alpha(topicCol, 0.3), main="Odourants detected by 2 or more ORs")
legend("bottomright", c("T1", "T2", "T3", "T4"), pch=19, col=alpha(c(2, 3, 4, 7), 0.3), cex=0.5)

#plotXYcorrelation2(x=odsDPT[ods_2ORsPlus], y=airMucusPC[ods_2ORsPlus], xlab = "avgDPT", 
#                   ylab="airMucusCoeff", col=odsMaxTopic[ods_2ORsPlus], main="Odourants detected by 2 or more ORs")
#legend("bottomleft", c("T1", "T2", "T3", "T4"), pch=19, col=1:4)

boxplot(airMucusPC[ods_2ORsPlus], main="air/mucus partition coeff.")

dev.off()

#tmp<-data.frame(mean3Dind=odsDPT[ods_2ORsPlus], airMucusPC=airMucusPC[ods_2ORsPlus])
#write.csv(tmp, "../supp_tables/Fig5C_CAS.csv")

airMucusData<-data.frame(meanDPT=odsDPT[ods_2ORsPlus], airMucusPCoeff=airMucusPC[ods_2ORsPlus], meanDPTrand=odsDPTrand[ods_2ORsPlus])
rownames(airMucusData)<-ods_2ORsPlus
write.csv(airMucusData, "../supp_tables/airMucusPC_odsDPT.csv")

airMucusData_pred<-data.frame(meanDPT=odsDPT_pred[ods_2ORsPlus], airMucusPCoeff=airMucusPC[ods_2ORsPlus])
rownames(airMucusData_pred)<-ods_2ORsPlus
airMucusData_pred<-airMucusData_pred[-which(is.na(airMucusData_pred$meanDPT)),]
write.csv(airMucusData_pred, "../supp_tables/airMucusPC_odsDPT_pred.csv")

airMucusData_NotPred<-data.frame(meanDPT=odsDPT_NotPred[ods_2ORsPlus], airMucusPCoeff=airMucusPC[ods_2ORsPlus])
rownames(airMucusData_NotPred)<-ods_2ORsPlus
airMucusData_NotPred<-airMucusData_NotPred[-which(is.na(airMucusData_NotPred$meanDPT)),]
write.csv(airMucusData_NotPred, "../supp_tables/airMucusPC_odsDPT_NotPred.csv")

pdf(file = "../figures_pdf2/F6_airMucusPC_rev.pdf", width=4, height = 4)

topicCol<-odsMaxTopic[ods_2ORsPlus]+1
topicCol[which(topicCol==5)]<-7
plotXYcorrelation(x=airMucusData$meanDPT, y=airMucusData$airMucusPCoeff, xlab = "avgDPT", 
                  ylab="airMucusCoeff", col=alpha(topicCol, 0.3), main="Predicted + not predicted")
legend("bottomright", c("T1", "T2", "T3", "T4", "T5"), pch=19, col=alpha(c(2, 3, 4, 7, 6), 0.3))

plotXYcorrelation(x=airMucusData_pred$meanDPT, y=airMucusData_pred$airMucusPCoeff, xlab = "avgDPT", 
                  ylab="airMucusCoeff", col=alpha(topicCol[rownames(airMucusData_pred)], 0.3), main="Predicted")
legend("bottomright", c("T1", "T2", "T3", "T4", "T5"), pch=19, col=alpha(c(2, 3, 4, 7, 6), 0.3))

plotXYcorrelation(x=airMucusData_NotPred$meanDPT, y=airMucusData_NotPred$airMucusPCoeff, xlab = "avgDPT", 
                  ylab="airMucusCoeff", col=alpha(topicCol[rownames(airMucusData_NotPred)], 0.3), main="Not predicted")
legend("bottomright", c("T1", "T2", "T3", "T4", "T5"), pch=19, col=alpha(c(2, 3, 4, 7, 6), 0.3))


dev.off()

pdf(file = "../figures_pdf2/FS5_airMucusPC_NotPred_rev.pdf", width=4, height = 4)

plotXYcorrelation(x=airMucusData_NotPred$meanDPT, y=airMucusData_NotPred$airMucusPCoeff, xlab = "avgDPT", 
                  ylab="airMucusCoeff", col=alpha(topicCol[rownames(airMucusData_NotPred)], 0.3), main="Not predicted")
legend("bottomright", c("T1", "T2", "T3", "T4", "T5"), pch=19, col=alpha(c(2, 3, 4, 7, 6), 0.3))

dev.off()
#### DPT of ORs detecting the same odourant vs others

connected<-matrix(0, length(unique(or.odourant$Olfr.Nomenclature)), 
                  length(unique(or.odourant$Olfr.Nomenclature)))
colnames(connected)<-unique(or.odourant$Olfr.Nomenclature)
rownames(connected)<-unique(or.odourant$Olfr.Nomenclature)

for(i in rownames(connected)){
  for(j in colnames(connected)){
    if(length(intersect(or.odourant$Odorant.CAS.Number[which(or.odourant$Olfr.Nomenclature==i)], 
                        or.odourant$Odorant.CAS.Number[which(or.odourant$Olfr.Nomenclature==j)]))>0){
      connected[i, j] <- 1
    }
  }
}


DPTdiffs<-matrix(NA, length(unique(or.odourant$Olfr.Nomenclature)), 
                  length(unique(or.odourant$Olfr.Nomenclature)))
colnames(DPTdiffs)<-unique(or.odourant$Olfr.Nomenclature)
rownames(DPTdiffs)<-unique(or.odourant$Olfr.Nomenclature)

for(i in rownames(DPTdiffs)){
  for(j in colnames(DPTdiffs)){
    DPTdiffs[i, j] <- abs(or.info[i, "DPTindex"] - or.info[j, "DPTindex"])
  }
}

connectedORsDPTdiffs<-DPTdiffs[lower.tri(DPTdiffs)][which(connected[lower.tri(connected)]==1)]
NotConnectedORsDPTdiffs<-DPTdiffs[lower.tri(DPTdiffs)][which(connected[lower.tri(connected)]==0)]

pdf(file = "../figures_pdf2/F6_DPTdiffs.pdf", width=7, height = 7)
boxplot(connectedORsDPTdiffs, NotConnectedORsDPTdiffs, names=c("connectedORs", "notConnectedORs"), ylab="absolute DPTindexes difference")
dev.off()
wilcox.test(connectedORsDPTdiffs, NotConnectedORsDPTdiffs)

### NO C1 ORs

connected_noC1<-matrix(0, length(unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[-which(rownames(or.info) %in% class1ORs)])])), 
                    length(unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[-which(rownames(or.info) %in% class1ORs)])])))
colnames(connected_noC1)<-unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[-which(rownames(or.info) %in% class1ORs)])])
rownames(connected_noC1)<-unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[-which(rownames(or.info) %in% class1ORs)])])

for(i in rownames(connected_noC1)){
  for(j in colnames(connected_noC1)){
    if(length(intersect(or.odourant$Odorant.CAS.Number[which(or.odourant$Olfr.Nomenclature==i)], 
                        or.odourant$Odorant.CAS.Number[which(or.odourant$Olfr.Nomenclature==j)]))>0){
      connected_noC1[i, j] <- 1
    }
  }
}

DPTdiffs_noC1<-matrix(NA, length(unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[-which(rownames(or.info) %in% class1ORs)])])), 
                   length(unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[-which(rownames(or.info) %in% class1ORs)])])))
colnames(DPTdiffs_noC1)<-unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[-which(rownames(or.info) %in% class1ORs)])])
rownames(DPTdiffs_noC1)<-unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[-which(rownames(or.info) %in% class1ORs)])])


for(i in rownames(DPTdiffs_noC1)){
  for(j in colnames(DPTdiffs_noC1)){
    DPTdiffs_noC1[i, j] <- abs(or.info[i, "DPTindex"] - or.info[j, "DPTindex"])
  }
}

connectedORsDPTdiffs_noC1<-DPTdiffs_noC1[lower.tri(DPTdiffs_noC1)][which(connected_noC1[lower.tri(connected_noC1)]==1)]
NotConnectedORsDPTdiffs_noC1<-DPTdiffs_noC1[lower.tri(DPTdiffs_noC1)][which(connected_noC1[lower.tri(connected_noC1)]==0)]


pdf(file = "../figures_pdf2/FS5_DPTdiffs_noC1_rev.pdf", width=7, height = 7)
boxplot(connectedORsDPTdiffs_noC1, NotConnectedORsDPTdiffs_noC1, names=c("connectedORs", "notConnectedORs"), ylab="absolute DPTindexes difference")
dev.off()
wilcox.test(connectedORsDPTdiffs_noC1, NotConnectedORsDPTdiffs_noC1)


### Predicted ORs

connected_p<-matrix(0, length(unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[which(!is.na(or.info$index_RFpredicted))])])), 
                  length(unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[which(!is.na(or.info$index_RFpredicted))])])))
colnames(connected_p)<-unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[which(!is.na(or.info$index_RFpredicted))])])
rownames(connected_p)<-unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[which(!is.na(or.info$index_RFpredicted))])])

for(i in rownames(connected_p)){
  for(j in colnames(connected_p)){
    if(length(intersect(or.odourant$Odorant.CAS.Number[which(or.odourant$Olfr.Nomenclature==i)], 
                        or.odourant$Odorant.CAS.Number[which(or.odourant$Olfr.Nomenclature==j)]))>0){
      connected_p[i, j] <- 1
    }
  }
}

DPTdiffs_p<-matrix(NA, length(unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[which(!is.na(or.info$index_RFpredicted))])])), 
                    length(unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[which(!is.na(or.info$index_RFpredicted))])])))
colnames(DPTdiffs_p)<-unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[which(!is.na(or.info$index_RFpredicted))])])
rownames(DPTdiffs_p)<-unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[which(!is.na(or.info$index_RFpredicted))])])


for(i in rownames(DPTdiffs_p)){
  for(j in colnames(DPTdiffs_p)){
    DPTdiffs_p[i, j] <- abs(or.info[i, "DPTindex"] - or.info[j, "DPTindex"])
  }
}

connectedORsDPTdiffs_p<-DPTdiffs_p[lower.tri(DPTdiffs_p)][which(connected_p[lower.tri(connected_p)]==1)]
NotConnectedORsDPTdiffs_p<-DPTdiffs_p[lower.tri(DPTdiffs_p)][which(connected_p[lower.tri(connected_p)]==0)]


pdf(file = "../figures_pdf2/F6_DPTdiffs_Pred_rev.pdf", width=7, height = 7)
boxplot(connectedORsDPTdiffs_p, NotConnectedORsDPTdiffs_p, names=c("connectedORs", "notConnectedORs"), ylab="absolute DPTindexes difference")
dev.off()
wilcox.test(connectedORsDPTdiffs_p, NotConnectedORsDPTdiffs_p)

### Not Predicted ORs

connected_np<-matrix(0, length(unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[which(is.na(or.info$index_RFpredicted))])])), 
                    length(unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[which(is.na(or.info$index_RFpredicted))])])))
colnames(connected_np)<-unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[which(is.na(or.info$index_RFpredicted))])])
rownames(connected_np)<-unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[which(is.na(or.info$index_RFpredicted))])])

for(i in rownames(connected_np)){
  for(j in colnames(connected_np)){
    if(length(intersect(or.odourant$Odorant.CAS.Number[which(or.odourant$Olfr.Nomenclature==i)], 
                        or.odourant$Odorant.CAS.Number[which(or.odourant$Olfr.Nomenclature==j)]))>0){
      connected_np[i, j] <- 1
    }
  }
}

DPTdiffs_np<-matrix(NA, length(unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[which(is.na(or.info$index_RFpredicted))])])), 
                   length(unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[which(is.na(or.info$index_RFpredicted))])])))
colnames(DPTdiffs_np)<-unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[which(is.na(or.info$index_RFpredicted))])])
rownames(DPTdiffs_np)<-unique(or.odourant$Olfr.Nomenclature[which(or.odourant$Olfr.Nomenclature %in% rownames(or.info)[which(is.na(or.info$index_RFpredicted))])])


for(i in rownames(DPTdiffs_np)){
  for(j in colnames(DPTdiffs_np)){
    DPTdiffs_np[i, j] <- abs(or.info[i, "DPTindex"] - or.info[j, "DPTindex"])
  }
}

connectedORsDPTdiffs_np<-DPTdiffs_np[lower.tri(DPTdiffs_np)][which(connected_np[lower.tri(connected_np)]==1)]
NotConnectedORsDPTdiffs_np<-DPTdiffs_np[lower.tri(DPTdiffs_np)][which(connected_np[lower.tri(connected_np)]==0)]


pdf(file = "../figures_pdf2/F6_DPTdiffs_notPred_rev.pdf", width=7, height = 7)
boxplot(connectedORsDPTdiffs_np, NotConnectedORsDPTdiffs_np, names=c("connectedORs", "notConnectedORs"), ylab="absolute DPTindexes difference")
dev.off()
wilcox.test(connectedORsDPTdiffs_np, NotConnectedORsDPTdiffs_np)

#write.csv(or.odourant)

#### Other way around

ORsairMucusPC<-data.frame(Olfr.Nomenclature=or.odourant$Olfr.Nomenclature, airMucusPC=airMucusPC[or.odourant$Odorant.CAS.Number])
ORsairMucusPC<-tapply(ORsairMucusPC$airMucusPC, ORsairMucusPC$Olfr.Nomenclature, function(x){mean(x, na.rm=T)})
ORsairMucusPC<-ORsairMucusPC[-which(is.na(ORsairMucusPC))]

pdf(file = "../figures_pdf2/F6_airMucusPC_ORs.pdf", width=7, height = 7)

plotXYcorrelation(x=or.info[names(ORsairMucusPC),]$DPTindex, y=ORsairMucusPC, xlab = "DPT", 
                  ylab="mean airMucusCoeff", col=1, main="all data")

plotXYcorrelation2(x=or.info[names(ORsairMucusPC),]$DPTindex, y=ORsairMucusPC, xlab = "DPT", 
                   ylab="mean airMucusCoeff", col=1, main="all data")

ORs_2OdsPlus<-table(or.odourant$Olfr.Nomenclature)
ORs_2OdsPlus<-names(ORs_2OdsPlus)[which(ORs_2OdsPlus>1)]

plotXYcorrelation(x=or.info[ORs_2OdsPlus,]$DPTindex, y=ORsairMucusPC[ORs_2OdsPlus], xlab = "DPT", 
                  ylab="mean airMucusCoeff", col=1, main="ORs detecting 2 or more odourants")

plotXYcorrelation2(x=or.info[ORs_2OdsPlus,]$DPTindex, y=ORsairMucusPC[ORs_2OdsPlus], xlab = "DPT", 
                  ylab="mean airMucusCoeff", col=1, main="ORs detecting 2 or more odourants")



dev.off()

## No averages

ORsairMucusDF<-data.frame(Olfr.Nomenclature=or.odourant$Olfr.Nomenclature, airMucusPC=airMucusPC[or.odourant$Odorant.CAS.Number], dpt=or.info[or.odourant$Olfr.Nomenclature,]$DPTindex)

pdf(file = "../figures_pdf2/FigS7_airMucusPC_noMeans.pdf", width=7, height = 7)

plotXYcorrelation(x=ORsairMucusDF$dpt, y=ORsairMucusDF$airMucusPC, xlab = "DPT", 
                  ylab="airMucusCoeff", col=1, main="all data")

ORs_2OdsPlus<-table(or.odourant$Olfr.Nomenclature)
ORs_2OdsPlus<-names(ORs_2OdsPlus)[which(ORs_2OdsPlus>1)]

plotXYcorrelation(x=ORsairMucusDF$dpt[which(ORsairMucusDF$Olfr.Nomenclature %in% ORs_2OdsPlus)], y=ORsairMucusDF$airMucusPC[which(ORsairMucusDF$Olfr.Nomenclature %in% ORs_2OdsPlus)], xlab = "DPT", 
                  ylab="mean airMucusCoeff", col=1, main="ORs detecting 2 or more odourants")


dev.off()

## Behavioural data
meanNoNA<-function(x){mean(x, na.rm=T)}

BD<-read.csv("../supp_tables/behaviour.csv")
BD<-cbind(BD, airMucusData[BD$CAS,])

pdf("../figures_pdf2/FS7_behaviours_noAvgs.pdf", 7, 7)

par(mfrow=c(3, 2))
for(i in 2:19){
  plotXYcorrelation(x=BD$meanDPT, y=BD[,i], xlab = "DPT", 
                    ylab=colnames(BD)[i], col=1, main="")
}

par(mfrow=c(3, 2))
for(i in 2:19){
  plotXYcorrelation(x=BD$airMucusPCoeff, y=BD[,i], xlab = "air/mucus PC", 
                    ylab=colnames(BD)[i], col=1, main="")
}
dev.off()

pdf("../figures_pdf2/FS7_behaviours.pdf", 7, 7)

par(mfrow=c(3, 2))
for(i in 2:19){
  plotXYcorrelation(x=tapply(BD$meanDPT, BD$CAS, mean), y=tapply(BD[,i], BD$CAS, mean), xlab = "DPT", 
                    ylab=colnames(BD)[i], col=1, main="")
}

par(mfrow=c(3, 2))
for(i in 2:19){
  plotXYcorrelation(x=tapply(BD$airMucusPCoeff, BD$CAS, mean), y=tapply(BD[,i], BD$CAS, mean), xlab = "air/mucus PC", 
                    ylab=colnames(BD)[i], col=1, main="")
}

dev.off()

## Pearson

pdf("../figures_pdf2/FS7_behaviours_noAvgs_Pearson.pdf", 7, 7)

par(mfrow=c(3, 2))
for(i in 2:19){
  plotXYcorrelation2(x=BD$meanDPT, y=BD[,i], xlab = "DPT", 
                    ylab=colnames(BD)[i], col=1, main="")
}

par(mfrow=c(3, 2))
for(i in 2:19){
  plotXYcorrelation2(x=BD$airMucusPCoeff, y=BD[,i], xlab = "air/mucus PC", 
                    ylab=colnames(BD)[i], col=1, main="")
}
dev.off()

pdf("../figures_pdf2/FS7_behaviours_Pearson.pdf", 7, 7)

par(mfrow=c(3, 2))
for(i in 2:19){
  plotXYcorrelation2(x=tapply(BD$meanDPT, BD$CAS, mean), y=tapply(BD[,i], BD$CAS, mean), xlab = "DPT", 
                    ylab=colnames(BD)[i], col=1, main="")
}

par(mfrow=c(3, 2))
for(i in 2:19){
  plotXYcorrelation2(x=tapply(BD$airMucusPCoeff, BD$CAS, mean), y=tapply(BD[,i], BD$CAS, mean), xlab = "air/mucus PC", 
                    ylab=colnames(BD)[i], col=1, main="")
}

dev.off()


library(gplots)

BD<-BD[-grep('NA', rownames(BD)),]
BD_scaled<-apply(BD[,2:19], 2, function(x) (x-mean(x))/sd(x))
#PCcor<-cor(t(physChemOds2ORs), method = "spearman", use="complete.obs")
BDcor_scaled<-cor(t(BD_scaled), method = "spearman", use="complete.obs")

pdf("../figures_pdf2/FS7_BehaviourData_Correlogram_scaled.pdf", 15, 15)

#heatmap.2(as.matrix(PCcor), trace='none', key.title=NA, key.ylab=NA, col=cm.colors(20), 
#          cexRow = 0.8, key.xlab = "coefficient", cexCol = 0.8, keysize = 1, main="Odourants properties correlogram", breaks=seq(-1, 1,0.1))

#heatmap.2(as.matrix(PCcor), trace='none', key.title=NA, key.ylab=NA, col=cm.colors(10), 
#          cexRow = 0.8, key.xlab = "coefficient", cexCol = 0.8, keysize = 1, main="Odourants properties correlogram", breaks=seq(0, 1,0.1))

heatmap.2(as.matrix(BDcor_scaled), trace='none', key.title=NA, key.ylab=NA, col=cm.colors(20), 
          cexRow = 0.5, key.xlab = "coefficient", cexCol = 0.5, keysize = 1, main="Behaviour data correlogram", breaks=seq(-1, 1,0.1))

dev.off()

#airMucusPC, DPT and behaviours correlogram

behaviour_airMucus<-sapply(2:19, function(i){cor(x=tapply(BD$airMucusPCoeff, BD$CAS, mean), y=tapply(BD[,i], BD$CAS, mean), method="spearman", use = "complete.obs")})
behaviour_DPT<-sapply(2:19, function(i){cor(x=tapply(BD$meanDPT, BD$CAS, mean), y=tapply(BD[,i], BD$CAS, mean), method="spearman", use = "complete.obs")})
BD_corDPT_AM<-data.frame(behaviour_airMucus, behaviour_DPT)
rownames(BD_corDPT_AM)<-colnames(BD)[2:19]

pdf("../figures_pdf2/FS7_BehaviourDataMeans_DPT_AM_Correlogram.pdf", 7, 7)

heatmap.2(as.matrix(BD_corDPT_AM), trace='none', key.title=NA, key.ylab=NA, col=cm.colors(20), 
          cexRow = 0.5, key.xlab = "coefficient", cexCol = 0.5, keysize = 1, main="Behaviour data correlogram", breaks=seq(-1, 1,0.1))

dev.off()

behaviour_airMucus<-sapply(2:19, function(i){cor(x=BD$airMucusPCoeff, y=BD[,i], method="spearman", use = "complete.obs")})
behaviour_DPT<-sapply(2:19, function(i){cor(x=BD$meanDPT, y=BD[,i], method="spearman", use = "complete.obs")})
BD_corDPT_AM<-data.frame(behaviour_airMucus, behaviour_DPT)
rownames(BD_corDPT_AM)<-colnames(BD)[2:19]

pdf("../figures_pdf2/FS7_BehaviourData_DPT_AM_Correlogram.pdf", 7, 7)

heatmap.2(as.matrix(BD_corDPT_AM), trace='none', key.title=NA, key.ylab=NA, col=cm.colors(20), 
          cexRow = 0.5, key.xlab = "coefficient", cexCol = 0.5, keysize = 1, main="Behaviour data correlogram", breaks=seq(-1, 1,0.1))

dev.off()

#### Bulk
rpm<-function(x){x*1000000/sum(x)} #get reads per million

bulk<-read.csv("../../bulk_ORs.csv", stringsAsFactors = F, header=F, sep=";") #from Luis's paper https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004593#s5
bulk<-bulk[-which(duplicated(bulk$V2)),]
#bulk.norm=bulk[-(1:2),14:19]
bulk.norm=bulk[-(1:2),8:13]
rownames(bulk.norm)<-bulk$V2[-c(1:2)]
bulk.norm<-apply(bulk.norm, 1, function(x){as.numeric(x)})
bulk.norm<-apply(bulk.norm, 1, rpm)
bulk.means<-rowMeans(bulk.norm)

bulkMeans4Ods<-bulk.means[match(or.odourant$Olfr.Nomenclature, names(bulk.means))]
OdsBulkMeans<-tapply(bulkMeans4Ods, or.odourant$Odorant.CAS.Number, meanNoNA)
Ods_BulkMeans<-OdsBulkMeans[match(rownames(airMucusData), names(OdsBulkMeans))]

bulkMeans4ORs<-bulk.means[match(or.odourant$Olfr.Nomenclature, names(bulk.means))]
names(bulkMeans4ORs)<-or.odourant$Olfr.Nomenclature

pdf("../figures_pdf2/FS7_BulkMeanExp_airMucusPC.pdf", 7, 7)

par(mfrow=c(1, 1))
plotXYcorrelation(x=Ods_BulkMeans, y=airMucusData$airMucusPCoeff, ylab="air/mucus PC", col=1, main="", xlab="Odourants Bulk Means")
dev.off()

pdf("../figures_pdf2/FS7_BulkMeanExp_DPT.pdf", 7, 7)
par(mfrow=c(1, 1))
plotXYcorrelation(x=Ods_BulkMeans, y=airMucusData$meanDPT, ylab="Odourants mean DPT", col=1, main="", xlab="Odourants Bulk Means")

dev.off()

pdf("../figures_pdf2/FS7_BulkMeanExp_airMucusPC_noMeans.pdf", 7, 7)

par(mfrow=c(1, 1))
plotXYcorrelation(x=bulkMeans4ORs[which(names(bulkMeans4ORs) %in% ORs_2OdsPlus)], 
                  y=ORsairMucusDF$airMucusPC[which(ORsairMucusDF$Olfr.Nomenclature %in% ORs_2OdsPlus)], ylab="air/mucus PC", col=1, main="", xlab="ORs Bulk Means")

dev.off()

pdf("../figures_pdf2/FS7_BulkMeanExp_DPT_noMeans.pdf", 7, 7)

par(mfrow=c(1, 1))
plotXYcorrelation(x=bulkMeans4ORs[which(names(bulkMeans4ORs) %in% ORs_2OdsPlus)], 
                  y=ORsairMucusDF$dpt[which(ORsairMucusDF$Olfr.Nomenclature %in% ORs_2OdsPlus)], ylab="DPT index", col=1, main="", xlab="ORs Bulk Means")

dev.off()

# Physicochemical properties correlogram

physChemOds2ORs<-physico.chem[rownames(airMucusData),]
physChemOds2ORs<-physChemOds2ORs[-which(is.na(physChemOds2ORs$MW)),]
physChemOds2ORs_scaled<-apply(physChemOds2ORs, 2, function(x) (x-mean(x))/sd(x))

physChemOds2ORs_pred<-physico.chem[rownames(airMucusData_pred),]
physChemOds2ORs_pred<-physChemOds2ORs_pred[-which(is.na(physChemOds2ORs_pred$MW)),]

physChemOds2ORs_NotPred<-physico.chem[rownames(airMucusData_NotPred),]
physChemOds2ORs_NotPred<-physChemOds2ORs_NotPred[-which(is.na(physChemOds2ORs_NotPred$MW)),]

write.csv(physChemOds2ORs, "../supp_tables/physChemOds2ORs.csv")
write.csv(physChemOds2ORs_pred, "../supp_tables/physChemOds2ORs_pred.csv")
write.csv(physChemOds2ORs_NotPred, "../supp_tables/physChemOds2ORs_NotPred.csv")

PCcor<-cor(t(physChemOds2ORs), method = "spearman", use="complete.obs")
PCcor_scaled<-cor(t(physChemOds2ORs_scaled), method = "spearman", use="complete.obs")

avDPT_physChemCors<-apply(physChemOds2ORs, 2, function(x){cor(airMucusData[rownames(physChemOds2ORs),]$meanDPT, x, use = "complete.obs", method="spearman")})
avDPT_physChemCorsPval<-apply(physChemOds2ORs, 2, 
                          function(x){cor.test(airMucusData[rownames(physChemOds2ORs),]$meanDPT, x, use = "complete.obs", method="spearman")$p.value})
avDPT_physChemRobPval<-apply(physChemOds2ORs, 2, 
                              function(x){rococo.test(airMucusData[rownames(physChemOds2ORs),]$meanDPT[which((!is.na(airMucusData[rownames(physChemOds2ORs),]$meanDPT))&(!is.na(x)))], x[which((!is.na(airMucusData[rownames(physChemOds2ORs),]$meanDPT))&(!is.na(x)))])@p.value})

avDPT_physChemCors_pred<-apply(physChemOds2ORs_pred, 2, function(x){cor(airMucusData_pred[rownames(physChemOds2ORs_pred),]$meanDPT, x, use = "complete.obs", method="spearman")})
avDPT_physChemCorsPval_pred<-apply(physChemOds2ORs_pred, 2, 
                              function(x){cor.test(airMucusData_pred[rownames(physChemOds2ORs_pred),]$meanDPT, x, use = "complete.obs", method="spearman")$p.value})

avDPT_physChemCors_NotPred<-apply(physChemOds2ORs_NotPred, 2, function(x){cor(airMucusData_NotPred[rownames(physChemOds2ORs_NotPred),]$meanDPT, x, use = "complete.obs", method="spearman")})
avDPT_physChemCorsPval_NotPred<-apply(physChemOds2ORs_NotPred, 2, 
                                   function(x){cor.test(airMucusData_NotPred[rownames(physChemOds2ORs_NotPred),]$meanDPT, x, use = "complete.obs", method="spearman")$p.value})


avDPT_physChemCors_rand<-apply(physChemOds2ORs, 2, function(x){cor(airMucusData[rownames(physChemOds2ORs),]$meanDPTrand, x, use = "complete.obs", method="spearman")})
avDPT_physChemCorsPval_rand<-apply(physChemOds2ORs, 2, 
                              function(x){cor.test(airMucusData[rownames(physChemOds2ORs),]$meanDPTrand, x, use = "complete.obs", method="spearman")$p.value})

physChemCV<-apply(physChemOds2ORs, 2, function(x){var(x, na.rm=T)/mean(x, na.rm=T)})

pdf("../figures_pdf2/FS7_PCpropsCorrelogramHM_scaled.pdf", 7, 7)

#heatmap.2(as.matrix(PCcor), trace='none', key.title=NA, key.ylab=NA, col=cm.colors(20), 
#          cexRow = 0.8, key.xlab = "coefficient", cexCol = 0.8, keysize = 1, main="Odourants properties correlogram", breaks=seq(-1, 1,0.1))

#heatmap.2(as.matrix(PCcor), trace='none', key.title=NA, key.ylab=NA, col=cm.colors(10), 
#          cexRow = 0.8, key.xlab = "coefficient", cexCol = 0.8, keysize = 1, main="Odourants properties correlogram", breaks=seq(0, 1,0.1))

heatmap.2(as.matrix(PCcor_scaled), trace='none', key.title=NA, key.ylab=NA, col=cm.colors(20), 
          cexRow = 0.8, key.xlab = "coefficient", cexCol = 0.8, keysize = 1, main="Odourants properties correlogram", breaks=seq(-1, 1,0.1))

dev.off()

pvalCol=rep(1, length(avDPT_physChemCors))
pvalCol[which(avDPT_physChemCorsPval<0.05)]<-2
pvalCol<-pvalCol[order(avDPT_physChemCors)]
names(pvalCol)<-names(avDPT_physChemCors)

pdf("../figures_pdf2/FS7_PCpropsCorrelogram.pdf", 10, 7)
plot(1:length(sort(avDPT_physChemCors)), sort(avDPT_physChemCors), ylab="Rho", ylim=c(-1, 1), lty=1, 
     type="b", pch=19, xlab="properties", main="Correlation between properties and mean DPT", col=alpha(pvalCol, 0.3))
dev.off()

top10<-names(sort(avDPT_physChemCors, decreasing = T))[1:10]
bottom10<-names(sort(avDPT_physChemCors))[1:10]

pdf("../figures_pdf2/FS7_PCpropsTop10Cors_1.pdf", 7, 7)

par(mfrow=c(3, 2))
for(prop in top10){
  #plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
  #                  col=1, main="")
  plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
                    col=1, main="")
}

par(mfrow=c(3, 2))
for(prop in bottom10){
  plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
                    col=1, main="")
}
dev.off()

## Continuous descriptors

#is.continuous<-apply(physChemOds2ORs, 2, function(x){length(unique(x))>70})
is.continuous<-apply(physChemOds2ORs, 2, function(x){max(table(x))<20})
continuous<-names(is.continuous)[which(is.continuous==TRUE)]

top10<-names(sort(avDPT_physChemCors[continuous], decreasing = T))[1:10]
bottom10<-names(sort(avDPT_physChemCors[continuous]))[1:10]

pvalCol2=rep(1, length(continuous))
pvalCol2[which(avDPT_physChemCorsPval[continuous]<0.05)]<-2
pvalCol2<-pvalCol2[order(avDPT_physChemCors[continuous])]
#names(pvalCol)<-names(avDPT_physChemCors)

pdf("../figures_pdf2/FS7_PCpropsCorrelogram_continuous.pdf", 10, 7)
plot(1:length(continuous), sort(avDPT_physChemCors[continuous]), ylab="Rho", ylim=c(-1, 1), lty=1, 
     type="b", pch=19, xlab="properties", main="Correlation between properties and mean DPT", col=alpha(pvalCol2, 0.3))
dev.off()

top10<-names(sort(avDPT_physChemCors_pred[continuous], decreasing = T))[1:10]
bottom10<-names(sort(avDPT_physChemCors_pred[continuous]))[1:10]

pdf("../figures_pdf2/FS7_PCpropsTop10Cors_continuous_v2_pred_rev.pdf", 7, 7)

par(mfrow=c(3, 2))
for(prop in top10){
  #plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
  #                  col=1, main="")
  plotXYcorrelation(x=airMucusData_pred[rownames(physChemOds2ORs_pred),]$meanDPT, y=physChemOds2ORs_pred[,prop], xlab="mean DPT", ylab=prop, 
                    col=topicCol[rownames(physChemOds2ORs_pred)], main="")
}

par(mfrow=c(3, 2))
for(prop in bottom10){
  #plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
  #                  col=1, main="")
  plotXYcorrelation(x=airMucusData_pred[rownames(physChemOds2ORs_pred),]$meanDPT, y=physChemOds2ORs_pred[,prop], xlab="mean DPT", ylab=prop, 
                    col=topicCol[rownames(physChemOds2ORs_pred)], main="")
}


par(mfrow=c(3, 2))
for(prop in top10){
  #plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
  #                  col=1, main="")
  plotXYcorrelation(x=airMucusData_pred[rownames(physChemOds2ORs_pred),]$airMucusPCoeff, y=physChemOds2ORs_pred[,prop], xlab="air/mucus PC", ylab=prop, 
                    col=1, main="")
}

par(mfrow=c(3, 2))
for(prop in bottom10){
  plotXYcorrelation(x=airMucusData_pred[rownames(physChemOds2ORs_pred),]$airMucusPCoeff, y=physChemOds2ORs_pred[,prop], xlab="air/mucus PC", ylab=prop, 
                    col=1, main="")
}
dev.off()

topSigNotPred<-names(avDPT_physChemCorsPval_NotPred[continuous])[which(avDPT_physChemCorsPval_NotPred[continuous]<0.01 
                                                                       & abs(avDPT_physChemCors_NotPred[continuous])>0.5)]
top10<-names(sort(avDPT_physChemCors_NotPred[continuous], decreasing = T))[1:10]
bottom10<-names(sort(avDPT_physChemCors_NotPred[continuous]))[1:10]

pdf("../figures_pdf2/FS7_PCpropsTop10Cors_continuous_v2_NotPred_rev.pdf", 7, 7)

par(mfrow=c(3, 2))
for(prop in top10){
  #plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
  #                  col=1, main="")
  plotXYcorrelation(x=airMucusData_NotPred[rownames(physChemOds2ORs_NotPred),]$meanDPT, y=physChemOds2ORs_NotPred[,prop], xlab="mean DPT", ylab=prop, 
                    col=topicCol[rownames(physChemOds2ORs_NotPred)], main="")
}

par(mfrow=c(3, 2))
for(prop in bottom10){
  #plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
  #                  col=1, main="")
  plotXYcorrelation(x=airMucusData_NotPred[rownames(physChemOds2ORs_NotPred),]$meanDPT, y=physChemOds2ORs_NotPred[,prop], xlab="mean DPT", ylab=prop, 
                    col=topicCol[rownames(physChemOds2ORs_NotPred)], main="")
}


par(mfrow=c(3, 2))
for(prop in top10){
  #plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
  #                  col=1, main="")
  plotXYcorrelation(x=airMucusData_NotPred[rownames(physChemOds2ORs_NotPred),]$airMucusPCoeff, y=physChemOds2ORs_NotPred[,prop], xlab="air/mucus PC", ylab=prop, 
                    col=1, main="")
}

par(mfrow=c(3, 2))
for(prop in bottom10){
  plotXYcorrelation(x=airMucusData_NotPred[rownames(physChemOds2ORs_NotPred),]$airMucusPCoeff, y=physChemOds2ORs_NotPred[,prop], xlab="air/mucus PC", ylab=prop, 
                    col=1, main="")
}
dev.off()

top10<-names(sort(avDPT_physChemCors[continuous], decreasing = T))[1:10]
bottom10<-names(sort(avDPT_physChemCors[continuous]))[1:10]

pdf("../figures_pdf2/FS7_PCpropsTop10Cors_continuous_v2.pdf", 7, 7)

par(mfrow=c(3, 2))
for(prop in top10){
  #plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
  #                  col=1, main="")
  plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
                    col=topicCol[rownames(physChemOds2ORs)], main="")
}

par(mfrow=c(3, 2))
for(prop in bottom10){
  plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
                    col=topicCol[rownames(physChemOds2ORs)], main="")
}

par(mfrow=c(3, 2))
for(prop in top10){
  #plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
  #                  col=1, main="")
  plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$airMucusPCoeff, y=physChemOds2ORs[,prop], xlab="air/mucus PC", ylab=prop, 
                    col=1, main="")
}

par(mfrow=c(3, 2))
for(prop in bottom10){
  plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$airMucusPCoeff, y=physChemOds2ORs[,prop], xlab="air/mucus PC", ylab=prop, 
                    col=1, main="")
}
dev.off()

pdf("../figures_pdf2/FS7_PCpropsTopSigCorsNotPred_continuous_v2_NotPred_rev.pdf", 7, 7)

par(mfrow=c(3, 2))
for(prop in topSigNotPred){
  #plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
  #                  col=1, main="")
  plotXYcorrelation(x=airMucusData_NotPred[rownames(physChemOds2ORs_NotPred),]$meanDPT, y=physChemOds2ORs_NotPred[,prop], xlab="mean DPT", ylab=prop, 
                    col=topicCol[rownames(physChemOds2ORs_NotPred)], main="")
  
  plotXYcorrelation(x=airMucusData_NotPred[rownames(physChemOds2ORs_NotPred),]$airMucusPCoeff, y=physChemOds2ORs_NotPred[,prop], xlab="air/mucus PC", ylab=prop, 
                    col=1, main="")
}

dev.off()

pdf("../figures_pdf2/FS7_PCpropsTopSigCorsNotPred_continuous_v2_pred_rev.pdf", 7, 7)

par(mfrow=c(3, 2))
for(prop in topSigNotPred){
  #plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
  #                  col=1, main="")
  plotXYcorrelation(x=airMucusData_pred[rownames(physChemOds2ORs_pred),]$meanDPT, y=physChemOds2ORs_pred[,prop], xlab="mean DPT", ylab=prop, 
                    col=topicCol[rownames(physChemOds2ORs_pred)], main="")
  
  plotXYcorrelation(x=airMucusData_pred[rownames(physChemOds2ORs_pred),]$airMucusPCoeff, y=physChemOds2ORs_pred[,prop], xlab="air/mucus PC", ylab=prop, 
                    col=1, main="")
}

dev.off()

pdf("../figures_pdf2/FS7_PCpropsTopSigCorsNotPred_continuous_v2_rev.pdf", 7, 7)

par(mfrow=c(3, 2))
for(prop in topSigNotPred){
  #plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
  #                  col=1, main="")
  plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
                    col=topicCol[rownames(physChemOds2ORs)], main="")
  
  plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$airMucusPCoeff, y=physChemOds2ORs[,prop], xlab="air/mucus PC", ylab=prop, 
                    col=1, main="")
}

dev.off()

pdf("../figures_pdf2/FS7_PCpropsTop10Cors_continuous_rand.pdf", 7, 7)

par(mfrow=c(3, 2))
for(prop in top10){
  #plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
  #                  col=1, main="")
  plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPTrand, y=physChemOds2ORs[,prop], xlab="mean DPT (randomized)", ylab=prop, 
                    col=1, main="")
}

par(mfrow=c(3, 2))
for(prop in bottom10){
  plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPTrand, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
                    col=1, main="")
}

plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPTrand, y=airMucusData[rownames(physChemOds2ORs),]$airMucusPCoeff, xlab = "avgDPT (randomized)", 
                  ylab="airMucusCoeff", col=1, main="")
dev.off()


### EXCLUDE CLASS 1

avDPT_physChemCorsNoC1<-apply(physChemOds2ORs[names(odsDPTnoC1),], 2, function(x){cor(airMucusData[names(odsDPTnoC1),]$meanDPT, x, use = "complete.obs", method="spearman")})
avDPT_physChemCorsPvalNoC1<-apply(physChemOds2ORs[names(odsDPTnoC1),], 2, 
                              function(x){cor.test(airMucusData[names(odsDPTnoC1),]$meanDPT, x, use = "complete.obs", method="spearman")$p.value})

#avDPT_physChemCorsNoC1_pred<-apply(physChemOds2ORs_pred[which(rownames(physChemOds2ORs_pred) %in% names(odsDPTnoC1)),], 2, function(x){cor(airMucusData_pred[which(rownames(airMucusData_pred) %in% names(odsDPTnoC1)),]$meanDPT, x, use = "complete.obs", method="spearman")})

avDPT_physChemCorsNoC1_pred<-apply(physChemOds2ORs_pred[names(odsDPTnoC1),], 2, function(x){cor(airMucusData_pred[names(odsDPTnoC1),]$meanDPT, x, use = "complete.obs", method="spearman")})
avDPT_physChemCorsPvalNoC1_pred<-apply(physChemOds2ORs_pred[names(odsDPTnoC1),], 2, 
                                  function(x){cor.test(airMucusData_pred[names(odsDPTnoC1),]$meanDPT, x, use = "complete.obs", method="spearman")$p.value})

avDPT_physChemCorsNoC1_NotPred<-apply(physChemOds2ORs_NotPred[names(odsDPTnoC1),], 2, function(x){cor(airMucusData_NotPred[names(odsDPTnoC1),]$meanDPT, x, use = "complete.obs", method="spearman")})
avDPT_physChemCorsPvalNoC1_NotPred<-apply(physChemOds2ORs_NotPred[names(odsDPTnoC1),], 2, 
                                       function(x){cor.test(airMucusData_NotPred[names(odsDPTnoC1),]$meanDPT, x, use = "complete.obs", method="spearman")$p.value})


top10<-names(sort(avDPT_physChemCorsNoC1, decreasing = T))[1:10]
bottom10<-names(sort(avDPT_physChemCorsNoC1))[1:10]

pvalCol3=rep(1, length(names(odsDPTnoC1)))
pvalCol3[which(avDPT_physChemCorsPvalNoC1<0.05)]<-2
pvalCol3<-pvalCol3[order(avDPT_physChemCorsNoC1)]

pdf("../figures_pdf2/FS7_PCpropsCorrelogram_NoC1.pdf", 10, 7)
plot(1:length(sort(avDPT_physChemCorsNoC1)), sort(avDPT_physChemCorsNoC1), ylab="Rho", ylim=c(-1, 1), lty=1, 
     type="b", pch=19, xlab="properties", main="Correlation between properties and mean DPT", col=alpha(pvalCol3, 0.3))
dev.off()


pdf("../figures_pdf2/FS7_PCpropsTop10Cors_NoC1.pdf", 7, 7)

par(mfrow=c(3, 2))
for(prop in top10){
  #plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
  #                  col=1, main="")
  plotXYcorrelation(x=airMucusData[names(odsDPTnoC1),]$meanDPT, y=physChemOds2ORs[names(odsDPTnoC1),prop], xlab="mean DPT", ylab=prop, 
                    col=1, main="")
}

par(mfrow=c(3, 2))
for(prop in bottom10){
  plotXYcorrelation(x=airMucusData[names(odsDPTnoC1),]$meanDPT, y=physChemOds2ORs[names(odsDPTnoC1),prop], xlab="mean DPT", ylab=prop, 
                    col=1, main="")
}
dev.off()

top10<-names(sort(avDPT_physChemCorsNoC1[continuous], decreasing = T))[1:10]
bottom10<-names(sort(avDPT_physChemCorsNoC1[continuous]))[1:10]

pvalCol4=rep(1, length(continuous))
pvalCol4[which(avDPT_physChemCorsPvalNoC1[continuous]<0.05)]<-2
pvalCol4<-pvalCol4[order(avDPT_physChemCorsNoC1[continuous])]

pdf("../figures_pdf2/FS7_PCpropsCorrelogram_NoC1continuous.pdf", 10, 7)
plot(1:length(continuous), sort(avDPT_physChemCorsNoC1[continuous]), ylab="Rho", ylim=c(-1, 1), lty=1, 
     type="b", pch=19, xlab="properties", main="Correlation between properties and mean DPT", col=alpha(pvalCol4, 0.3))
dev.off()

pdf("../figures_pdf2/FS7_PCpropsTop10Cors_NoC1continuous.pdf", 7, 7)

par(mfrow=c(3, 2))
for(prop in top10){
  #plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
  #                  col=1, main="")
  plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
                    col=1, main="")
}

par(mfrow=c(3, 2))
for(prop in bottom10){
  plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
                    col=1, main="")
}
dev.off()

avDPT_physChemRobPval<-c(avDPT_physChemRobPval, rococo.test(airMucusData$meanDPT[which((!is.na(airMucusData$meanDPT)) & (!is.na(airMucusData$airMucusPCoeff)))], airMucusData$airMucusPCoeff[which((!is.na(airMucusData$meanDPT)) & (!is.na(airMucusData$airMucusPCoeff)))])@p.value)
#avDPT_physChemCorsPval<-c(avDPT_physChemCorsPval, cor.test(airMucusData[rownames(physChemOds2ORs),]$meanDPT, airMucusData[rownames(physChemOds2ORs),]$airMucusPCoeff, use = "complete.obs", method="spearman")$p.value)
#avDPT_physChemCors<-c(avDPT_physChemCors, cor(airMucusData[rownames(physChemOds2ORs),]$meanDPT, airMucusData[rownames(physChemOds2ORs),]$airMucusPCoeff, use = "complete.obs", method="spearman"))
avDPT_physChemCorsPval<-c(avDPT_physChemCorsPval, cor.test(airMucusData$meanDPT, airMucusData$airMucusPCoeff, use = "complete.obs", method="spearman")$p.value)
avDPT_physChemCors<-c(avDPT_physChemCors, cor(airMucusData$meanDPT, airMucusData$airMucusPCoeff, use = "complete.obs", method="spearman"))

avDPT_physChemCorsPval_pred<-c(avDPT_physChemCorsPval_pred, cor.test(airMucusData_pred$meanDPT, airMucusData_pred$airMucusPCoeff, use = "complete.obs", method="spearman")$p.value)
avDPT_physChemCors_pred<-c(avDPT_physChemCors_pred, cor(airMucusData_pred$meanDPT, airMucusData_pred$airMucusPCoeff, use = "complete.obs", method="spearman"))

avDPT_physChemCorsPval_NotPred<-c(avDPT_physChemCorsPval_NotPred, cor.test(airMucusData_NotPred$meanDPT, airMucusData_NotPred$airMucusPCoeff, use = "complete.obs", method="spearman")$p.value)
avDPT_physChemCors_NotPred<-c(avDPT_physChemCors_NotPred, cor(airMucusData_NotPred$meanDPT, airMucusData_NotPred$airMucusPCoeff, use = "complete.obs", method="spearman"))

df_p<-data.frame(avDPT_physChemCors_pred, avDPT_physChemCorsPval_pred, avDPT_physChemCorsPvalAdj=p.adjust(avDPT_physChemCorsPval_pred, method="fdr"))
#write.csv(df_p, "./avDPT_physChemCors_pred.csv")

df_np<-data.frame(avDPT_physChemCors_NotPred, avDPT_physChemCorsPval_NotPred, avDPT_physChemCorsPvalAdj=p.adjust(avDPT_physChemCorsPval_NotPred, method="fdr"))
#write.csv(df_np, "./avDPT_physChemCors_NotPred.csv")

df<-data.frame(avDPT_physChemCors, avDPT_physChemCorsPval, avDPT_physChemCorsPvalAdj=p.adjust(avDPT_physChemCorsPval, method="fdr"))
#write.csv(df, "./avDPT_physChemCors.csv")
avDPT_physChemCorsPval_rand<-c(avDPT_physChemCorsPval_rand, cor.test(airMucusData$meanDPTrand, airMucusData$airMucusPCoeff, use = "complete.obs", method="spearman")$p.value)
avDPT_physChemCors_rand<-c(avDPT_physChemCors_rand, cor(airMucusData$meanDPTrand, airMucusData$airMucusPCoeff, use = "complete.obs", method="spearman"))


physChemCV<-c(physChemCV, var(airMucusData$airMucusPCoeff, na.rm=T)/abs(mean(airMucusData$airMucusPCoeff, na.rm=T)))

shapes<-rep(15, length(avDPT_physChemCors))
shapes[which(names(avDPT_physChemCors) %in% continuous)]<-16
shapes[length(shapes)]<-16

pvalCol5<-rep(1, length(avDPT_physChemCors))
pvalCol5[which(p.adjust(avDPT_physChemCorsPval, method = "fdr")<0.05)]<-4
pvalCol5[which(names(avDPT_physChemCors) %in% names(which(p.adjust(avDPT_physChemCorsPvalNoC1, method = "fdr")<0.05)))]<-2
pvalCol5[length(pvalCol5)]<-7

pvalCol5_pred<-rep(1, length(avDPT_physChemCors_pred))
pvalCol5_pred[which(p.adjust(avDPT_physChemCorsPval_pred, method = "fdr")<0.05)]<-4
pvalCol5_pred[which(names(avDPT_physChemCors_pred) %in% names(which(p.adjust(avDPT_physChemCorsPvalNoC1_pred, method = "fdr")<0.05)))]<-2
pvalCol5_pred[length(pvalCol5_pred)]<-7

pvalCol5_NotPred<-rep(1, length(avDPT_physChemCors_NotPred))
pvalCol5_NotPred[which(p.adjust(avDPT_physChemCorsPval_NotPred, method = "fdr")<0.05)]<-4
pvalCol5_NotPred[which(names(avDPT_physChemCors_NotPred) %in% names(which(p.adjust(avDPT_physChemCorsPvalNoC1_NotPred, method = "fdr")<0.05)))]<-2
pvalCol5_NotPred[length(pvalCol5_NotPred)]<-7

pdf("../figures_pdf2/FS7_PCprops_coeffs.pdf", 7, 7)
plot(avDPT_physChemCors, -log10(p.adjust(avDPT_physChemCorsPval, method="fdr")), pch=shapes, col=pvalCol5, xlab="correlation coefficient", ylab="-log10 adj p. val", xlim=c(-1, 1))
legend("bottomright", c("not significant", "significant taking all ORs", "significant excluding C1 ORs", "air/mucus PC", "< 20 repeated values", "> 20 repeated values"), col=c(1, 4, 2, 7, 1, 1), pch=c(16, 16, 16, 16, 1, 0), cex=0.5)

shapes2<-rep(NA, length(avDPT_physChemCors))
shapes2[which(names(avDPT_physChemCors) %in% continuous)]<-16
shapes2[length(shapes)]<-16

shapes2_pred<-rep(NA, length(avDPT_physChemCors_pred))
shapes2_pred[which(names(avDPT_physChemCors_pred) %in% continuous)]<-16
shapes2_pred[length(shapes2_pred)]<-16

shapes2_NotPred<-rep(NA, length(avDPT_physChemCors_NotPred))
shapes2_NotPred[which(names(avDPT_physChemCors_NotPred) %in% continuous)]<-16
shapes2_NotPred[length(shapes2_NotPred)]<-16

labs<-rep(NA, length(avDPT_physChemCors))
labs[which(((avDPT_physChemCors<(-0.5))) | (avDPT_physChemCors>0.5))]<-names(avDPT_physChemCors)[which(((avDPT_physChemCors<(-0.5))) | (avDPT_physChemCors>0.5))]
labs[length(labs)]<-"air/mucus PC"
labs[which(is.na(shapes2))]<-NA

labs_pred<-rep(NA, length(avDPT_physChemCors_pred))
labs_pred[which(((avDPT_physChemCors_pred<(-0.5))) | (avDPT_physChemCors_pred>0.5))]<-names(avDPT_physChemCors_pred)[which(((avDPT_physChemCors_pred<(-0.5))) | (avDPT_physChemCors_pred>0.5))]
labs_pred[length(labs_pred)]<-"air/mucus PC"
labs_pred[which(is.na(shapes2_pred))]<-NA

labs_NotPred<-rep(NA, length(avDPT_physChemCors_NotPred))
labs_NotPred[which(((avDPT_physChemCors_NotPred<(-0.5))) | (avDPT_physChemCors_NotPred>0.5))]<-names(avDPT_physChemCors_NotPred)[which(((avDPT_physChemCors_NotPred<(-0.5))) | (avDPT_physChemCors_NotPred>0.5))]
labs_NotPred[length(labs_NotPred)]<-"air/mucus PC"
labs_NotPred[which(is.na(shapes2_NotPred))]<-NA


plot(avDPT_physChemCors, -log10(p.adjust(avDPT_physChemCorsPval, method="fdr")), pch=shapes2, col=pvalCol5, xlab="correlation coefficient", ylab="-log10 adj p. val", xlim=c(-1, 1))
abline(v = min(avDPT_physChemCors_rand[which(!is.na(avDPT_physChemCors_rand))]))
abline(v = max(avDPT_physChemCors_rand[which(!is.na(avDPT_physChemCors_rand))]))

abline(v = 0.5, col="red")
abline(v = -0.5, col="red")

text(-log10(p.adjust(avDPT_physChemCorsPval, method="fdr"))~avDPT_physChemCors, labels=labs,cex=0.4, font=2)
legend("bottomright", c("not significant", "significant taking all ORs", "significant excluding C1 ORs", "air/mucus PC"), col=c(1, 4, 2, 7), pch=16, cex=0.5)

plot(avDPT_physChemCors_rand, -log10(p.adjust(avDPT_physChemCorsPval_rand, method="fdr")), pch=shapes2, col=pvalCol5, xlab="correlation coefficient", ylab="-log10 adj p. val", xlim=c(-1, 1), main="Randomization")
legend("bottomright", c("not significant", "significant taking all ORs", "significant excluding C1 ORs", "air/mucus PC"), col=c(1, 4, 2, 7), pch=16, cex=0.5)

shapes3<-rep(NA, length(avDPT_physChemCors))
shapes3[which((names(avDPT_physChemCors) %in% names(which(p.adjust(avDPT_physChemRobPval, method="fdr")<0.05))) & (names(avDPT_physChemCors) %in% continuous ))]<-16
shapes3[length(shapes)]<-16
plot(avDPT_physChemCors, -log10(p.adjust(avDPT_physChemCorsPval, method="fdr")), pch=shapes3, col=pvalCol5, xlab="correlation coefficient", ylab="-log10 adj p. val", xlim=c(-1, 1))
abline(v = min(avDPT_physChemCors_rand[which(!is.na(avDPT_physChemCors_rand))]))
abline(v = max(avDPT_physChemCors_rand[which(!is.na(avDPT_physChemCors_rand))]))
legend("bottomright", c("not significant", "significant taking all ORs", "significant excluding C1 ORs", "air/mucus PC"), col=c(1, 4, 2, 7), pch=16, cex=0.5)

plot(avDPT_physChemCors_rand, -log10(p.adjust(avDPT_physChemCorsPval_rand, method="fdr")), pch=shapes2, col=pvalCol5, xlab="correlation coefficient", ylab="-log10 adj p. val", xlim=c(-1, 1), main="Randomization")
legend("bottomright", c("not significant", "significant taking all ORs", "significant excluding C1 ORs", "air/mucus PC"), col=c(1, 4, 2, 7), pch=16, cex=0.5)


plot(log10(physChemCV), avDPT_physChemCors, pch=shapes2, col=pvalCol5, xlab="coefficient of variation", ylab="correlation coefficient")
abline(h = min(avDPT_physChemCors_rand[which(!is.na(avDPT_physChemCors_rand))]))
abline(h = max(avDPT_physChemCors_rand[which(!is.na(avDPT_physChemCors_rand))]))
legend("topright", c("not significant", "significant taking all ORs", "significant excluding C1 ORs", "air/mucus PC"), col=c(1, 4, 2, 7), pch=16, cex=0.5)
#points(log10(physChemCV)[length(physChemCV)], avDPT_physChemCors[length(avDPT_physChemCors)], col=7)

plot(log10(physChemCV), avDPT_physChemCors_rand, pch=shapes2, col=pvalCol5, xlab="coefficient of variation", ylab="correlation coefficient", main="Randomization")
legend("topright", c("not significant", "significant taking all ORs", "significant excluding C1 ORs", "air/mucus PC"), col=c(1, 4, 2, 7), pch=16, cex=0.5)

dev.off()


pdf("../figures_pdf2/F5x_PCprops_coeffs.pdf", 7, 7)

plot(avDPT_physChemCors, -log10(p.adjust(avDPT_physChemCorsPval, method="fdr")), pch=shapes2, col=pvalCol5, xlab="correlation coefficient", ylab="-log10 adj p. val", xlim=c(-1, 1))
abline(v = min(avDPT_physChemCors_rand[which(!is.na(avDPT_physChemCors_rand))]))
abline(v = max(avDPT_physChemCors_rand[which(!is.na(avDPT_physChemCors_rand))]))

abline(v = 0.5, col="red")
abline(v = -0.5, col="red")

text(-log10(p.adjust(avDPT_physChemCorsPval, method="fdr"))~avDPT_physChemCors, labels=labs,cex=0.4, font=2)
legend("bottomright", c("not significant", "significant taking all ORs", "significant excluding C1 ORs", "air/mucus PC"), col=c(1, 4, 2, 7), pch=16, cex=0.5)

dev.off()

pdf("../figures_pdf2/F6supp_PCprops_coeffs_rev.pdf", 7, 7)

## Predicted ORs
plot(avDPT_physChemCors_pred, -log10(p.adjust(avDPT_physChemCorsPval_pred, method="fdr")), pch=shapes2_pred, col=pvalCol5_pred, xlab="correlation coefficient", ylab="-log10 adj p. val", xlim=c(-1, 1))
abline(v = 0.5, col="red")
abline(v = -0.5, col="red")

text(-log10(p.adjust(avDPT_physChemCorsPval_pred, method="fdr"))~avDPT_physChemCors_pred, labels=labs_pred,cex=0.4, font=2)
legend("bottomright", c("not significant", "significant taking all ORs", "significant excluding C1 ORs", "air/mucus PC"), col=c(1, 4, 2, 7), pch=16, cex=0.5)

## Not predicted ORs
plot(avDPT_physChemCors_NotPred, -log10(p.adjust(avDPT_physChemCorsPval_NotPred, method="fdr")), pch=shapes2_NotPred, col=pvalCol5_NotPred, xlab="correlation coefficient", ylab="-log10 adj p. val", xlim=c(-1, 1))
abline(v = 0.5, col="red")
abline(v = -0.5, col="red")

text(-log10(p.adjust(avDPT_physChemCorsPval_NotPred, method="fdr"))~avDPT_physChemCors_NotPred, labels=labs_NotPred,cex=0.4, font=2)
legend("bottomright", c("not significant", "significant taking all ORs", "significant excluding C1 ORs", "air/mucus PC"), col=c(1, 4, 2, 7), pch=16, cex=0.5)

dev.off()

df$continuous<-(!is.na(shapes2))
rownames(df)[length(rownames(df))]<-"air/mucusPC"
write.csv(df, "../supp_tables/avDPT_physChemCors.csv")

df_p$continuous<-(!is.na(shapes2_pred))
rownames(df_p)[length(rownames(df_p))]<-"air/mucusPC"
write.csv(df_p, "../supp_tables/avDPT_physChemCors_pred.csv")

df_np$continuous<-(!is.na(shapes2_NotPred))
rownames(df_np)[length(rownames(df_np))]<-"air/mucusPC"
write.csv(df_np, "../supp_tables/avDPT_physChemCors_NotPred.csv")
### NEW
physChemOds2ORs<-physico.chem[unique(or.odourant$Odorant.CAS.Number),]
meanDPT4phys<-tapply(or.info[or.odourant$Olfr.Nomenclature,]$DPTindex, or.odourant$Odorant.CAS.Number, function(x){mean(x, na.rm=T)})
meanDPT4phys<-meanDPT4phys[-which(is.na(physChemOds2ORs$MW))]
physChemOds2ORs<-physChemOds2ORs[-which(is.na(physChemOds2ORs$MW)),]

avDPT_physChemCors<-apply(physChemOds2ORs, 2, function(x){cor(meanDPT4phys, x, use = "complete.obs", method="spearman")})
avDPT_physChemCorsPval<-apply(physChemOds2ORs, 2, 
                              function(x){cor.test(meanDPT4phys, x, use = "complete.obs", method="spearman")$p.value})
top5<-names(sort(avDPT_physChemCors, decreasing = T))[1:5]
bottom5<-names(sort(avDPT_physChemCors))[1:5]

pvalCol=rep(1, length(avDPT_physChemCors))
pvalCol[which(avDPT_physChemCorsPval<0.05)]<-2
pvalCol<-pvalCol[order(avDPT_physChemCors)]

pdf("../figures_pdf2/FS7_PCpropsCorrelogram_allOds.pdf", 10, 7)
plot(1:length(sort(avDPT_physChemCors)), sort(avDPT_physChemCors), ylab="Rho", ylim=c(-1, 1), lty=1, 
     type="b", pch=19, xlab="properties", main="Correlation between properties and mean DPT", col=alpha(pvalCol, 0.3))
dev.off()


#physChemOds2ORs<-physChemOds2ORs[-which(is.na(physChemOds2ORs$MW)),]


pdf("../figures_pdf2/FS7_PCpropsTop10Cors.pdf", 7, 7)

par(mfrow=c(3, 2))
for(prop in top5){
  #plotXYcorrelation(x=airMucusData[rownames(physChemOds2ORs),]$meanDPT, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
  #                  col=1, main="")
  plotXYcorrelation(x=meanDPT4phys, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
                                     col=1, main="")
}

par(mfrow=c(3, 2))
for(prop in bottom5){
  plotXYcorrelation(x=meanDPT4phys, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
                    col=1, main="")
}
dev.off()

pdf("../figures_pdf2/FS7_solubilityRelatedPCpropsCors.pdf", 7, 7)
par(mfrow=c(3, 2))
for(prop in c("Hy", "MW", "AMW", "MLOGP", "MLOGP2", "ALOGP", "ALOGP2")){
  plotXYcorrelation(x=meanDPT4phys, y=physChemOds2ORs[,prop], xlab="mean DPT", ylab=prop, 
                    col=1, main="")
}
dev.off()

tmp<-data.frame(cor=avDPT_physChemCors, pval=avDPT_physChemCorsPval)
tmp<-tmp[order(tmp$pval),]
write.csv(tmp, "../supp_tables/PCprops_meanDPT_cor.csv")
#### Air mucus PC of odourants binding the same OR vs others

connected<-matrix(0, length(unique(or.odourant$Odorant.CAS.Number)), 
                  length(unique(or.odourant$Odorant.CAS.Number)))
colnames(connected)<-unique(or.odourant$Odorant.CAS.Number)
rownames(connected)<-unique(or.odourant$Odorant.CAS.Number)

for(i in rownames(connected)){
  for(j in colnames(connected)){
    if(length(intersect(or.odourant$Olfr.Nomenclature[which(or.odourant$Odorant.CAS.Number==i)], 
                        or.odourant$Olfr.Nomenclature[which(or.odourant$Odorant.CAS.Number==j)]))>0){
      connected[i, j] <- 1
    }
  }
}


airMucusPCdiffs<-matrix(NA, length(unique(or.odourant$Odorant.CAS.Number)), 
                 length(unique(or.odourant$Odorant.CAS.Number)))
colnames(airMucusPCdiffs)<-unique(or.odourant$Odorant.CAS.Number)
rownames(airMucusPCdiffs)<-unique(or.odourant$Odorant.CAS.Number)

for(i in rownames(airMucusPCdiffs)){
  for(j in colnames(airMucusPCdiffs)){
    airMucusPCdiffs[i, j] <- abs(airMucusPC[i] - airMucusPC[j])
  }
}

connectedOdsairMucusPCdiffs<-airMucusPCdiffs[lower.tri(airMucusPCdiffs)][which(connected[lower.tri(connected)]==1)]
NotConnectedOdsairMucusPCdiffs<-airMucusPCdiffs[lower.tri(airMucusPCdiffs)][which(connected[lower.tri(connected)]==0)]

pdf(file = "../figures_pdf2/F6_airMucusPCdiffs.pdf", width=7, height = 7)
boxplot(connectedOdsairMucusPCdiffs, NotConnectedOdsairMucusPCdiffs, names=c("connectedOds", "notConnectedOds"), ylab="absolute airMucusPC difference")
dev.off()
wilcox.test(connectedOdsairMucusPCdiffs, NotConnectedOdsairMucusPCdiffs)


######OR TOPICS VS CLUSTER OF ODOURANTS (ONLY DETECTED GENES) #######

#take all the OR in a given topic, and see which odourant they can detect 
temp<-or.odourant

#keep only the ORs and odourants for which we have a topic estimation and CAS number
topic.est<-row.names(or.info)[which(!is.na(or.info$maxTopic) & or.info$Predicted=="FALSE")]
temp<-temp[temp$Olfr.Nomenclature%in%topic.est,]
temp<-temp[temp$Odorant.CAS.Number%in%names(clust.colour),]

or.used<-unique(temp$Olfr.Nomenclature)
od.used<-unique(temp$Odorant.CAS.Number)

length(or.used)
table(or.info[or.used,"maxTopic"])
table(or.info[or.used,"maxTopic"], or.info[or.used, "Class"])

length(od.used)


#out of the OR used, how many were previously studied?
length(or.used)#OR used
length(which(!is.na(or.info[or.used, "indexMiyamichiReal"])))#by Myamichi
length(which(!is.na(or.info[or.used, "ZolfrReal"])))#by Mombaerts
or.info[or.used, "indexMiyamichiChemS"]
plot(or.info[or.used, "maxTopic"],or.info[or.used, "indexMiyamichiChemS"])#agreement between chem senses indexes and ours

#replace OR with topic and odourant with cluster
temp$Olfr.Nomenclature<-or.info[temp$Olfr.Nomenclature, "maxTopic"]
temp$Odorant.CAS.Number<-clust.colour[temp$Odorant.CAS.Number]



topic.vs.clusters<-table(temp[,2], temp[,3])
chisq.test(topic.vs.clusters)
#norm.topic.vs.clusters<-t(apply(topic.vs.clusters, 1, function(x) round(100*(x/sum(x)), 1)))


norm.topic.vs.clusters<-topic.vs.clusters
for(i in 1:nrow(norm.topic.vs.clusters)){
  for(j in 1:ncol(norm.topic.vs.clusters)){
    #ni<-table(or.info[or.used,"maxTopic"])[i]#number of ORs in topic 1
    nj<-table(clust.colour[od.used])[colnames(topic.vs.clusters)[j]]#number of odourants in cluster j  
    #nj<-table(clust.colour[od.used])[rownames(topic.vs.clusters)[j]]#number of odourants in cluster j  
    
    norm.topic.vs.clusters[i,j]= round(norm.topic.vs.clusters[i,j]/(nj),2)
    
  }
}

norm.topic.vs.clusters<-round(t(apply(norm.topic.vs.clusters, 1, function(x) x/max(x))),2)
#normalized by dividing by tot number of odourant in a cluster and by max in each row

norm.topic.vs.clusters
row.names(norm.topic.vs.clusters)<-c("Topic 1(42)", "Topic 2(20)",
                                     "Topic 3(29)", "Topic 4(6)",
                                     "Topic 5 (1)")
colnames(norm.topic.vs.clusters)<-paste0(colnames(norm.topic.vs.clusters),
                                         "(",
                                         table(clust.colour[od.used])[colnames(norm.topic.vs.clusters)],
                                         ")")
#norm.topic.vs.clustersHM<-norm.topic.vs.clusters[c("yellow", "blue", "turquoise", "brown", "green"),]
#norm.topic.vs.clustersHM<-norm.topic.vs.clusters[,c("blue", "brown", "turquoise", "yellow")]
#norm.topic.vs.clustersHM<-norm.topic.vs.clusters[,c("yellow", "turquoise", "brown", "blue")]

#norm.topic.vs.clustersHM<-norm.topic.vs.clusters[,c("pink", "black", "blue", "red", "yellow", "turquoise", "green", "brown")]
#norm.topic.vs.clustersHM<-norm.topic.vs.clusters[,c("brown", "blue", "turquoise", "yellow")]
#norm.topic.vs.clustersHM<-norm.topic.vs.clusters[,c("black", "red", "magenta", "yellow", "purple", "greenyellow",
#                                                    "brown", "salmon", "pink", "turquoise", "tan", "green", "blue")]

## New
#norm.topic.vs.clustersHM<-norm.topic.vs.clusters[,c("green", "brown", "yellow", "blue", "pink", "magenta", "turquoise", 
#                                                    "black", "red")]
#norm.topic.vs.clustersHM<-norm.topic.vs.clusters[,c("yellow", "blue", "turquoise", "green", 
#                                                   "brown", "red")]

norm.topic.vs.clustersHM<-norm.topic.vs.clusters[c(1,2,3,4,5),
                                                 c("black(17)", "red(13)", "green(18)",
                                                   "turquoise(38)","yellow(19)","pink(12)",
                                                   "blue(22)","brown(17)")]

norm.topic.vs.clustersHM2<-norm.topic.vs.clusters[c(1,2,3,4),
                       c("black(17)", "red(13)", "green(18)",
                         "turquoise(38)","yellow(19)","pink(12)",
                         "blue(22)", "brown(17)")]

library(gtools)
library(viridis)
library(gplots)
heatmap.2(t(as.matrix(norm.topic.vs.clusters)), xlab="maxTopic", ylab="Cluster", trace='none', key.title=NA, key.ylab=NA, col=viridis(10), cexRow = 0.8, cexCol = 0.8, key.xlab = "normIntersection", Colv = F, Rowv=F,
          breaks=seq(0,1,0.1))
#pdf(file = "../figures_pdf2/F6C.2V3.pdf", width=7, height = 7) #min clust size = 5
pdf(file = "../figures_pdf2/F6C.2V4.pdf", width=7, height = 7) #min clust size = 10, clust method = ward
#pdf(file = "../figures_pdf2/F6C.2.pdf", width=7, height = 7) #min clust size = 5
#pdf(file = "../figures_pdf2/F6C.3.pdf", width=7, height = 7) #min clust size = 10
#pdf(file = "../figures_pdf2/F6C.4.pdf", width=7, height = 7) #min clust size = 3
heatmap.2(t(as.matrix(norm.topic.vs.clustersHM)), xlab="maxTopic", ylab="Cluster", trace='none', key.title=NA, key.ylab=NA, col=cividis(10), cexRow = 0.8, cexCol = 0.8, key.xlab = "normIntersection", Colv = F, Rowv=F,
          breaks=seq(0,1,0.1))
heatmap.2(t(as.matrix(norm.topic.vs.clustersHM2)), xlab="maxTopic", ylab="Cluster", trace='none', key.title=NA, key.ylab=NA, col=cividis(10), cexRow = 0.8, cexCol = 0.8, key.xlab = "normIntersection", Colv = F, Rowv=F,
          breaks=seq(0,1,0.1))
dev.off()

write.csv(data.frame(odourant=names(clust.colour), cluster=clust.colour), "../supp_tables/odourantsClustersDragon.csv")
write.csv(topic.vs.clusters, "../supp_tables/topicsClustersAssociation_Dragon_DEORs.csv")

######OR TOPICS VS CLUSTER OF ODOURANTS (ONLY GENES ANALYSED IN JOEL'S PAPER) #######

#take all the OR in a given topic, and see which odourant they can detect 
temp<-or.odourant

#keep only the ORs and odourants for which we have a topic estimation and CAS number
topic.est<-row.names(or.info)[which(!is.na(or.info$maxTopic) & or.info$Predicted=="FALSE")]
temp<-temp[temp$Olfr.Nomenclature%in%topic.est,]
temp<-temp[temp$Odorant.CAS.Number%in%names(clust.colour),]


or.joel<-c("Olfr429",
           "Olfr1352",
           "Olfr796",
           "Olfr1395",
           "Olfr1377",
           "Olfr323",
           "Olfr311",
           "Olfr749",
           "Olfr221",
           "Olfr15",
           "Olfr19",
           "Olfr167",
           "Olfr168",
           "Olfr171",
           "Olfr202",
           "Olfr109",
           "Olfr340",
           "Olfr1104",
           "Olfr1264",
           "Olfr1019",
           "Olfr1062",
           "Olfr1079",
           "Olfr1341",
           "Olfr447",
           "Olfr638",
           "Olfr554",
           "Olfr653",
           "Olfr556",
           "Olfr683",
           "Olfr558",
           "Olfr685",
           "Olfr569",
           "Olfr599",
           "Olfr609",
           "Olfr611",
           "Olfr620",
           "Olfr67",
           "Olfr64",
           "Olfr65",
           "Olfr715",
           "Olfr508",
           "Olfr514",
           "Olfr532",
           "Olfr979",
           "Olfr983",
           "Olfr895",
           "Olfr876",
           "Olfr1324")

keep<-which(temp$Olfr.Nomenclature %in% or.joel)

temp<-temp[keep,]

or.used<-unique(temp$Olfr.Nomenclature)
od.used<-unique(temp$Odorant.CAS.Number)

length(or.used)
table(or.info[or.used,"maxTopic"])
table(or.info[or.used,"maxTopic"], or.info[or.used, "Class"])

length(od.used)


#out of the OR used, how many were previously studied?
length(or.used)#OR used
length(which(!is.na(or.info[or.used, "indexMiyamichiReal"])))#by Myamichi
length(which(!is.na(or.info[or.used, "ZolfrReal"])))#by Mombaerts
or.info[or.used, "indexMiyamichiChemS"]
plot(or.info[or.used, "maxTopic"],or.info[or.used, "indexMiyamichiChemS"])#agreement between chem senses indexes and ours

#replace OR with topic and odourant with cluster
temp$Olfr.Nomenclature<-or.info[temp$Olfr.Nomenclature, "maxTopic"]
temp$Odorant.CAS.Number<-clust.colour[temp$Odorant.CAS.Number]



topic.vs.clusters<-table(temp[,2], temp[,3])
chisq.test(topic.vs.clusters)
#norm.topic.vs.clusters<-t(apply(topic.vs.clusters, 1, function(x) round(100*(x/sum(x)), 1)))


norm.topic.vs.clusters<-topic.vs.clusters
for(i in 1:nrow(norm.topic.vs.clusters)){
  for(j in 1:ncol(norm.topic.vs.clusters)){
    #ni<-table(or.info[or.used,"maxTopic"])[i]#number of ORs in topic 1
    nj<-table(clust.colour[od.used])[colnames(topic.vs.clusters)[j]]#number of odourants in cluster j  
    #nj<-table(clust.colour[od.used])[rownames(topic.vs.clusters)[j]]#number of odourants in cluster j  
    
    norm.topic.vs.clusters[i,j]= round(norm.topic.vs.clusters[i,j]/(nj),2)
    
  }
}

norm.topic.vs.clusters<-round(t(apply(norm.topic.vs.clusters, 1, function(x) x/max(x))),2)
#normalized by dividing by tot number of odourant in a cluster and by max in each row

norm.topic.vs.clusters

row.names(norm.topic.vs.clusters)<-c("Topic 1(14)", "Topic 2(8)",
                                     "Topic 3(6)", "Topic 4(1)",
                                     "Topic 5 (1)")
colnames(norm.topic.vs.clusters)<-paste0(colnames(norm.topic.vs.clusters),
                                         "(",
                                         table(clust.colour[od.used])[colnames(norm.topic.vs.clusters)],
                                         ")")
#norm.topic.vs.clustersHM<-norm.topic.vs.clusters[,c("blue", "brown", "turquoise", "yellow")]
#norm.topic.vs.clustersHM<-norm.topic.vs.clusters[,c("yellow", "turquoise", "brown", "blue")]
#norm.topic.vs.clustersHM<-norm.topic.vs.clusters[,c("pink", "red", "blue", "turquoise", "yellow", "brown")]

## New
#norm.topic.vs.clustersHM<-norm.topic.vs.clusters[,c("green", "brown", "yellow", "blue", "pink", "magenta", "turquoise", 
#                                                    "black", "red")]
#norm.topic.vs.clustersHM<-norm.topic.vs.clusters[,c("yellow", "blue", "turquoise", "green", 
#                                                    "brown", "red")]
norm.topic.vs.clustersHM<-norm.topic.vs.clusters[c(1,2,3,4,5),
                                                 c("black(13)", "red(3)", "green(13)",
                                                   "turquoise(28)","yellow(17)","pink(11)",
                                                   "blue(7)","brown(7)")]

norm.topic.vs.clustersHM2<-norm.topic.vs.clusters[c(1,2,3),
                                                  c("black(13)", "red(3)", "green(13)",
                                                    "turquoise(28)","yellow(17)","pink(11)",
                                                    "blue(7)", "brown(7)")]

#norm.topic.vs.clustersHM<-norm.topic.vs.clusters[,c("grey", "blue", "brown", "turquoise", "yellow")]

#pdf(file = "../figures_pdf2/F6x.2.pdf", width=7, height = 7)
#pdf(file = "../figures_pdf2/F6x.2V3.pdf", width=7, height = 7) #min clust size = 5
pdf(file = "../figures_pdf2/F6x.2V4.pdf", width=7, height = 7) #min clust size = 10, clustering method = ward
heatmap.2(t(as.matrix(norm.topic.vs.clustersHM)), xlab="maxTopic", ylab="Cluster", trace='none', key.title=NA, key.ylab=NA, col=cividis(10), cexRow = 0.8, cexCol = 0.8, key.xlab = "normIntersection", Colv = F, Rowv=F,
          breaks=seq(0,1,0.1))
heatmap.2(t(as.matrix(norm.topic.vs.clustersHM2)), xlab="maxTopic", ylab="Cluster", trace='none', key.title=NA, key.ylab=NA, col=cividis(10), cexRow = 0.8, cexCol = 0.8, key.xlab = "normIntersection", Colv = F, Rowv=F,
          breaks=seq(0,1,0.1))
dev.off()

write.csv(topic.vs.clusters, "../supp_tables/topicsClustersAssociation_Dragon_JoelORs.csv")

######OR TOPICS VS CLUSTER OF ODOURANTS (ALSO PREDICTED OR GENES) #######

#take all the OR in a given topic, and see which odourant they can detect 
temp<-or.odourant

#keep only the ORs and odourants for which we have a topic estimation and CAS number
topic.est<-row.names(or.info)[which(!is.na(or.info$maxTopic))]
temp<-temp[temp$Olfr.Nomenclature%in%topic.est,]
temp<-temp[temp$Odorant.CAS.Number%in%names(clust.colour),]

or.used<-unique(temp$Olfr.Nomenclature)
od.used<-unique(temp$Odorant.CAS.Number)

length(or.used)
table(or.info[or.used,"maxTopic"])
table(or.info[or.used,"maxTopic"], or.info[or.used, "Class"])

length(od.used)


#out of the OR used, how many were previously studied?
length(or.used)#OR used
length(which(!is.na(or.info[or.used, "indexMiyamichiReal"])))#by Myamichi
length(which(!is.na(or.info[or.used, "ZolfrReal"])))#by Mombaerts
or.info[or.used, "indexMiyamichiChemS"]
plot(or.info[or.used, "maxTopic"],or.info[or.used, "indexMiyamichiChemS"])#agreement between chem senses indexes and ours

#replace OR with topic and odourant with cluster
temp$Olfr.Nomenclature<-or.info[temp$Olfr.Nomenclature, "maxTopic"]
temp$Odorant.CAS.Number<-clust.colour[temp$Odorant.CAS.Number]



topic.vs.clusters<-table(temp[,2], temp[,3])
chisq.test(topic.vs.clusters)
#norm.topic.vs.clusters<-t(apply(topic.vs.clusters, 1, function(x) round(100*(x/sum(x)), 1)))


norm.topic.vs.clusters<-topic.vs.clusters
for(i in 1:nrow(norm.topic.vs.clusters)){
  for(j in 1:ncol(norm.topic.vs.clusters)){
    #ni<-table(or.info[or.used,"maxTopic"])[i]#number of ORs in topic 1
    nj<-table(clust.colour[od.used])[colnames(topic.vs.clusters)[j]]#number of odourants in cluster j  
    #nj<-table(clust.colour[od.used])[rownames(topic.vs.clusters)[j]]#number of odourants in cluster j  
    
    norm.topic.vs.clusters[i,j]= round(norm.topic.vs.clusters[i,j]/(nj),2)
    
  }
}

norm.topic.vs.clusters<-round(t(apply(norm.topic.vs.clusters, 1, function(x) x/max(x))),2)
#normalized by dividing by tot number of odourant in a cluster and by max in each row

norm.topic.vs.clusters

row.names(norm.topic.vs.clusters)<-c("Topic 1(61)", "Topic 2(36)",
                                     "Topic 3(39)", "Topic 4(6)",
                                     "Topic 5 (1)")
colnames(norm.topic.vs.clusters)<-paste0(colnames(norm.topic.vs.clusters),
                                         "(",
                                         table(clust.colour[od.used])[colnames(norm.topic.vs.clusters)],
                                         ")")

#norm.topic.vs.clustersHM<-norm.topic.vs.clusters[,c("blue", "brown", "turquoise", "yellow")]
#norm.topic.vs.clustersHM<-norm.topic.vs.clusters[,c("yellow", "turquoise", "brown", "blue")]
#norm.topic.vs.clustersHM<-norm.topic.vs.clusters[,c("pink", "black", "blue", "red", "yellow", "turquoise", "green", "brown")]

## New
#norm.topic.vs.clustersHM<-norm.topic.vs.clusters[,c("green", "brown", "yellow", "blue", "pink", "magenta", "turquoise", 
#                                                    "black", "red")]
#norm.topic.vs.clustersHM<-norm.topic.vs.clusters[,c("yellow", "blue", "turquoise", "green", 
#                                                    "brown", "red")]

norm.topic.vs.clustersHM<-norm.topic.vs.clusters[c(1,2,3,4,5),
                                                 c("black(19)", "red(21)", "green(22)",
                                                   "turquoise(41)","yellow(23)","pink(17)",
                                                   "blue(26)","brown(24)")]

norm.topic.vs.clustersHM2<-norm.topic.vs.clusters[c(1,2,3,4),
                                                  c("black(19)", "red(21)", "green(22)",
                                                    "turquoise(41)","yellow(23)","pink(17)",
                                                    "blue(26)","brown(24)")]

#pdf(file = "../figures_pdf2/F6D.2.pdf", width=7, height = 7)
#pdf(file = "../figures_pdf2/F6D.2V3.pdf", width=7, height = 7) #min clust size = 5
pdf(file = "../figures_pdf2/F6D.2V4.pdf", width=7, height = 7) #min clust size = 10, clust method = ward
#heatmap.2(t(as.matrix(norm.topic.vs.clustersHM)), xlab="maxTopic", ylab="Cluster", trace='none', key.title=NA, key.ylab=NA, col=viridis(10), cexRow = 0.8, key.xlab = "normIntersection", Colv = F, Rowv=F,
#          breaks=seq(0,1,0.1))
heatmap.2(t(as.matrix(norm.topic.vs.clustersHM)), xlab="maxTopic", ylab="Cluster", trace='none', key.title=NA, key.ylab=NA, col=cividis(10), cexRow = 0.8, cexCol = 0.8, key.xlab = "normIntersection", Colv = F, Rowv=F,
          breaks=seq(0,1,0.1))
heatmap.2(t(as.matrix(norm.topic.vs.clustersHM2)), xlab="maxTopic", ylab="Cluster", trace='none', key.title=NA, key.ylab=NA, col=cividis(10), cexRow = 0.8, cexCol = 0.8, key.xlab = "normIntersection", Colv = F, Rowv=F,
          breaks=seq(0,1,0.1))
dev.off()

write.csv(topic.vs.clusters, "../supp_tables/topicsClustersAssociation_Dragon_allORs.csv")

######OR DPT INDEXES VS CLUSTER OF ODOURANTS (ALSO PREDICTED OR GENES) ######

#take all the OR in a given topic, and see which odourant they can detect 
temp<-or.odourant

#keep only the ORs and odourants for which we have a topic estimation and CAS number
topic.est<-row.names(or.info)[which(!is.na(or.info$maxTopic))]
temp<-temp[temp$Receptor.Name.Olfr.Nomenclature%in%topic.est,]
temp<-temp[temp$Detected.Odorant.CAS.Number%in%names(clust.colour),]

or.used<-unique(temp$Receptor.Name.Olfr.Nomenclature)
od.used<-unique(temp$Detected.Odorant.CAS.Number)

length(or.used)
table(or.info[or.used,"maxTopic"])
table(or.info[or.used,"maxTopic"], or.info[or.used, "Class"])

length(od.used)

id<-temp$Receptor.Name.Olfr.Nomenclature
#replace OR with DPT and odourant with cluster
temp$Receptor.Name.Olfr.Nomenclature<-as.numeric(or.info[temp$Receptor.Name.Olfr.Nomenclature, "DPTindex"])
temp$Detected.Odorant.CAS.Number<-clust.colour[temp$Detected.Odorant.CAS.Number]



boxplot(temp$Receptor.Name.Olfr.Nomenclature~as.factor(temp$Detected.Odorant.CAS.Number), 
        las=2, ylab="DPT index", xlab="odourant clusters")




topic.vs.clusters<-table(temp[,1], temp[,2])

#norm.topic.vs.clusters<-t(apply(topic.vs.clusters, 1, function(x) round(100*(x/sum(x)), 1)))


norm.topic.vs.clusters<-topic.vs.clusters
for(i in 1:nrow(norm.topic.vs.clusters)){
  for(j in 1:ncol(norm.topic.vs.clusters)){
    #ni<-table(or.info[or.used,"maxTopic"])[i]#number of ORs in topic 1
    nj<-table(clust.colour[od.used])[colnames(topic.vs.clusters)[j]]#number of odourants in cluster j  
    
    norm.topic.vs.clusters[i,j]= round(norm.topic.vs.clusters[i,j]/(nj),2)
    
  }
}

norm.topic.vs.clusters<-round(t(apply(norm.topic.vs.clusters, 1, function(x) x/max(x))),2)
#normalized by dividing by tot number of odourant in a cluster and by max in each row

norm.topic.vs.clusters



######OR CLASS VS CLUSTER OF ODOURANTS (ONLY DETECTED OR GENES) #######

#take all the OR in a given topic, and see which odourant they can detect 
temp<-or.odourant

#keep only the ORs and odourants for which we have a topic estimation and CAS number
topic.est<-row.names(or.info)[which(!is.na(or.info$maxTopic) & or.info$Predicted=="FALSE")]
temp<-temp[temp$Receptor.Name.Olfr.Nomenclature%in%topic.est,]
temp<-temp[temp$Detected.Odorant.CAS.Number%in%names(clust.colour),]

or.used<-unique(temp$Receptor.Name.Olfr.Nomenclature)
od.used<-unique(temp$Detected.Odorant.CAS.Number)

length(or.used)
table(or.info[or.used,"maxTopic"])
table(or.info[or.used,"maxTopic"], or.info[or.used, "Class"])

length(od.used)

#consider only ORs in topic 1
select<-unique(temp$Receptor.Name.Olfr.Nomenclature[or.info[temp$Receptor.Name.Olfr.Nomenclature,"maxTopic"]==1])
temp<-temp[temp$Receptor.Name.Olfr.Nomenclature%in%select,]

or.info[temp$Receptor.Name.Olfr.Nomenclature, "Class"]

#replace OR with class and odourant with cluster
temp$Receptor.Name.Olfr.Nomenclature<-or.info[temp$Receptor.Name.Olfr.Nomenclature, "Class"]
temp$Detected.Odorant.CAS.Number<-clust.colour[temp$Detected.Odorant.CAS.Number]



topic.vs.clusters<-table(temp[,1], temp[,2])

norm.topic.vs.clusters<-topic.vs.clusters
for(i in 1:nrow(norm.topic.vs.clusters)){
  for(j in 1:ncol(norm.topic.vs.clusters)){
    #ni<-table(or.info[or.used,"maxTopic"])[i]#number of ORs in topic 1
    nj<-table(clust.colour[od.used])[colnames(topic.vs.clusters)[j]]#number of odourants in cluster j  
    
    norm.topic.vs.clusters[i,j]= round(norm.topic.vs.clusters[i,j]/(nj),2)
    
  }
}

norm.topic.vs.clusters<-round(t(apply(norm.topic.vs.clusters, 1, function(x) x/max(x))),2)
norm.topic.vs.clusters


######OR CLASS VS CLUSTER OF ODOURANTS (ALSO PREDICTED OR GENES) #######

#take all the OR in a given topic, and see which odourant they can detect 
temp<-or.odourant

#keep only the ORs and odourants for which we have a topic estimation and CAS number
topic.est<-row.names(or.info)[which(!is.na(or.info$maxTopic))]
temp<-temp[temp$Receptor.Name.Olfr.Nomenclature%in%topic.est,]
temp<-temp[temp$Detected.Odorant.CAS.Number%in%names(clust.colour),]

or.used<-unique(temp$Receptor.Name.Olfr.Nomenclature)
od.used<-unique(temp$Detected.Odorant.CAS.Number)

length(or.used)
table(or.info[or.used,"maxTopic"])
table(or.info[or.used,"maxTopic"], or.info[or.used, "Class"])

length(od.used)

#consider only ORs in topic 1
select<-unique(temp$Receptor.Name.Olfr.Nomenclature[or.info[temp$Receptor.Name.Olfr.Nomenclature,"maxTopic"]==1])
temp<-temp[temp$Receptor.Name.Olfr.Nomenclature%in%select,]

or.info[temp$Receptor.Name.Olfr.Nomenclature, "Class"]

#replace OR with class and odourant with cluster
temp$Receptor.Name.Olfr.Nomenclature<-or.info[temp$Receptor.Name.Olfr.Nomenclature, "Class"]
temp$Detected.Odorant.CAS.Number<-clust.colour[temp$Detected.Odorant.CAS.Number]



topic.vs.clusters<-table(temp[,1], temp[,2])

norm.topic.vs.clusters<-topic.vs.clusters
for(i in 1:nrow(norm.topic.vs.clusters)){
  for(j in 1:ncol(norm.topic.vs.clusters)){
    #ni<-table(or.info[or.used,"maxTopic"])[i]#number of ORs in topic 1
    nj<-table(clust.colour[od.used])[colnames(topic.vs.clusters)[j]]#number of odourants in cluster j  
    
    norm.topic.vs.clusters[i,j]= round(norm.topic.vs.clusters[i,j]/(nj),2)
    
  }
}

norm.topic.vs.clusters<-round(t(apply(norm.topic.vs.clusters, 1, function(x) x/max(x))),2)
norm.topic.vs.clusters

######OR CLASS VS CLUSTER OF ODOURANTS (ONLY GENES ANALYSED IN JOEL'S PAPER) #######

#take all the OR in a given topic, and see which odourant they can detect 
temp<-or.odourant

#keep only the ORs and odourants for which we have a topic estimation and CAS number
topic.est<-row.names(or.info)[which(!is.na(or.info$maxTopic) & or.info$Predicted=="FALSE")]
temp<-temp[temp$Receptor.Name.Olfr.Nomenclature%in%topic.est,]
temp<-temp[temp$Detected.Odorant.CAS.Number%in%names(clust.colour),]

or.joel<-c("Olfr429",
           "Olfr1352",
           "Olfr796",
           "Olfr1395",
           "Olfr1377",
           "Olfr323",
           "Olfr311",
           "Olfr749",
           "Olfr221",
           "Olfr15",
           "Olfr19",
           "Olfr167",
           "Olfr168",
           "Olfr171",
           "Olfr202",
           "Olfr109",
           "Olfr340",
           "Olfr1104",
           "Olfr1264",
           "Olfr1019",
           "Olfr1062",
           "Olfr1079",
           "Olfr1341",
           "Olfr447",
           "Olfr638",
           "Olfr554",
           "Olfr653",
           "Olfr556",
           "Olfr683",
           "Olfr558",
           "Olfr685",
           "Olfr569",
           "Olfr599",
           "Olfr609",
           "Olfr611",
           "Olfr620",
           "Olfr67",
           "Olfr64",
           "Olfr65",
           "Olfr715",
           "Olfr508",
           "Olfr514",
           "Olfr532",
           "Olfr979",
           "Olfr983",
           "Olfr895",
           "Olfr876",
           "Olfr1324")

keep<-which(temp$Receptor.Name.Olfr.Nomenclature %in% or.joel)

temp<-temp[keep,]

or.used<-unique(temp$Receptor.Name.Olfr.Nomenclature)
od.used<-unique(temp$Detected.Odorant.CAS.Number)




length(or.used)
table(or.info[or.used,"maxTopic"])
table(or.info[or.used,"maxTopic"], or.info[or.used, "Class"])

length(od.used)

#replace OR with class and odourant with cluster
temp$Receptor.Name.Olfr.Nomenclature<-or.info[temp$Receptor.Name.Olfr.Nomenclature, "Class"]
temp$Detected.Odorant.CAS.Number<-clust.colour[temp$Detected.Odorant.CAS.Number]



topic.vs.clusters<-table(temp[,1], temp[,2])

norm.topic.vs.clusters<-topic.vs.clusters
for(i in 1:nrow(norm.topic.vs.clusters)){
  for(j in 1:ncol(norm.topic.vs.clusters)){
    #ni<-table(or.info[or.used,"maxTopic"])[i]#number of ORs in topic 1
    nj<-table(clust.colour[od.used])[colnames(topic.vs.clusters)[j]]#number of odourants in cluster j  
    
    norm.topic.vs.clusters[i,j]= round(norm.topic.vs.clusters[i,j]/(nj),2)
    
  }
}

norm.topic.vs.clusters<-round(t(apply(norm.topic.vs.clusters, 1, function(x) x/max(x))),2)
norm.topic.vs.clusters



#####CONVERT CAS-NUMBERS TO CHEMICAL NAMES IN THE CLUSTERS ######

cas.chem<-read.csv(file="../supp_tables/Dsstox_CAS_number_name.csv", sep=";")
row.names(cas.chem)<-cas.chem[,1]

cas.name<-cas.chem$preferred_name
names(cas.name)<-cas.chem$casrn

clust.name.colour<-clust.colour
names(clust.name.colour)<-cas.name[names(clust.colour)]

names(clust.name.colour[clust.name.colour=="green"])
names(clust.name.colour[clust.name.colour=="blue"])
names(clust.name.colour[clust.name.colour=="brown"])
names(clust.name.colour[clust.name.colour=="turquoise"])
names(clust.name.colour[clust.name.colour=="yellow"])


# PLOT ODOURANT PROPERTIES AS FUNCTION OF THEIR CLUSTER ########

df<-data.frame(cluster=clust.colour, property=physico.chem[names(clust.colour), "ALOGP"])
boxplot(df$property~df$cluster, las=2)

#pdf(file = "../figures_pdf2/F6B.2.pdf", width=7, height = 7)
#pdf(file = "../figures_pdf2/F6B.2V3.pdf", width=7, height = 7) #min clust size = 5
pdf(file = "../figures_pdf2/F6B.2V4.pdf", width=7, height = 7) #min clust size = 10, clust method = ward

for(property in c("MW", "Hy", "ALOGP", "ALOGP2", "MLOGP", "MLOGP2", "nRCOOH", "nArCOOH", "Sp",
                  "Mp", "Pol", "UNIP", "TPSA.NO.", "TPSA.Tot.")){
  df<-data.frame(cluster=clust.colour, property=physico.chem[names(clust.colour), property])
  boxplot(df$property~df$cluster, las=2, ylab=property, xlab="cluster")
}

for(property in c("MW", "Hy", "ALOGP", "ALOGP2", "MLOGP", "MLOGP2", "nRCOOH", "nArCOOH", "Sp",
                  "Mp", "Pol", "UNIP", "TPSA.NO.", "TPSA.Tot.")){
  df<-data.frame(or.odourant, maxTopic=or.info[or.odourant$data.Receptor.Name.Olfr.Nomenclature, "maxTopic"], property=physico.chem[or.odourant$data.Detected.Odorant.CAS.Number, property])
  boxplot(df$property~df$maxTopic, las=2, ylab=property, xlab="Topic")
}

#for(property in c("MW", "nH", "nO", "nS", "nAromBond", "nHBAcc", "nHBDon", "SLogP")){
#  df<-data.frame(cluster=clust.colour, property=physico.chem[names(clust.colour), property])
#  boxplot(df$property~df$cluster, las=2, ylab=property, xlab="cluster")
#}

#for(property in c("MW", "nH", "nO", "nS", "nAromBond", "nHBAcc", "nHBDon", "SLogP")){
#  df<-data.frame(or.odourant, maxTopic=or.info[or.odourant$Olfr.Nomenclature, "maxTopic"], property=physico.chem[or.odourant$Odorant.CAS.Number, property])
#  boxplot(df$property~df$maxTopic, las=2, ylab=property, xlab="Topic")
#}

dev.off()

#Find physico-chem properties that are different between clusters #####
g1<-"green"
g2<-"red"
#collect properties
prop.g1<-physico.chem[names(clust.colour)[clust.colour==g1],]
prop.g2<-physico.chem[names(clust.colour)[clust.colour==g2],]
tests<-vector(length=ncol(physico.chem))
names(tests)<-colnames(physico.chem)
for(i in 1:ncol(prop.g1)){
  tests[i]<-wilcox.test(prop.g1[,i],
                        prop.g2[,i])$p.value
}
tests<-sort(tests, decreasing = F)
tests[1:50]
boxplot(prop.g1[,"ALOGP"],
        prop.g2[,"ALOGP"])

# par(mfrow=c(1,2))
# 
# df<-data.frame(cluster=clust.colour[intersect(names(clust.colour), row.names(solub))],
#                property=solub[intersect(names(clust.colour), row.names(solub)), "V3"])
# boxplot(df$property~df$cluster,
#         ylim=c(-12,2),
#         las=2)



boxplot(solub$V3)



####CHECK SOLUBILITY OF THE ODOURANTS DETECTED BY CLASS I ODOURANTS ########

od.classI<-unique(or.odourant[or.info[or.odourant$Receptor.Name.Olfr.Nomenclature,"Class"]==1,"Detected.Odorant.CAS.Number"])
od.classI<-od.classI[!is.na(od.classI)]#find the odourants that are detected by Class I ORs

od.classII<-unique(or.odourant[or.info[or.odourant$Receptor.Name.Olfr.Nomenclature,"Class"]==2,"Detected.Odorant.CAS.Number"])
od.classII<-od.classII[!is.na(od.classII)]#find the odourants that are detected by Class II ORs


par(mfrow=c(1,1))
boxplot(solub[setdiff(od.classI, od.classII),3], 
        solub[setdiff(od.classII, od.classI),3], 
        ylab="Solubility", names=c("Only Class I", "Only Class II"))

par(mfrow=c(1,1))
boxplot(physico.chem[setdiff(od.classI, od.classII),"Hy"], 
        physico.chem[setdiff(od.classII, od.classI),"Hy"], 
        ylab="Hydrophilic factor", names=c("Only Class I", "Only Class II"))


#check correlation between solubility and ALOGP
plot(solub[intersect(row.names(solub), row.names(physico.chem)), "V3"], 
     physico.chem[intersect(row.names(solub), row.names(physico.chem)), "ALOGP"], 
     ylab="ALOGP", xlab="Solubility")

#check correlation between solubility and Hy
plot(solub[intersect(row.names(solub), row.names(physico.chem)), "V3"], 
     physico.chem[intersect(row.names(solub), row.names(physico.chem)), "Hy"], 
     ylab="Hy", xlab="Solubility", col=clust.colour[intersect(row.names(solub), row.names(physico.chem))])




#####HOW MANY OR/ODOURANTS ARE WE REALLY USING? #########
length(unique(or.odourant$Receptor.Name.Olfr.Nomenclature))#OR in the odourant/or table
length(unique(or.odourant$Detected.Odorant.CAS.Number))#odourant in the odourant/or table
dim(physico.chem)#Odourants for which we have physico-chem properties (and CAS number)
length(unique(or.odourant[or.odourant$Detected.Odorant.CAS.Number%in%names(clust.colour),"Receptor.Name.Olfr.Nomenclature"]))
# ORs corresponding to the odourants with physico-chem properties

sum(table(or.info[unique(or.odourant[or.odourant$Detected.Odorant.CAS.Number%in%names(clust.colour),"Receptor.Name.Olfr.Nomenclature"]), "maxTopic"]))
#...out of these, we have that many with estimated topic membership

#Finally, we are doing the analysis on that many ORs and odourants:
#52 ORsx 74odourants
#Out of the 52, we have 27 in T1, 11 in T2, 10 in T3, 3 in T4, and only 1 in T5


#That becomes 81 if we consider also predicted ORs:
#41 in T1, 23 in T2, 13 in T3, 3 in T4 and 1 in T5
#we haven't gained anything in 4 and 5...










##################*************BITS & PIECES #################


id<-names(predicted.topic)[predicted.topic==2 & 
                         or.info[names(predicted.topic), "Class"]==1]
or.info[id, "ZolfrReal"]


id<-(detected)[or.info[(detected), "Class"]==1]
table(or.info[id,"maxTopic"])

or.info["Olfr653",]


temp$Receptor.Name.Olfr.Nomenclature<-or.info[temp$Receptor.Name.Olfr.Nomenclature, "Class"]

id1<-unique(temp[or.info[temp$Receptor.Name.Olfr.Nomenclature, "Class"]==1,"Detected.Odorant.CAS.Number"])
id2<-unique(temp[or.info[temp$Receptor.Name.Olfr.Nomenclature, "Class"]==2,"Detected.Odorant.CAS.Number"])


boxplot(solub[id1,3], solub[id2,3])





(sort(pca$rotation[,1], decreasing=F))[1:50]

(sort(pca$rotation[,2], decreasing=F))[1:50]

pca$rotation[,2]

#mod<-glmnet(x = pca$x[,1:20], y = y, family = "mgaussian")
mod<-lm(y[,2]~pca$x[,1:20])
summary(mod)

mod.shuffled<-lm(y.shuffled[,2]~pca$x[,1:20])
summary(mod.shuffled)

#plot(pca$x[,8],y[,1])
df<-data.frame(x=pca$x[,4], y=y[,4], z=as.factor(all.dob[row.names(physico.chem),"max_zone"]))

plot.std(df, "PC", "DOB", legend="", title="")



##############

sort(pca$rotation[,1], decreasing=T)[1:50]

#############

estimate<-(predict(mod, as.data.frame(pca$x[-remove,1:20])))
plot(estimate, y[-remove,2], ylim=c(0,1), xlim=c(0,1))
abline(a=0,b=1, col="red")

cor(estimate, y[-remove,2])



estimate.shuffled<-(predict(mod.shuffled, as.data.frame(pca$x[,1:20])))
plot(estimate.shuffled, y[,2], ylim=c(0,1), xlim=c(0,1))
abline(a=0,b=1, col="red")

cor(estimate.shuffled, y[,2])



# install.packages("igraph") 
# install.packages("network") 
# install.packages("sna") 
# install.packages("ndtv")

require(igraph)
require(network)     
require(sna)
require(ndtv)

nodes<-c(unique(or.odourant[,1]), unique(or.odourant[,2]))
edges<-or.odourant


net <- graph_from_data_frame(d=edges, vertices=nodes) 

V(net)$size <- 3
V(net)$label <- names(V(net))
E(net)$width <- 1
E(net)$arrow.size <- .1

V(net)[1:length(unique(or.odourant[,1]))]$type="Olfr"
V(net)[(length(unique(or.odourant[,1]))+1):(length(unique(or.odourant[,1]))+length(unique(or.odourant[,2])))]$type="Odourant"

colrs <- c("gray50", "tomato")
names(colrs)<-c("Olfr","Odourant")
V(net)$color <- colrs[V(net)$type]

plot(net)



##### Analysis of distance between ORs detecting the same odourants ######

#select only odourants for which we have more than 1 associated OR
od<-names(which(table(or.odourant[,2])>1))

#for each odourant, find the corresponding ORs and calculate the distance between them #########

#using DPT
cv<-vector(length=length(od))
for(id in 1:length(od)){
  #od<-id[2]
  or<-or.odourant[which(or.odourant[,"Detected.Odorant.CAS.Number"]==od[id]),"Receptor.Name.Olfr.Nomenclature"]
  
  or<-or[which(or%in%row.names(or.info))]
  
  or.dpt<-or.info[or,"DPTindex"]
  #or.dpt<-or.dpt[!is.na(or.dpt)]
  
  cv[id]<-sd(or.dpt)/mean(or.dpt)
  
  
}


#calculate the same quantity, but after shuffling the matrix
rand.or.odourant<-or.odourant[which(or.odourant$Detected.Odorant.CAS.Number%in%od),]
rand.or.odourant[,2]<-rand.or.odourant[sample(1:nrow(rand.or.odourant)),2]

table(rand.or.odourant$Detected.Odorant.CAS.Number)

#using DPT
rand.cv<-vector(length=length(od))
for(id in 1:length(od)){
  #od<-id[2]
  or<-rand.or.odourant[which(rand.or.odourant[,"Detected.Odorant.CAS.Number"]==od[id]),"Receptor.Name.Olfr.Nomenclature"]
  
  or<-or[which(or%in%row.names(or.info))]
  
  or.dpt<-or.info[or,"DPTindex"]
  #or.dpt<-or.dpt[!is.na(or.dpt)]
  
  rand.cv[id]<-sd(or.dpt)/mean(or.dpt)
  
  
}

boxplot(cv, rand.cv)
wilcox.test(cv, rand.cv)


#The difference doesn't look significant, maybe there are certain odourants that are more localized
#than expected by chance, but overall the signal is weak

#COMPARE TOPIC DISTANCE AND ODOURANT DISTANCE BETWEEN OLF REC GENES ######

#for each pair of olfactory receptor genes (at least those for which we have odourant info), compute the dpt distance
id<-unique(or.odourant$Receptor.Name.Olfr.Nomenclature)

id<-id[id%in%row.names(or.info)]

id<-id[!is.na(or.info[id,"T1"])]


pairs.or<-t(combn(id, 2))

dim(pairs.or)

euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

topic.distance<-vector(length=nrow(pairs.or))
for(i in 1:length(topic.distance)){
  
 
  g1<-pairs.or[i,1]
  g2<-pairs.or[i,2]
  
  t1<-or.info[g1, c("T1","T2","T3","T4","T5")]
  t2<-or.info[g2, c("T1","T2","T3","T4","T5")]
  
  topic.distance[i]<-euc.dist(t1,t2)
}


################
# Look at the physico-chem properties of all the odourants detected in T1, T2, ... #####


#select the ORs that are included in the OR/odourant table
or.included<-intersect(unique(or.odourant[,1]), row.names(or.info)[!is.na(or.info$maxTopic)])

table(or.info[or.included,"maxTopic"])


#select all the odourant in each topic 


odourant.topic<-list()

topics<-1:5
for(t in topics){
  id<-or.included[or.info[or.included,"maxTopic"]==t]
  odourant.topic[[t]]<-unique(or.odourant[which(or.odourant$Receptor.Name.Olfr.Nomenclature%in%id),"Detected.Odorant.CAS.Number"])
  
}
 

## plot the odourant in the PCA or diff map colored according to whether they are detected in topic 1-5 or not


data.scaled<-apply(physico.chem, 2, function(x) (x-mean(x))/sd(x))

dm<-DiffusionMap(data=data.scaled)

plot(dm$DC1, dm$DC2)
dpt<-DPT(dm)

par(mfrow=c(3,2))
for(topic in 1:5){
label=vector(length=nrow(physico.chem))

names(label)<-row.names(physico.chem)
label[intersect(names(label),odourant.topic[[topic]])]<-"TRUE"

df<-data.frame(x=dm$DC1, y=dm$DC2, z=label)


p1<-plot.std(df, "DC1", "DC2", legend="", title=paste0("Topic ", topic))
print(p1)

}

dpt<-DPT(dm)

dpt.odourant<-dpt$DPT2
names(dpt.odourant)<-row.names(physico.chem)

########

boxplot(dpt.odourant[odourant.topic[[1]]], 
        dpt.odourant[odourant.topic[[2]]],
        dpt.odourant[odourant.topic[[3]]],
        dpt.odourant[odourant.topic[[4]]],
        dpt.odourant[odourant.topic[[5]]],
        dpt.odourant)

wilcox.test(dpt.odourant[odourant.topic[[5]]], dpt.odourant)



