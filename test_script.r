## load packages
library("phytools")
library("nlme")
library("multcomp")
library("evomap")
library("geiger")
library("phylolm")

phy<-read.tree(file = "tree_64_sp.tree")

table<-read.table("Complete_Table_VF.tsv", header = TRUE)
table

table <- table[match(phy$tip.label,table$Genus_species),]

source("Prehensile-tail-evolution-main_Alexis/functions.R")

extract.group<-function(tree, data_frame, group1){
  vec_gr1<-c()
  vec_gr2<-c()
  for (i in 1:nrow(data_frame)){
    if(data_frame[i,3]==group1){
      vec_gr1<-c(vec_gr1, which(tree$tip.label==data_frame[i,1]))
    }
    else{
      vec_gr2<-c(vec_gr2, which(tree$tip.label==data_frame[i,1]))
    }
  }
  return(list(vec_gr1, vec_gr2))
}

AICc<-function(model, tree){
  AICc_value <- -2 * model$logLik + (2 * as.numeric(model$p) * (length(tree$tip.label)/(length(tree$tip.label) - as.numeric(model$p) - 1)))
  return(AICc_value)
}

spp<-table[,1]

pancova<-function(temp_data, phy, df, iterator){
  
  print(colnames(temp_data)[1])
  
  temp_data[,1]<-log(temp_data[,1])
  
  rownames(temp_data)<-df[,1]
  
  colnames(temp_data)<-c("Dependent","Independent")
  
  ####################################################
  #Case study 1: Platyrrhines versus Catarrhines
  
  #Set groups to be compared
  Non_p<-unlist(extract.group(phy, df, "N")[1])
  P<-unlist(extract.group(phy, df, "N")[2])
  
  #Plot the two groups being compared
  plot(temp_data$Dependent~temp_data$Independent,col="white",pch=19,xlab="", ylab="",asp=1,cex.lab=2)
  pGLS.plotGrade("Dependent","Independent",temp_data,phy,model="BM",group=Non_p,col="#D2B48C",lwd=5,cex=1,pch=19)
  pGLS.plotGrade("Dependent","Independent",temp_data,phy,model="BM",group=P,col="#40826D",lwd=5,cex=1,pch=19)
  #Highlight them in the tree to double check
  tipCol<-rep("#40826D",length(phy$tip.label))
  tipCol[Non_p]<-"#D2B48C"
  plot(phy,tip.col=tipCol,cex=0.6)
  df$Prehensility<-as.factor(df$Prehensility)
  df_temp<-df[,c(1,3,iterator,10)]
  colnames(df_temp)<-c("Genus_species", "Group", "Y", "X")
  print(df_temp)
  corBM<-corBrownian(phy=phy,form=~spp)
  primate.ancova<-gls(log(Y)~log(X)*Group -1, correlation=corBM, data = df_temp)
  post.hoc<-glht(primate.ancova,linfct=mcp(Group="Tukey"))
  print(summary(post.hoc))
  df<-primate.ancova$dims$N-primate.ancova$dims$p
  alpha<-0.05
  ciBeta_Cathemeral<-coef(primate.ancova)["log(X)"]+
  qt(c(alpha/2,1-alpha/2),df)*
  sqrt(primate.ancova$varBeta["log(X)","log(X)"])
  ciBeta_Cathemeral

  ciBeta_Nocturnal<-coef(primate.ancova)["log(X)"]+
  coef(primate.ancova)["log(X):GroupP"]+
  qt(c(alpha/2,1-alpha/2),df)*
  sqrt(primate.ancova$varBeta[ "log(X)","log(X)"]+
  primate.ancova$varBeta["log(X):GroupP","log(X):GroupP"]+2*primate.ancova$varBeta["log(X)","log(X):GroupP"])
  CIs<-data.frame(parameter= paste("beta(",levels(table$Prehensility),")",sep=""),lower=c(ciBeta_Cathemeral[1],ciBeta_Nocturnal[1]), estimate=coef(primate.ancova)["log(X)"]+c(0, coef(primate.ancova)[c("log(X):GroupP")]), upper=c(ciBeta_Cathemeral[2],ciBeta_Nocturnal[2]))
  rownames(CIs)<-NULL
  CIs
  print(intervals(primate.ancova))
}

for (i in c(4:9)){
  pancova(table[,c(i,10)], phy, table, i)
}  


#Test PGLS only with climbing species
table
table_preh <- table%>%filter(Prehensility=="P")
table_np <- table%>%filter(Prehensility=="N")
table_climbing <- table_np%>%filter(Climbing=="CLI")
table_climbing_VF <- rbind(table_preh, table_climbing)

phy_cl <- phy
setdiff(table_climbing_VF$Genus_species,phy_cl$tip.label)
todrop2 <- setdiff(phy_cl$tip.label,table_climbing_VF$Genus_species)
phy_cl<-drop.tip(phy_cl,todrop2)
plot(phy_cl)

table_climbing_VF <- table_climbing_VF[match(phy_cl$tip.label,table_climbing_VF$Genus_species),]


for (i in c(4:9)){
  pancova(table_climbing_VF[,c(i,10)], phy_cl, table_climbing_VF, i, table_climbing_VF[,1])
}  


##PGLS all specimens 1 group
#tuto pgls phytools
#creer modele brownien
bm<-corBrownian(1,phy)
bm
#modele de pgls
#RI all vertebrae
model1<-gls(log(RImean)~log(sizelsr), data=table, correlation=bm)
plot(log(RImean)~log(sizelsr),data=table,cex=1.5,pch=19,bg="grey", col=c("#D2B48C","#40826D")[as.factor(table$Prehensility)])
text(log(RImean)~log(sizelsr),data=table, labels=table$Nomenclature, cex = 0.6, pos = 1)
abline(model1,lwd=2,col="darkgrey",lty="dashed")
summary(model1)
residuals(model1)

#RI distal
model1<-gls(log(RIdistmean)~log(sizelsr), data=table, correlation=bm)
plot(log(RIdistmean)~log(sizelsr),data=table,cex=1.5,pch=19,bg="grey", col=c("#D2B48C","#40826D")[as.factor(table$Prehensility)])
text(log(RIdistmean)~log(sizelsr),data=table, labels=table$Nomenclature, cex = 0.6, pos = 1)
abline(model1,lwd=2,col="darkgrey",lty="dashed")
summary(model1)
residuals(model1)

#RI last quarter
model1<-gls(log(RIlqmean)~log(sizelsr), data=table, correlation=bm)
plot(log(RIlqmean)~log(sizelsr),data=table,cex=1.5,pch=19,bg="grey", col=c("#D2B48C","#40826D")[as.factor(table$Prehensility)])
text(log(RIlqmean)~log(sizelsr),data=table, labels=table$Nomenclature, cex = 0.6, pos = 1)
abline(model1,lwd=2,col="darkgrey",lty="dashed")
summary(model1)
residuals(model1)

#TPEI all vertebrae
model1<-gls(log(TPEImean)~log(sizelsr), data=table, correlation=bm)
plot(log(TPEImean)~log(sizelsr),data=table,cex=1.5,pch=19,bg="grey", col=c("#D2B48C","#40826D")[as.factor(table$Prehensility)])
text(log(TPEImean)~log(sizelsr),data=table, labels=table$Nomenclature, cex = 0.6, pos = 1)
abline(model1,lwd=2,col="darkgrey",lty="dashed")
summary(model1)
residuals(model1)

#TPEI distal
model1<-gls(log(TPEIdistmean)~log(sizelsr), data=table, correlation=bm)
plot(log(TPEIdistmean)~log(sizelsr),data=table,cex=1.5,pch=19,bg="grey", col=c("#D2B48C","#40826D")[as.factor(table$Prehensility)])
text(log(TPEIdistmean)~log(sizelsr),data=table, labels=table$Nomenclature, cex = 0.6, pos = 1)
abline(model1,lwd=2,col="darkgrey",lty="dashed")
summary(model1)
residuals(model1)

#TPEI last quarter
model1<-gls(log(TPEIlqmean)~log(sizelsr), data=table, correlation=bm)
plot(log(TPEIlqmean)~log(sizelsr),data=table,cex=1.5,pch=19,bg="grey", col=c("#D2B48C","#40826D")[as.factor(table$Prehensility)])
text(log(TPEIlqmean)~log(sizelsr),data=table, labels=table$Nomenclature, cex = 0.6, pos = 1)
abline(model1,lwd=2,col="darkgrey",lty="dashed")
summary(model1)
residuals(model1)

#ANOVA phylo Alexis (13/05/24)
#RI all
vec <- table$Prehensility
names(vec) <- table$Genus_species
vecy <- table$RImean
names(vecy) <- table$Genus_species
phylANOVA(phy,x = vec, y = vecy, nsim = 1000)
#RI dist
vec <- table$Prehensility
names(vec) <- table$Genus_species
vecy <- table$RIdistmean
names(vecy) <- table$Genus_species
phylANOVA(phy,x = vec, y = vecy, nsim = 1000)
#RI lastq
vec <- table$Prehensility
names(vec) <- table$Genus_species
vecy <- table$RIlqmean
names(vecy) <- table$Genus_species
phylANOVA(phy,x = vec, y = vecy, nsim = 1000)
#TPEI all
vec <- table$Prehensility
names(vec) <- table$Genus_species
vecy <- table$TPEImean
names(vecy) <- table$Genus_species
phylANOVA(phy,x = vec, y = vecy, nsim = 1000)
#TPEI dist
vec <- table$Prehensility
names(vec) <- table$Genus_species
vecy <- table$TPEIdistmean
names(vecy) <- table$Genus_species
phylANOVA(phy,x = vec, y = vecy, nsim = 1000)
#TPEI lastq
vec <- table$Prehensility
names(vec) <- table$Genus_species
vecy <- table$TPEIlqmean
names(vecy) <- table$Genus_species
phylANOVA(phy,x = vec, y = vecy, nsim = 1000)

i<-6

df_temp<-table[,c(2,3,i,10)]
colnames(df_temp)<-c("Genus_species", "Group", "Y", "X")

corBM<-corBrownian(phy=phy,form=~spp)
primate.ancova<-gls(log(Y)~log(X)*
                      Group -1,
                    correlation=corBM,
                    data = df_temp)
primate.ancova
post.hoc<-glht(primate.ancova,linfct=mcp(Group="Tukey"))
summary(post.hoc)