library(vegan)
library(ggplot2)
library(grid)
library(zoo)
library(reshape2)
library(reshape)
library(diversityR)
library(dplyr)
library(lme4)
library(rgl)
library(ggpubr)
library(lmerTest)
library(nlme)
#make the OTU tables for repeated samples, caculate the similarity of repeated samples
dts<-read.csv("stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample_2.shared",header=T,sep = "\t")
rS<-dts[,4:ncol(dts)]
dS<- as.matrix(vegdist(rS,method="bray"))
dS<-1-dS
reS<-c()
n=1
while(n<ncol(dS)){
   m<-dS[n,n+1]
   reS<-c(reS,m)
   n=n+2
} 
fn=paste("Repeated_Similarity",".csv",sep="")
write.csv(reS,fn,row.names=TRUE)

#Fig1a Caculate averaged pairwised similarity within each baseline/within generous donor
dts1<-read.csv("stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample_1.shared",header = T,sep = "\t")
AVE<-c()
dts<-dts1[c(237:257),]#Subject4
dts<-dts1[c(258:315),]#Subject5A
dts<-dts1[c(316:415),]#Subject5B
dts<-dts1[c(416:430),]#Subject5C
dts<-dts1[c(431:485),]#Subject5D
dts<-dts1[c(486:527),]#Subject9
dts<-dts1[c(1:21),]#Subject15
dts<-dts1[c(22:96),]#Subject17A
dts<-dts1[c(97:193),]#Subject17B
dts<-dts1[c(194:236),]#Subject23

#caculation for different generous donors
dts<-dts1[c(538:539),]#GD2
dts<-dts1[c(540:545),]#GD3
dts<-dts1[c(546:555),]#GD


rS<-dts[,4:ncol(dts1)]
dS<- as.matrix(vegdist(rS,method="bray"))
dS<-1-dS
ave<-(rowSums(dS)-1)/(ncol(dS)-1)
AVE<-c(AVE,ave)

fn=paste("Averaged_SimilarityToOthers_Raw",".csv",sep="")
write.csv(AVE,fn,row.names=TRUE)

#Fig.1_plot Averaged BC-similairty of each baseline and replicate controls
md<-read.csv("Averaged_SimilarityToOthers_meta.csv",header = T,sep = ",")
myColors <-c("gold1","dodgerblue","purple","red","#C39953","#568203","gray50","black","dark red")
names(myColors) <- levels(md$Subject)
colScale <- scale_colour_manual(name = "Subject",values = myColors)
#version1
ggplot(md,aes(x=md$Baseline,y=md$Similarity))+geom_rect(xmin = 11-0.5,xmax = 12+0.5,ymin=-Inf,ymax=Inf,fill="gray88",alpha=0.2)+
  geom_boxplot(outlier.size = 0,outlier.shape=NA)+colScale+
  geom_jitter(aes(color=as.factor(md$Subject)),shape=16,size=1.6,position=position_jitter(0.3),alpha=0.65)+
  labs(x="                               Baseline periods                                        Controls       ", y="Bray-Curtis \n similarity")+
  theme(axis.title.y = element_text(angle=0,vjust=0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),legend.position ="none",title = element_text(size = rel(1.8)),legend.text = element_text(size = 18),legend.key=element_blank(),legend.key.size=unit(1.0,"cm"),
        axis.text.x= element_text(size=20,margin=margin(t=5)),axis.text.y = element_text(size=20,margin=margin(r=5)),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))

#Linear mixed model based on Fig.1 
md<-read.csv("Averaged_SimilarityToOthers_meta.csv",header = T,sep = ",")
fit<-lmer(md$Similarity~md$Group+(1|md$Subject),data=md)
summary(fit)

#Fig.2 Bacterial load of each baseline
md<-read.csv("Averaged_SimilarityToOthers_meta.csv",header = T,sep = ",")
md<-md[1:527,]
myColors <-c("gold1","dodgerblue","purple","red","#C39953","#568203","gray50","black","dark red")
names(myColors) <- levels(md$Subject)
colScale <- scale_colour_manual(name = "Subject",values = myColors)
#Fig.2_plot
ggplot(md,aes(x=md$Baseline,y=md$BacterialLoad))+geom_boxplot(outlier.size = 0,outlier.shape=NA)+colScale+geom_jitter(aes(color=as.factor(md$Subject)),shape=16,size=1.6,position=position_jitter(0.3),alpha=0.65)+labs(x="Baseline periods", y="Log(Copies/mL)")+theme(axis.title.y = element_text(angle=0,vjust=0.5))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),legend.position ="none",title = element_text(size = rel(1.8)),legend.text = element_text(size = 18),legend.key=element_blank(),legend.key.size=unit(1.0,"cm"),
axis.text.x= element_text(size=20,margin=margin(t=5)),axis.text.y = element_text(size=20,margin=margin(r=5)), axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))

#Fig.3 Antibiotic associated the bacterial community change
dts1<-read.csv("stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample_1.shared",header = T,sep = "\t")
dts1<-dts1[1:527,]
myColors <-c("dodgerblue","gold1","purple","red","#C39953","#568203")
dts1$Abx<-factor(dts1$Abx,levels = c("None","Inhaled colistin","Azithromycin","Azithromycin + inhaled aztreonam","Azithromycin + inhaled tobramycin"),labels = c("None","Inhaled colistin","Azithromycin","Azithromycin + inhaled aztreonam","Azithromycin + inhaled tobramycin"))
names(myColors) <- levels(dts1$Abx)
colScale <- scale_colour_manual(name = "Maintenance antibiotics",values = myColors)

dts<-dts1[c(237:257),]#Subject4
dts<-dts1[c(258:485),]#Subject5
dts<-dts1[c(258:315),]#Subject5A_Abx
dts<-dts1[c(316:415),]#Subject5B_Abx
dts<-dts1[c(416:430),]#Subject5C
dts<-dts1[c(431:485),]#Subject5D_Abx
dts<-dts1[c(486:527),]#Subject9_Abx
dts<-dts1[c(1:21),]#Subject15
dts<-dts1[c(22:193),]#Subject17
dts<-dts1[c(22:96),]#Subject17A_Abx
dts<-dts1[c(97:193),]#Subject17B_Abx
dts<-dts1[c(194:236),]#Subject23

rs<-dts[,c(6:ncol(dts))]
set.seed(6)
ord<-metaMDS(rs,distance = "bray")
ord
Dismatrix<-ord$diss
nmdsdf<-data.frame(axis1 = ord$points[,1], axis2 = ord$points[,2],group=dts$Abx)
#centroid
gg<-data.frame(merge(nmdsdf,aggregate(cbind(mean.x=axis1 ,mean.y=axis2)~group,nmdsdf,mean),by="group"))
ggplot(data=gg, aes(gg$axis1,gg$axis2,color=gg$group))+geom_point(size=3.0)+xlim(-1.3,1)+ylim(-0.8,1)+geom_segment(aes(x=gg$mean.x, y=gg$mean.y,xend=gg$axis1,yend=gg$axis2))+geom_point(aes(x=gg$mean.x,y=gg$mean.y),size=1.5,stroke=2,color="Black")+labs(x="nMDS1", y="nMDS2",color="Baseline")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),title = element_text(size = rel(1.8)),legend.position ="none",legend.text = element_text(size = 18),legend.key=element_blank(),legend.key.size=unit(1.0,"cm"),axis.text = element_text(size=16))+colScale                                                

#Legend
dts<-dts1[c(22:193,316:415),]#Subject with Abx_change
sdts<-dts[,c(6:ncol(dts))]
rs<-sdts
set.seed(6)
ord<-metaMDS(rs,distance = "bray")
ord
Dismatrix<-ord$diss
nmdsdf<-data.frame(axis1 = ord$points[,1], axis2 = ord$points[,2],group=dts$Abx)
p1<-ggplot(data = nmdsdf, aes(nmdsdf$axis1, nmdsdf$axis2)) + geom_point(aes(color = nmdsdf$group),size=5.5)+labs(x="nMDS1", y="nMDS2")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=2),title = element_text(size = rel(2.3)),legend.text = element_text(size = 22),legend.key=element_blank(),legend.key.size=unit(1.1,"cm"),axis.text = element_text(size=16))+colScale+guides(col = guide_legend(ncol = 2)) 
library(gridExtra)
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 
legend <- g_legend(p1)
grid.draw(legend)

#Fig.4A_subject2 the microbial community of different baselines within subject2
dts<-read.csv("Subject2_4Baseline_Replicates.shared",header = T,sep = "\t")
rs<-dts[,c(7:ncol(dts))]
set.seed(5)
ord<-metaMDS(rs,distance = "bray")
ord 
Dismatrix<-ord$diss
dts$Baseline <- factor(dts$Baseline,levels = c("2a","2b","2c","2d"),labels = c("2a","2b","2c","2d")) 
nmdsdf<-data.frame(axis1 = ord$points[,1], axis2 = ord$points[,2],group=dts$Baseline,pair=dts$Pair)
myColors <-c("gold1","dodgerblue","red","purple")
names(myColors) <- levels(dts$Baseline)
colScale <- scale_colour_manual(name = "Baseline",values = myColors)
nmdsdf1<-nmdsdf[c(1:228),]
#plot
gg1<-data.frame(merge(nmdsdf1,aggregate(cbind(mean.x=axis1 ,mean.y=axis2)~group,nmdsdf1,mean),by="group"))
gg1$group1<-factor(gg1$group,levels = c("2a","2b","2c","2d"),labels = c("a","b","c","d"))
ggplot(data=gg1, aes(gg1$axis1,gg1$axis2,color=gg1$group))+xlim(-1.3,1.2)+ylim(-0.7,0.8)+geom_point(size=2.5,alpha=0.6)+geom_segment(aes(x=gg1$mean.x, y=gg1$mean.y,xend=gg1$axis1,yend=gg1$axis2),alpha=0.6)+geom_point(aes(x=gg1$mean.x,y=gg1$mean.y),size=7.3,shape=15,color="white")+geom_text(aes(x=gg1$mean.x,y=gg1$mean.y),label=gg1$group1,size=8,color="Black")+labs(x="nMDS1", y="nMDS2",color="Baseline")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),title = element_text(size = rel(1.8)),legend.key=element_blank(),legend.text = element_text(size = 20),legend.key.size=unit(1.0,"cm"),axis.text = element_text(size=16))+colScale                                              
#Fig.4B_subject5 the microbial community of different baselines within subject5
dts<-read.csv("Subject5_2Baseline_Replicates.shared",header = T,sep = "\t")
rs<-dts[,c(7:ncol(dts))]
set.seed(5)
ord<-metaMDS(rs,distance = "bray")
ord 
Dismatrix<-ord$diss
dts$Baseline <- factor(dts$Baseline,levels = c("5a","5b"),labels = c("5a","5b")) 
nmdsdf<-data.frame(axis1 = ord$points[,1], axis2 = ord$points[,2],group=dts$Baseline,pair=dts$Pair)
myColors <-c("gold1","dodgerblue","red","purple")
names(myColors) <- levels(dts$Baseline)
colScale <- scale_colour_manual(name = "Baseline",values = myColors)
nmdsdf1<-nmdsdf[c(1:172),]
#plot
gg1<-data.frame(merge(nmdsdf1,aggregate(cbind(mean.x=axis1 ,mean.y=axis2)~group,nmdsdf1,mean),by="group"))
gg1$group1<-factor(gg1$group,levels = c("5a","5b"),labels = c("a","b"))
ggplot(data=gg1, aes(gg1$axis1,gg1$axis2,color=gg1$group))+xlim(-1.3,1.2)+ylim(-0.7,0.8)+geom_point(size=2.5,alpha=0.6)+geom_segment(aes(x=gg1$mean.x, y=gg1$mean.y,xend=gg1$axis1,yend=gg1$axis2),alpha=0.6)+geom_point(aes(x=gg1$mean.x,y=gg1$mean.y),size=7.3,shape=15,color="white")+geom_text(aes(x=gg1$mean.x,y=gg1$mean.y),label=gg1$group1,size=8,color="Black")+labs(x="nMDS1", y="nMDS2",color="Baseline")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),title = element_text(size = rel(1.8)),legend.key=element_blank(),legend.text = element_text(size = 20),legend.key.size=unit(1.0,"cm"),axis.text = element_text(size=16))+colScale                                              

#Fig.5 Interval days
#Interval days_BC similarity
dts1<-read.csv("Subject6.shared",header=T,sep = "\t")
rS<-dts1[,4:ncol(dts1)]
dS<- as.matrix(vegdist(rS,method="bray"))
M<-NULL
#i indicate the interval days
for (i in 1:nrow(rS)){
  da<-c()
  n=1
  b=nrow(rS)-i
  while(n<=b){
    m<-dS[n,n+i]
    da<-c(da,m)
    n<-n+1
  }
  M<-cbind(M,da)
}
for (i in 1: nrow(M)){
  for (j in 1: ncol(M)){
    if (i>nrow(rS)-j){
      M[i,j]=0
    }
  }
}
TransM<-t(M)
rownames(TransM)<-c(1:ncol(TransM))
md<-melt(TransM,id.vars=rownames(TransM))
##remove rows with value of 0,1
md<-md[!(md$value==0|md$value==1|md$value=="NA"),]
md<-md[complete.cases(md), ]
colnames(md) <- c("Interval Day", "Start Day","Dissimilarity")
md$Similarity<-1-md$Dissimilarity
fn=paste("Subject6_Interval",".csv",sep="")
write.csv(md,fn,row.names=TRUE)
#Combine all subjects together and filter the target interval days for plot

#Fig.5_plot
md<-read.csv("AllSubject_Interval_5Days.csv",header = T,sep = ",")
myColors <-c("#dba404","gray45")
names(myColors) <- levels(md$Abx)
colScale <- scale_colour_manual(name = "Maintenance \n  antibiotics",values = myColors)
md$Interval.Day<-as.factor(md$Interval.Day)
ggplot(md,aes(x=md$Interval.Day,y=md$Similarity))+geom_boxplot(outlier.size = 0,outlier.shape=NA,lwd=1)+ylim(0.1,1)+
  colScale+geom_jitter(aes(color=as.factor(md$Abx)),shape=16,size=1.2,position=position_jitter(0.31),alpha=0.75)+labs(x="Interval days", y="Bray-Curtis \n similarity")+theme(axis.title.y = element_text(angle=0,vjust=0.5),legend.position ="none")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),title = element_text(size = rel(1.8)),legend.key=element_blank(),legend.text = element_text(size = 18),legend.key.size=unit(1.0,"cm"),axis.text.y = element_text(size=16,margin=margin(r=5)),axis.text.x = element_text(size=16,margin=margin(t=5)))
#Legend
p1<-ggplot(md,aes(x=md$Interval.Day,y=md$Similarity))+geom_boxplot(outlier.size = 0,outlier.shape=NA,lwd=1)+ylim(0.1,1)+
  colScale+geom_jitter(aes(color=as.factor(md$Abx)),shape=16,size=3.5,position=position_jitter(0.31),alpha=0.75)+labs(x="Interval days", y="Bray-Curtis \n similarity")+theme(axis.title.y = element_text(angle=0,vjust=0.5))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),title = element_text(size = rel(1.8)),legend.key=element_blank(),legend.text = element_text(size = 18),legend.key.size=unit(1.0,"cm"),axis.text.y = element_text(size=16,margin=margin(r=5)),axis.text.x = element_text(size=16,margin=margin(t=5)))
library(gridExtra)
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 
legend <- g_legend(p1)
grid.draw(legend)

##Model based on Fig.5
md<-read.csv("AllSubject_Interval.csv",header = T,sep = ",")
fit<-lmer(md$Similarity~md$Interval.Day+md$Abx+(1|md$Subject)+(0+md$Interval.Day|md$Subject)+(1|md$Baseline),data=md)
summary(fit)


#Fig.E1
# ranked group barchart_by water 
dds<-read.table("Sample_control_water.shared",header=T,sep = "\t")
dd<- subset(dds, Water>0.01)
dd$Rank[order((dd[,4]))] <- 1:nrow(dd)
dd<-dd[,c(1,2,3,4,5,6)]
ddt<-melt(dd,id.vars = c("Group","Taxa","Rank"))
colnames(ddt)[4] <- "Type"
colnames(ddt)[5] <- "Relative_abundance"
ddt_sort<- ddt[order(ddt$Rank), ]
myorder <- ddt_sort %>% filter(Type=="Water") %>% arrange(desc(Relative_abundance)) %>% .$Group %>% as.character
mytaxa<-ddt_sort %>% filter(Type=="Water") %>% arrange(desc(Relative_abundance)) %>% .$Taxa %>% as.character
ddt_sort<- ddt[order(ddt$Rank), ]
ggplot(ddt_sort, aes(x=Group, y=Relative_abundance, fill = Type)) + labs(x="Taxa", y="Relative \n abundance")+
  geom_bar(stat="identity", position = "dodge") + scale_fill_manual(values=c("grey","Black","#FFBF00"))+
  scale_x_discrete(limits=myorder,labels=mytaxa)+ylim(0,0.2)+
  theme(title = element_text(size = rel(1.5)),axis.title.y = element_text(angle=0,vjust=0.5,size=18,margin = margin(r=15)),axis.title.x = element_text(size=18,margin = margin(t=15)),legend.position ="none",legend.text=element_text(size = 17),legend.key.size=unit(0.7,"cm"),axis.text.y = element_text(size=16,margin=margin(r=5)),axis.text.x = element_text(angle=25,size=16,hjust=1,vjust=1,margin=margin(t=5)),panel.background = element_blank(),axis.line = element_line(colour = "black", size=1))
# ranked group barchart_by reagent control 
dds<-read.table("Sample_control_water.shared",header=T,sep = "\t")
colnames(dds)[5] <- "Reagent control"
dd<- subset(dds, dds$`Reagent control`>0.01)
dd$Rank[order((dd[,5]))] <- 1:nrow(dd)
dd<-dd[,c(1,2,3,4,5,6)]
ddt<-melt(dd,id.vars = c("Group","Taxa","Rank"))
colnames(ddt)[4] <- "Type"
colnames(ddt)[5] <- "Relative_abundance"
ddt_sort<- ddt[order(ddt$Rank), ]
myorder <- ddt_sort %>% filter(Type=="Reagent control") %>% arrange(desc(Relative_abundance)) %>% .$Group %>% as.character
mytaxa<-ddt_sort %>% filter(Type=="Reagent control") %>% arrange(desc(Relative_abundance)) %>% .$Taxa %>% as.character
ddt_sort<- ddt[order(ddt$Rank), ]
ggplot(ddt_sort, aes(x=Group, y=Relative_abundance, fill = Type)) + labs(x="Taxa", y="Relative \n abundance")+
  geom_bar(stat="identity", position = "dodge") + scale_fill_manual(values=c("grey","Black","#FFBF00"))+
  scale_x_discrete(limits=myorder,labels=mytaxa)+ylim(0,0.2)+
  theme(title = element_text(size = rel(1.5)),axis.title.y = element_text(angle=0,vjust=0.5,size=18,margin = margin(r=15)),axis.title.x = element_text(size=18,margin = margin(t=15)),legend.position ="none",legend.text=element_text(size = 17),legend.key.size=unit(0.7,"cm"),axis.text.y = element_text(size=16,margin=margin(r=5)),axis.text.x = element_text(angle=25,size=16,hjust=1,vjust=1,margin=margin(t=5)),panel.background = element_blank(),axis.line = element_line(colour = "black", size=1))
# ranked group barchart_by Samples 
dds<-read.table("Sample_control_water.shared",header=T,sep = "\t")
dd<- subset(dds, dds$Samples>0.01)
dd$Rank[order((dd[,6]))] <- 1:nrow(dd)
dd<-dd[,c(1,2,3,4,5,6)]
ddt<-melt(dd,id.vars = c("Group","Taxa","Rank"))
colnames(ddt)[4] <- "Type"
colnames(ddt)[5] <- "Relative_abundance"
ddt_sort<- ddt[order(ddt$Rank), ]
myorder <- ddt_sort %>% filter(Type=="Samples") %>% arrange(desc(Relative_abundance)) %>% .$Group %>% as.character
mytaxa<-ddt_sort %>% filter(Type=="Samples") %>% arrange(desc(Relative_abundance)) %>% .$Taxa %>% as.character
ddt_sort<- ddt[order(ddt$Rank), ]
ggplot(ddt_sort, aes(x=Group, y=Relative_abundance, fill = Type)) + labs(x="Taxa", y="Relative \n abundance")+
  geom_bar(stat="identity", position = "dodge") + scale_fill_manual(values=c("grey","Black","#FFBF00"))+
  scale_x_discrete(limits=myorder,labels=mytaxa)+ylim(0,0.2)+
  theme(title = element_text(size = rel(1.5)),axis.title.y = element_text(angle=0,vjust=0.5,size=18,margin = margin(r=15)),axis.title.x = element_text(size=18,margin = margin(t=15)),legend.position ="none",legend.text=element_text(size = 17),legend.key.size=unit(0.7,"cm"),axis.text.y = element_text(size=16,margin=margin(r=5)),axis.text.x = element_text(angle=25,size=16,hjust=1,vjust=1,margin=margin(t=5)),panel.background = element_blank(),axis.line = element_line(colour = "black", size=1))
##Legend
p1<-ggplot(ddt_sort, aes(x=Group, y=Relative_abundance, fill = Type)) + labs(x="OTUs", y="Relative \n abundance")+
  geom_bar(stat="identity", position = "dodge") + scale_fill_manual(values=c("grey","Black","#FFBF00"))+
  scale_x_discrete(limits=myorder,labels=mytaxa)+ylim(0,0.2)+
  theme(title = element_text(size = rel(1.5)),axis.title.y = element_text(angle=0,vjust=0.5,size=18,margin = margin(r=15)),axis.title.x = element_text(size=18,margin = margin(t=15)),legend.text=element_text(size = 17),legend.key.size=unit(0.7,"cm"),axis.text.y = element_text(size=16,margin=margin(r=5)),axis.text.x = element_text(angle=25,size=16,hjust=1,vjust=1,margin=margin(t=5)),panel.background = element_blank(),axis.line = element_line(colour = "black", size=1))
library(gridExtra)
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 
legend <- g_legend(p1)
grid.draw(legend)

##Fig.E2 stackedBar plot with missing day dashed linetype (Average day to all other day version)
col<-c("grey","#A52A2A","#C39953","#DBE9F4","#89CFF0","#FFAA1D","#CC0000","#006A4E","#873260","#536872","#B5A642","#CB4154","#1DACD6","#BF94E4","#66FF00","#F4C2C2","#009E73","#CD7F32","#737000","#007FFF","#E7FEFF","#7BB661","#F0DC82","#CC5500","#FFC1CC","#5F9EA0","#ED872D","#702963","#FFF600","#A67B5B","#4B3621","#A3C1AD","#A9A9A9","#3B7A57","#FFBF00","#A52A2A","#B284BE","#0048BA")
dts1<-read.table("stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample_1_taxa.shared",header=T,sep="\t")
dts<-dts1[c(237:257),]#Subject4
dts<-dts1[c(258:315),]#Subject5A
dts<-dts1[c(316:415),]#Subject5B
dts<-dts1[c(416:430),]#Subject5C
dts<-dts1[c(431:485),]#Subject5D
dts<-dts1[c(486:527),]#Subject9
dts<-dts1[c(1:21),]#Subject15
dts<-dts1[c(22:96),]#Subject17A
dts<-dts1[c(97:193),]#Subject17B
dts<-dts1[c(194:236),]#Subject23

minsub1<-cbind(dts[,2],dts[,6:ncol(dts)])
colnames(minsub1)[1] <- "Group"
minsub1<-melt(minsub1,id.vars = "Group")
time.x <- as.numeric(minsub1[,1])
begin <- time.x[1]
time.x <- as.numeric(minsub1$Group)-begin+1
colnames(minsub1)[2] <- "OTUs"

#Day to the other day similarity
dts2<-read.csv("Averaged_SimilarityToOthers_meta_1.csv",header = T,sep = ",")
#for baseline with no missing days
rS<-dts2[c(237:257),]#Subject4
rS<-dts2[c(416:430),]#Subject5C
rS<-dts2[c(431:485),]#Subject5D

time.x1 <- as.numeric(rS[,1])
begin <- time.x1[1]
time.x1<- as.numeric(rS$Sample)-begin+1
dtd<-data.frame(cbind(time.x1,rS$Similarity))

#plot with label day interval 10
ggplot() + 
  geom_bar(data = minsub1,aes(x=time.x,y=value, fill=OTUs), stat = "identity") + scale_fill_manual(values = col)+labs(x="Daily samples",y="Relative abundance")+theme(axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)),title = element_text(size = rel(1.8)),legend.text=element_text(size = 18),legend.key.size=unit(0.6,"cm"),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18),legend.position ="none",panel.background = element_blank(),axis.line = element_line(colour = "gray"))+scale_x_continuous("Daily samples",breaks=c(0,10,20,40,60,80,100,120),expand = c(0.001,0.001))+scale_y_continuous(expand = c(0,0))+
  geom_point(data =dtd, aes(time.x1, dtd$V2),size=2.5)+geom_line(data = dtd, aes(time.x1, dtd$V2))
#plot with label day interval 20 
ggplot() + 
  geom_bar(data = minsub1,aes(x=time.x,y=value, fill=OTUs), stat = "identity") + scale_fill_manual(values = col)+labs(x="Daily samples",y="Relative abundance")+theme(axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)),title = element_text(size = rel(1.8)),legend.text=element_text(size = 18),legend.key.size=unit(0.6,"cm"),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18),legend.position ="none",panel.background = element_blank(),axis.line = element_line(colour = "gray"))+scale_x_continuous("Daily samples",breaks=c(0,20,40,60,80,100,120),expand = c(0.001,0.001))+scale_y_continuous(expand = c(0,0))+
  geom_point(data =dtd, aes(time.x1, dtd$V2),size=2.5)+geom_line(data = dtd, aes(time.x1, dtd$V2))

#for baseline with missing days
rS<-dts2[c(258:315),]#Subject5A
rS<-dts2[c(316:415),]#Subject5B
rS<-dts2[c(486:527),]#Subject9
rS<-dts2[c(1:21),]#Subject15
rS<-dts2[c(22:96),]#Subject17A
rS<-dts2[c(97:193),]#Subject17B
rS<-dts2[c(194:236),]#Subject23

time.x1 <- as.numeric(rS[,1])
begin <- time.x1[1]
time.x1<- as.numeric(rS$Sample)-begin+1
dtd<-data.frame(cbind(time.x1,rS$Similarity))
dtd$Segment<- rep(1,nrow(dtd))
dtd$linet<- rep(1,nrow(dtd))
dtd<-within(dtd,linet[dtd$Segment==2] <-"dashed")
fn=paste("dtd_6_prep",".csv",sep="")
write.csv(dtd,fn,row.names=F)
#edit in .cvs file with dummy days 
dtd<-read.csv("dtd_6_prep.csv",header = T,sep = ",")
gg<-ggplot() + geom_bar(data = minsub1,aes(x=time.x,y=value, fill=OTUs), stat = "identity") + 
  scale_fill_manual(values = col)+
  labs(x="Daily samples",y="Relative abundance")+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)),title = element_text(size = rel(1.8)),legend.text=element_text(size = 18),legend.key.size=unit(0.6,"cm"),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18),legend.position ="none",panel.background = element_blank(),axis.line = element_line(colour = "gray"))+scale_x_continuous("Daily samples",breaks=c(0,20,40,60,80,100,120),expand = c(0.001,0.001))+scale_y_continuous(expand = c(0,0))+
  geom_point(data =dtd, aes(time.x1, dtd$V2),size=2.5)
for (i in 1:7) gg <- gg + 
  geom_line(data=subset(dtd, Segment==i), 
            aes(x=time.x1, y=V2,linetype=as.character(linet), group=i))
gg
##Legend
p1<-ggplot() + 
  geom_bar(data = minsub1,aes(x=time.x,y=value, fill=OTUs), stat = "identity") + scale_fill_manual(values = col)+labs(x="Daily samples",y="Relative abundance")+theme(axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)),title = element_text(size = rel(1.8)),legend.text=element_text(size = 18),legend.key.size=unit(0.6,"cm"),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18),panel.background = element_blank(),axis.line = element_line(colour = "gray"))+scale_x_continuous("Daily samples",breaks=c(0,20,40,60,80,100,120),expand = c(0.001,0.001))+scale_y_continuous(expand = c(0,0))+
  geom_point(data =dtd, aes(time.x1, dtd$V2),size=2.5)+geom_line(data = dtd, aes(time.x1, dtd$V2))
library(gridExtra)
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 
legend <- g_legend(p1)
grid.draw(legend)

#Fig.E3
#define the outliers 
md<-read.csv("Averaged_SimilarityToOthers_meta.csv",header = T,sep = ",")
#define outliers by baseline
md1<-md[c(237:257),]#Subject4
md1<-md[c(258:315),]#Subject5A
md1<-md[c(316:415),]#Subject5B
md1<-md[c(416:430),]#Subject5C
md1<-md[c(431:485),]#Subject5D
md1<-md[c(486:527),]#Subject9
md1<-md[c(1:21),]#Subject15
md1<-md[c(22:96),]#Subject17A
md1<-md[c(97:193),]#Subject17B
md1<-md[c(194:236),]#Subject23
OutVals = boxplot(md1$Similarity)$out
length(which(md1$Similarity %in% OutVals))
M1<-cbind(which(md1$Similarity %in% OutVals),OutVals)
#Fig.E3_plot
md<-read.csv("Averaged_SimilarityToOthers_meta.csv",header = T,sep = ",")
md<-md[c(1:527),]
ggplot(md,aes(x=md$Outlier_defined,y=md$BacterialLoad))+ylim(7,10.3)+geom_boxplot(outlier.size = 0,outlier.shape=NA)+colScale+geom_jitter(shape=16,size=1.6,position=position_jitter(0.3),alpha=0.65)+labs(x=" ", y="Log(Copies/mL)")+theme(axis.title.y = element_text(angle=0,vjust=0.5))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),legend.position ="none",title = element_text(size = rel(1.8)),legend.text = element_text(size = 18),legend.key=element_blank(),legend.key.size=unit(1.0,"cm"),
axis.text.x= element_text(size=20,margin=margin(t=5)),axis.text.y = element_text(size=20,margin=margin(r=5)), axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))

#model of Fig.E3
fit<-lmer(md$BacterialLoad~md$Outlier_defined+(1|md$Subject)+(1|md$Baseline),data=md)
summary(fit)


#Fig.E4. TBL change of Abx effect
md<-read.csv("Averaged_SimilarityToOthers_meta.csv",header = T,sep = ",")
md<-md[1:527,]
myColors<-c("dodgerblue","gold1","purple","red","#C39953","#568203")
md$Abx<-factor(md$Abx,levels = c("None","Inhaled colistin","Azithromycin","Azithromycin + inhaled aztreonam","Azithromycin + inhaled tobramycin"),labels = c("None","Inhaled colistin","Azithromycin","Azithromycin + inhaled aztreonam","Azithromycin + inhaled tobramycin"))
names(myColors) <- levels(md$Abx)
colScale <- scale_colour_manual(name = "Maintenance antibiotics",values = myColors)

md1<-md[c(258:315),]#Subject5A_Abx
md1<-md[c(316:415),]#Subject5B_Abx
md1<-md[c(431:485),]#Subject5D_Abx
md1<-md[c(486:527),]#Subject9_Abx
md1<-md[c(22:96),]#Subject17A_Abx
md1<-md[c(97:193),]#Subject17B_Abx

ggplot(md1,aes(x=md1$Abx,y=md1$BacterialLoad))+geom_boxplot(outlier.size = 0,outlier.shape=NA)+ylim(6.5,10.6)+
  colScale+geom_jitter(aes(color=as.factor(md1$Abx)),shape=16,size=2,position=position_jitter(0.3),alpha=0.65)+
  labs(x="", y="Log(Copies/mL)")+theme(axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)),panel.grid.major = element_blank(),legend.position ="none",axis.text.x = element_text(angle=25,size=18,hjust=1,vjust=1),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),title = element_text(size = rel(1.8)),legend.text = element_text(size = 18),legend.key=element_blank(),legend.key.size=unit(1.0,"cm"),axis.text = element_text(size=16))

#model of Fig.E4
md<-read.csv("Averaged_SimilarityToOthers_meta0.csv",header = T,sep = ",")
md<-md[1:527,]
md1<-md[md$Baseline == '2b',] #used baseline 2a as example
m1 = gls(BacterialLoad~Abx+Day, correlation = corAR1(form=~Day),data=md1)
summary(m1)


##Fig.E5-Subject2
dts<-read.csv("Subject2_4Baseline_Replicates.shared",header = T,sep = "\t")
rs<-dts[,c(7:ncol(dts))]
set.seed(5)
ord<-metaMDS(rs,distance = "bray")
ord 
Dismatrix<-ord$diss
dts$Baseline <- factor(dts$Baseline,levels = c("2a","2b","2c","2d"),labels = c("2a","2b","2c","2d")) 
nmdsdf<-data.frame(axis1 = ord$points[,1], axis2 = ord$points[,2],group=dts$Baseline,pair=dts$Pair)
myColors <-c("gold1","dodgerblue","red","purple")
names(myColors) <- levels(dts$Baseline)
colScale <- scale_colour_manual(name = "Baseline",values = myColors)
nmdsdf1<-nmdsdf[c(1:228),]
gg1<-data.frame(merge(nmdsdf1,aggregate(cbind(mean.x=axis1 ,mean.y=axis2)~group,nmdsdf1,mean),by="group"))
fn=paste("gg1_output",".csv",sep="")
write.csv(gg1,fn,row.names=TRUE)
fn=paste("nmdsdf_output",".csv",sep="")
write.csv(nmdsdf,fn,row.names=TRUE)
#insert repeated day pionts for the dashed line 
gg<-read.csv("gg_input.csv",header = T,sep = ",")
myshapes <-c(16,1)
names(myshapes) <- levels(gg$Type)
shapeScale <- scale_shape_manual(name = "Type",values = myshapes)
gg1<-gg[c(1:228),]
gg2<-gg[c(229:243),]
ggplot()+geom_point(data=gg1, aes(gg1$axis1,gg1$axis2,color=gg1$group,shape=gg1$Type),size=3,color="grey",alpha=0.2)+xlim(-1.3,1.2)+ylim(-0.7,0.8)+
  geom_segment(aes(x=gg1$mean.x, y=gg1$mean.y,xend=gg1$axis1,yend=gg1$axis2),alpha=0.2)+
  geom_point(data=gg2, aes(gg2$axis1,gg2$axis2,color=gg2$group,shape=gg2$Type),size=3)+
  geom_point(data=gg2, aes(gg2$r1,gg2$r2,color=gg2$group),size=3)+
  geom_segment(aes(x=gg2$r1, y=gg2$r2,xend=gg2$axis1,yend=gg2$axis2),color="black")+shapeScale+
  guides(color = guide_legend(order=1),shape = guide_legend(order=2))+
  labs(x="nMDS1", y="nMDS2",color="Baseline")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),title = element_text(size = rel(1.8)),legend.key=element_blank(),legend.text = element_text(size = 18),legend.key.size=unit(1.0,"cm"),axis.text = element_text(size=16))+colScale                                              
#Fig.E5_subject5
dts<-read.csv("Subject5_2Baseline_Replicates.shared",header = T,sep = "\t")
rs<-dts[,c(7:ncol(dts))]
set.seed(5)
ord<-metaMDS(rs,distance = "bray")
ord 
Dismatrix<-ord$diss
dts$Baseline <- factor(dts$Baseline,levels = c("5a","5b"),labels = c("5a","5b")) 
nmdsdf<-data.frame(axis1 = ord$points[,1], axis2 = ord$points[,2],group=dts$Baseline,pair=dts$Pair)
myColors <-c("gold1","dodgerblue","red","purple")
names(myColors) <- levels(dts$Baseline)
colScale <- scale_colour_manual(name = "Baseline",values = myColors)
nmdsdf1<-nmdsdf[c(1:172),]
gg1<-data.frame(merge(nmdsdf1,aggregate(cbind(mean.x=axis1 ,mean.y=axis2)~group,nmdsdf1,mean),by="group"))
fn=paste("gg1_s5_output",".csv",sep="")
write.csv(gg1,fn,row.names=TRUE)
fn=paste("nmdsdf_s5_output",".csv",sep="")
write.csv(nmdsdf,fn,row.names=TRUE)
gg<-read.csv("gg_s5_input.csv",header = T,sep = ",")
myshapes <-c(16,1)
names(myshapes) <- levels(gg$Type)
shapeScale <- scale_shape_manual(name = "Type",values = myshapes)
#version2_seperate replicates only
gg1<-gg[c(1:172),]
gg2<-gg[c(173:179),]
ggplot()+geom_point(data=gg1, aes(gg1$axis1,gg1$axis2,color=gg1$group,shape=gg1$Type),size=3,color="grey",alpha=0.2)+xlim(-1.1,0.8)+ylim(-0.7,0.8)+
  geom_segment(aes(x=gg1$mean.x, y=gg1$mean.y,xend=gg1$axis1,yend=gg1$axis2),alpha=0.2)+
  geom_point(data=gg2, aes(gg2$axis1,gg2$axis2,color=gg2$group,shape=gg2$Type),size=3)+
  geom_point(data=gg2, aes(gg2$r1,gg2$r2,color=gg2$group),size=3)+
  geom_segment(aes(x=gg2$r1, y=gg2$r2,xend=gg2$axis1,yend=gg2$axis2),color="black")+shapeScale+
  guides(color = guide_legend(order=1),shape = guide_legend(order=2))+
  labs(x="nMDS1", y="nMDS2",color="Baseline")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),title = element_text(size = rel(1.8)),legend.key=element_blank(),legend.text = element_text(size = 18),legend.key.size=unit(1.0,"cm"),axis.text = element_text(size=16))+colScale                                              

#Fig.E6
##Same Abx
md<-read.csv("Subject6_Interval_5days.csv",header = T,sep = ",")
md$Interval.Day<-as.factor(md$Interval.Day)
ggplot(md,aes(x=md$Interval.Day, y=md$Similarity, fill=md$Abx))+
  geom_boxplot(alpha = 0.80)+
  scale_fill_manual(values=c("gray"))+
  geom_point(size=1.2, shape=NA, position = position_jitterdodge())+ylim(0.1,1)+
  labs(x="Interval Days", y="Bray-Curtis \n similarity",fill="Maintenance \n antibiotics")+theme(axis.title.y = element_text(angle=0,vjust=0.5))+
  theme(panel.grid.major = element_blank(), legend.position ="none",panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),title = element_text(size = rel(1.8)),legend.text = element_text(size = 18),legend.key=element_blank(),legend.key.size=unit(1.0,"cm"),axis.text.y = element_text(size=16,margin=margin(r=5)),axis.text.x = element_text(size=16,margin=margin(t=5)))

##Seperate by Abx
md<-read.csv("Subject5_Interval_5days.csv",header = T,sep = ",")
md$Interval.Day<-as.factor(md$Interval.Day)
ggplot(md,aes(x=md$Interval.Day, y=md$Similarity, fill=md$Abx))+
  geom_boxplot(alpha = 0.80)+
  scale_fill_manual(values=c("#FFBF00","gray"))+
  geom_point(size=1.2, shape=NA, position = position_jitterdodge())+ylim(0.1,1)+
  labs(x="Interval Days", y="Bray-Curtis \n similarity",fill="Maintenance \n antibiotics")+theme(axis.title.y = element_text(angle=0,vjust=0.5))+
  theme(panel.grid.major = element_blank(), legend.position ="none",panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),title = element_text(size = rel(1.8)),legend.text = element_text(size = 18),legend.key=element_blank(),legend.key.size=unit(1.0,"cm"),axis.text.y = element_text(size=16,margin=margin(r=5)),axis.text.x = element_text(size=16,margin=margin(t=5)))

#Fig.E6_legend
p1<-ggplot(md,aes(x=md$Interval.Day, y=md$Similarity, fill=md$Abx))+
  geom_boxplot(alpha = 0.80)+
  scale_fill_manual(values=c("#FFBF00","gray"))+
  geom_point(size=2.0, shape=NA, position = position_jitterdodge())+
  labs(x="Interval Days", y="Bray-Curtis \n similarity",fill="Maintenance \n antibiotics")+theme(axis.title.y = element_text(angle=0,vjust=0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=2),title = element_text(size = rel(1.8)),legend.text = element_text(size = 18),legend.key=element_blank(),legend.key.size=unit(1.0,"cm"),axis.text = element_text(size=16))
library(gridExtra)
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 
legend <- g_legend(p1)
grid.draw(legend)







