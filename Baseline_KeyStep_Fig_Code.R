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

###All code of the key step for plot and analysis, with attched data set. For repeated plot/analysis, only chose one as example.
##Fig.1
##caclulate the average of day to all other days similarity in each baselines
dts1<-read.csv("stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample_1.shared",header = T,sep = "\t")
AVE<-c()
dts<-dts1[c(237:257),]# use Baseline 1 as example
rS<-dts[,6:ncol(dts)]
dS<- as.matrix(vegdist(rS,method="bray"))
dS<-1-dS
ave<-(rowSums(dS)-1)/(ncol(dS)-1)
AVE<-c(AVE,ave)
#output file 
fn=paste("Averaged_SimilarityToOthers_Baseline1",".csv",sep="")
write.csv(AVE,fn,row.names=TRUE)
#
md<-read.csv("Averaged_SimilarityToOthers.csv",header = T,sep = ",")
md$Subject <- factor(md$Subject,levels = c("4","5","9","15","17","23","Replicate","GD"),labels = c("1","2","3","4","5","6","Repeat samples","Generous donor"))
md$Baseline <- factor(md$Baseline,levels = c("4A","5A","5B","5C","5D","9A","15A","17A","17B","23A","Replicate","GD"),labels = c("1","2a","2b","2c","2d","3","4","5a","5b","6","TR1","TR2"))
myColors <-c("gold1","dodgerblue","purple","red","#C39953","#568203","gray50","black","dark red")
names(myColors) <- levels(md$Subject)
colScale <- scale_colour_manual(name = "Subject",values = myColors)
#plot
ggplot(md,aes(x=md$Baseline,y=md$Similarity.to.others))+geom_rect(xmin = 11-0.5,xmax = 12+0.5,ymin=-Inf,ymax=Inf,fill="gray88",alpha=0.2)+
  geom_boxplot(outlier.size = 0,outlier.shape=NA)+colScale+
  geom_jitter(aes(color=as.factor(md$Subject)),shape=16,size=1.6,position=position_jitter(0.3),alpha=0.65)+
  labs(x="                               Baseline periods                                        Controls       ", y="Bray-Curtis \n similarity")+
  theme(axis.title.y = element_text(angle=0,vjust=0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),legend.position ="none",title = element_text(size = rel(1.8)),legend.text = element_text(size = 18),legend.key=element_blank(),legend.key.size=unit(1.0,"cm"),
        axis.text.x= element_text(size=20,margin=margin(t=5)),axis.text.y = element_text(size=20,margin=margin(r=5)),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))

#Linear mixed model based on Fig.1 
md<-read.csv("Averaged_SimilarityToOthers.csv",header = T,sep = ",")
fit<-lmer(md$Similarity.to.others~md$Sample_type+(1|md$Subject),data=md)
summary(fit)
confint(fit)#confidence interval
#Model for the inter-and intra-subjects difference
md<-read.csv("AllSubject_InterIntraSubject.csv",header=T,sep = ",")
fit<-lm(md$Similarity~Subject,data=md) #same or different subject group 
summary(fit)
confint(fit)
#Model for the association of similarity with baseline period length
fit<-lmer(md$Similarity~md$Length+(1|md$Subject),data=md)
summary(fit)
confint(fit)

##Fig.2
md<-read.csv("All_TBL.csv",header = T,sep = ",")
md$Subject <- factor(md$Subject,levels = c("4","5","9","15","17","23"),labels = c("1","2","3","4","5","6"))
myColors <-c("gold1","dodgerblue","purple","red","#C39953","#568203","gray50","black","dark red")
names(myColors) <- levels(md$Subject)
colScale <- scale_colour_manual(name = "Subject",values = myColors)
md$Sub_Baseline <- factor(md$Sub_Baseline,levels = c("4A","5A","5B","5C","5D","9A","15A","17A","17B","23A"),labels = c("1","2a","2b","2c","2d","3","4","5a","5b","6"))
#plot
ggplot(md,aes(x=md$Sub_Baseline,y=md$LOG))+geom_boxplot(outlier.size = 0,outlier.shape=NA)+colScale+geom_jitter(aes(color=as.factor(md$Subject)),shape=16,size=1.6,position=position_jitter(0.3),alpha=0.65)+labs(x="Baseline periods", y="Log(Copies/mL)")+theme(axis.title.y = element_text(angle=0,vjust=0.5))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),legend.position ="none",title = element_text(size = rel(1.8)),legend.text = element_text(size = 18),legend.key=element_blank(),legend.key.size=unit(1.0,"cm"),
                                                                                                                                                                                                                                                                                                                        axis.text.x= element_text(size=20,margin=margin(t=5)),axis.text.y = element_text(size=20,margin=margin(r=5)), axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))
##Fig.3
dts1<-read.csv("stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample_1_1.shared",header = T,sep = "\t")
myColors <-c("dodgerblue","gold1","purple","red","#C39953","#568203")
dts1$Abx<-factor(dts1$Abx,levels = c("None","InhaledColi","OralAzi","OralAzi+InhaledAztre","OralAzi+InhaledTobra"),labels = c("None","Inhaled colistin","Azithromycin","Azithromycin + inhaled aztreonam","Azithromycin + inhaled tobramycin"))
names(myColors) <- levels(dts1$Abx)
colScale <- scale_colour_manual(name = "Maintenance antibiotics",values = myColors)
dts<-dts1[c(258:315),] #Use baseline 2a as example
sdts<-dts[,c(6:ncol(dts))]
rs<-sdts
set.seed(6)
ord<-metaMDS(rs,distance = "bray")
ord #stress value
Dismatrix<-ord$diss
nmdsdf<-data.frame(axis1 = ord$points[,1], axis2 = ord$points[,2],group=dts$Abx)
#plot
gg<-data.frame(merge(nmdsdf,aggregate(cbind(mean.x=axis1 ,mean.y=axis2)~group,nmdsdf,mean),by="group"))
ggplot(data=gg, aes(gg$axis1,gg$axis2,color=gg$group))+geom_point(size=3.0)+xlim(-1.3,0.8)+ylim(-0.8,0.8)+geom_segment(aes(x=gg$mean.x, y=gg$mean.y,xend=gg$axis1,yend=gg$axis2))+geom_point(aes(x=gg$mean.x,y=gg$mean.y),size=1.5,stroke=2,color="Black")+labs(x="nMDS1", y="nMDS2",color="Baseline")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),title = element_text(size = rel(1.8)),legend.position ="none",legend.text = element_text(size = 18),legend.key=element_blank(),legend.key.size=unit(1.0,"cm"),axis.text = element_text(size=16))+colScale                                                


##Fig.4
dts<-read.csv("P5_4Baseline_Replicates.shared",header = T,sep = "\t") #used subject 2 as example
sdts<-dts[,c(7:ncol(dts))]
rs<-sdts
set.seed(5)
ord<-metaMDS(rs,distance = "bray")
ord 
Dismatrix<-ord$diss
dts$Baseline <- factor(dts$Baseline,levels = c("A","B","C","D"),labels = c("2a","2b","2c","2d")) #Subject5
nmdsdf<-data.frame(axis1 = ord$points[,1], axis2 = ord$points[,2],group=dts$Baseline,pair=dts$Pair)
myColors <-c("gold1","dodgerblue","red","purple")
names(myColors) <- levels(dts$Baseline)
colScale <- scale_colour_manual(name = "Baseline",values = myColors)
nmdsdf1<-nmdsdf[c(1:228),]
#plot
gg1<-data.frame(merge(nmdsdf1,aggregate(cbind(mean.x=axis1 ,mean.y=axis2)~group,nmdsdf1,mean),by="group"))
gg1$group1<-factor(gg1$group,levels = c("2a","2b","2c","2d"),labels = c("a","b","c","d"))
ggplot(data=gg1, aes(gg1$axis1,gg1$axis2,color=gg1$group))+geom_point(size=2.5,alpha=0.6)+xlim(-1.1,0.8)+ylim(-0.7,0.8)+geom_segment(aes(x=gg1$mean.x, y=gg1$mean.y,xend=gg1$axis1,yend=gg1$axis2),alpha=0.6)+geom_point(aes(x=gg1$mean.x,y=gg1$mean.y),size=7.3,shape=15,color="white")+geom_text(aes(x=gg1$mean.x,y=gg1$mean.y),label=gg1$group1,size=8,color="Black")+labs(x="nMDS1", y="nMDS2",color="Baseline")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),title = element_text(size = rel(1.8)),legend.key=element_blank(),legend.text = element_text(size = 20),legend.key.size=unit(1.0,"cm"),axis.text = element_text(size=16))+colScale                                              

##Fig.5
#Caculate the BC similarity of all possible interval days within each baseline 
dts1<-read.csv("P17B.shared",header=T,sep = "\t") #use baseline 5a as example
rS<-dts1[,6:ncol(dts1)]
dS<- as.matrix(vegdist(rS,method="bray"))

M<-NULL
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
md<-md[!(md$value==0|md$value==1|md$value=="NA"),]#remove rows with value of 0,1,NA (missing days)
md<-md[complete.cases(md), ]
colnames(md) <- c("Interval Day", "Start Day","Dissimilarity")
md$Similarity<-1-md$Dissimilarity
fn=paste("P17B_AllIntervalDays",".csv",sep="")
write.csv(md,fn,row.names=TRUE)

#Caculate the antibiotics (same,different) of all possible interval days within each baseline
dS<-read.csv("P17B_abx.shared",header=F,sep = "\t")#use baseline 5a as example
M<-NULL
for (i in 1:nrow(dS)){
  da<-c()
  n=1
  b=nrow(dS)-i
  while(n<=b){
    m<-dS[n,n+i]
    da<-c(da,m)
    n<-n+1
  }
  M<-cbind(M,da)
}
for (i in 1: nrow(M)){
  for (j in 1: ncol(M)){
    if (i>nrow(dS)-j){
      M[i,j]=0
    }
  }
}
TransM<-t(M)
rownames(TransM)<-c(1:ncol(TransM))
md<-melt(TransM,id.vars=rownames(TransM))
md<-md[!(md$value==0|md$value==1|md$value=="NA"),] #remove rows with value of 0,1,NA (missing days)
md<-md[complete.cases(md), ]
colnames(md) <- c("Interval Day", "Start Day","Abx_span")
fn=paste("P17B_Abx_IntervalDays",".csv",sep="")
write.csv(md,fn,row.names=TRUE)
#merge the interval similarity and antibiotics to file, filter only 1,3,7,14,30 interval days for Fig.
md<-read.csv("All_IntervalDays_baseline_Abxspan.csv",header = T,sep = ",") #only 1,3,7,14,30 days included
myColors <-c("#dba404","gray45")
names(myColors) <- levels(md$Abx)
colScale <- scale_colour_manual(name = "Maintenance \n  antibiotics",values = myColors)
md$Interval.Day<-as.factor(md$Interval.Day)
ggplot(md,aes(x=md$Interval.Day,y=md$Similarity))+geom_boxplot(outlier.size = 0,outlier.shape=NA,lwd=1)+ylim(0.1,1)+
  colScale+geom_jitter(aes(color=as.factor(md$Abx)),shape=16,size=1.2,position=position_jitter(0.31),alpha=0.75)+labs(x="Interval days", y="Bray-Curtis \n similarity")+theme(axis.title.y = element_text(angle=0,vjust=0.5),legend.position ="none")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),title = element_text(size = rel(1.8)),legend.key=element_blank(),legend.text = element_text(size = 18),legend.key.size=unit(1.0,"cm"),axis.text.y = element_text(size=16,margin=margin(r=5)),axis.text.x = element_text(size=16,margin=margin(t=5)))

##Model based on Fig.5
md<-read.csv("All_AllIntervalDays_Abx_shortterm.csv",header = T,sep = ",")
fit<-lmer(md$Similarity~md$Interval.Day+md$Abx+(1|md$Subject)+(0+md$Interval.Day|md$Subject)+(1|md$Baseline),data=md)
summary(fit)
confint(fit)

##Suppl Fig.1_ reagent control 
dds<-read.table("SampleControl_All.csv",header=T,sep = ",")
dds<-cbind(dds[,c(1:3)], dds[,c(4:6)]/1068)
colnames(dds)[4] <- "Samples"
colnames(dds)[5] <- "Reagent Control"
colnames(dds)[6] <- "Water"
dd<- subset(dds, dds$`Reagent Control`>0.01) # use reagent controlas example, set cut off >0.01 
dd$Rank[order((dd[,5]))] <- 1:nrow(dd) #ranked by reagent control
dd<-dd[,c(1,2,3,4,5,6)]
ddt<-melt(dd,id.vars = c("Group","Taxa","Rank"))
colnames(ddt)[4] <- "Type"
colnames(ddt)[5] <- "Relative_abundance"
ddt_sort<- ddt[order(ddt$Rank), ]
myorder <- ddt_sort %>% filter(Type=="Reagent Control") %>% arrange(desc(Relative_abundance)) %>% .$Group %>% as.character
mytaxa<-ddt_sort %>% filter(Type=="Reagent Control") %>% arrange(desc(Relative_abundance)) %>% .$Taxa %>% as.character
ddt_sort<- ddt[order(ddt$Rank), ]
ggplot(ddt_sort, aes(x=Group, y=Relative_abundance, fill = Type)) + labs(x="Taxa", y="Relative \n abundance")+
  geom_bar(stat="identity", position = "dodge") + scale_fill_manual(values=c("#FFBF00","Black","grey"))+
  scale_x_discrete(limits=myorder,labels=mytaxa)+ylim(0,0.2)+
  theme(title = element_text(size = rel(1.5)),axis.title.y = element_text(angle=0,vjust=0.5,size=18,margin = margin(r=15)),axis.title.x = element_text(size=18,margin = margin(t=15)),legend.position ="none",legend.text=element_text(size = 17),legend.key.size=unit(0.7,"cm"),axis.text.y = element_text(size=16,margin=margin(r=5)),axis.text.x = element_text(angle=25,size=16,hjust=1,vjust=1,margin=margin(t=5)),panel.background = element_blank(),axis.line = element_line(colour = "black", size=1))

##Suppl Fig.2
col<-c("grey","#A52A2A","#C39953","#DBE9F4","#89CFF0","#FFAA1D","#CC0000","#006A4E","#873260","#536872","#B5A642","#CB4154","#1DACD6","#66FF00","#BF94E4","#F4C2C2","#009E73","#CD7F32","#737000","#007FFF","#E7FEFF","#7BB661","#F0DC82","#CC5500","#FFC1CC","#5F9EA0","#ED872D","#702963","#FFF600","#A67B5B","#4B3621","#A3C1AD","#A9A9A9","#3B7A57","#FFBF00","#A52A2A","#B284BE","#0048BA")
dts1<-read.table("stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample_2_taxa.shared",header=T,sep="\t")
dts<-dts1[c(258:315),]#use baseline 2a as example
minsub1<-cbind(dts[,2],dts[,6:32])
colnames(minsub1)[1] <- "Group"
minsub1<-melt(minsub1,id.vars = "Group")
time.x <- as.numeric(minsub1[,1])
begin <- time.x[1]
time.x <- as.numeric(minsub1$Group)-begin+1
colnames(minsub1)[2] <- "OTUs"
dts2<-read.csv("Averaged_SimilarityToOthers.csv",header = T,sep = ",")#averaged to all other similarity
rS<-dts2[c(258:315),]# baseline 2a
time.x1 <- as.numeric(rS[,1])
begin <- time.x1[1]
time.x1<- as.numeric(rS$Sample)-begin+1
dtd<-data.frame(cbind(time.x1,rS$Similarity.to.others))
#edit the dtd file for dummy time points of missing days
dtd<-read.csv("dtd_2a.csv",header = T,sep = ",")
gg<-ggplot() + geom_bar(data = minsub1,aes(x=time.x,y=value, fill=OTUs), stat = "identity") + 
  scale_fill_manual(values = col)+
  labs(x="Daily samples",y="Relative abundance")+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)),title = element_text(size = rel(1.8)),legend.text=element_text(size = 18),legend.key.size=unit(0.6,"cm"),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18),legend.position ="none",panel.background = element_blank(),axis.line = element_line(colour = "gray"))+scale_x_continuous("Daily samples",breaks=c(0,20,40,60,80,100,120),expand = c(0.001,0.001))+scale_y_continuous(expand = c(0,0))+
  geom_point(data =dtd, aes(time.x1, dtd$V2),size=2.5)
for (i in 1:7) gg <- gg + 
  geom_line(data=subset(dtd, Segment==i), 
            aes(x=time.x1, y=V2,linetype=as.character(linet), group=i))
gg

##Suppl Fig.3
md<-read.csv("Averaged_SimilarityToOthers_outlier.csv",header = T,sep = ",")
md<-md[c(1:527),]
#plot for total bacterial load 
ggplot(md,aes(x=md$Outlier,y=md$TBL_LOG))+ylim(7,10.2)+geom_boxplot(outlier.size = 0,outlier.shape=NA)+colScale+geom_jitter(shape=16,size=1.6,position=position_jitter(0.3),alpha=0.65)+labs(x=" ", y="Log(Copies/mL)")+theme(axis.title.y = element_text(angle=0,vjust=0.5))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),legend.position ="none",title = element_text(size = rel(1.8)),legend.text = element_text(size = 18),legend.key=element_blank(),legend.key.size=unit(1.0,"cm"), axis.text.x= element_text(size=20,margin=margin(t=5)),
 axis.text.y = element_text(size=20,margin=margin(r=5)), axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))
#model
fit<-lmer(md$TBL_LOG~md$Outlier+(1|md$Subject)+(1|md$Baseline),data=md)
summary(fit)
confint(fit)

##Suppl Fig.4
dts1<-read.csv("stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample_1_taxa2.shared",header = T,sep = "\t")
dts<-dts1[c(258:485),]#Subject 2 as example
myshapes <-c(16,1,8)
names(myshapes) <- levels(dts$Outlier)
shapeScale <- scale_shape_manual(name = "Type",values = myshapes)
rs<-dts[,c(6:ncol(dts))]
set.seed(6)
ord<-metaMDS(rs,distance = "bray")
ord
Dismatrix<-ord$diss
nmdsdf<-data.frame(axis1 = ord$points[,1], axis2 = ord$points[,2],Outlier=dts$Outlier)
vf <- envfit(ord,rs, perm =100)
spp.scrs <- as.data.frame(scores(vf, display = "vectors")[vf$vectors$pval <=0.01&vf$vectors$r>0.35, ])
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
# Biplot
ggplot(data=nmdsdf, aes(nmdsdf$axis1,nmdsdf$axis2))+geom_point(aes(shape=Outlier),size=2.8)+labs(x="nMDS1", y="nMDS2",color="Subject")+
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), colour = "grey38") +
  geom_text_repel(data = spp.scrs, aes(x = NMDS1, y = NMDS2, label = Species),size =4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),legend.position ="none",axis.line = element_line(color = 'black',size=1.2),title = element_text(size = rel(1.8)),legend.text = element_text(size = 18),legend.key=element_blank(),legend.key.size=unit(1.0,"cm"),axis.text = element_text(size=16))+shapeScale

##Suppl Fig.5
md<-read.csv("Averaged_SimilarityToOthers.csv",header = T,sep = ",")
md<-md[c(1:527),]
ggplot(md,aes(x=md$Outlier_defined,y=md$Days.at.4))+geom_boxplot(outlier.size = 0,outlier.shape=NA)+colScale+geom_jitter(shape=16,size=1.6,position=position_jitter(0.3),alpha=0.65)+labs(x=" ", y="Sample refrigerated days")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),legend.position ="none",title = element_text(size = rel(1.8)),legend.text = element_text(size = 18),legend.key=element_blank(),legend.key.size=unit(1.0,"cm"),
axis.text.x= element_text(size=20,margin=margin(t=5)),axis.text.y = element_text(size=20,margin=margin(r=5)), axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))
#Model 
fit<-lmer(md$Days.at.4~md$Outlier_defined+(1|md$Subject)+(1|md$Baseline),data=md)
summary(fit)
confint(fit)

##Suppl Fig.6
#Function that returns two lists (N and Power) for plotting. Alpha, p0, p1 and sample size may be adjusted.
binom.exact <- function(alpha=0.2, p0=0.5, p1=0.95,start=3,stop=25) {
  N <- start:stop
  #Find critical values for specified alpha, N and p0
  critical.values <- qbinom(p = 1 - alpha, size = N, prob = p0)
  #Calculate power using critical values from previous step and p1
  Power <- 1 - pbinom(critical.values, N, p1)    
  invisible(list(N=N, Power=Power))
}
#Retrive power vs N data for three different p1 values
o85 <- binom.exact(0.2,0.5,0.85,2,13)
o9 <- binom.exact(0.2,0.5,0.9,2,13)
o95 <- binom.exact(0.2,0.5,0.95,2,13)
#Initial plot to get dimensions for legend
plot(o85$N,o85$Power,type="b",ylim=c(0,1),xlim=c(2,13),col="blue",ylab="Power",xlab="Sample Size",pch=1,xaxt="n",main=expression(paste(alpha,"=0.2")))
pos <- legend("bottomright", legend=c("0.95", "0.90", "0.85"),
              col=c("red", "green","blue"), pch = 3:1, lty=1, cex=1,bg="white",title="Population proportion of \nnon-outlier samples")
xleft <- pos$rect[["left"]]
ytop <- pos$rect[["top"]]
ybottom <- ytop - pos$rect[["h"]]
xright <- xleft + pos$rect[["w"]]
dev.off()

#Plotting power vs N as a line with points
pdf("alpha02.pdf")
plot(o85$N,o85$Power,type="b",ylim=c(0,1),xlim=c(2,13),col="blue",ylab="Power",xlab="Sample Size",pch=1,xaxt="n",main=expression(paste(alpha,"=0.2")))
lines(o9$N,o9$Power,type="b",pch=2,col="green")
lines(o95$N,o95$Power,type="b",pch=3,col="red")
axis(1, at = seq(1, 13, by = 1), las=1)
grid(NA,NULL)
abline(v=o85$N,col='gray',lty=3)
rect(xleft, ybottom-.1, xright, ytop+0.05,col="white")
legend("bottomright",bty="n", legend=c("0.95", "0.90", "0.85"),
       col=c("red", "green","blue"), pch = 3:1, lty=1, cex=1,title="Population proportion of \nnon-outlier samples")
dev.off()

##Suppl Fig.7
md<-read.csv("All_TBL.csv",header = T,sep = ",")
myColors <-c("dodgerblue","gold1","purple","red","#C39953","#568203")
md$Abx<-factor(md$Abx,levels = c("None","InhaledColi","OralAzi","OralAzi+InhaledAztre","OralAzi+InhaledTobra"),labels = c("None","Inhaled colistin","Azithromycin","Azithromycin + inhaled aztreonam","Azithromycin + inhaled tobramycin"))
names(myColors) <- levels(md$Abx)
colScale <- scale_colour_manual(name = "Maintenance antibiotics",values = myColors)
md1<-md[c(22:79),]#Baseline 2a as example
ggplot(md1,aes(x=md1$Abx,y=md1$LOG))+geom_boxplot(outlier.size = 0,outlier.shape=NA)+ylim(6.5,10.6)+
  colScale+geom_jitter(aes(color=as.factor(md1$Abx)),shape=16,size=2,position=position_jitter(0.3),alpha=0.65)+
  labs(x="", y="Log(Copies/mL)")+theme(axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)),panel.grid.major = element_blank(),legend.position ="none",axis.text.x = element_text(angle=25,size=18,hjust=1,vjust=1),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),title = element_text(size = rel(1.8)),legend.text = element_text(size = 18),legend.key=element_blank(),legend.key.size=unit(1.0,"cm"),axis.text = element_text(size=16))
#model 
md<-read.csv("All_TBL.csv",header = T,sep = ",")
md1<-md[md$Sub_Baseline == '5A',] #used baseline 2a as example
m1 = gls(LOG ~Abx+Day, correlation = corAR1(form=~Day),data=md1)
summary(m1)


##Suppl Fig.8
dts<-read.csv("P5_4Baseline_Replicates.shared",header = T,sep = "\t") #use subject 2 as example
sdts<-dts[,c(7:ncol(dts))]
rs<-sdts
set.seed(5)
ord<-metaMDS(rs,distance = "bray")
ord 
Dismatrix<-ord$diss
dts$Baseline <- factor(dts$Baseline,levels = c("A","B","C","D"),labels = c("2a","2b","2c","2d")) #Subject5
nmdsdf<-data.frame(axis1 = ord$points[,1], axis2 = ord$points[,2],group=dts$Baseline,pair=dts$Pair)
myColors <-c("gold1","dodgerblue","red","purple")
names(myColors) <- levels(dts$Baseline)
colScale <- scale_colour_manual(name = "Baseline",values = myColors)
nmdsdf1<-nmdsdf[c(1:228),]
gg1<-data.frame(merge(nmdsdf1,aggregate(cbind(mean.x=axis1 ,mean.y=axis2)~group,nmdsdf1,mean),by="group"))
fn=paste("gg_plot_file1_output",".csv",sep="")
write.csv(gg1,fn,row.names=TRUE)
fn=paste("nmdsdf_file1",".csv",sep="")
write.csv(nmdsdf,fn,row.names=TRUE)
#merge two output file and edit the data format
gg<-read.csv("gg_plot_file1.csv",header = T,sep = ",")
myshapes <-c(1,16)
names(myshapes) <- levels(gg$Type)
shapeScale <- scale_shape_manual(name = "Type",values = myshapes)
gg1<-gg[c(1:228),]
gg2<-gg[c(229:243),]
ggplot()+geom_point(data=gg1, aes(gg1$axis1,gg1$axis2,color=gg1$group,shape=gg1$Type),size=3,color="grey",alpha=0.2)+xlim(-1.1,0.8)+ylim(-0.7,0.8)+
  geom_segment(aes(x=gg1$mean.x, y=gg1$mean.y,xend=gg1$axis1,yend=gg1$axis2),alpha=0.2)+
  geom_point(data=gg2, aes(gg2$axis1,gg2$axis2,color=gg2$group,shape=gg2$Type),size=3)+
  geom_point(data=gg2, aes(gg2$r1,gg2$r2,color=gg2$group),size=3)+
  geom_segment(aes(x=gg2$r1, y=gg2$r2,xend=gg2$axis1,yend=gg2$axis2),color="black")+shapeScale+
  guides(color = guide_legend(order=1),shape = guide_legend(order=2))+
  labs(x="nMDS1", y="nMDS2",color="Baseline")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),title = element_text(size = rel(1.8)),legend.key=element_blank(),legend.text = element_text(size = 18),legend.key.size=unit(1.0,"cm"),axis.text = element_text(size=16))+colScale                                              

##Suppl Fig.9
#Same Abx
md<-read.csv("P4_IntervalDays.csv",header = T,sep = ",")
md$Interval.Day<-as.factor(md$Interval.Day)
ggplot(md,aes(x=md$Interval.Day, y=md$Similarity, fill=md$Abx_Change))+
  geom_boxplot(alpha = 0.80)+
  scale_fill_manual(values=c("gray"))+
  geom_point(size=1.2, shape=NA, position = position_jitterdodge())+ylim(0.1,1)+
  labs(x="Interval Days", y="Bray-Curtis \n similarity",fill="Maintenance \n antibiotics")+theme(axis.title.y = element_text(angle=0,vjust=0.5))+
  theme(panel.grid.major = element_blank(), legend.position ="none",panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),title = element_text(size = rel(1.8)),legend.text = element_text(size = 18),legend.key=element_blank(),legend.key.size=unit(1.0,"cm"),axis.text.y = element_text(size=16,margin=margin(r=5)),axis.text.x = element_text(size=16,margin=margin(t=5)))

#Different Abx
md<-read.csv("P17_IntervalDays_Abxspan.csv",header = T,sep = ",")
md$Interval.Day<-as.factor(md$Interval.Day)
ggplot(md,aes(x=md$Interval.Day, y=md$Similarity, fill=md$Abx_Change))+
  geom_boxplot(alpha = 0.80)+
  scale_fill_manual(values=c("#FFBF00","gray"))+
  geom_point(size=1.2, shape=NA, position = position_jitterdodge())+ylim(0.1,1)+
  labs(x="Interval Days", y="Bray-Curtis \n similarity",fill="Maintenance \n antibiotics")+theme(axis.title.y = element_text(angle=0,vjust=0.5))+
  theme(panel.grid.major = element_blank(), legend.position ="none",panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(color = 'black',size=1.2),title = element_text(size = rel(1.8)),legend.text = element_text(size = 18),legend.key=element_blank(),legend.key.size=unit(1.0,"cm"),axis.text.y = element_text(size=16,margin=margin(r=5)),axis.text.x = element_text(size=16,margin=margin(t=5)))


