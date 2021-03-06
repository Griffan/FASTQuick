##################################################
#####
require(ggplot2)
require(scales)
require(knitr)
require(dplyr)
require(rmarkdown)
library(ggplot2)
library(scales)
library(dplyr)


#multiplot function
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
findBump<-function(m,x,cntCol,numLimits)
{
  cnt=0
  DepthCntMax=max(m[,cntCol])
  pivot=which.max(m[,cntCol])
  prev=DepthCntMax
  for(i in pivot:1)
  {
    if(cnt==numLimits)
    {
     break
    }

    if(m[i,2]>prev*1.2)#turning points
    {
      cnt=cnt+1
    }

    prev=m[i,2]
    min=i
  }

  cnt=0
  prev=DepthCntMax
  #print(dim(m))
  for(i in pivot:length(m[,1]))
  {
    if(cnt==numLimits)
    {
      break
    }

    if(m[i,2]>prev*1.2)
     {
      cnt=cnt+1
    }

      prev=m[i,2]
      max=i
  }

  return(list(MAX=max,MIN=min))
}
create.DenDist=function(InsertTable,positionCol,countCol){
  start_pos=InsertTable[1,positionCol]
  count=InsertTable[1,countCol]
  NewT=c(-1,0)
  for(i in 1:length(InsertTable[,positionCol]))
  {
    if(InsertTable[i,positionCol]<start_pos+10)
    {
      count=count+InsertTable[i,countCol]
    }
    else
    {
      NewT=rbind(c(start_pos,count),NewT)
      start_pos=InsertTable[i,positionCol]
      count=InsertTable[i,countCol]
    }
  }
  return(data.frame(NewT,row.names=NULL))
}

##################################################
#usage:
#Rscript RPlotScript.R <output_prefix> <SVDPrefix> <FASTQuickInstallDir>


print("Usage:Rscript RPlotScript.R <output_prefix> <SVDPrefix> <FASTQuickInstallDir>")

args=commandArgs(trailingOnly = TRUE)
input=args[1]
SVDPrefix=args[2]
FASTQuickInstallDir=args[3]
#print(args)

pdf(file=paste(input,".pdf",sep=""))
par(mfrow=c(2,2))

mydata= read.table(paste(input,".DepthDist",sep=""),header=FALSE)
mydata=mydata[2:150,]
colnames(mydata)=c("Depth","SiteCount")
lmt=findBump(mydata,1,2,3)
q1=ggplot(mydata,aes(x=Depth,y=SiteCount))+geom_line(color="#00BFC4")+
  ggtitle("Depth Distribution")+coord_cartesian(xlim=c(mydata[lmt$MIN,1],mydata[lmt$MAX,1]))+ theme(plot.title = element_text(hjust = 0.5))


mydata= read.table(paste(input,".EmpCycleDist",sep=""),header=FALSE)
maxCycle = 0
for (i in 1:150)
{
  maxCycle = i
  if(mydata[i,3] == 0) break
}
mydata=mydata[1:maxCycle,]
colnames(mydata)=c("SequencingCycle","VariantCount","BaseCount","EmpiricalQuality","ReadCount")
q2=ggplot(mydata)+geom_line(aes(y=EmpiricalQuality,x=SequencingCycle),color="#00BFC4")+
  ggtitle("Sequencing Cycle V.S. Empirical Quality")+coord_cartesian(xlim=c(0,maxCycle))+ylim(0,45)+
  theme(plot.title = element_text(hjust = 0.5,size=10))


mydata=read.table(paste(input,".EmpRepDist",sep=""),header=FALSE)
mydata=mydata[1:40,]
colnames(mydata)=c("SequencingQuality","VariantCount","BaseCount","EmpiricalQuality")
mydata$BaseCount=mydata$BaseCount
q3=ggplot(mydata)+geom_line(aes(y=EmpiricalQuality,x=SequencingQuality),color="#00BFC4")+
  ggtitle("Sequencing Quality V.S. Empirical Quality")+geom_abline(intercept=0, slope=1,color="purple",linetype="dotted")+
  xlim(0,40)+ylim(0,40)+ theme(plot.title = element_text(hjust = 0.5,size=10))

q4=ggplot(mydata)+geom_line(aes(y=BaseCount,x=SequencingQuality),color="red",linetype="dotted")+
  ggtitle("Base Count Distribution")+xlim(0,40)+ theme(plot.title = element_text(hjust = 0.5))

mydata= read.table(paste(input,".GCDist",sep=""),header=FALSE)
mydata=mydata[2:101,c(1,3,4)]
colnames(mydata)=c("GC","SiteCount","NormalizedMeanDepth")
total = sum(mydata$SiteCount)
#interpolate
test2=mydata %>% mutate(CumCount=cumsum(SiteCount)/total*100)
PercentileVSGC <- data.frame(
  with(test2,
       approx(CumCount, GC, xout = seq(0, 100, by = 0.05), method = "linear")
  )
)
PercentileVSGC$Depth=apply(PercentileVSGC, 1, function(x) { sum(test2[test2$GC<=x[2],]$SiteCount*test2[test2$GC<=x[2],]$NormalizedMeanDepth)/sum(test2[test2$GC<=x[2],]$SiteCount)})
colnames(PercentileVSGC)=c("GCPercentile","GCPercentage","NormalizedMeanDepth")
apf<- approxfun(test2$CumCount, test2$GC)
f <- Vectorize(function(x) {
  if (x == 0) return(1/1e10)
  apf(x)
})

q5=ggplot(PercentileVSGC,aes(x=GCPercentile,y=NormalizedMeanDepth))+geom_line(color="#00BFC4")+
  #ggtitle("Normalized Mean Depth V.S. GC Content")+
  #coord_cartesian(xlim = c(0,100), ylim = c(0,1.5), expand = TRUE)+
  ylim(0,1.5)+
  theme(plot.title = element_text(hjust = 0.5), axis.title.x =element_text(hjust = 0.5))+
  scale_x_continuous(limits=c(0,100),expand = c(0, 0), sec.axis = sec_axis(tran=~f(.), name = "GCPercentage", breaks = seq(20,60,5)))+
  geom_abline(intercept=1, slope=0,color="red",linetype="dotted")+ theme(plot.title = element_text(hjust = 0.5))


mydata= read.table(paste(input,".AdjustedInsertSizeDist",sep="",row.names=NULL),colClasses=c("numeric","numeric"),header=FALSE)
mydata=mydata[-1,]
NewTable=create.DenDist(mydata,1,2)
Adjust.Table=NewTable[order(NewTable[,1]),]
sumNewTable=sum(Adjust.Table[,2])
Adjust.Table[,2]=Adjust.Table[,2]/sumNewTable
Adjust.Table=data.frame(Adjust.Table,rep("AdjustedInsertSize",length(Adjust.Table[,1])),row.names=NULL)
colnames(Adjust.Table)=c("InsertSize","Frequency","Category")

mydata= read.table(paste(input,".RawInsertSizeDist",sep="",row.names=NULL),colClasses=c("numeric","numeric"),header=FALSE)
mydata=mydata[-1,]
NewTable=create.DenDist(mydata,1,2)
Raw.Table=NewTable[order(NewTable[,1]),]
sumNewTable=sum(Raw.Table[,2])
Raw.Table[,2]=Raw.Table[,2]/sumNewTable
Raw.Table=data.frame(Raw.Table,rep("RawInsertSize",length(Raw.Table[,1])),row.names=NULL)
colnames(Raw.Table)=c("InsertSize","Frequency","Category")

Combined.Table=rbind(Raw.Table,Adjust.Table)
Combined.Table$Category=as.factor(Combined.Table$Category)

lmt=findBump(Adjust.Table,1,2,3)
q6=ggplot(Combined.Table,aes(x=InsertSize,y=Frequency,colour=Category))+geom_line()+
  ggtitle("InsertSize Distribution")+coord_cartesian(xlim=c(ifelse(Adjust.Table[lmt$MIN,1]>100,100,Adjust.Table[lmt$MIN,1]),ifelse(Adjust.Table[lmt$MAX,1]<1000,1000,Adjust.Table[lmt$MAX,1])))+
  theme(plot.title = element_text(hjust = 0.5),legend.position=c(0.85,0.9),legend.key.size = unit(0.2, "cm"),legend.text =element_text(size=5), legend.title = element_text(size=5) )


mydata <- scan(paste(input,".Summary",sep=""), what="", sep="\n")
fileLen=length(mydata)-1#minus contamination line
ExpectedDepth=as.numeric(strsplit(strsplit(mydata[fileLen-13]," ")[[1]][5], "\\[")[[1]][1])
EstimatedDepth=as.numeric(strsplit(strsplit(mydata[fileLen-12]," ")[[1]][5], "\\[")[[1]][1])

AccessibleFraction=as.numeric(strsplit(strsplit(mydata[fileLen-11]," ")[[1]][8],"%")[[1]][1])
EstimatedQ20Depth=as.numeric(strsplit(mydata[fileLen-3],"\\s+|:")[[1]][8])
EstimatedQ30Depth=as.numeric(strsplit(mydata[fileLen-2],"\\s+|:")[[1]][8])

Q20BaseFraction=as.numeric(strsplit(mydata[fileLen-5],":")[[1]][2])
Q30BaseFraction=as.numeric(strsplit(mydata[fileLen-4],":")[[1]][2])
Depth1=as.numeric(strsplit(mydata[fileLen-9],":")[[1]][2])
Depth2=as.numeric(strsplit(mydata[fileLen-8],":")[[1]][2])
Depth5=as.numeric(strsplit(mydata[fileLen-7],":")[[1]][2])
Depth10=as.numeric(strsplit(mydata[fileLen-6],":")[[1]][2])

plotdata=c(EstimatedQ30Depth,EstimatedQ20Depth,EstimatedDepth,ExpectedDepth)
#print(plotdata)
plotname=c("EstimatedQ30Depth","EstimatedQ20Depth","EstimatedDepth","ExpectedDepth")
mydate=cbind(read.table(text=plotname),plotdata)
q7=ggplot(mydate,aes(x=mydate[,1],y=mydate[,2]))+ggtitle("Depth")+
  scale_x_discrete(limits=c("EstimatedQ30Depth","EstimatedQ20Depth","EstimatedDepth","ExpectedDepth"))+
  theme(axis.text.x = element_text(size=8,angle = 310, hjust = 0),axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5))+
  geom_bar(stat="identity",fill="#00BFC4",alpha=0.5)+ylab('Depth')#+xlim(0,length(plotdata)+1)+ylim(0,max(plotdata))

plotdata2=c(Q20BaseFraction,Q30BaseFraction,Depth1,Depth2,Depth5,Depth10)
plotname2=c("Q20","Q30","Depth 1","Depth 2","Depth 5","Depth 10")
#mydata2=cbind(read.table(text=plotname2),plotdata2)
mydata2=data.frame(plotname2,plotdata2)
q8=ggplot(mydata2,aes(x=mydata2[,1],y=mydata2[,2]))+ggtitle("Summary")+
  scale_x_discrete(limits=c("Q20","Q30","Depth 1","Depth 2","Depth 5","Depth 10"))+
  geom_bar(stat="identity",fill="#00BFC4",alpha=0.5)+ylab('Fraction')+ 
  theme(axis.text.x = element_text(angle = 310, hjust = 0.3),axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5))

#ancestry


alphaScale=scale_alpha_discrete(range = c(0.9,0.3),guide=FALSE)
sizeScale=scale_size(range=c(1.5,1),guide=FALSE)
cat("Background data points:1000g phase3\n")


#set 1000g color scale
colScale = scale_color_manual(values=c('ESN'='#FFCD00','GWD'='#FFB900','LWK'='#CC9933','MSL'='#E1B919','YRI'='#FFB933','ACB'='#FF9900','ASW'='#FF6600',
                                       'CLM'='#CC3333','MXL'='#E10033','PEL'='#FF0000','PUR'='#CC3300','CDX'='#339900','CHB'='#ADCD00','CHS'='#00FF00',
                                       'JPT'='#008B00','KHV'='#00CC33','CEU'='#0000FF','FIN'='#00C5CD','GBR'='#00EBFF','IBS'='#6495ED','TSI'='#00008B',
                                       'BEB'='#8B008B','GIH'='#9400D3','ITU'='#B03060','PJL'='#E11289','STU'='#FF00FF','AFR'='#FFCD33','AFR/AMR'='#FF9900',
                                       'AMR'='#FF3D3D','EAS'='#ADFF33','EUR'='#64EBFF','SAS'='#FF30FF','UserSample'='#000000'),
                              breaks=c('ESN','GWD','LWK','MSL','YRI','ACB','ASW','CLM','MXL','PEL','PUR','CDX','CHB','CHS','JPT','KHV','CEU','FIN','GBR',
                                       'IBS','TSI','BEB','GIH','ITU','PJL','STU','AFR','AFR/AMR','AMR','EAS','EUR','SAS','UserSample'))
#set 1000g coordinates
POP=read.table(paste0(FASTQuickInstallDir,"/resource/1000g.pop"),header = F)
RefCoord.1kg=read.table(paste0(SVDPrefix,".V"),header = F)

PCdim=dim(RefCoord.1kg)[2] - 1

if(PCdim >= 4) 
{
  print("plot using 4 PCs")
  RefCoord.1kg=RefCoord.1kg[,1:5]
  RefCoord.1kg['POP'] <- POP$V2[match(RefCoord.1kg$V1, POP$V1)]
  RefCoord.1kg['Size'] <- rep(0.5,length(RefCoord.1kg$POP))
  colnames(RefCoord.1kg)=c("ID","PC1","PC2","PC3","PC4","POP","DotSize")
}else if(PCdim >=2)
{
  print("plot using 2 PCs")
  RefCoord.1kg=RefCoord.1kg[,1:3]
  RefCoord.1kg['POP'] <- POP$V2[match(RefCoord.1kg$V1, POP$V1)]
  RefCoord.1kg['Size'] <- rep(0.5,length(RefCoord.1kg$POP))
  colnames(RefCoord.1kg)=c("ID","PC1","PC2","POP","DotSize")
}
RefCoord=RefCoord.1kg


#assuming the input target sample has format ID PC1 PC2 PC3 PC4
TargetSample=read.table(file=paste0(input,".Ancestry"),header=T)
TargetSample=data.frame(ID=c("IntendedSample"),PC1=c(TargetSample[1,3]), PC2=c(TargetSample[2,3]), PC3=c(TargetSample[3,3]), PC4=c(TargetSample[4,3]), POP=c("UserSample"), DotSize=c(1))
colnames(TargetSample)=c("ID","PC1","PC2","PC3","PC4","POP","DotSize")
#plot RefCoord as grey
#CombinedData=TargetSample
#q9=ggplot(data=RefCoord,aes(PC1,PC2))+geom_point(color="grey")+geom_point(data=CombinedData,aes(PC1,PC2,color=POP))+
#  colScale+alphaScale+sizeScale
#prepare dataset for plot
if(PCdim >= 4) 
{
  CombinedData=rbind(RefCoord,TargetSample)
}else if(PCdim >=2)
{
  CombinedData=rbind(RefCoord,TargetSample[,c(1,2,3,6,7)])
}
#plot RefCoord as colorful
q10=ggplot()+geom_point(data=CombinedData,aes(PC1,PC2,color=POP, size=DotSize), alpha=0.5)+#geom_text(data=CombinedData,aes(PC1,PC2,label=POP),size=1)+
  colScale+alphaScale+scale_size(guide = 'none')#+sizeScale

if(PCdim >= 4) 
{
q11=ggplot()+geom_point(data=CombinedData,aes(PC3,PC4,color=POP, size=DotSize), alpha=0.5)+#geom_text(data=CombinedData,aes(PC1,PC2,label=POP),size=1)+
  colScale+alphaScale+scale_size(guide = 'none')#+sizeScale
}

#output file
multiplot(q1,q2,q3,q4, cols=2)
multiplot(q5,q7,q6,q8,cols=2)
multiplot(q10,cols=1)
if(PCdim >= 4) 
{
multiplot(q11,cols=1)
}
dev.off()



