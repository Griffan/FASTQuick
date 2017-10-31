library(ggplot2) 
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
findBump<-function(m,x,y,z)
{
  cnt=0
  DCMax=max(m[,y])
  pivot=which.max(m[,y])
  prev=DCMax
  for(i in pivot:1)
  {
    if(cnt==z)
    {
     break
    }
  
    if(m[i,2]>prev)
    {
      cnt=cnt+1
    }
    
    prev=m[i,2]
    min=i
  }
 
  
  cnt=0
  prev=DCMax
  print(dim(m))
  for(i in pivot:length(m[,1]))
  {
    if(cnt==z)
    { 
      break
    }
    
    if(m[i,2]>prev)
     {
      cnt=cnt+1
    }
      
      prev=m[i,2]
      max=i
  }
  
  return(list(MAX=max,MIN=min))
}
create.DenDist=function(InsertTable,col1,col2){
  start_pos=InsertTable[1,col1]
  count=InsertTable[1,col2]
  NewT=c(-1,0)
  for(i in 1:length(InsertTable[,col1]))
  {
    if(InsertTable[i,col1]<start_pos+5)
    {
      count=count+InsertTable[i,col2]
    }
    else
    {
      NewT=rbind(c(start_pos,count),NewT)
      start_pos=InsertTable[i,col1]
      count=InsertTable[i,col2]
    }
  }
  return(NewT)
}




args=commandArgs(trailingOnly = TRUE)
input=args[1]
print(args)

pdf(file=paste(input,".pdf",sep=""))
par(mfrow=c(2,2))

mydata= read.table(paste(input,".DepthDist",sep=""),header=FALSE)
mydata=mydata[1:100,]
colnames(mydata)=c("Depth","SiteCount")
lmt=findBump(mydata,1,2,3)
q1=ggplot(mydata,aes(x=Depth,y=SiteCount))+geom_line(color="#00BFC4")+
  ggtitle("Depth Distribution")+coord_cartesian(xlim=c(mydata[lmt$MIN,1],mydata[lmt$MAX,1]))+ theme(plot.title = element_text(hjust = 0.5))


mydata= read.table(paste(input,".EmpCycleDist",sep=""),header=FALSE)
mydata=mydata[1:100,]
colnames(mydata)=c("Cycle","VariantCount","BaseCount","EmpericalQuality","ReadCount")
q2=ggplot(mydata)+geom_line(aes(y=EmpericalQuality,x=Cycle),color="#00BFC4")+
  ggtitle("Sequencing Cycle V.S. Emperical Quality")+coord_cartesian(xlim=c(0,100))+ylim(0,45)+
  theme(plot.title = element_text(hjust = 0.5,size=10))


mydata=read.table(paste(input,".EmpRepDist",sep=""),header=FALSE)
mydata=mydata[1:40,]
colnames(mydata)=c("SequencingQuality","VariantCount","BaseCount","EmpericalQuality")
mydata$BaseCount=mydata$BaseCount
q3=ggplot(mydata)+geom_line(aes(y=EmpericalQuality,x=SequencingQuality),color="#00BFC4")+
  ggtitle("Sequencing Quality V.S. Emperical Quality")+geom_abline(intercept=0, slope=1,color="purple",linetype="dotted")+
  xlim(0,40)+ylim(0,40)+ theme(plot.title = element_text(hjust = 0.5,size=10))

q4=ggplot(mydata)+geom_line(aes(y=BaseCount,x=SequencingQuality),color="red",linetype="dotted")+
  ggtitle("Base Count Distribution")+xlim(0,40)+ theme(plot.title = element_text(hjust = 0.5))

mydata= read.table(paste(input,".GCDist",sep=""),header=FALSE)
mydata=mydata[2:100,]
colnames(mydata)=c("GC","dummy1","dummy2","AvgDepth")
q5=ggplot(mydata,aes(x=GC,y=AvgDepth))+geom_line(color="#00BFC4")+
  ggtitle("GC Distribution")+coord_cartesian(xlim=c(0,100))+
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
ExpectedDepth=as.numeric(strsplit(mydata[length(mydata)-12]," ")[[1]][5])
EstimatedDepth=as.numeric(strsplit(mydata[length(mydata)-11]," ")[[1]][4])
AccessibleFraction=as.numeric(strsplit(strsplit(mydata[length(mydata)-10]," ")[[1]][8],"/")[[1]][1])
EstimatedQ20Depth=as.numeric(strsplit(mydata[length(mydata)-9]," ")[[1]][7])
EstimatedQ30Depth=as.numeric(strsplit(mydata[length(mydata)-8]," ")[[1]][7])

Q20BaseFraction=as.numeric(strsplit(mydata[length(mydata)-5],":")[[1]][2])
Q30BaseFraction=as.numeric(strsplit(mydata[length(mydata)-4],":")[[1]][2])
Depth1=as.numeric(strsplit(mydata[length(mydata)-3],":")[[1]][2])
Depth2=as.numeric(strsplit(mydata[length(mydata)-2],":")[[1]][2])
Depth5=as.numeric(strsplit(mydata[length(mydata)-1],":")[[1]][2])
Depth10=as.numeric(strsplit(mydata[length(mydata)],":")[[1]][2])

plotdata=c(EstimatedQ30Depth,EstimatedQ20Depth,EstimatedDepth,ExpectedDepth)
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

multiplot(q1,q2,q3,q4, cols=2)
multiplot(q5,q7,q6,q8,cols=2)
dev.off()

