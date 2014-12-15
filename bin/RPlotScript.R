library(ggplot2) 

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
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
args=commandArgs(trailingOnly = TRUE);
input=args[1]
input="C:\\Users\\PC\\Dropbox\\workingspace\\FastqA\\FastPopCon\\FastqA\\bin\\HG00553"
print(args)

pdf(file=paste(input,".pdf",sep=""))
par(mfrow=c(2,2));

mydata= read.table(paste(input,".DepthDist",sep=""),header=FALSE);
Depth=mydata[,1];
DepCount=mydata[,2];
DCMax=max(mydata[,2]);
q1=ggplot(mydata,aes(x=Depth,y=DepCount))+geom_line()+ggtitle("DepthDist")+xlim(0,100)+ylim(0,DCMax);


mydata= read.table(paste(input,".EmpCycleDist",sep=""),header=FALSE);
mydata=mydata[1:100,]
EmpericalQuality1=mydata[,4];
Cycle=mydata[,1];
#ReadCount=mydata[,5];
#CyMax=max(mydata[,1]);
q2=ggplot(mydata)+geom_line(aes(y=EmpericalQuality1,x=Cycle))+ggtitle("EmpCycleDist")+xlim(0,100)+ylim(0,45)+ylab('EmpericalQuality')


mydata= read.table(paste(input,".EmpRepDist",sep=""),header=FALSE);
mydata=mydata[1:40,];
EmpericalQuality2=mydata[,4];
ReportQuality=mydata[,1];
#Scale=max(mydata[,3])/40;
BaseCount=mydata[,3]/1000000;
q3=ggplot(mydata)+geom_line(aes(y=EmpericalQuality2,x=ReportQuality))+geom_line(y=BaseCount,x=ReportQuality,color="red",linetype="dotted")+ggtitle("EmpRepDist")+geom_abline(intercept=0, slope=1,color="purple")+xlim(0,40)+ylim(0,40)+ylab('EmpericalQuality');


mydata= read.table(paste(input,".GCDist",sep=""),header=FALSE);
GC=mydata[,1];
AvgDepth=mydata[,4];
#MaxDepth=max(AvgDepth);
q4=ggplot(mydata,aes(x=GC,y=AvgDepth))+geom_line()+ggtitle("GCDist")+xlim(0,100)+ylim(0,max(AvgDepth))+geom_abline(intercept=1, slope=0,color="red",linetype="dotted");


mydata= read.table(paste(input,".InsertSizeDist",sep=""),header=FALSE);
mydata=mydata[mydata[,2]<(mean(mydata[,2])+3*sd(mydata[,2])),]
InsertSize=mydata[,1];
IS_Count=mydata[,2];
q5=ggplot(mydata,aes(x=InsertSize,y=IS_Count))+geom_line()+ggtitle("InsertSizeDist")+xlim(0,1000)+ylim(0,max(IS_Count[1:length(IS_Count)]));


mydata <- scan(paste(input,".summary",sep=""), what="", sep="\n")
ExpectedDepth=as.numeric(strsplit(mydata[14]," ")[[1]][5])
EstimatedDepth=as.numeric(strsplit(mydata[14]," ")[[1]][5])
AccessibleFraction=as.numeric(strsplit(strsplit(mydata[16]," ")[[1]][8],"/")[[1]][1]);
EstimatedQ20Depth=as.numeric(strsplit(mydata[17]," ")[[1]][7]);
EstimatedQ30Depth=as.numeric(strsplit(mydata[18]," ")[[1]][7]);

Q20BaseFraction=as.numeric(strsplit(mydata[21],":")[[1]][2]);
Q30BaseFraction=as.numeric(strsplit(mydata[22],":")[[1]][2]);
Depth1=as.numeric(strsplit(mydata[23],":")[[1]][2]);
Depth2=as.numeric(strsplit(mydata[24],":")[[1]][2]);
Depth5=as.numeric(strsplit(mydata[25],":")[[1]][2]);
Depth10=as.numeric(strsplit(mydata[26],":")[[1]][2]);

plotdata=c(ExpectedDepth,EstimatedDepth,EstimatedQ20Depth,EstimatedQ30Depth)
plotname=c("ExpectedDepth","EstimatedDepth","EstimatedQ20Depth","EstimatedQ30Depth")
mydate=cbind(read.table(text=plotname),plotdata)
q6=ggplot(mydate,aes(x=mydate[,1],y=mydate[,2]))+ggtitle("Depth")+scale_x_discrete(limits=c("ExpectedDepth","EstimatedDepth","EstimatedQ20Depth","EstimatedQ30Depth"))+geom_bar(stat="identity",fill="grey")+xlab('Stat')+ylab('Depth')#+xlim(0,length(plotdata)+1)+ylim(0,max(plotdata));

plotdata2=c(Q20BaseFraction,Q30BaseFraction,Depth1,Depth2,Depth5,Depth10)
plotname2=c("Q20","Q30","Depth 1","Depth 2","Depth 5","Depth 10")
#mydata2=cbind(read.table(text=plotname2),plotdata2)
mydata2=data.frame(plotname2,plotdata2)
q7=ggplot(mydata2,aes(x=mydata2[,1],y=mydata2[,2]))+ggtitle("Summary")+scale_x_discrete(limits=c("Q20","Q30","Depth 1","Depth 2","Depth 5","Depth 10"))+geom_bar(stat="identity",fill="grey")+xlab('Stat')+ylab('Fraction')

multiplot(q1,q2,q3,q4, cols=2);
multiplot( q6, q7,q5,cols=1);
dev.off();

