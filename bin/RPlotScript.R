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
input="E:\\Dropbox\\workingspace\\FastqA\\FastPopCon\\FastqA\\bin\\HG00553"
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
#CyMax=max(mydata[,1]);
q2=ggplot(mydata)+geom_point(aes(y=EmpericalQuality1,x=Cycle))+ggtitle("EmpCycleDist")+xlim(0,100)+ylim(0,80)+ylab('EmpericalQuality')


mydata= read.table(paste(input,".EmpRepDist",sep=""),header=FALSE);
mydata=mydata[1:40,];
EmpericalQuality2=mydata[,4];
ReportQuality=mydata[,1]-33;
q3=ggplot(mydata)+geom_point(aes(y=EmpericalQuality2,x=ReportQuality))+ggtitle("EmpRepDist")+xlim(0,45)+ylim(0,45)+ylab('EmpericalQuality');


mydata= read.table(paste(input,".GCDist",sep=""),header=FALSE);
GC=mydata[,1];
AvgDepth=mydata[,4];
#MaxDepth=max(AvgDepth);
q4=ggplot(mydata,aes(x=GC,y=AvgDepth))+geom_line()+ggtitle("GCDist")+xlim(0,100)+ylim(0,max(AvgDepth));


mydata= read.table(paste(input,".InsertSizeDist",sep=""),header=FALSE);
mydata=mydata[mydata[,2]<(mean(mydata[,2])+3*sd(mydata[,2])),]
InsertSize=mydata[,1];
IS_Count=mydata[,2];
q5=ggplot(mydata,aes(x=InsertSize,y=IS_Count))+geom_point()+ggtitle("InsertSizeDist")+xlim(0,1000)+ylim(0,max(IS_Count[1:length(IS_Count)]));

mydata <- scan(paste(input,".summary",sep=""), what="", sep="\n")
ExpectedDepth=as.numeric(strsplit(mydata[14]," ")[[1]][5])
EstimatedDepth=as.numeric(strsplit(mydata[14]," ")[[1]][5])
AccessibleFraction=as.numeric(strsplit(strsplit(mydata[16]," ")[[1]][8],"/")[[1]][1]);
EstimatedQ20Depth=as.numeric(strsplit(mydata[17]," ")[[1]][7]);
EstimatedQ30Depth=as.numeric(strsplit(mydata[18]," ")[[1]][7]);
plotdata=c(ExpectedDepth,EstimatedDepth,EstimatedQ20Depth,EstimatedQ30Depth)
plotname=c("ExpectedDepth","EstimatedDepth","EstimatedQ20Depth","EstimatedQ30Depth")
q6=ggplot(mydata,aes(x=plotname,y=plotdata))+geom_bar()+ggtitle("Summary")+xlim(0,length(plotdata))+ylim(0,max(plotdata));
multiplot(q1,q2,q3,q4, cols=2);
multiplot(q5, Q6cols=2);
dev.off();

