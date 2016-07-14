require(ggplot2)
require(scales)
##require(cowplot)

## input parameters
args<-commandArgs(TRUE)

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

cbPalette <- c("#999999", "#E69F00", "#56B4E9","#009E73", "mediumpurple3","#F0E442", "#0072B2",  "#CC79A7")

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

#                 V1       V2
#1       Total PETs: 0.500000
#2      Expect PETs: 0.045756
#3 Expect both PETs: 0.045756
#4        Chim PETs: 0.000041
#5      1Empty PETs: 0.065329
#6      2Empty PETs: 0.388321
#7       Valid PETs: 0.036749
QCplotA <-function(trimstatfile){
  trim=read.table(trimstatfile,sep="\t")
  trim$V2 = trim$V2/1000000
  type = c("Both","1Empty","2Empty","Chimeric","Valid")
  x = c(rep("All",4), "Valid")
  count = c(trim$V2[3], trim$V2[5], trim$V2[6], trim$V2[4], trim$V2[7])
  trimdf <- data.frame(type,x,count)
  
  trimdf$labelpos <- c(trimdf$count[1]/2, 
                       sum(trimdf$count[1:2])-trimdf$count[2]/2,
                       sum(trimdf$count[1:3])-trimdf$count[3]/2,
                       sum(trimdf$count[1:4])-trimdf$count[4]/2,
                       trimdf$count[5]/2)
  trimdf$label <- trimdf$count/trim$V2[1]*100
  figa <- ggplot(trimdf,aes(x=x,y=count,fill=type))+
    geom_bar(stat="identity")
  figa <- figa +geom_text(aes(x=x, y=labelpos, label=sprintf("%1.2f%%", label)),fontface="bold", size=2.5)+
    ylab("PET Count (Millions)")+xlab(NULL)+
    ggtitle("Linker Trimming")+
    scale_fill_manual(breaks=c("Chimeric","2Empty","1Empty","Both","Valid"),values=c("lightskyblue3",  "lightblue2","mediumpurple1", "gray","mediumpurple3" ))
}

#                V1      V2
#1              All 1892134
#2  Reads1 Low MAPQ  401037
#3  Reads1 Unmapped  235968
#4  Reads2 Low MAPQ  404841
#5  Reads2 Unmapped  238269
#6      Output PETs 1572621
#7       Duplicates   23706
#8          Uniques 1548915
#9            Intra  565352
#10           Inter  983563
#11          OneEnd  283563
QCplotB <-function(statfile){
  sams=read.table(statfile,sep="\t")
  if (dim(sams)[1]<11){stop("Input file has not enough rows")}
  samsNs = sams$V2/1000000
  type = c("High MAPQ","Low MAPQ","Unmapped","High MAPQ","Low MAPQ","Unmapped",
           "Deduplicated","Duplicates")
  x = c(rep("R1",3), rep("R2",3),rep("PET",2) )
  count = c(samsNs[1]-samsNs[2]-samsNs[3], samsNs[2], samsNs[3],
            samsNs[1]-samsNs[4]-samsNs[5], samsNs[4], samsNs[5],
            samsNs[8], samsNs[7])
  
  df <- data.frame(type,x,count)
  
  df$labelpos <- c( count[1]/2, count[1]+count[2]/2, count[1]+count[2]+count[3]/2,
                    count[4]/2, count[4]+count[5]/2, count[4]+count[5]+count[6]/2,
                    count[7]/2, count[7]+count[8]/2)
  
  df$label <- c( count[1:6]/samsNs[1], count[7:8]/samsNs[6] )*100
  df$x <- as.character(df$x)
  df$x <- factor(df$x, levels=unique(df$x))
  
  figb <- ggplot(df,aes(x=x,y=count,fill=type))+
    geom_bar(stat="identity")
  figb <- figb +geom_text(aes(x=x, y=labelpos, label=sprintf("%1.2f%%", label)),fontface="bold", size=2.5)+
    ylab("Count (Millions)")+xlab(NULL)+
    ggtitle("Read Alignment")+scale_fill_brewer(palette=1,breaks=c("Unmapped","Low MAPQ","High MAPQ","Duplicates","Deduplicated"))
  figb
}



#                V1      V2
#All  106695748
#Reads1 Low MAPQ	49282906
#Reads1 Unmapped	5496441
#Reads2 Low MAPQ	46326471
#Reads2 Unmapped	10592972
#Output PETs	49296777
#Duplicates	46792665
#Uniques	2504112
#Intra	1310707
#Inter	490330
#OneEndMapped	703075
#PETs in the same peak	1242867
#Intra PETs bewteen two peaks	27128
#Inter PETs bewteen two peaks	455182

QCplotC <-function(statfile){
  sams=read.table(statfile,sep="\t")
  if (dim(sams)[1]<14) {stop("Input file has not enough rows")}
  samsNs = sams$V2/1000000
  type = c("Intra_Chrom","Inter_Chrom", "OneEndMapped",  "Intra_BetweenPeaks","Inter_BetweenPeaks","InSamePeak","Others")
  x = c( rep("All",3), rep("Filted",4) )
  count = c(samsNs[9], samsNs[10], samsNs[11],
            samsNs[13], samsNs[14], samsNs[12], samsNs[9]+samsNs[10]-sum(samsNs[12:14]))
  
  df <- data.frame(type,x,count)  
  df$labelpos <- c( count[1]/2, count[1]+count[2]/2, count[1]+count[2]+count[3]/2,
                    count[4]/2, count[4]+count[5]/2, count[4]+count[5]+count[6]/2, sum(count[4:6])+count[7]/2 )
  
  df$label <- c( count[1:3]/samsNs[8], count[4:7]/sum(count[4:7]) )*100
  df$x <- as.character(df$x)
  df$x <- factor(df$x, levels=unique(df$x))
  
  figc <- ggplot(df,aes(x=x,y=count,fill=type))+geom_bar(stat="identity")
  figc <- figc +geom_text(aes(x=x, y=labelpos, label=sprintf("%1.2f%%", label)),fontface="bold", size=2.5)+
    ylab("Count (Millions)")+xlab(NULL)+ggtitle("Read Pairs")
  figc+scale_fill_manual(values=cbPalette[1:7],breaks=c("OneEndMapped","Inter_Chrom","Intra_Chrom","Others","InSamePeak","Inter_BetweenPeaks","Intra_BetweenPeaks"))
}


#####
#QCplotD 
# args[1]: directory
# args[2]: name
prefix = paste(args[1],"/",args[2], sep='')
bedpefile1 = paste(prefix,".interactions.intra.bedpe", sep='')
bedpefile2 = paste(prefix,".interactions.inter.bedpe", sep='')
pdfout = paste(prefix,".QCplot.pdf", sep='')

#bedpefile1 = "test.interactions.inter.bedpe"
#bedpefile2 = "test.interactions.intra.bedpe"
#pdfout = "ChIA-PET2-QCplot.pdf"

dfinter <- read.table(bedpefile2, colClasses = c(rep("NULL", 10), "integer"),head=F)
interN <- length(dfinter$V11)
interN4 <- sum(dfinter$V11>=4)
dfintra <- read.table(bedpefile1, colClasses = c(rep("NULL", 10), "integer"),head=F)

intraN <- length(dfintra$V11)
intraN4 <- sum(dfintra$V11>=4)
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank())
#    plot.title=element_text(size=10, face="bold"))

#dfa <- data.frame(group = c("Inter", "Intra"), value = c(interN, intraN))
#bp<- ggplot(dfa, aes(x="", y=value, fill=group))+geom_bar(width = 1, stat = "identity")
#pie <- bp + coord_polar("y", start=0)+ blank_theme +theme(axis.text.x=element_blank()) +
#  geom_text(aes(y = value/2 + c(0, cumsum(value)[-length(value)]), label = percent(value/sum(value))), size=3)
#figD1a <- pie + scale_fill_manual(values=c("#9999FF", "#996699"))+theme(legend.position="right")

dfb <- data.frame(group = c("Interchromosomal", "Intrachromosomal"), value = c(interN4, intraN4))
bp<- ggplot(dfb, aes(x="", y=value, fill=group))+geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0)+ blank_theme +theme(axis.text.x=element_blank()) +
  geom_text(aes(y = value/2 + c(0, cumsum(value)[-length(value)]), label = sprintf("%1.2f%%", 100*value/sum(value)) ), size=2.5)
figD1b <- pie + scale_fill_manual("PETs 4+",values=c("#9999FF", "#996699"))+ggtitle("Chromatin Interactions")+theme(legend.position="right")

#########################
dat <- table(dfintra)
df <- as.data.frame(dat)
ge9id <-as.numeric(as.character(df$dfintra))>=9
ge9count <- sum(df$Freq[ge9id])
newdf <- df[!ge9id,]
adddf <- data.frame(dfintra='9+',Freq=ge9count)
newdf <- rbind(newdf,adddf)
colnames(newdf) <- c("PET","Number")
figD2 <- ggplot( newdf, aes(x=PET, y = Number, fill=PET)) +
  geom_bar(stat="identity")+ geom_text(aes(label=Number),vjust=-0.25,size=2.5)+
  scale_fill_brewer(palette=3)+
  theme(legend.position="none") + xlab('Intrachromosomal PET count')

#######################
dat <- table(dfinter)
df <- as.data.frame(dat)
ge9id <-as.numeric(as.character(df$dfinter))>=9
ge9count <- sum(df$Freq[ge9id])
#print(ge9count)
newdf <- df[!ge9id,]
adddf <- data.frame(dfinter='9+',Freq=ge9count)
newdf <- rbind(newdf,adddf)
colnames(newdf) <- c("PET","Number")
figD3 <- ggplot( newdf, aes(x=PET, y = Number, fill=PET)) +
  geom_bar(stat="identity")+ geom_text(aes(label=Number),vjust=-0.25,size=2.5)+
  scale_fill_brewer(palette=1)+
  theme(legend.position="none") + xlab('Interchromosomal PET count')


figA <- QCplotA(paste(prefix,'.trim.stat', sep=''))
figB <- QCplotB(paste(prefix,'.bedpe.stat', sep=''))

figC <- QCplotC(paste(prefix,'.bedpe.stat', sep=''))


##pp <- plot_grid(figA,figB, figC,figD1b, figD2, figD3, labels = c("A","B","C","D","E","F"), ncol=2,label_size = 12)
pdf(pdfout, width=7.5,height=8)
multiplot(figA,figB, figC,figD1b, figD2, figD3, layout=matrix(c(1,2,3,4,5,6), ncol = 2, byrow = TRUE))
##pp
dev.off()
#bedpef1 = "ctcf2.interactions.intra.bedpe"
#bedpef2 = "ctcf2.interactions.intra.bedpe"
#QCplotD(bedpef1,bedpef2,'QCplot-D.pdf')


