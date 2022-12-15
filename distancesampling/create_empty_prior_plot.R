library(jpeg)
jpeg("priorplot.jpg", width=6.5,height=6.5,res=300,units="in")

par(mar=c(5,5,1,1),mgp=c(3.25,1,0))
plot(c(0.01,100),c(1,100),type="n",log="xy",
     xlab="groups per sq km",
     ylab="individuals per group",
     xaxt="n",
     las=1,
     cex.axis=2,
     cex.lab=2)
axis(1,at=(c(0.01,0.1,1,10,100)),labels=as.character(c(0.01,0.1,1,10,100)),cex.axis=2)
dev.off()
