x <- read.table("hetAtr.ibc",header=T)
pdf("hetAtr.ibc.pdf",height=5,width=5)
plot(x$Fhat1,x$Fhat2,xlab="Fhat1",ylab="Fhat2")
dev.off()
