x <- read.table("hetAtr.CV")
pdf("hetAtr.admix.pdf",height=10,width=5)
par(mfrow=c(6,1),mar=c(4,4,2,2))
plot(x$V1,x$V2,xlab="K",ylab="CV")
Q1 <- as.matrix(read.table("hetAtr.ld_pruned.1.Q"))
Q2 <- as.matrix(read.table("hetAtr.ld_pruned.2.Q"))
Q3 <- as.matrix(read.table("hetAtr.ld_pruned.3.Q"))
Q4 <- as.matrix(read.table("hetAtr.ld_pruned.4.Q"))
Q5 <- as.matrix(read.table("hetAtr.ld_pruned.5.Q"))
barplot(t(Q1),col=rainbow(6))
barplot(t(Q2),col=rainbow(2))
barplot(t(Q3),col=rainbow(3))
barplot(t(Q4),col=rainbow(4))
barplot(t(Q5),col=rainbow(5))
dev.off()