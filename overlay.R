plot(c(),c(),ylim=c(0,3000),xlim = c(0,5e6),type = "n",frame.plot = FALSE)
overlay<-function(contig,col="gray"){
  contig<-contig[order(contig$start),]
  lines(contig$start,contig$scores,col=col,lwd=2,type="l")
  contig[contig$scores>2000,]
}
