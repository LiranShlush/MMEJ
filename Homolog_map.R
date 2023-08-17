library(karyoploteR)
library(gridExtra)
chr_list = c('1','2','3','4', '5', '6', '7', '8', '9','10', '11', '12', 
             '13', '14', '15', '16', '17', '18','19','20', '21',  '22', 'X', 'Y')

kp <- plotKaryotype(genome="hg19") #<- use the genome you need, if not human
for(i in 1:24){
  dsf = read.table(input_file,  sep = '\t',header = TRUE)
  dsf$Chr = paste("chr", chr_list[i], sep = "")
  dsf <- dsf[, c('Chr', 'Start')]
 

  snps.gr <- toGRanges(dsf[,c(1,2,2)])

  kp <- kpPlotDensity(kp, data=snps.gr)
}


