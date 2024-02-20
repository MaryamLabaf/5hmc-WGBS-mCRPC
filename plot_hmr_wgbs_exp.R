#plot the function 
#######################################################################
############################# PLOT TRACKS #############################
#######################################################################
colsglobal <- c("Benign"="steelblue3", "Loc-PCa"="darkgoldenrod3",
                "N-Adjacent-tissue"="green4", "mCRPC-tSCNC"="gray2",
                "mCRPC-adeno"="firebrick3",
                "mCRPC (CMP)"="tomato", "mCRPC-cfDNA-befor1stARSi"="maroon2")
psam_col <- c("G23"="darkred","I24"="coral3","G24"="deepskyblue4","I23"="blueviolet")
numplots <- 7
pdf(file=filename, onefile=T, width=10, height=numplots/1.2)
#pdf(file=filename, onefile=T, width=10, height=numplots+0.5)
#matrows <- c(1,1,2,2,2,2,3,3,3,3,4,5,5,6:numplots)
matrows <- c(1,1,1,2,2,2,3,3,3,4,4,5,5,6,6,7)
layout(matrix(matrows, length(matrows), 1, byrow=T))

## Plot genes
par(mai=c(0.4,1,0.05,0.5))
exonrows <- bedgene$chr==chr & overlap(bedgene$start, bedgene$end, start, end)
if(sum(exonrows)>0) {
  plotGenes(bedgene[exonrows,],chr,start,end)
  labelgenome(chr,start,end,n=20,scale="Mb")
} else {
  plot.new()
}

par(mai=c(0.05,1,0.05,0.5))

## Plot 5hmc

plotBed(beddata=fivehmcslice.df,chrom=chr,chromstart = start,chromend = end,
        row='supplied',
        rownumber=roworder,
        color=fivehmcslice.df$color, border = T)
axis(side=2,las=2,tcl=.2,lwd = 0.5)
mtext('5hmc',side=2,line=4,cex=0.5)
abline(v = c(tss_start,tss_end), col="gray80", lwd=1, lty=1)

## Plot HMR
plotBed(beddata=hmrslice.df,chrom=chr,chromstart = start,chromend = end,
        row='supplied',
        rownumber=hmr_roworder,
        color=hmrslice.df$color)

axis(side=2,las=2,tcl=.2)
mtext('HMR',side=2,line=4,cex=0.5)
abline(v = c(tss_start,tss_end), col="gray80", lwd=1, lty=1)

## Plot HMR-PSMA samples
plotBed(beddata=hmrslice_psma.df,chrom=chr,chromstart = start,chromend = end,
        row='supplied',
        rownumber=hmr_psma_roworder,
        color=hmrslice_psma.df$color)

axis(side=2,las=2,tcl=.2)
mtext('HMR-PSMA-data',side=2,line=4,cex=0.5)
abline(v = c(tss_start,tss_end), col="gray80", lwd=1, lty=1)

## plot legend-PSMA data
#abline(h=0, col='black')
#legend("topright",inset=0.01,
       #legend=c('met', "met(CPM)", "met(be1stARSi)",'tSCNC','loc','NT','benign'),
#       legend=c('23G+', "24I+", "24D-",'23I-'),
#       fill=opaque(psam_col),border=psam_col,cex=0.7)


## Plot rHMR
'plotBed(beddata=rhmrslice.df,chrom=chr,chromstart = start,chromend = end,
        row="supplied",
        rownumber=rhmr_roworder,
        color=rhmrslice.df$color)

axis(side=2,las=2,tcl=.2)
mtext("rHMR",side=2,line=4,cex=0.5)'

## plot CpGs

plotBed(cpg,chr,start,end,type="region",row='supplied',rownumber=1)
axis(side=2,las=2,tcl=.2)
mtext('CPG',side=2,line=4,cex=0.5)
abline(v = c(tss_start,tss_end), col="gray80", lwd=1, lty=1)

## plot legend
#abline(h=0, col='black')
#legend("topright",inset=0.01,
#legend=c('met', "met(CPM)", "met(be1stARSi)",'tSCNC','loc','NT','benign'),
#       legend=c('Benign', "Loc", "NT",'tSCNC','met','met-CMP',
#                'met-bef1stARSi'),
#       fill=opaque(colsglobal),border=colsglobal,cex=0.7)

## Plot correlation with expression
cor_range <- c(-1,1)
colscor <- c('darkred','dodgerblue3')
plotBedgraph(peaks[,c('chr','start','end','hmc_cor')],chr,start,end,color=colscor[2],range=cor_range,transparency=0.6)
plotBedgraph(peaks[,c('chr','start','end','wgbs_cor')],chr,start,end,color=colscor[1],range=cor_range,transparency=0.6,overlay=T)
legend("topleft",inset=0,legend=c('WGBS','5hMC'),
       fill=opaque(colscor),border=colscor,cex=0.5)
axis(side=2,las=2,tcl=.2)
mtext('Rho',side=2,line=4,cex=0.5)
abline(h=0, col='black')


## Plot correlation with expression of mCRPC-tSCNC samples
#cor_range <- c(-1,1)
#colscor <- c('darkred','dodgerblue3')
#plotBedgraph(peaks_tscnc[,c('chr','start','end','hmc_cor')],chr,start,end,color=colscor[2],range=cor_range,transparency=0.6)
#plotBedgraph(peaks_tscnc[,c('chr','start','end','wgbs_cor')],chr,start,end,color=colscor[1],range=cor_range,transparency=0.6,overlay=T)
#legend("topleft",inset=0,legend=c('WGBS','5hMC'),
#       fill=opaque(colscor),border=colscor,cex=0.5)
#axis(side=2,las=2,tcl=.2)
#mtext('Rho-t-scnc',side=2,line=4,cex=0.5)
#abline(h=0, col='black')

dev.off()
