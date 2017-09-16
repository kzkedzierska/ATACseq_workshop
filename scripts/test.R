###  checks for the installed packages 
packages.needed <- c("ATACseqQC", "Diffbind", "MotifDb", "TxDb.Mmusculus.UCSC.mm10.knownGene")
which.ones <- packaged.needed %in% installed.packages()
if (sum(which.ones) > 0) {
		source("https://bioconductor.org/biocLite.R")
}
for (f in pakgaes.needed[]) {
		biocLite(f)
}

### check again
packages.needed <- c("ATACseqQC", "Diffbind", "MotifDb", "TxDb.Mmusculus.UCSC.mm10.knownGene")
which.ones <- packaged.needed %in% installed.packages()
if (sum(which.ones) > 0) {
		print("Houston we have a problem!")
	print(paste0(paste(packages.needed[which.ones], sep=', '), " were not installed."))
}
