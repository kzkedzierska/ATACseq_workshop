###  checks for the installed packages 
packages_needed <- c("ATACseqQC", "DiffBind", "MotifDb", "TxDb.Mmusculus.UCSC.mm10.knownGene")
which_ones <- ! packages_needed %in% installed.packages()
if (sum(which_ones) > 0) {
		source("https://bioconductor.org/biocLite.R")
	biocLite(packages_needed[which_ones])
}

### check again
which_ones <- ! packages_needed %in% installed.packages()
if (sum(which_ones) > 0) {
	print("Houston we have a problem!")
	print(paste0(paste(packages_needed[which_ones], sep=', '), " were not installed."))
} else {
	print("All's good! :)")
}
