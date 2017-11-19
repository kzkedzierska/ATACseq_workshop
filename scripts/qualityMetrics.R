require(ATACseqQC)
require(grid)
sampleSheet <- read.csv("sampleSheet_up.csv", 
                        sep = "\t", 
                        header = FALSE, 
                        stringsAsFactors = FALSE)
bamfiles <- sampleSheet$V6
bamfiles.labels <- as.character(mapply(FUN = gsub, 
                                       basename(bamfiles), 
                                       MoreArgs = list(pattern = "_part.bam", replacement = "")))

# fragSize <- fragSizeDist(bamfiles, bamfiles.labels) 
# Helper loop because of the issue with the fragSizeDist function
for (a in 1:length(bamfiles)) {
	grid.newpage()
	fragSize <- fragSizeDist(bamfiles[a], bamfiles.labels[a])
}
