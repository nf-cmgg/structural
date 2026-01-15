#!/usr/bin/env Rscript

bamfile<-"$bam"
name <- "${task.ext.prefix ?: meta.id}"
tax_id <- 12345 #TODO fix this
annotations<-"$annotations"

#Functions
log2offset <- function(offset=.Machine\$double.xmin) {
    # Get offset?
    if (missing(offset)) {
        offset <- getOption("QDNAseq::log2offset", .Machine\$double.xmin)
        offset <- as.double(offset);
        stopifnot(is.finite(offset));
        return(offset);
    }

    # Reset offset?
    if (is.null(offset)) offset <- eval(formals(log2offset)\$offset);

    # Set offset
    stopifnot(length(offset) == 1L);
    offset <- as.double(offset);
    stopifnot(is.finite(offset));
    options("QDNAseq::log2offset"=offset);

    offset;
}
log2adhoc <- function(x, offset=log2offset(), inv=FALSE) {
    if (!inv) {
        x[x < 0] <- 0
        x <- x + offset
        log2(x)
    } else {
        x <- 2^x
        x - offset
    }
}
makeSegments <- function(data,chrdata) {
    previous    <- 2000
    chrpr       <- -100
    values      <- c()
    start       <- c()
    end         <- c()
    chrs       <- c()
    el          <- nrow(data)
    data <- rbind(data,c(-10000,-10000,-10000)) #add value to allow data[i+1]
    for (i in 1:el) {
        if ((data[i,3] != previous & previous != data[i+1,3]) | chrdata[i] != chrpr) { #bug repaired 12/06/09
            start   <- c(start, data[i,1])
            last    <- i - 1
            chrs <- c(chrs, chrdata[i])
            if (last > 0) end <- c(end, data[last,2])
            values  <- c(values, data[i,3])
        }
        previous    <- data[i,3]
        chrpr <- chrdata[i]
    }
    end     <- c(end, data[el,2])
    result  <- cbind(chrs, start, end, values)
    result
}
#Load Libraries
library('QDNAseq')
library('CGHcall')
library('lsr')
#Calculations

#QDNAseq
attach(annotations, name="bins")
readCounts <- binReadCounts(bins, bamfile, minMapq=40)
if(grepl('noresidual',annotations,ignore.case=TRUE) | grepl('noblacklist',annotations,ignore.case=TRUE)){
    if( grepl('noresidual',annotations,ignore.case=TRUE) && !grepl('noblacklist',annotations,ignore.case=TRUE)){
        readCountsFiltered <- applyFilters(readCounts, residual=NA)
    }else if(!grepl('noresidual',annotations,ignore.case=TRUE) && grepl('noblacklist',annotations,ignore.case=TRUE)){
        readCountsFiltered <- applyFilters(readCounts, blacklist=NA)
    }else if( grepl('noresidual',annotations,ignore.case=TRUE) && grepl('noblacklist',annotations,ignore.case=TRUE)){
        readCountsFiltered <- applyFilters(readCounts, residual=NA, blacklist=NA)
    }
}else{
    readCountsFiltered <- applyFilters(readCounts)
}
readCountsFiltered <- estimateCorrection(readCountsFiltered)
if ( tax_id != "10090"){
    readCountsFiltered <- applyFilters(readCountsFiltered,chromosomes =c())
}
copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
copyNumbersSegmented <- segmentBins(copyNumbersSmooth,transformFun=QDNAseq:::sqrtadhoc, alpha=0.005 )
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
exportBins(copyNumbersSegmented, file=paste(name,".bed", sep = ""), format="bed")
exportBins(copyNumbersSmooth, file=paste(name,"cna", sep="."),format="tsv")

#segmentation
segmented <- assayDataElement(copyNumbersSegmented, "segmented")
segmented <- log2adhoc(segmented)
segmented_starts <- bpstart(copyNumbersSegmented)
segmented_ends <- bpend(copyNumbersSegmented)
segmented_chr <- chromosomes(copyNumbersSegmented)
dSegs <- data.frame(segmented_chr,segmented_starts,segmented_ends,c(segmented))
dSegsNew <- dSegs[!is.na(dSegs[,'c.segmented.']),]
dSegsNew <- dSegsNew[complete.cases(dSegsNew[,'segmented_chr']),]
segment <- makeSegments(data.matrix(dSegsNew[,c('segmented_starts','segmented_ends','c.segmented.')]), data.matrix(dSegsNew[,'segmented_chr']))
write.table(segment, paste(name,"_segments.txt", sep = ""), sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE)

#calculate avg/chr
cghdata <- read.table(paste(name,"cna", sep="."), header=TRUE)
colnames(cghdata)[5]<-"value"
averages <- tapply(cghdata\$value, cghdata\$chromosome, mean)
medians <- tapply(cghdata\$value, cghdata\$chromosome, median)
names(averages) <- paste("chr",names(averages),sep="")
names(medians) <- paste("median_chr",names(medians),sep="")

#calculate stdev
std <- sd(cghdata[cghdata\$chromosome != c("X","Y"),]\$value)
mad <- mad(cghdata[cghdata\$chromosome != c("X","Y"),]\$value)
aad <- aad(cghdata[cghdata\$chromosome != c("X","Y"),]\$value)
write.table(averages, file="statistics.out", col.names = FALSE, row.names=names(averages), dec = "," )
write.table(medians, file="statistics.out", append = TRUE, col.names = FALSE, row.names=names(medians), dec = "," )
write.table(std, file="statistics.out", append = TRUE, col.names = FALSE, row.names=c("stdev"), dec = ",")
write.table(mad, file="statistics.out", append = TRUE, col.names = FALSE, row.names=c("mad"), dec = ",")
write.table(aad, file="statistics.out", append = TRUE, col.names = FALSE, row.names=c("aad"), dec = ",")

## VERSIONS FILE
writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
        paste('    r-lsr:', as.character(packageVersion('lsr'))),
        paste('    bioconductor-qdnaseq:', as.character(packageVersion('QDNAseq')))
    ),
'versions.yml')
