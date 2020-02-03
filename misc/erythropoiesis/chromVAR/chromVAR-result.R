# chromVAR-result.R
#------------------------------------------------------------------------------------------------------------------------
# from first use of chromVAR, on brand cord atac:
#    khaleesi ~/github/TrenaProjectErythropoiesis/misc/chromVAR/quickLook/success.R
#  mtx.counts[regions.with.high.variance,]
#       d04.1 d04.2 d08.1 d10.1 d10.2 d11.1 d11.2 d12.1 d12.2 d16.1 d16.2
#  [1,]     0     0     0     0     1     4     3     2     4     1     4
#  [2,]     1     2     6     3     3     4     2     6     5     1     2
# regions.with.high.variance
#  [1] 33236 37245
#  sd(mtx.counts[regions.with.high.variance[1],]) #
#  [1] 1.737292
# rowRanges(counts_filtered)[33236,]
#  GRanges object with 1 range and 7 metadata columns:
#        seqnames            ranges strand |                               name     score signalValue    pValue    qValue      peak              bias
#    [1]    chr19 50376261-50377475      * | ATAC_Cord_d16_rep1_hg38_peak_21024       262     9.48854    29.559  26.21702        81 0.631275720164609
#
# showGenomicRegion(igv, "chr19:50376261-50377475")
#------------------------------------------------------------------------------------------------------------------------
track <- DataFrameAnnotationTrack("roi", data.frame(chrom="chr19", start=50376261, end=50377475, stringsAsFactors=FALSE),
                                  color="darkRed", trackHeight=25)
displayTrack(igv, track)
source("~/github/regulatoryGenomePaper/demos/common.R")
tbl.atac <- getATACseq("chr19", 50376000, 50378000)
tbl.atac.2 <- rbind(data.frame(chrom=rep("chr19",2), start=rep(50376000,2), end=rep(50378000,2),
                               name=c("d04_rep1", "d04_rep2"), c5=rep(0,2), stringsAsFactors=FALSE),
                  tbl.atac)


for(r in seq_len(nrow(tbl.atac))){
   track.name <- sub("_hg38_macs2_default_peak.*$", "", sub("ATAC_Cord_", "", tbl.atac$name[r]))
   track <- DataFrameQuantitativeTrack(track.name, tbl.atac[r, c(1,2,3,5)], autoscale=FALSE, min=0, max=500, color="random")
   displayTrack(igv, track)
   }

tbl.gh <- getEnhancers(tp, targetGene="NR1H2")
track <- DataFrameQuantitativeTrack("GH", tbl.gh[, c("chrom", "start", "end", "combinedscore")],
                                    autoscale=FALSE, min=0, max=50, color="darkGreen")
displayTrack(igv, track)

roi <- getGenomicRegion(igv)
tbl.chip <- getChipSeq(tp, roi$chrom, roi$start, roi$end)
