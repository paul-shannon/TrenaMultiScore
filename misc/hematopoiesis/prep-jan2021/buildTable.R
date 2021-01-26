#------------------------------------------------------------------------------------------------------------------------
library(httr)
library(jsonlite)
library(raster) # for cv
library(RUnit)
library(igvR)
library(TrenaMultiScore)
library(TrenaProjectErythropoiesis)
library(org.Hs.eg.db)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tpe"))
    tpe <- TrenaProjectErythropoiesis(quiet=FALSE)
#------------------------------------------------------------------------------------------------------------------------
gois <- c("CHD4",
          "AFP",
          "AKAP4",
          "ATP6V1B1-AS1",
          "C10orf107",
          "C1orf194",
          "C8orf48",
          "CALR3",
          "CAMK2N2",
          "CCL1",
          "CLCA2",
          "CYP4F11",
          "DHRS2",
          "DRD2",
          "ENTHD1",
          "EPB41L4A-AS2",
          "HMMR-AS1",
          "KRT16",
          "LINC01602",
          "LUM",
          "MIP",
          "MIR5581",
          "MMP12",
          "NLRP11",
          "SLC25A48",
          "SVEP1",
          "TEX19",
          "TMEM26",
          "ZNF560")
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_getRegion()
   test_CHD4()
   test_build.one()
   test_goi()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
getRegion <- function(goi, shoulder)
{
   uri <- sprintf("http://localhost:8000/geneLoc")
   body.jsonString <- sprintf('%s', toJSON(list(gene=goi, genome="hg38", shoulder=shoulder)))
   r <- POST(uri, body=body.jsonString)
   tbl.geneLoc <- fromJSON(httr::content(r)[[1]])

   list(chrom=tbl.geneLoc$chrom, start=tbl.geneLoc$start-shoulder, end=tbl.geneLoc$end+shoulder)

} # getRegion
#------------------------------------------------------------------------------------------------------------------------
test_getRegion <- function()
{
    message(sprintf("--- test_getRegion"))
    x <- getRegion("CHD4", 0)
    checkEquals(x, list(chrom="chr12", start=6556886, end=6608430))

      # entrez: chr12:6,570,082-6,614,524(GRCh38/hg38),  Size:44,443 basesOrientation:Minus strand
      # igv: chr12:6,569,082-6,608,379

} # test_getRegion
#------------------------------------------------------------------------------------------------------------------------
test_CHD4 <- function()
{
   tms.chd4 <- TrenaMultiScore(tpe, "CHD4");
   getGeneHancerRegion(tms.chd4)
   findOpenChromatin(tms.chd4)
   findFimoTFBS(tms.chd4, fimo.threshold=1e-3)
   scoreMotifHitsForConservation(tms.chd4)
   scoreMotifHitsForGeneHancer(tms.chd4)
   addDistanceToTSS(tms.chd4)

   mtx <- getExpressionMatrix(tpe, "brandLabDifferentiationTimeCourse-27171x28-namesCorrected")

   addGeneExpressionCorrelations(tms.chd4, mtx)
   addGenicAnnotations(tms.chd4)
   tbl <- getMultiScoreTable(tms.chd4)
   tbl$targetGene <- "CHD4"
   tbl.strongHits <- as.data.frame(sort(table(subset(tbl, abs(cor) > 0.75 & p.value < 1e-6)$tf)))
   tbl.strongHits <- subset(tbl.strongHits, Freq >= 3)
   checkTrue(all(c("E2F6", "ZNF410", "KLF16", "PATZ1", "ZNF263", "MAZ") %in% tbl.strongHits$Var1))

} # test_CHD4
#------------------------------------------------------------------------------------------------------------------------
get.goi <- function()
{

    oxidative.phosphorylation <-  c("ATP5PB","ATP5PF", "ATP5PO", "COX4I1", "COX5B", "COX7A2", "COX7C", "COX15", "CYC1",
                                    "DLD", "MECP2", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5", "NDUFA7", "NDUFA8", "NDUFA9",
                                    "NDUFAB1", "NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB7",
                                    "NDUFB8", "NDUFB9", "NDUFB10", "NDUFC2", "NDUFS1", "NDUFS2", "NDUFS3", "NDUFV1",
                                    "NDUFS4", "NDUFS5", "NDUFS6", "NDUFS8", "NDUFV2", "NDUFV3", "UQCRB", "UQCRC1",
                                    "UQCRC2", "UQCRH", "COX5A", "ATP5PD", "ATP5MG", "UQCRQ", "UQCR10", "STOML2",
                                    "NDUFA12", "NDUFS7", "COX6C", "COX8A", "ECSIT", "NDUFB11", "SLC25A6", "VDAC1",
                                    "SLC25A13", "NRF1", "TFAM", "ATP2A3", "CEBPA", "ATF2", "FXN", "SSBP1", "SUPV3L1",
                                    "TOMM20", "TOMM40", "MTX2", "IMMT", "TIMM10", "TIMM9", "COA3", "TIMM21", "PAM16",
                                    "GDAP1", "TOMM7", "CHCHD3", "ATAD3A", "RHOT1", "DNAJC11", "AGK", "TOMM22", "APOO",
                                    "DNAJC19", "APOOL", "ROMO1", "TOMM5", "CISD1", "NME6", "MYC", "HBZ", "PRDX1", "PNKD",
                                    "PSAT1", "KDM3A", "PM20D2")

    translation <- c("MRPL49", "EEF1A2", "RPL10A", "RPL5", "RPL6", "RPL11", "RPL17", "RPL23A", "RPL37A", "RPL38", "RPLP0",
                     "RPLP1", "RPLP2", "RPS3A", "RPS5", "RPS10", "RPS12", "RPS13", "RPS14", "RPS18", "RPS19", "RPS21",
                     "RPS27", "RPS28", "SRP9", "SRP14", "SSR4", "EIF5B", "MRPS30", "SEC61B", "RPL35", "RPL36", "MRPS18B",
                     "MRPL13", "MRPL15", "MRPS16", "MRPS23", "MRPS33", "SARS2", "MRPL16", "PTCD3", "MRPS18A", "MRPL47",
                     "MRPS25", "MRPS9", "MRPS6", "MRPL41", "MRPL40", "MRPL38", "MRPL11", "MRPS34", "MRPL43", "MRPS36",
                     "MRPL54", "CHCHD1", "MRPL52", "ABCF1", "DDX3X", "GLE1", "ILF3", "NDUFA7", "TARBP2", "TCOF1", "PUM1",
                     "METAP1", "PUM2", "GEMIN5", "COA3", "HSPA14", "TRMT10C", "SRBD1", "UPF3B", "EIF2A", "SUPV3L1",
                     "TFAM", "MAP1A", "CLASP2", "TMOD2", "MYC", "CHMP1B", "VDAC1", "TOMM7")

    ribosome <- c("RPL10A", "RPL5", "RPL6", "RPL11", "RPL17", "RPL23A", "RPL37A", "RPL38", "RPLP0", "RPLP1", "RPLP2",
                  "RPS3A", "RPS5", "RPS10", "RPS12", "RPS13", "RPS14", "RPS18", "RPS19", "RPS21", "RPS27", "RPS28",
                  "RPL35", "RPL36", "MRPL13", "MRPL15", "MRPS16", "MRPL16", "MRPS18A", "MRPS9", "MRPS6", "MRPL11",
                  "SRP9", "SRP14", "SSR4", "SEC61B", "GPAA1", "SEC16A", "UBAC2", "EEF1A2", "SLC25A6", "UQCRC2",
                  "TOMM20", "TOMM40", "MTX2", "CHP1", "TIMM10", "TIMM9", "TIMM21", "PAM16", "GDAP1", "TOMM7", "AGK",
                  "TOMM22", "DNAJC19", "ROMO1", "TOMM5", "HIST1H1D", "HNRNPU", "ILF3", "HNRNPM", "YBX1", "TCOF1",
                  "H1FX", "EIF5B", "DKC1", "HIST1H1B", "DNAJA1", "H2AFY", "UPF3B", "DDX6", "SUPV3L1", "PUM1", "PUM2",
                  "EDC4", "MYEF2", "ZC3H14", "DDX3X", "GLE1", "EIF2A", "TRMT10C", "MRM3", "NOP10", "DHX30", "CLASP2",
                  "SNRNP35", "PUF60", "GEMIN5", "RBMX", "TRMT6", "FYTTD1", "PM20D2", "CRIPT", "TM9SF3", "B2M", "TAF15",
                  "ELMO1", "DLD", "NDUFAB1", "PSAT1", "TARBP2", "HSPA14", "YES1", "PRKRA", "SARS2")


    cholesterol <- c("EBP", "MVD", "SQLE", "HMGCR", "CYP51A1", "DHCR7", "HMGCS1", "ACLY", "MVK", "HSD17B7")

    small.molecule.metabolism <- c("NDUFAF4", "B3GALT6", "PNMT", "HMGCR", "CSPG4", "ACOT1", "ACOT4", "CMBL",
                                   "HMOX2", "MCCC2", "ELOVL2", "CHST14", "SPR", "NQO1", "DDAH1", "DCTPP1",
                                   "DHCR24", "LPCAT2", "TCN1", "SQLE", "PANK1", "NME1", "ALDH1B1", "PUDP",
                                   "MVK", "PLA2G3", "SLC27A2", "XYLB", "HSD3B2", "CYP1B1", "MVD", "CYP51A1",
                                   "LUM", "GNG13", "HMGCS1", "PLPP1", "SLC19A1", "AGMAT", "PDSS1", "PFAS",
                                   "GPD1L", "UGT1A6", "CBR1", "ALB", "DHCR7", "ADRA2A", "FASN", "GCSH",
                                   "GALE", "HS6ST1", "BDH1", "HSD17B7", "HPGDS", "EBP", "SCD", "MAOB", "ACLY",
                                   "CYP4F11", "KHK", "LIPG", "SLC5A6", "AACS", "HPGD", "FABP5", "MGST1", "MGST2")


    others <- c("TAL1", "LYL1", "GATA2", "GFI1B", "CIAPIN1", "ADNP", "PFKP", "ATP2A3", "ATP2B4", "VDAC1", "GNA14", "MYC")

    original.set <- c("AFP", "AKAP4", "ATP6V1B1-AS1", "C10orf107", "C1orf194", "C8orf48", "CALR3", "CAMK2N2", "CCL1", "CLCA2",
                      "CYP4F11", "DHRS2", "DRD2", "ENTHD1", "EPB41L4A-AS2", "HMMR-AS1", "KRT16", "LINC01602", "LUM", "MIP",
                      "MIR5581", "MMP12", "NLRP11", "SLC25A48", "SVEP1", "TEX19", "TMEM26", "ZNF560")
    new.genes <- c("CHD4")

    length(oxidative.phosphorylation)   #  98
    length(translation)                 #  82
    length(ribosome)                    # 105
    length(cholesterol)                 #  10
    length(small.molecule.metabolism)   #  66
    length(others)                      #  12
    length(original.set)                #  28
    length(new.genes)                   #   1

    goi <- c(oxidative.phosphorylation, translation, ribosome, cholesterol,
             small.molecule.metabolism, others, original.set, new.genes)
    length(goi)          # 402
    length(unique(goi))  # 310
    x <- toupper(unique(goi))

    name.fixes <- list("C10ORF107"="CABCOCO1",
                       "EPB41L4A-AS2"="EPB41L4A-DT",
                       #"HIST1H1B"="H1-5",
                       # "HIST1H1D"="H1-3",
                       #"H1FX"="H1-10",
                       # "H2AFY"="MACROH2A1",
                       "C1ORF194"="C1orf194",
                       "C8ORF48"="C8orf48")
    matched.indices <- match(names(name.fixes), x)
    x[matched.indices] <- as.character(name.fixes)

    sort(x)

} # get.goi
#------------------------------------------------------------------------------------------------------------------------
# make sure that all of our genes have standardized gene symbol names, which can be used
# to look up hg38 location, and in genehancer, to get promoter and enhancer regions,
# and which match the expression matrix use to calculate tf/targetGene correlations
test_goi <- function()
{
   message(sprintf("--- test_goi"))

   goi <- get.goi()
   tbl.syms <- select(org.Hs.eg.db, keytype="SYMBOL", keys=goi, columns=c("ENTREZID", "SYMBOL")) # "ENSEMBL",
   checkEquals(nrow(tbl.syms), length(goi))
   mtx <- getExpressionMatrix(tpe, "brandLabDifferentiationTimeCourse-27171x28-namesCorrected")
   all(goi %in% rownames(mtx))

   # missing.genes <- setdiff(goi, rownames(mtx))
   # missing.genes.original.names.in.matrix <- intersect(names(name.fixes), rownames(mtx))
   #
   #                    #   mtx    genehancer & current hg38 annotations
   #  matrix.fixes <- list("ATP5L"="ATP5MG",
   #                      "ATP5F1"="ATP5PB",
   #                      "ATP5H"="ATP5PD",
   #                      "ATP5J"="ATP5PF",
   #                      "ATP5O"="ATP5PO",
   #                      "C10orf107"="CABCOCO1",
   #                      "EPB41L4A-AS2"="EPB41L4A-DT")
   # mtx.new <- mtx
   # rownames(mtx.new)[match(names(matrix.fixes), rownames(mtx.new))] <- as.character(matrix.fixes)
   # all(goi %in% rownames(mtx.new))
   # mtx <- mtx.new
   # save(mtx, file="~/github/TrenaProjectErythropoiesis/inst/extdata/expression/brandLabDifferentiationTimeCourse-27171x28-namesCorrected.RData")

} # test_goi
#------------------------------------------------------------------------------------------------------------------------
build.one <- function(gene)
{
   tms.gene <- TrenaMultiScore(tpe, gene);
   getGeneHancerRegion(tms.gene)
   findOpenChromatin(tms.gene)
   findFimoTFBS(tms.gene, fimo.threshold=1e-3)
   scoreMotifHitsForConservation(tms.gene)
   scoreMotifHitsForGeneHancer(tms.gene)
   addDistanceToTSS(tms.gene)

   mtx <- getExpressionMatrix(tpe, "brandLabDifferentiationTimeCourse-27171x28-namesCorrected")
   addGeneExpressionCorrelations(tms.gene, mtx)
   addGenicAnnotations(tms.gene)
   tbl <- getMultiScoreTable(tms.gene)
   tbl$targetGene <- gene

   invisible(tbl)

} # build.one
#------------------------------------------------------------------------------------------------------------------------
test_build.one <- function()
{
   tbl.cab <- build.one("CABCOCO1")
   tbl.sub <- subset(tbl.cab, gh > 10 & abs(cor) > 0.7)
   tbl.freq <- as.data.frame(table(tbl.sub$tf))
   checkTrue(all(c("ETV4", "GATA2", "MEF2C", "TEAD4") %in% tbl.freq$Var1))

} # test_build.one
#------------------------------------------------------------------------------------------------------------------------
