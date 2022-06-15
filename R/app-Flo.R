###################################################################
# Flo App, created by Florian Vetsch for population genetic analysis
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodFlo <- function(input = NA, output = NA, param = NA,
                           htmlFile = "00index.html") {
  # #setwdNew(basename(output$getColumn("Report")))
  dataset <- input$meta
  ans4Report <- list() # a list of results for rmarkdown report
  ans4Report[["dataset"]] <- dataset

  output_dir <- basename(output$getColumn("Report"))
  prefix <- file.path(output_dir, "vcf_stats")

  # For Rmd
  # SNP counts
  snp_counts <- file.path(output_dir, "vcf_stats.snps")

  # InDel counts
  # ToDo

  # Private SNP counts
  private_snp_counts <- file.path(output_dir, "vcf_stats.private")

  # Shared SNP counts
  shared_snp_counts <- file.path(output_dir, "vcf_stats.shared")

  # Transions/Transversions
  tstv <- file.path(output_dir, "vcf_stats.samples-tstv")

  # run vcf-stats
  cmd <- paste("vcf-stats", file.path("/srv/gstore/projects", input$getColumn("Filtered VCF")), "-p", prefix)
  result <- ezSystem(cmd)
  gc()

  
  # convert vcf to gds format
  library(SNPRelate)
  library(gdsfmt)
  
  snp_pa <- file.path("/srv/gstore/projects", input$getColumn("Filtered VCF"))
  # Reformat
  snpgdsVCF2GDS(snp_pa, "snp.gds", method="biallelic.only")
  # open a GDS file
  snp <- snpgdsOpen("snp.gds")
  # Run PCA
  # algorithm, num.thread, bayesian
  pca <- snpgdsPCA(snp, autosome.only=F, remove.monosnp=F)
  vars <- pca$varprop[is.nan(pca$varprop) == F]
  vars_sum <- cumsum(vars)[cumsum(vars) < 0.8 ] * 100
  # make a data.frame
  df <- data.frame(sample.id = pca$sample.id,
                    EV1 = pca$eigenvect[,1],    # the first eigenvector
                    EV2 = pca$eigenvect[,2],    # the second eigenvector
                    EV3 = pca$eigenvect[,3],
                    stringsAsFactors = T)
 
  
  # Phylogenetic tree
  library("fastreeR", lib.loc = "/misc/GT/analysis/florian/R_LIBS") # without the lib.loc it doesn't find library
  library(ape)

  # exmp <- vcf2dist("~/ragi_highcov_sa0001_1k.vcf.gz")
  # exmp <- dist2tree(exmp)
  tree_newick <- vcf2tree(snp_pa)
  tr <- ape::read.tree(text=tree_newick)
  
  ## Copy the style files and templates
  styleFiles <- file.path(
    system.file("templates", package = "ezRun"),
    c(
      "fgcz.css", "Flo.Rmd",
      "fgcz_header.html", "banner.png"
    )
  )
  file.copy(from = styleFiles, to = ".", overwrite = TRUE)

  ### generate the main reports
  rmarkdown::render(
    input = "Flo.Rmd", envir = new.env(),
    output_dir = ".", output_file = htmlFile, quiet = TRUE
  )

  html_files <- c("00index.html",  "banner.png",  "fgcz.css",  "fgcz_header.html")
  file.copy(from = html_files, to = "vcf_stats")
  if(file.exists("00index_files")){
    cmd <- "mv rmarkdownLib vcf_stats; mv 00index_files vcf_stats"
  }else{
    cmd <- "mv rmarkdownLib vcf_stats"
  }
  ezSystem(cmd)

  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodMinimal(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run
##' @section Functions:
##' \itemize{
##'   \item{\code{plotReadCountToLibConc(dataset, colname): }}
##'   {Plots \code{colname} from \code{dataset} against read counts in millions.}
##'   \item{\code{getQualityMatrix(inputFile): }}
##'   {Gets a quality count matrix from a fastq or gziped fastq.gz file with dimensions read quality and read length.}
##'   \item{\code{plotQualityMatrixAsHeatmap(qualMatrixList, isR2=FALSE, xScale=1, yScale=1): }}
##'   {Returns a png table of quality matrices interpreted as heatmaps.}
##'   \item{\code{plotQualityHeatmap(result, name=NULL, colorRange=c(0,sqrt(40)), colors=gray((1:256)/256), main=NULL, pngFileName=NULL, xScale=1, yScale=1): }}
##'   {Creates and returns the images used by \code{plotQualityMatrixAsHeatmap()}.}
##' }

EzAppFlo <-
  setRefClass("EzAppFlo",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodFlo
        name <<- "EzAppFlo"
        appDefaults <<- rbind(perLibrary = ezFrame(Type = "logical", DefaultValue = TRUE, Description = "Flo brabra"))
      }
    )
  )

