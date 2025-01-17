ezMethodMageckTest = function(input=NA, output=NA, param=NA){
  require(Herper)
  require(stringr)

  # Loading the variables
  dataset <- input$meta
  dir.create(param$comparison, showWarnings=FALSE)
  outputPrefix <- file.path(param$comparison, output$getNames())
  
  mergedCountFileName <- paste0(output$getNames(), ".merged.count.tsv")
  mergedCountFileLoc <- file.path(param$comparison, mergedCountFileName)
  sampleNames <- rownames(dataset)
  
  # Combining the count files into a single count file
  countFilePaths <- file.path(param$dataRoot, dataset$`Count [File]`)
  counts <- lapply(countFilePaths, data.table::fread)
  mergeCounts <- counts %>% reduce(inner_join, by="sgRNA")
  mergeCounts <- mergeCounts %>%
    rename("Gene"="Gene.x") %>%
    select(!starts_with("Gene."))
  
  # Rename the sample columns to the actual sample names
  sampleColumns <- !(colnames(mergeCounts) %in% c("sgRNA", "Gene"))
  colnames(mergeCounts)[sampleColumns] <- sampleNames
  
  ezWrite.table(mergeCounts, file=mergedCountFileLoc, row.names=FALSE)
  
  # We give the design of the experiment as indices corresponding to the
  # columns (skipping the first 2 positions) which are sample vs ref groups
  fullColumnName <- paste(param$grouping, "[Factor]")
  #which(input$getColumn(param$grouping) == param$refGroup)-1
  rId <- paste(which(dataset[[fullColumnName]] == param$refGroup) - 1, collapse=",")
  sId <- paste(which(dataset[[fullColumnName]] == param$sampleGroup) - 1, collapse=",")
  
  # Load the conda environment
  local_CondaEnv("mageckenv", pathToMiniConda = "/usr/local/ngseq/miniconda3")

  opt <- c(
    "test",
    "-k",
    mergedCountFileLoc,
    "-t",
    sId,
    "-c",
    rId,
    "-n",
    outputPrefix,
    as.vector(str_split(param$specialOptions, "\ +", simplify = TRUE))
  )
  
  # Execute the command
  system2("mageck", args=opt)
  
  # We convert the raw outputs to xlsx files
  lapply(c(".sgrna_summary", ".gene_summary"), function(fileComp) {
    dat <- ezRead.table(paste0(outputPrefix, fileComp, ".txt"), row.names=NULL)
    writexl::write_xlsx(dat, paste0(outputPrefix, fileComp, ".xlsx"))
  })
  
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodMageckTest(input=NA, output=NA, param=NA)
##' @description Use this reference class to run Mageck Test
##' @author Falko Noé
EzAppMageckTest <-
  setRefClass("EzAppMageckTest",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMageckTest
                  name <<- "EzAppMageckTest"
                  appDefaults <<- rbind(
                    outputDir = ezFrame(
                      Type = "character",
                      DefaultValue = ".",
                      Description = "Output directory"
                    )
                  )
                }
              )
  )
