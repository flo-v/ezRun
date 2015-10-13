###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Gets the last value of an object
##' @param x a vector, object or a list.
##' @template roxygen-template
##' @examples
##' lastVal(1:10)
##' obj1 <- 1:5
##' obj2 <- "a string"
##' lastVal(c(obj1,obj2))
lastVal = function(x){
  if (is.list(x)){
    x[[length(x)]]
  } else {
    x[length(x)]
  }
}

##' @title Switches to a working directory
##' @description If the directory does not exist, the function will create it recursively.
##' @param dir a character string specifying the desired directory.
##' @template roxygen-template
##' @examples
##' setwdNew("newDirectory")
setwdNew = function(dir){
  if (!file.exists(dir)){
    dir.create(dir, recursive=TRUE)
  }
  setwd(dir)
}

##' @title Creates a Venn diagram based on the overlapping elements in vectors
##' @description This is a convenience function that converts the list of vectors first into a logical matrix and subsequently applies the venn diagramm function from the limma package
##' @param setList a named list of length two or three. 
##' @template roxygen-template
##' @examples
##' aList = list(a=1:5,b=3:6)
##' vennFromSets(aList)
vennFromSets = function(setList){
  stopifnot(!is.null(names(setList)) && length(setList) %in% 2:3)
  library(limma, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  x = ezMatrix(FALSE, rows=unique(unlist(setList)), cols=names(setList))
  for (i in 1:length(setList)){
    x[match(setList[[i]], rownames(x)), i] = TRUE
  }
  vc = vennCounts(x)
  vennDiagram(vc)
}

##' @title Creates a contingency table from a list of vectors.
##' @description The list should contain two or three elements.
##' @param setList a named list of two or three elements.
##' @template roxygen-template
##' @examples
##' aList = list(a=1:5,b=3:6)
##' tableFromSets(aList)
tableFromSets= function(setList){
  stopifnot(!is.null(names(setList)) && length(setList) %in% 2:3)
  x = data.frame(row.names=unique(unlist(setList)))
  for (i in 1:length(setList)){
    x[names(setList)[i]] = rownames(x) %in% setList[[i]]
  }
	table(x)
} 

##' @title Shrink the dynamic range of a numeric vector, matrix, or data frame
##' @description Values outside the range will be set as the lower or upper boundary of the range.
##' @param x the values to modify.
##' @param theRange two values specifying the lower and upper boundary of the range.
##' @return Returns the shrinked object.
##' @template roxygen-template
##' @examples
##' shrinkToRange(1:10,c(2,6))
shrinkToRange = function(x, theRange){
  x[ x > theRange[2]] = theRange[2]
  x[ x < theRange[1]] = theRange[1]
  return(x)
}

##' @title Combines two matrices by their columns
##' @description The matrices need to have the same dimensions and column names must be equal. The first matrix needs row names as well. This is used to merge, e.g. matrices for expression signal and present flag by sample.
##' @param x first matrix containing both row and column names.
##' @param y second matrix containing the same column names as the first one.
##' @return Returns the combined matrix.
##' @template roxygen-template
##' @examples
##' m1 = matrix(1:10,2)
##' m2 = matrix(11:20,2)
##' colnames(m1) = as.character(1:5)
##' rownames(m1) = c("a","b")
##' colnames(m2) = as.character(1:5)
##' interleaveMatricesByColumn(m1,m2)
interleaveMatricesByColumn = function(x, y, suffixes=c("[Signal]", "[Present]")){

  if (any(colnames(x) != colnames(y))){
    stop("incompatible matrices\n", colnames(x), "\n", colnames(y))
  }
  combined = ezMatrix(NA, rows=rownames(x), cols=1:(ncol(x)*2))
  combined[ , (1:ncol(x))*2 -1] = x
  combined[ , (1:ncol(x))*2 ] = y
  colnames(combined) = paste(rep(colnames(x), each=2), suffixes)
  return(combined)
}

##' @title Normalization
##' @description Normalizes a matrix according to the provided method; normalization is performed column-wise.
##' @param x the matrix to normalize.
##' @param method which normalization method to use. The default does not normalize. Possible methods:
##' \itemize{
##'  \item{"none"}
##'  \item{"quantile"}
##'  \item{"logMean"}
##'  \item{"median"}
##'  \item{"vsn"}
##' }
##' @param presentFlag a binary matrix with the same size as \code{x} which indicates if a values is considered as measured correctly. Default: NULL
##' @return Returns the modified value.
##' @template roxygen-template
##' @examples 
##' ezNorm(4)
ezNorm = function(x, method="none", presentFlag=NULL){

  switch(method, "none"=x,
    "quantile"=ezQuantileNorm(x),
    "logMean"=ezLogmeanNorm(x, presentFlag=presentFlag),
		"median"=ezMedianNorm(x, presentFlag=presentFlag),
    "vsn"=ezVsnNorm(x),
    stop("Unsupported normalization method: ", method)
  )
}

##' @title Quantile Normalization
##' @description A convenience call to the corresponding method from the preprocessCore package that keeps row and column names of the input matrix.
##' @param x the matrix to normalize.
##' @return Returns the normalized matrix.
##' @seealso \code{\link[preprocessCore]{normalize.quantiles}}
##' @template roxygen-template
##' @examples 
##' m1 = matrix(1:20,4)
##' ezQuantileNorm(m1)
ezQuantileNorm = function(x){
  library(preprocessCore, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  norm = normalize.quantiles(x)
  colnames(norm) = colnames(x)
  rownames(norm) = rownames(x)
  norm
}

##' @title Variance Stabilizing Normalization
##' @description Normalizes the matrix with the Variance Stabilizing Normalization and keeps row and column names.
##' @param x the matrix to normalize.
##' @return Returns the normalized matrix.
##' @seealso \code{\link[vsn]{justvsn}}
##' @template roxygen-template
##' @examples 
##' m1 = matrix(1:200,50)
##' ezVsnNorm(m1)
ezVsnNorm = function(x, lts.quantile=0.6){
  require(vsn)
  return(2^justvsn(x, lts.quantile=lts.quantile))
}

##' @title Convert numeric to factor
##' @description Cuts a numeric vector into factors using defined breaks.
##' @param x a numeric vector to be cut.
##' @param breaks The values by which the vector should be cut. 2 or more are needed.
##' @param prefix adds a prefix to the level names.
##' @param labels use labels for the factors instead of displaying the ranges provided by breaks.
##' @return Returns the numeric input as factors.
##' @template roxygen-template
##' @examples
##' ezCut(1:10,breaks=c(2,5,7),prefix=letters[1:4])
ezCut = function(x, breaks, prefix=NULL, labels=NULL){
	
  if (is.null(labels)){
    labels = paste("<=", breaks[1])
    for (i in 2:length(breaks)){
      labels = c(labels, paste("(", breaks[i-1], " - ", breaks[i], "]", sep=""))
    }
    labels = c(labels, paste(">", breaks[length(breaks)]))
    if (!is.null(prefix)){
      labels = paste(prefix, labels)
    }
  }
  classes = cut(x, breaks = c(-Inf, breaks, Inf), labels = labels)
	return(classes)
}

##' @title Tests if x contains an error
##' @description Returns TRUE if x is a list with at least one element called error.
##' @param x any R object, but only an error in a list can be found.
##' @return Returns FALSE or TRUE.
##' @template roxygen-template
##' @examples
##' isError("error")
##' isError(list(a=3:5,error=3))
##' isError(list(errrrror=3))
isError = function(x){
	if (is.list(x)){
		if (!is.null(x$error)){
			return(TRUE)
		}
	}
	return(FALSE)
}

##' @title Matches patterns and returns a logical vector
##' @description Searches for entries in \code{x} that match \code{patterns}.
##' @param patterns A value, vector or list to match the provided values with.
##' @param x the values to match.
##' @param combine whether just one value or all from \code{patterns} need to match to return TRUE:
##' \itemize{
##'  \item{"or"}{ only one value from \code{patterns} needs to match.}
##'  \item{"and"}{ all values from \code{patterns} need to match.}
##' }
##' @seealso \code{\link[base]{grepl}}
##' @return Returns a logical vector containing the results of the tested pattern.
##' @template roxygen-template
##' @template addargs-template
##' @templateVar fun \code{grepl()}
##' @examples
##' ezGrepl(3,1:40)
##' ezGrepl(c(2,4),1:100)
##' ezGrepl(c(2,4),1:100,combine="and")
ezGrepl <- function(patterns, x, combine="or", ...){

  combine = match.arg(combine, c("or", "and"), several.ok = FALSE)
  result = grepl(patterns[1], x, ...)
  if (length(patterns) == 1){
    return(result)
  }
  for (pt in patterns[-1]){
    if (combine == "or"){
      result = result | grepl(pt, x, ...)
    } else {
      result = result & grepl(pt, x, ...)
    }
  }
  return(result)
}

##' @title Seperates a character vector into a matrix by splitting it.
##' @description The \code{split} vector needs to divide the \code{x} vector into pieces of equal lengths. The split is removed from the original input. 
##' @param x the character vector to be split element-wise.
##' @param split a character vector defining with what to split.
##' @seealso \code{\link[base]{strsplit}}
##' @seealso \code{\link[limma]{strsplit2}}
##' @return Returns a matrix containing the split vector and \code{x} as the rownames.
##' @template roxygen-template
##' @template addargs-template
##' @templateVar fun \code{strsplit()}
##' @examples
##' ezSplit(letters[1:5],"b")
##' ezSplit(rep("abcde",4),letters[1:4])
ezSplit = function(x, split, ...){

  splitList <- strsplit(x, split, ...)
  lengths = sapply(splitList, length)
  idx = which(lengths != lengths[1])
  if (length(idx) >0){
    stop(paste("Row ", idx[1], " length ", lengths[idx[1]], " but expected was ", lengths[1]))    
  }
  result <- matrix(unlist(splitList), nrow=length(splitList), ncol=lengths[1], byrow=TRUE)
  rownames(result) <- x
  result
}

##' @title Trims white space
##' @description  Removes white space at the end or in the front of a character.
##' @param x the character to trim.
##' @return Returns the trimmed character.
##' @template roxygen-template
##' @examples
##' trimWhiteSpace("    bla    ")
trimWhiteSpace = function (x){
    sub("[ \t\n\r]*$", "", sub("^[ \t\n\r]*", "", x))
}

##' @title Creates a matrix
##' @description Either use the rows and columns to define the matrix or the dimensions.
##' @param x a vector containing the matrix elements
##' @param rows a vector, whose length defines the number of rows.
##' @param cols a vector, whose length defines the number of columns.
##' @param dim a vector of length 2 defining the dimensions of the matrix.
##' @return Returns a matrix with \code{rows} and \code{cols} defining the names or using generic ones, if \code{dim} is used instead.
##' @template roxygen-template
##' @examples
##' ezMatrix(1,rows=1:4,cols=1:3)
##' ezMatrix(3:6,dim=c(4,6))
ezMatrix <- function(x, rows=NULL, cols=NULL, dim=NULL){
  if(is.null(rows) && is.null(cols)){
    return(ezMatrix(x, rows=1:dim[1], cols=1:dim[2]))
  }
  matrix(x, nrow=length(rows), ncol=length(cols), dimnames=list(rows, cols))
}

##' @title Scales columns of a matrix
##' @description The columns will be scaled by multiplying them with the \code{scaling}. 
##' @param x the matrix to scale.
##' @param scaling a vector containing the scale for each column of \code{x}.
##' @return Returns a matrix with scaled columns.
##' @template roxygen-template
##' @examples
##' m1 = matrix(1:20,5)
##' ezScaleColumns(m1,1:4)
ezScaleColumns = function(x, scaling){
  for (i in 1:length(scaling)){
    x[,i] <- x[ ,i] * scaling[i]
  }
  return(x)
}

##' @title Scales columns of a matrix to median
##' @description Columns of the matrix will be scaled to an overall median or to a defined target.
##' @param x the matrix to scale.
##' @param use a logical vector defining which rows to use.
##' @param target a value or vector defining the median(s) by which columns should be scaled. The default will use the overall median of the matrix for each column.
##' @param presentFlag a binary matrix with the same size as \code{x} which indicates if a values is considered as measured correctly.
##' @return Returns a matrix with columns normalized to a median/medians.
##' @template roxygen-template
##' @examples
##' m1 = matrix(1:20,5)
##' ezMedianNorm(m1)
##' ezMedianNorm(m1,target=10)
##' ezMedianNorm(m1,use=c(TRUE,FALSE))
ezMedianNorm = function(x, use=NULL, target=NULL, presentFlag=NULL){
  
  sf = getMedianScalingFactor(x, use=use, target=target, presentFlag=presentFlag)
  return(ezScaleColumns(x, sf))
}

# or use @rdname as a simpler approach to add a function to the help file of another.
##' @describeIn ezMedianNorm Calculates the scaling factor for the main function.
getMedianScalingFactor = function(x, use=NULL, target=NULL, presentFlag=NULL){
  
  if (is.null(use)){
    use = rep(TRUE, nrow(x))
  }
  if (!is.null(presentFlag)){
    isAllPresent = apply(presentFlag, 1, all)
    use = use & isAllPresent
  }
  medians <- apply(x[use, ], 2, median, na.rm=TRUE)
  if (is.null(target)){
    target = median(medians)
  }
  sf = target / medians
  names(sf) = colnames(x)
  return(sf)
}

##' @title Scales columns of a matrix to logarithmic mean
##' @description Columns of the matrix will be scaled to an overall logarithmic mean or to a defined target.
##' @param target a value or vector defining the means by which columns should be scaled. The default will use the overall logarithmic mean of the matrix for each column.
##' @inheritParams ezMedianNorm
##' @return Returns a matrix with columns normalized to a logarithmic mean/means.
##' @template roxygen-template
##' @examples
##' m1 = matrix(1:20,5)
##' ezLogmeanNorm(m1)
##' ezLogmeanNorm(m1,target=10)
##' ezLogmeanNorm(m1,use=c(TRUE,FALSE))
ezLogmeanNorm = function(x, use=NULL, target=NULL, presentFlag=NULL){

	sf = ezLogmeanScalingFactor(x, use=use, target=target, presentFlag=presentFlag)
  return(ezScaleColumns(x, sf))
}

##' @describeIn ezLogmeanNorm Calculates the scaling factor for the main function.
ezLogmeanScalingFactor = function(x, use=NULL, target=NULL, presentFlag=NULL){

  x.log <- log(x)
  if (is.null(use)){
    use = rep(TRUE, nrow(x))
  }
  if (!is.null(presentFlag)){
    isAllPresent = apply(presentFlag, 1, all)
    use = use & isAllPresent
  }
  means <- apply(x.log[use, ], 2, mean, na.rm=TRUE)
	if (is.null(target)){
		target = exp(mean(means))
	}
	sf = target / exp(means)
	names(sf) = colnames(x)
	return(sf)
}

##' @title Geometric mean
##' @description Calculates the geometric mean of a numeric argument.
##' @param x a value, vector or matrix for which to calculate the geometric mean.
##' @return Returns the geometric mean of \code{x}.
##' @template roxygen-template
##' @template addargs-template
##' @templateVar fun \code{mean()}
##' @examples
##' ezGeomean(1:10)
ezGeomean <- function(x, ...){
  exp(mean(log(x), ...))
}

##' @title Averages columns together
##' @description Rearranges and averages columns according to \code{by}.
##' @param x the matrix whose columns should be averaged.
##' @param by an integer or vector by which to average and/or rearrange columns.
##' @param func the function to apply to the result. Default: mean with removing NA's.
##' @return Returns a vector or matrix of averaged columns.
##' @template roxygen-template
##' @examples
##' m1 = matrix(1:20,5)
##' rownames(m1) = 1:5
##' averageColumns(m1,1)
##' averageColumns(m1,c(4,2,3,1))
##' averageColumns(m1,c(1,1,2,2))
averageColumns = function(x, by=NULL, func=function(x){mean(x, na.rm=TRUE)}){

  cols = sort(unique(by))
  result = ezMatrix(NA, rows=rownames(x), cols=cols)
  for (c in cols){
    result[ ,c] = apply(x[ , c == by, drop=FALSE], 1, func)
  }
  result
}

## --> is horribly slow
# averageColumns = function(x, by=NULL, func=mean, ...){
#   return(t(averageRows(t(as.matrix(x)), by=by, func, ...)))
# }

##' @title Averages rows together
##' @description Rearranges and averages rows according to \code{by}.
##' @param data the matrix whose rows should be averaged.
##' @param by an integer or vector by which to average and/or rearrange rows.
##' @param func the function to apply to the result. Default: mean.
##' @return Returns a vector or matrix of averaged rows.
##' @template roxygen-template
##' @template addargs-template
##' @templateVar fun \code{aggregate()}
##' @examples
##' m1 = matrix(1:20,5)
##' averageRows(m1,c(1,1,2,2,3))
# why by=labels here and before by=NULL? Neither seem to work if not defined by the user.
averageRows = function(data, by=labels, func=mean, ...){

  by <- list(AveragingID=by)
  data <- aggregate(data, by=by, func, ...)
  rownames(data) = data$AveragingID
  data$AveragingID = NULL
  return(data)
}

##' @title Inverse mapping
##' @description Strips the names and values of a list apart and rematches all the names that apply to each value.
##' @param xList a named list.
##' @return Returns a list for the input values containing all the names that belong to each value in the input list.
##' @template roxygen-template
##' @examples
##' l1 = list(a=1:3, b=c(2,5), c=4:8)
##' inverseMapping(l1)
inverseMapping = function(xList){
  
	mm = makeMultiMapping(xList)
  invMap = tapply(mm$source, mm$target, function(x){list(x)})
  return(invMap)
}

##' @describeIn inverseMapping Unlists the input and returns a data.frame with seperated names and values.
makeMultiMapping = function(xList){

	target = unlist(xList, use.names=FALSE)
	counts = sapply(xList, length)
  source = rep(names(xList), times=counts)
	data.frame(source=source, target=target, stringsAsFactors = FALSE)
}

# TODO(Rsge not supported anymore, still waiting for a reply of the authors to use their source code.)
.ezSgelapply = function(jobList, FUN, param, queue="GT", cores=4, ram=10, scratch=50, mailto=NULL,
                        saveGlobal=TRUE, removeFiles=TRUE){
  library(Rsge, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  sge.options("sge.qsub.options"=paste("-cwd -q ", queue,
                                       " -pe smp ", cores,
                                       ##" -l C=1",
                                       " -l R=", round(ram/cores, digits=2),
                                       " -l S=", round(scratch/cores, digits=2),
                                       sep=""))
  sge.options("sge.save.global"=saveGlobal)
  sge.options("sge.remove.files"=removeFiles)
  ## TODO jobs must inspect NSLOTS in order to know how many threads to use!!
  
  sgeFunc = function(x, FUN=NULL, param=NULL){
    cwd = getwd()
    wd = paste("/scratch/rjob", ezTime(), Sys.getpid(), sep="_")
    stopifnot(!file.exists(wd))
    setwdNew(wd)
    result = FUN(x, param)
    setwd(cwd)
    unlink(wd, recursive=TRUE)
    return(result)
    
  }
  sge.parLapply(jobList, sgeFunc, FUN, param, njobs=length(jobList))
}

##' @title Parallel version of \code{lapply()}
##' @description This function is a modified version of \code{mclapply()} of the parallel package and allows a parallelized usage of lapply.
##' @param x a list to apply the function to.
##' @param FUN the function to apply to each list element.
##' @return Returns a list of the same length as \code{x} with \code{FUN} applied to its elements.
##' @template roxygen-template
##' @template addargs-template
##' @templateVar fun \code{lapply()} and \code{mclapply()}
##' @seealso \code{\link[parallel]{mclapply}}
##' @examples
##' l1 = list(a=1:3, b=c(2,5), c=4:8)
##' ezMclapply(l1,sum)
## TODOP: errorLogFile is currently unused.
ezMclapply = function(x, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE, mc.silent = FALSE, mc.cores=min(length(x), ezThreads()), errorLogFile=""){
  mc.cores = min(mc.cores, length(x))
  if (mc.cores == 1){
    return(lapply(x, FUN, ...))
  }
  result = mclapply(x, FUN, ..., mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed,
                    mc.silent = mc.silent, mc.cores = mc.cores)
  gc()
  isError = sapply(result, function(x){any(grepl("error", class(x), ignore.case=TRUE))})
  if (any(isError)){
    sapply(result[isError], print)
    stop("mclapply failed")
  }
  isNull = sapply(result, is.null) 
  if (any(isNull)){
    stop("mclapply returned NULL results: ", sum(isNull), " / ", length(isNull))
  }
  return(result)
}

##' @title Which values are duplicated?
##' @description Keeps duplicated values by setting them to TRUE.
##' @param x a vector, matrix or list.
##' @param mode specifies how to keep duplicated values. Possible modes:
##' \itemize{
##'  \item{"keepFirst"}{ uses \code{duplicated()} to keep duplicates in the original order.}
##'  \item{"keepLast"}{ keeps duplicates in a reversed order.}
##'  \item{"random"}{ keeps duplicates in a randomized order.}
##'  \item{"all"}{ keeps all duplicated values independent of the order.}
##' }
##' @return Returns a logical indicating which values are duplicated according to the specified \code{mode}.
##' @template roxygen-template
##' @examples
##' v1 = c(1,2,3,4,5,4,3,2,3,4,5,6,7,8,7)
##' ezDuplicated(v1)
##' ezDuplicated(v1,"all")
ezDuplicated = function(x, mode="keepFirst"){
	switch(mode,
				 "keepFirst"=duplicated(x),
				 "keepLast" = rev(duplicated(rev(x))),
				 "random"={
						n = length(x);
						idx = sample(1:n, n, replace = FALSE);
						isDup = rep(FALSE, n);
						isDup[idx] = duplicated(x[idx]);
						return(isDup)
				 },
				 "all"={
						dups = unique(x[duplicated(x)])
						return(x %in% dups)
				 })
}

##' @title Which values occur at least \code{n} times?
##' @description Keeps multiplicated values by setting them to TRUE.
##' @param x a vector, matrix or list.
##' @param n a positive integer specifying how many times a value needs to occur to return true.
##' @param mode specifies how to keep multiplicated values. Possible modes:
##' \itemize{
##'  \item{"keepFirst"}{ keeps values occuring \code{n} times in the original order.}
##'  \item{"keepLast"}{ keeps values occuring \code{n} times in a reversed order.}
##'  \item{"random"}{ keeps values occuring \code{n} times in a randomized order.}
##'  \item{"all"}{ keeps all values occuring \code{n} times independent of the order.}
##' }
##' @return Returns a logical indicating which values are multiplicated according to the specified \code{mode} and \code{n}.
##' @template roxygen-template
##' @examples
##' v1 = c(1,2,3,4,5,4,3,2,3,4,5,6,7,8,7)
##' ezMultiplicated(v1)
##' ezMultiplicated(v1,3)
##' ezMultiplicated(v1,2,"all")
ezMultiplicated = function(x, n=2, mode="keepFirst"){
	stopifnot(n >= 1)
	if (n == 1){
	  return(rep(TRUE, length(x)))
	}
	if (mode == "all"){
    return (x %in% unique(x[ezMultiplicated(x, n=n, mode="keepFirst")]))
  }
  idx = switch(mode,
             keepFirst=1:length(x),
             random=sample(1:length(x), length(x), replace=FALSE),
             keepLast=length(x):1)
  x = x[idx]
  isMulti = ezGetMulti(x) >= n
  isMulti[idx] = isMulti
	return(isMulti)
}

##' @describeIn ezMultiplicated Returns values incrementing seperately for each occurance of each value in the input.
ezGetMulti = function(x){
  #x = c(1,1,3,3,3,3, 2,2)
  idx = unsplit(tapply(x, x, function(y){1:length(y)}), x)
}

##' @title Collapses a set of values
##' @description Collapses a set of values into a character and pastes them. Useful for documentation purposes.
##' @param x a vector, matrix or list.
##' @param sep the separator to use between values.
##' @param na.rm a logical specifying whether to remove \code{NA}'s.
##' @param empty.rm a logical specifying whether to remove empty values.
##' @param uniqueOnly a logical specifying whether to keep only unique values.
##' @return Returns the values collapsed into one character.
##' @template roxygen-template
##' @examples
##' l1 = list(a=c(1,"",6),c=c("rsrg","yjrt",NA,6))
##' collapse(l1,sep="_")
##' collapse(l1,na.rm=T,empty.rm=T,uniqueOnly=T)
collapse = function(x, sep="; ", na.rm=FALSE, empty.rm=FALSE, uniqueOnly=FALSE){
  if (length(x) == 0){
    return("")
  }
  x = unlist(x)
  if (na.rm){
    x = x[!is.na(x)]
  }
  if (empty.rm){
    x = x[x != ""];    
  }
  if (uniqueOnly){
    x = unique(x)
  }
  paste(x, collapse=sep)
}


