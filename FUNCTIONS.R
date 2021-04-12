##############################################################################
## R Functions for preparing data the Metabolon GWAS
##
## By: Laura Corbin, Caroline Bull, David Hughes
## Date begun: Feb 10 2020
##
##############################################################################

###############################
## Perform winzorization on
## a matrix of data. 
## Samples in rows.
## feautres in columns. 
###############################
winsorize_x = function(x, nsd=5) {
	cut_point_top = mean(x,na.rm=T) + (nsd*(sd(x,na.rm=T)))
	cut_point_bottom = mean(x,na.rm=T) - (nsd*(sd(x,na.rm=T)))
	i = which(x > cut_point_top)
	x[i] = cut_point_top
	j = which(x < cut_point_bottom)
	x[j] = cut_point_bottom
	return(x)
}


##################################
## Count the number of samples
## that had to be winzorized
## for each feature
##################################
winsorize_count = function(x, nsd=5) {
	cut_point_top = mean(x,na.rm=T) + (nsd*(sd(x,na.rm=T)))
	cut_point_bottom = mean(x,na.rm=T) - (nsd*(sd(x,na.rm=T)))
	i = which(x > cut_point_top)
	x[i] = cut_point_top
	j = which(x < cut_point_bottom)
	x[j] = cut_point_bottom
	return(length(c(i,j)))
}


################################
##
##  Z-transformation from the 
##	GENABEL package
##
################################
ztransform = function (formula, data, family = gaussian) 
{
  if (missing(data)) {
    if (is(formula, "formula")) 
      data <- environment(formula)
    else data <- environment()
  }
  else {
    if (is(data, "gwaa.data")) {
      data <- data@phdata
    }
    else if (!is(data, "data.frame")) {
      stop("data argument should be of gwaa.data or data.frame class")
    }
  }
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (is(try(formula, silent = TRUE), "try-error")) {
    formula <- data[[as(match.call()[["formula"]], "character")]]
  }
  if (is(formula, "formula")) {
    mf <- model.frame(formula, data, na.action = na.pass, 
                      drop.unused.levels = TRUE)
    mids <- complete.cases(mf)
    mf <- mf[mids, ]
    y <- model.response(mf)
    desmat <- model.matrix(formula, mf)
    lmf <- glm.fit(desmat, y, family = family)
    resid <- lmf$resid
  }
  else if (is(formula, "numeric") || is(formula, "integer") || 
           is(formula, "double")) {
    y <- formula
    mids <- (!is.na(y))
    y <- y[mids]
    resid <- y
    if (length(unique(resid)) == 1) 
      stop("trait is monomorphic")
    if (length(unique(resid)) == 2) 
      stop("trait is binary")
  }
  else {
    stop("formula argument must be a formula or one of (numeric, integer, double)")
  }
  y <- (resid - mean(resid))/sd(resid)
  tmeas <- as.logical(mids)
  out <- rep(NA, length(mids))
  out[tmeas] <- y
  out
}

################################
##
##  RNT-transformation from the 
##	GENABEL package, modified 
##  to randomly split tied values
##
################################
rntransform = function (formula, data, family = gaussian, split_ties = TRUE) 
{
  if (is(try(formula, silent = TRUE), "try-error")) {
    if (is(data, "gwaa.data")) 
      data1 <- phdata(data)
    else if (is(data, "data.frame")) 
      data1 <- data
    else stop("'data' must have 'gwaa.data' or 'data.frame' class")
    formula <- data1[[as(match.call()[["formula"]], "character")]]
  }
  var <- ztransform(formula, data, family)
  if(split_ties == TRUE){
    out <- rank(var, ties.method = "random") - 0.5
  } else {
    out <- rank(var) - 0.5
  }
  out[is.na(var)] <- NA
  mP <- 0.5/max(out, na.rm = T)
  out <- out/(max(out, na.rm = T) + 0.5)
  out <- qnorm(out)
  out
}


################################
##
##  Hypergeometric Test
##	
##  
##
################################
hgtest = function(allfeatures, selectedfeatures){
	## unique list of categories in the selected feature list
	cats = na.omit( unique(selectedfeatures) )
	## iterate over those cat and perform test
	
	HG_F_test = sapply(cats, function(i){
		g = length(allfeatures)
		## total number of tested features with some annotation
		N = length(allfeatures) - ( sum(is.na(allfeatures) | allfeatures == "") )
		## number of features in cat i
		m = sum(allfeatures == i, na.rm = TRUE ) 
		## number of features NOT in cat i
		n = N - m
		## number of features selected by X criteria, that have some annotation
		k = length(selectedfeatures) - ( sum(is.na(selectedfeatures) | selectedfeatures == "") )
		## number fo selected features in cat i
		x = sum(selectedfeatures == i , na.rm = TRUE) 

		
		## estiamte fold enrichment and p.value
		fold.enrichment <-  (x / k ) / (m / N)
		p.value <-  phyper(q=x -1, m=m, n=n, k=k, lower.tail=FALSE)

		## FISHER EXACT TEST
		## BUILD 2 x 2 contigency table
		dmat = matrix( c(x,k-x,m-x,n-(k-x)),2,2, byrow  = TRUE, dimnames = list(c("Selected","NotSelected"), c("FocusedCat","AllOtherCats") )) 
		ftest <- fisher.test(x=dmat, alternative="greater")


		## data out
		out = c(fold.enrichment, p.value, ftest$estimate, ftest$p.value); names(out) = c("fold.enrichment", "HG_pval", "Ftest_OR", "Ftest_pval")
		return(out)
		})
	## return to user
	return(t(HG_F_test))	

} 

