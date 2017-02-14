#' Genetic Correlation Heterogeneity for Causality (GCHC)
#' 
#' This function tests the difference between genetic correlation estimates 
#' reported by LD Score regression (LDSC), to infer causality.
#' 
#' @param rg1.log LDSC log file for the 1st genetic correlation estimation between the exposure and outcome phenotypes.
#' @param rg2.log LDSC log file for the 2nd genetic correlation estimation between the exposure and outcome phenotypes.
#' @param exposure.sumstats The gzipped file of exposure phenotype GWAS summary statistics, munged by LDSC.
#' @param outcome.sumstats1 The gzipped file of outcome phenotype GWAS (1st population) summary statistics, munged by LDSC.
#' @param outcome.sumstats2 The gzipped file of outcome phenotype GWAS (2nd population) summary statistics, munged by LDSC.
#' @param exposure Name of the exposure phenotype.
#' @param outcome Name of the outcome phenotype.
#' @param heading logical value that specified whether the software intro heading is printed.
#' 
#' @note GWAS of the exposure phenotype is required to be done in only one population, to keep its heritability fixed. 
#' Two different populations/sources are required for GWAS of the outcome phenotype.
#' 
#' @return A list is returned with:
#' \itemize{
#' \item{rg.diff }{The estimated difference between two genetic correlation estimates.}
#' \item{se }{The standard error of \code{rg.diff}.}
#' \item{p.value }{The p-value testing the null hypothesis of \code{rg.diff} = 0.}
#' \item{r.rg }{The estimated correlation between two genetic correlation estimates.}
#' }
#' 
#' @author Xia Shen
#' 
#' @references 
#' Shen X, Ning Z, Joshi PK, Lee Y, Wilson JF, Pawitan Y (2017). Genetic correlation heterogeneity detects causal factors
#' for complex traits. \emph{Submitted}.
#' 
#' @seealso 
#' GCHC homepage: http://gchc.shen.se
#' 
#' @examples 
#' \dontrun{
#' gchc(rg1.log = 'EDU_BMI1.log', 
#'                 rg2.log = 'EDU_BMI2.log', 
#'                 exposure.sumstats = 'EDU.sumstats.gz', 
#'                 outcome.sumstats1 = 'BMI1.sumstats.gz',
#'                 outcome.sumstats2 = 'BMI2.sumstats.gz',
#'                 exposure = 'EA', 
#'                 outcome = 'BMI')
#' }
#' @aliases gchc
#' @keywords causal inference, genetic correlation
#' 
gchc <-
function(rg1.log = NULL, 
                rg2.log = NULL, 
                exposure.sumstats = NULL, 
                outcome.sumstats1 = NULL,
                outcome.sumstats2 = NULL,
                exposure = 'Exposure', 
                outcome = 'Outcome',
				heading = TRUE)
{
    ## startup
	if (heading) {
	    print.heading()
    	cc <- match.call()
    	cat('Call:\n')
    	print(cc)
    	cat('\n')
	}
    t0 <- as.numeric(proc.time()[3])
    ## check input
    if (any(is.null(c(rg1.log, rg2.log)))) {
        stop('Insufficient log files from LDSC!')
    }
    if (any(is.null(c(exposure.sumstats, outcome.sumstats1, outcome.sumstats2)))) {
        indep <- TRUE
    } else {
        indep <- FALSE
    }
    ## load rg
    if (file.exists(rg1.log)) {
		out <- readLines(rg1.log)
		rgidx <- grep(pattern = 'Genetic Correlation:', out)
		if (length(rgidx) > 0) {
			rgline <- out[rgidx]
			spl <- unlist(strsplit(rgline, ' '))
			if (spl[3] != 'nan') {
				rg1 <- as.numeric(spl[3])
				len <- length(unlist(strsplit(spl[4], '')))
				serg1 <- as.numeric(substr(spl[4], 2, len - 1))
			}
		}
    } else {
        stop('rg1.log file not found!')
    }
    if (file.exists(rg2.log)) {
		out <- readLines(rg2.log)
		rgidx <- grep(pattern = 'Genetic Correlation:', out)
		if (length(rgidx) > 0) {
			rgline <- out[rgidx]
			spl <- unlist(strsplit(rgline, ' '))
			if (spl[3] != 'nan') {
				rg2 <- as.numeric(spl[3])
				len <- length(unlist(strsplit(spl[4], '')))
				serg2 <- as.numeric(substr(spl[4], 2, len - 1))
			}
		}
    } else {
        stop('rg2.log file not found!')
    }
    cat('Genetic correlation estimates loaded.\n')
    ## calculate parameters
    b <- rg1 - rg2
    if (indep) {
        s <- sqrt(serg1**2 + serg2**2)
        rzz <- 0
        cat('No munged summary statistics given, genetic correlations are assumed independent.\n')
    } else {
        c11 <- read.sumstats(outcome.sumstats1)
        c21 <- read.sumstats(outcome.sumstats2)
        ss <- read.sumstats(exposure.sumstats)
        snps <- ss$SNP[ss$SNP %in% c11$SNP & ss$SNP %in% c21$SNP]
        cat(length(snps), 'variants to estimate the correlation between genetic correlation estimates.\n')
        zz1 <- ss[snps,'Z']*c11[snps,'Z']*(1 - 2*(ss[snps,'A1'] != c11[snps,'A1']))
        zz2 <- ss[snps,'Z']*c21[snps,'Z']*(1 - 2*(ss[snps,'A1'] != c21[snps,'A1']))
        rzz <- cor(zz1, zz2, use = 'pairwise.complete.obs')
        cat('Estimated correlation between genetic correlation estimates:', rzz, '\n')
        s <- sqrt((serg1**2 + serg2**2 - 2*rzz*serg1*serg2)) 
    }
    ## Wald test
    wald <- b/s
    pv <- pchisq(wald**2, 1, lower.tail = FALSE)
    t1 <- as.numeric(proc.time()[3])
    cat('\n')
    cat('Results\n')
    cat('-------\n')
    cat('Difference in genetic correlation:', round(b, digits = 3), '\n')
    cat('Standard error:', round(s, digits = 3), '\n')
    cat('P:', pv, '*\n')
    cat('* Small P suggests causality exists from', exposure, 'to', outcome, '\n\n')
    cat('Analysis finished at', date(), '\n')
    cat('Total time elapsed:', t1 - t0, 's\n')
    return(list(rg.diff = b, se = s, p.value = pv, r.rg = rzz))
}
