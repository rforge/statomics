#' Bounds for sufficient-cause interactions
#' 
#' The function estimates bounds for sufficient-cause interaction.
#' 
#' @param data A string as file name for the input data, or an R data frame giving the data.
#' @param response A string that tells the name of the response variable in the data. Default value is 
#' "outcome". The corresponding column of this variable needs to be coded as 0 and 1 numeric values.
#' @param factor1 A string that tells the name of the first exposure variable in the data. Default value is 
#' "exposure1". The corresponding column of this variable needs to be factors. See \code{as.factor}.
#' @param level1 The level of \code{factor1} at which to estimate the interaction.
#' @param factor2 A string that tells the name of the second exposure variable in the data. Default value is 
#' "exposure2". The corresponding column of this variable needs to be factors. See \code{as.factor}.
#' @param level2 The level of \code{factor2} at which to estimate the interaction.
#' @param covariates A vector of strings that give the names of the covariates in the data. 
#' @param population.prevalence A numeric value that gives the population prevalence of response incidence.
#' If \code{NULL}, no weighting for population prevalence.
#' @param monotonicity Logical. If \code{FALSE}, then assumption-free bounds are returned. If \code{TRUE},
#' then bounds under the assumption of monotone exposure effects are returned.
#' @param weak Logical. If \code{FALSE}, then strong sufficient-cause interaction is considered. 
#' If \code{TRUE}, then weak sufficient-cause interaction is considered.
#' @param bootstrap Optional. The number of bootstrap replicates. If missing, 
#' then no estimated standard errors are returned.
#' @param cc Logical. Specify \code{TRUE} if the response is from a case-control study, and \code{population.prevalence}
#' is required.
#' 
#' @return The function returns returns a 4 × 4 matrix. The upper two elements in first column contain 
#' the estimated bounds on the proportion of individuals that the causal interaction exists, named as lower and upper, 
#' respectively. The lower two elements in the first column contain the estimated bounds on 
#' the proportion of effects that is due to interaction, named as p.lower and p.upper, respectively, which are calculated when \code{monotonicity = TRUE}. 
#' If \code{bootstrap} is specified, then the second column contains the bootstrap standard errors for the estimates, 
#' and the third and fourth columns contain the lower and upper limits of the 95% bootstrap percentile interval, respectively.
#' The last column contains the p-values from Wald test.
#' 
#' @author Arvid Sjolander, Xia Shen
#' 
#' @references 
#' Sjolander A, Lee W, Kallberg H, Pawitan Y. (2014). Bounds on sufficient-cause interactions. \emph{European Journal of Epidemiology} \bold{29}(11), 813-820. 
#' 
#' @seealso 
#' \code{\link{bounds}}
#' 
#' @examples 
#' \dontrun{
#' data(ex1)
#' bounds.suff(ex1, level1 = 1, level2 = 1, covariates = c('C1', 'C2'), bootstrap = 500)
#' bounds.suff(ex1, level1 = 1, level2 = 1, covariates = c('C1', 'C2'), 
#'             monotonicity = TRUE, bootstrap = 500)
#' bounds.suff(ex1, level1 = 1, level2 = 1, covariates = c('C1', 'C2'), 
#'             population.prevalence = 0.001, monotonicity = TRUE, bootstrap = 500)
#' }
#' @aliases bounds.suff
#' @keywords bounds, sufficient-cause

`bounds.suff` <- function(data, response = 'outcome', factor1 = 'exposure1', level1, factor2 = 'exposure1', level2, covariates = NULL, 
		population.prevalence = NULL, monotonicity = FALSE, weak = FALSE, bootstrap, cc = FALSE) {
	if (cc & is.null(population.prevalence)) {
		stop('population prevalence has to be specified for a case-control study!')
	}
	if (is.character(data)) {
		dd <- read.table(data, header = TRUE)
	} else {
		dd <- data
	}
	dd[,factor1] <- as.factor(dd[,factor1])
	dd[,factor2] <- as.factor(dd[,factor2])
	form <- paste(response, '~', factor1, '+', factor2, sep = '')
	if (!is.null(covariates)) {
		for (i in 1:length(covariates)) {
			form <- paste(form, '+', covariates[i], sep = '')
		}
	}
	w <- NULL
	if (!is.null(population.prevalence)) {
		pstar <- mean(dd[,response])
		w1 <- rep(population.prevalence/pstar, nrow(data))
		w0 <- rep((1 - population.prevalence)/(1 - pstar), nrow(data))
		w <- dd[,response]*w1 + (1 - dd[,response])*w0
	}
	arvid <- bounds.suff.arvid(as.formula(form), X = factor1, x = level1, Z = factor2, z = level2, data = dd, 
			monotonicity = monotonicity, weights = w, weak = weak, R = bootstrap, cc = cc)
	arvid[1,4] <- arvid[2,3] <- arvid[3,4] <- arvid[4,3] <- NA
	return(arvid)
}



