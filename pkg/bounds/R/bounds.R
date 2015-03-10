#' Bounds for causal interactions
#' 
#' The function estimates bounds for causal interaction.
#' 
#' @param data A string as file name for the input data, or an R data frame giving the data.
#' @param response A string that tells the name of the response variable in the data. Default value is 
#' "response". The corresponding column of this variable needs to be coded as 0 and 1 numeric values.
#' @param factor1 A string that tells the name of the first exposure variable in the data. Default value is 
#' "exposure1". The corresponding column of this variable needs to be factors. See \code{as.factor}.
#' @param factor2 A string that tells the name of the second exposure variable in the data. Default value is 
#' "exposure2". The corresponding column of this variable needs to be factors. See \code{as.factor}.
#' @param covariates A vector of strings that give the names of the covariates in the data. 
#' @param population.prevalence A numeric value that gives the population prevalence of response incidence.
#' If \code{NULL}, no weighting for population prevalence.
#' @param monotonicity Logical. If \code{FALSE}, then assumption-free bounds are returned. If \code{TRUE},
#' then bounds under the assumption of monotone exposure effects are returned.
#' @param bootstrap Optional. The number of bootstrap replicates. If missing, 
#' then no estimated standard errors are returned.
#' @param cc Logical. Specify \code{TRUE} if the response is from a case-control study, and \code{population.prevalence}
#' is required.
#' 
#' @return The function returns returns a 4 × 5 matrix. The upper two elements in first column contain 
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
#' Sjolander A, Lee W, Kallberg H, Pawitan Y. (2014). Bounds on causal interactions for binary outcomes. \emph{Biometrics} \bold{70}(3), 500-505.
#' 
#' @seealso 
#' \code{\link{bounds.suff}}
#' 
#' @examples 
#' \dontrun{
#' data(ex1)
#' bounds(ex1, covariates = c('C1', 'C2'), bootstrap = 500)
#' bounds(ex1, covariates = c('C1', 'C2'), monotonicity = TRUE, bootstrap = 500)
#' bounds(ex1, covariates = c('C1', 'C2'), population.prevalence = 0.001, 
#'        monotonicity = TRUE, bootstrap = 500)
#' }
#' @aliases bounds
#' @keywords bounds

`bounds` <- function(data, response = 'outcome', factor1 = 'exposure1', factor2 = 'exposure2', covariates = NULL, 
		population.prevalence = NULL, monotonicity = FALSE, bootstrap, cc = FALSE) {
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
	arvid <- bounds.arvid(as.formula(form), X = factor1, Z = factor2, data = dd, 
			monotonicity = monotonicity, weights = w, R = bootstrap, cc = cc)
	arvid[1,4] <- arvid[2,3] <- arvid[3,4] <- arvid[4,3] <- NA
	return(arvid)
}



