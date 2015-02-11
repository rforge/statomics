#' Exact p-values for genome-wide association analysis in inbred populations
#' 
#' The function imports GenABEL (gwaa.data class) data format
#' and calculates the exact dataset-specific p-values for each variant
#' or a given effect size and allele frequency.
#' 
#' @param gwaa.object An (optional) object of \code{\link{gwaa.data-class}}.
#' @param n An (optional) integer gives the sample size, only used when \code{gwaa.object = NULL}.
#' @param maf An (optional) vector gives minor allele frequencies across the genome, 
#' only used when \code{gwaa.object = NULL}.
#' @param beta An vector gives the effect sizes to be tested, and \code{beta.maf} must
#' also be given if not testing for the whole genome. 
#' @param beta.maf An (optional) vector gives the minor allele frequencies in accordance with
#' \code{beta}, must be used together with \code{beta}. If \code{beta.maf} 
#' is missing, exact p-values will be calculated for the whole genome.
#' @param low.maf A numeric value that givens the cut-off of the lowest minor allele frequency
#' allowed in the analysis.
#' @param type A string tells the statistical test type, can be \code{'one-sided'} or \code{'two-sided'}.
#' 
#' @note None.
#' 
#' @return The function returns a data frame of exact p-values (\code{$p.exact}) and corresponding effect
#' sizes \code{$beta} and MAFs \code{$MAF}.
#' 
#' @author Xia Shen
#' 
#' @references 
#' Xia Shen (2015). Beyond permutation test: calculating exact dataset-specific p-values 
#' for genome-wide association studies in inbred populations. \emph{Submitted}.
#' 
#' @seealso 
#' \code{\link{qtscore}}, \code{\link{t.test}}
#' 
#' @examples 
#' \dontrun{
#' ## loading example gwaa.data of data from Atwell et al. (2010) Nature
#' data(arab)
#' 
#' ## running a regular GWA analysis for height
#' y <- qnorm((rank(arab@phdata$X1_LD, na.last = "keep") - 0.5)/sum(!is.na(arab@phdata$X1_LD)))
#' qt1 <- qtscore(y, data = arab)
#' 
#' ## running the multivariate GWAS again
#' exact <- p.exact.gaussian(arab, beta = qt1[,'effB'])
#' }
#' @aliases p.exact.gaussian
#' @keywords p.exact, exact p-value, genome-wide false discovery rate
#' 
`p.exact.gaussian` <- function(gwaa.object = NULL, n = NULL, maf = NULL, beta, beta.maf = NULL, low.maf = .05, type = 'two-sided') {
	if (all(c(is.null(gwaa.object), is.null(n), is.null(maf)))) {
		stop('insufficient input!')
	} else if (is.null(gwaa.object) & any(is.null(c(n, maf)))) {
		stop('insufficient input!')
	} else if (all(c(!is.null(gwaa.object), !is.null(n), !is.null(maf)))) {
		cat('over-sufficient input provided, only n & maf in use.\n')
	}
	if (!is.null(gwaa.object)) {
		cat('processing gwaa.object ...\n')
		Z <- as.double(gwaa.object)
		n <- nrow(Z)
		n2 <- colSums(Z == 2)
		n0 <- colSums(Z == 0)
		minor <- as.numeric(n2 < n0)
		allfreq <- minor*n2/n + (1 - minor)*n0/n
	}
	if (is.null(gwaa.object) & !is.null(n) & !is.null(maf)) {
		allfreq <- maf
	}
	k <- table(allfreq)
	freq <- as.numeric(names(k))
	cat('MAF distribution obtained:\n')
	cat('lowest MAF =', freq[which(freq >= low.maf)[1]], 'count =', k[which(freq >= low.maf)[1]], '\n')
	cat('highest MAF =', freq[length(k)], 'count =', k[length(k)], '\n')
	if (!is.null(beta) & is.null(beta.maf) & length(beta) == length(allfreq)) {
		beta.maf <- allfreq
	}
	if (is.null(gwaa.object) & (is.null(beta) | is.null(beta.maf))) {
		stop('either beta or beta.maf is missing!')
	}
	if (is.null(gwaa.object) & length(beta) != length(beta.maf)) {
		stop('beta and beta.maf have different lengths!')
	}
	if (!is.null(beta) & !is.null(beta.maf) & length(beta) == length(beta.maf)) {
		cat('calculating exact p-values ...\n')
		p <- q <- rep(NA, length(beta))
		for (i in 1:length(beta)) {
			if (beta.maf[i] >= low.maf) {
				b <- abs(beta[i])
				freqidx <- which(freq >= beta.maf[i] - 1e-8)
				freqi <- freq[freqidx]
				if (type == 'one-sided') {
					logpj <- log(pnorm(b, 0, sqrt(1/n/freqi + 1/n/(1 - freqi))))
				} else if (type == 'two-sided') {
					logpj <- log(pchisq(b**2/(1/n/freqi + 1/n/(1 - freqi)), 1))
				} else {
					stop('wrong type of test specified!')
				}
				q[i] <- exp(sum(logpj*k[freqidx]))
				p[i] <- 1 - q[i]
			} else {
				q[i] <- 0
				p[i] <- 1
			}
			progress(i/length(beta)*100)
		}
	}
	res <- data.frame(beta = beta, MAF = beta.maf, p.exact = p)
	if (nrow(res) == length(allfreq) & !is.null(gwaa.object)) {
		rownames(res) <- gtdata(gwaa.object)@snpnames
	}
	return(res)
}
