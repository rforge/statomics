#' Exact p-values for genome-wide association analysis of case-control data in inbred populations
#' 
#' The function imports GenABEL (gwaa.data class) data format
#' and calculates the exact dataset-specific p-values of a case-control
#' phenotype for each variant or a given odds ratio and allele frequency.
#' 
#' @param pheno A string that gives the binary phenotype name in \code{gwaa.object} or
#' a vector that gives phenotypic values match the order in \code{gwaa.object}.
#' @param gwaa.object An object of \code{\link{gwaa.data-class}} to be analyzed.
#' @param or An (optional) vector gives the odds ratios to be tested, and \code{or.maf} must
#' also be given if not testing for the whole genome. 
#' @param or.maf An (optional) vector gives the minor allele frequencies in accordance with
#' \code{or}, must be used together with \code{or}. If \code{or.maf} or \code{or} 
#' is missing, exact p-values will be calculated for the whole genome.
#' @param low.maf A numeric value that gives the cut-off of the lowest minor allele frequency
#' allowed in the analysis.
#' @param high.ld A numeric value that gives the cut-off of the highest linkage disequilibrium R-square
#' allowed in the analysis, i.e. LD pruning.
#' @param method A string tells the method used, currently only "logOR" is available.
#' @param type A string tells the statistical test type, can be \code{'one-sided'} or \code{'two-sided'}.
#' @param con.table Genome-wide contingency tables. If \code{NULL}, to be calculated based on data.
#' If a string tells the file name, load from the file or to be calculated and saved into the file. If an R
#' matrix, regard as genome-wide contingency tables (see the saved R object for the format).
#' 
#' @note None.
#' 
#' @return The function returns a data frame of exact p-values (\code{$p.exact}) and corresponding odds
#' ratios \code{$or} and MAFs \code{$MAF}.
#' 
#' @author Xia Shen
#' 
#' @references 
#' Xia Shen (2015). Beyond permutation test: calculating exact dataset-specific p-values 
#' for genome-wide association studies in inbred populations. \emph{Submitted}.
#' 
#' @seealso 
#' \code{\link{ccfast}}, \code{\link{glm}}
#' 
#' @examples 
#' \dontrun{
#' ## loading example gwaa.data of data from Atwell et al. (2010) Nature
#' data(arab)
#' 
#' ## running a regular GWA analysis for AvrRPM1
#' cc1 <- ccfast('X33_.i.avrRpm1..i.', data = arab)
#' 
#' ## check the top finding using the exact p-value
#' top <- which.min(cc1[,'P1df'])
#' ctab <- table(arab@phdata$X33_.i.avrRpm1..i., as.double(arab[,top]))
#' ctab[ctab == 0] <- ctab[ctab == 0] + .5
#' or <- ctab[1,1]*ctab[2,2]/ctab[2,1]/ctab[1,2]
#' f <- summary(arab[,top])$Q.2
#' maf <- min(f, 1 - f)
#' exact <- p.exact.binary(pheno = 'X33_.i.avrRpm1..i.', gwaa.object = arab, 
#'          or = or, or.maf = maf, con.table = 'tab.AvrRPM1.RData')
#' }
#' @aliases p.exact.binary
#' @keywords p.exact, exact p-value, genome-wide false discovery rate
#' 

`p.exact.binary` <- function(pheno, gwaa.object, or = NULL, or.maf = NULL, low.maf = .05, high.ld = .9, method = 'logOR', type = 'two-sided', con.table = NULL) {
	if (is.character(pheno)) {
		y <- gwaa.object@phdata[,pheno]
	} else {
		y <- pheno
	}
	naidx <- c()
	if (any(is.na(y))) {
		naidx <- which(is.na(y))
		y <- y[-naidx]
	}
	if (!all(names(table(y)) == c('0', '1'))) {
		stop('binary phenotype needed and coded as 0, 1!\n')
	}
	prev <- mean(y)
	if (!is.null(gwaa.object)) {
		cat('processing gwaa.object ...\n')
		if (length(naidx) > 0) {
			Z <- as.double(gwaa.object)[-naidx,]
		} else {
			Z <- as.double(gwaa.object)
		}
		if (high.ld < 1) {
			cat('LD pruning ...\n')
			r2 <- rep(NA, ncol(Z) - 1)
			for (i in 1:length(r2)) {
				r2[i] <- cor(Z[,i], Z[,i + 1])**2
				progress(i/length(r2)*100)
			}
			cat('\n')
			bad <- which(r2 > high.ld)
			good <- which(!(1:ncol(Z) %in% bad))
			Z <- Z[,good]
			gwaa.object <- gwaa.object[,good]
		}
		n <- nrow(Z)
		n2 <- colSums(Z == 2)
		n0 <- colSums(Z == 0)
		minor <- as.numeric(n2 < n0)
		allfreq <- minor*n2/n + (1 - minor)*n0/n
	}
	k <- table(allfreq)
	freq <- as.numeric(names(k))
	cat('MAF distribution obtained:\n')
	cat('lowest MAF =', freq[which(freq >= low.maf)[1]], 'count =', k[which(freq >= low.maf)[1]], '\n')
	cat('highest MAF =', freq[length(k)], 'count =', k[length(k)], '\n')
	cat('calculating genome-wide odds ratios ...\n')
	cc1 <- ccfast(pheno, gwaa.object, quiet = TRUE)
	allor <- cc1[,'effB']
	idx <- which(allor %in% c(0, 10000))
	for (j in idx) {
		tabj <- c(table(y, Z[,j]))
		tabj[tabj == 0] <- tabj[tabj == 0] + .5
		allor[j] <- tabj[1]*tabj[4] - tabj[2]*tabj[3]
	}
	if (is.null(or) | is.null(or.maf)) {
		or <- allor
		or.maf <- allfreq
	}
	n1 <- round(n*prev)
	n2 <- n - n1
	if (!is.null(or) & !is.null(or.maf) & length(or) == length(or.maf)) {
		p <- q <- rep(NA, length(or))
		if (method != 'logOR') {
			stop('only logOR method is supported currently!')
			cat('calculating exact p-values ...\n')
			for (i in 1:length(or)) {
				if (or.maf[i] >= low.maf) {
					smallor <- or[i] < 1
					freqidx <- which(freq >= or.maf[i] - 1e-8)
					freqi <- freq[freqidx]
					nchoose <- round(n*freqi)
					logpj <- rep(NA, length(freqi))
					for (fi in 1:length(freqi)) {
						orfi <- rep(NA, nchoose[fi])
						for (a in 1:nchoose[fi]) {
							orfi[a] <- a*(n1 - round(n*freqi[fi]) + a)/(round(n*freqi[fi]) - a)/(n2 - a)
						}
						if (smallor) {
							ai <- max((1:nchoose[fi])[orfi < or[i]])
						} else {
							ai <- min((1:nchoose[fi])[orfi > or[i]])
						}
						logpj[fi] <- log(1 - min(1, phyper(ai + .01*(smallor*2 - 1), n1, n2, nchoose[fi], lower.tail = smallor)*(1 + two.sided)))
					}
					q[i] <- exp(sum(logpj*k[freqidx]))
					p[i] <- 1 - q[i]
				} else {
					q[i] <- 0
					p[i] <- 1
				}
				progress(i/length(or)*100)
			}
		} else {
			cat('processing contingency tables ...\n')
			if (is.null(con.table)) {
				con.table <- matrix(NA, length(allor), 4)
				for (j in 1:length(allor)) {
					con.table[j,] <- c(table(y, Z[,j]))
					progress(j/length(allor)*100)
				}
				cat('\n')
			}
			if (is.character(con.table)) {
				fn <- con.table
				if (file.exists(fn)) {
					load(fn)
				} else {
					con.table <- matrix(NA, length(allor), 4)
					for (j in 1:length(allor)) {
						con.table[j,] <- c(table(y, Z[,j]))
						progress(j/length(allor)*100)
					}
					save(con.table, file = fn)
					cat('\n')
				}
			}
			con.table[con.table == 0] <- con.table[con.table == 0] + .5
			vars <- rowSums(1/con.table)
			cat('calculating exact p-values ...\n')
			for (i in 1:length(or)) {
				if (or.maf[i] >= low.maf) {
					b <- abs(log(or[i]))
					freqidx <- which(freq >= or.maf[i] - 1e-8)
					freqi <- freq[freqidx]
					if (type == 'one-sided') {
						logp <- log(pnorm(b, 0, sqrt(vars[allfreq >= min(freqi) - 1e-8])))
					} else if (type == 'two-sided') {
						logp <- log(1 - pchisq(b**2/vars[allfreq >= min(freqi) - 1e-8], 1, lower.tail = FALSE))
					} else {
						stop('wrong type of test specified!')
					}
					q[i] <- exp(sum(logp))
					p[i] <- 1 - q[i]
				} else {
					q[i] <- 0
					p[i] <- 1
				}
				progress(i/length(or)*100)
			}
		}
	}
	cat('\n')
	res <- data.frame(or = or, MAF = or.maf, p.exact = p)
	rownames(res) <- gtdata(gwaa.object)@snpnames
	return(res)
}
