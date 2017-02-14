gchc <-
function(rg1.log = NULL, 
                rg2.log = NULL, 
                exposure.sumstats = NULL, 
                outcome.sumstats1 = NULL,
                outcome.sumstats2 = NULL,
                exposure = 'Exposure', 
                outcome = 'Outcome')
{
    ## startup
    cat('*********************************************************************\n')
    cat('* Genetic Correlation Heterogeneity for Causality (GCHC)\n')
    cat('* Version 1.0.0\n')
    cat('* (C) 2017 Xia Shen\n')
    cat('* Usher Institute of Univ. Edinburgh / MEB of Karolinska Institutet\n')
    cat('* GNU General Public License v3\n')
    cat('*********************************************************************\n\n')
    cc <- match.call()
    cat('Call:\n')
    print(cc)
    cat('\n')
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
