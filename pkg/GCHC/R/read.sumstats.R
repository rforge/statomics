read.sumstats <-
function(sumstats) {
    file <- gzfile(sumstats)
    x <- readLines(file)
    x <- gsub('\t\t\t\t', '\tNA\tNA\tNA\tNA', x)
    rnd <- as.character(round(runif(1)*1e6))
    writeLines(x, rnd)
    x <- read.table(rnd, header = TRUE)
    rownames(x) <- x$SNP
    file.remove(rnd)
    close(file)
    return(x)
}
