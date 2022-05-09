CM.pvalue <-
function (rw)
{
    a <- rw[1]
    e <- rw[2]
    f <- rw[3]
    g <- rw[4]
    h <- rw[5]
    J <- sum(rw) - a
    tmp <- c()
    m1 <- min(e + f, e + g)
    m2 <- max(e + f, e + g)
    for (astar in 0:J) {
        e <- astar:m1
        tmp <- c(tmp, sum((0.5)^e * choose(e, astar) * choose(m1,
            e) * choose(J - m1, m2 - e)/choose(J, m2)))
    }
    out <- NULL
    out <- sum(tmp[(0:J) >= a])
    return(out)
}

