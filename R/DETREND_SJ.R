DETREND_SJ <-
function(x, tt = 'linear', bp = c()) {

        
        if (!is.numeric(x) && !is.complex(x))
                stop("'x' must be a numeric or complex vector or matrix.")
        trendType <- pmatch(tt, c('constant', 'linear'), nomatch = 0)
        
        if (is.vector(x))
                x <- as.matrix(x)
        n <- nrow(x)
        if (length(bp) > 0 && !all(bp %in% 1:n))
                stop("Breakpoints 'bp' must elements of 1:length(x).")
        
        if (trendType == 1) {  # 'constant'
                if (!is.null(bp))
                        warning("Breakpoints not used for 'constant' trend type.")
                y <- x - matrix(1, n, 1) %*% apply(x, 2, mean)
                
        } else if (trendType == 2) {  # 'linear'
                
                if(length(bp)==0) {
                        bp <- 1
                }
                
                bp <- sort(unique(c(1, c(bp), n))) ## both ends
                lbp <- length(bp)
                
                lb <- length(bp) - 1
                
                a <- matrix(0, n, lbp)
                a[1:n,1] <- as.matrix(1:n)/n
                if(lbp>2) {
                        ktmp <- 2
                        
                        for(kt in 2:(lbp-1)) {
                                m <- n - bp[ktmp]
                                a[(bp[ktmp]+1):n,ktmp] <- as.matrix(1:m)/m
                                ktmp <- ktmp +1
                        }
                
                }
                a[1:n,lbp] <- 1
                y <- x - a %*% mldivide(a, x)
                
                A <- matrix(0,n,lbp)
                
                if(lbp>1) {
                        bp[length(bp)] <- n-1 
                        for(k in 1:(lbp-1)) {
                                A[(bp[k]+1):bp[k+1],k] <- 1:(bp[k+1]-bp[k])
                                A[(bp[k+1]+1):n,k] <- A[bp[k+1],k]
                        }  
                }

                
                A[1:n,lbp] <- 1
                
                z <- qr(A)
                Q <- qr.Q(z)
                R <- qr.R(z)
                
                ##p <- qr.solve(R, (t(Q)%*%x) ) # changed last
                b2 <- (t(Q)%*%x)
                p <- mldivide(R, b2, pinv = TRUE)
                
                #r <- x - A%*%p
                param_no <- length(p)
                
                #H <- A %*% solve( t(A) %*% A) %*% t(A)
                H <- A %*% pinv( t(A) %*% A) %*% t(A)
                
                dfSJ <- length(x) - 1.25 * sum(diag(H)) + 0.5
                normr <- max(svd(y)$d)
                
        } else {
                stop("Trend type 'tt' must be 'constant' or 'linear'.")
        }
        
        
        DTSJ.values <- list(
                "y" = y,                
                "param_no" = param_no,
                "p" = p,
                "R" = R,
                "dfSJ" = dfSJ,
                "normr" = normr
        )
        
        class(DTSJ.values) <- "DTSJ"
        return(DTSJ.values)
}
