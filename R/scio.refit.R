scio.refit <- function(S, Omega, thr = 1e-04, pkg=c("QUIC", "glasso"), ...) {
    ## By default, Omega from the smallest lambda to the largest 
    pkg <- match.arg(pkg, c("QUIC", "glasso"))
    
    if(pkg == "glasso" && !requireNamespace("glasso",  quietly = TRUE))  {
        warning("Refitting did not run because package glasso is not available!")
        return(list(w=NULL))
    }
    if(pkg == "QUIC" && !requireNamespace("QUIC",  quietly = TRUE))  {
        warning("Refitting did not run because package QUIC is not available!")
        return(list(w=NULL))
    }
    
    HUGE <- 1e8
    ss <- dim(Omega)

    if (length(ss)>2) {
        w <- 0*Omega; o <- 0*Omega
        for (jj in ss[3]:1) {
            rho <- HUGE*(abs(Omega[,,jj])<thr) + thr
            diag(rho) <-  thr           # not penalize the diagonal
            ## Warm start 
            if (pkg=="glasso") {
                if (jj < ss[3] ) {
                    tmp <-  glasso::glasso(S, rho, start="warm", w.init=o[,,jj+1], wi.init=w[,,jj+1], ...)
                    w[,,jj] <- tmp$wi
                    o[,,jj] <- tmp$w
                } else {
                    tmp <-  glasso::glasso(S, rho, ...)
                    w[,,jj] <- tmp$wi
                    o[,,jj] <- tmp$w
                }
            } else {
                if (jj < ss[3]) {
                    tmp <-  QUIC::QUIC(S, rho, W.init=o[,,jj+1], X.init=w[,,jj+1], ...)
                    w[,,jj] <- tmp$X
                    o[,,jj] <- tmp$W
                } else {
                    tmp <-  QUIC::QUIC(S, rho, ...)
                    w[,,jj] <- tmp$X
                    o[,,jj] <- tmp$W
                }
            }
        }
    } else {
        rho <- HUGE*(abs(Omega)<thr) + thr
        diag(rho) <-  thr           # not penalize the diagonal
        if (pkg=="glasso") {
            w <- glasso::glasso(S, rho, ...)$wi
        } else {
            w <- QUIC::QUIC(S, rho, ...)$X
        } 
    }
    return(list(w=w))
}
