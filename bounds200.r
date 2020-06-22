library(crossrun)

# temporary inclusion of functions bestbox and cutbox,
# before inclusion in the package crossrun
# taken directly from the script for 
# https://doi.org/10.1371/journal.pone.0233920
# with inclusion of Rmpfr:: before Rmpfr functions

bestbox <- function(pt0,
                    pts,
                    target = 0.925,
                    n1     = 100,
                    mult   = 2,
                    prec   = 120) {
  nill    <- Rmpfr::mpfr(0, prec)
  one     <- Rmpfr::mpfr(1, prec)
  two     <- Rmpfr::mpfr(2, prec)
  multm   <- Rmpfr::mpfr(mult, prec)
  targetm <- Rmpfr::mpfr(target, prec)
  targt   <- targetm * (multm ^ (n1 - 1)) # target on "times" scale
  pt0n    <- pt0[[n1]]
  ptsn    <- pts[[n1]]
  bpt0    <- boxprobt(pt0n)  # box probabilities for no shift
  bpttarg <- boxprobt(ptsn)  # box probabilities for target shift
  boxprt  <-
    two * (multm ^ (n1 - 1)) # initialize to impossible high value
  for (cc in 0:(n1 - 1))
    for (ll in 1:n1) {
      if (pt0n[cc + 1, ll] > nill &
          bpt0[cc + 1, ll] >= targt &
          bpttarg[cc + 1, ll] < boxprt) {
        c1 <- cc
        l1 <- ll
        boxprt <- bpttarg[cc + 1, ll]
      }
    }
  return(c(c1, l1))
} # end function bestbox

cutbox <- function(pt0,
                   pts,
                   target = 0.925,
                   n1     = 100,
                   c1     = 41,
                   l1     = 10,
                   mult   = 2,
                   prec   = 120) {
  nill      <- Rmpfr::mpfr(0, prec)
  one       <- Rmpfr::mpfr(1, prec)
  two       <- Rmpfr::mpfr(2, prec)
  multm     <- Rmpfr::mpfr(mult, prec)
  targetm   <- Rmpfr::mpfr(target, prec)
  targt     <- targetm * (multm ^ (n1 - 1)) # target on "times" scale
  pt0n      <- pt0[[n1]]
  ptsn      <- pts[[n1]]
  bpt0      <-
    boxprobt(pt0n)   # box probabilities for no shift, pt scale
  boxpt0    <-
    bpt0[c1 + 1, l1] # no shift probability of actual box, pt scale
  cornerpt0 <-
    pt0n[c1 + 1, l1] # no shift corner probability, pt scale
  finished  <- FALSE
  cbord     <- NA
  lbord     <- NA
  if (boxpt0 - cornerpt0 >= targt) {
    cutboxpt0 <-
      boxpt0 - cornerpt0 # pt of cutted box after removed corner
    cbord     <- c1 + 1
    lbord     <- l1 - 1
    while (!finished) {
      pt0n.directionc <- pt0n[cbord + 1, l1]
      pt0n.directionl <- pt0n[c1 + 1, lbord]
      ptsn.directionc <- ptsn[cbord + 1, l1]
      ptsn.directionl <- ptsn[c1 + 1, lbord]
      if ((cutboxpt0 - pt0n.directionc < targt |
           pt0n.directionc == 0) &
          (cutboxpt0 - pt0n.directionl < targt |
           pt0n.directionl == 0)) {
        finished <- TRUE
      } else if (cutboxpt0 - pt0n.directionc < targt |
                 pt0n.directionc == 0) {
        lstrip    <- pt0n[c1 + 1, lbord:1]
        nlstrip   <- length(lstrip)
        maxlstrip <- max((1:nlstrip)[lstrip > 0])
        lstrip    <- lstrip[(1:nlstrip) <= maxlstrip]
        lstripcum <- cumsum(lstrip)
        if (cutboxpt0 - max(lstripcum) >= targt) {
          lbord <- 0
        } else {
          # 0 cannot occurr
          lbord <-
            lbord + 1 - min((1:nlstrip)[cutboxpt0 - lstripcum < targt])
        }
        finished <- TRUE
      } else if (cutboxpt0 - pt0n.directionl < targt |
                 pt0n.directionl == 0) {
        cstrip    <- pt0n[(cbord + 1):n1, l1]
        ncstrip   <- length(cstrip)
        maxcstrip <- max((1:ncstrip)[cstrip > 0])
        cstrip    <- cstrip[(1:ncstrip) <= maxcstrip]
        cstripcum <- cumsum(cstrip)
        if (cutboxpt0 - max(cstripcum) >= targt) {
          cbord <- n1
        } else {
          # n1 cannot occurr
          cbord <-
            cbord + min((1:ncstrip)[cutboxpt0 - cstripcum < targt]) - 1
        }
        finished <- TRUE
      } else if (ptsn.directionc >= ptsn.directionl) {
        cbord     <- cbord + 1
        cutboxpt0 <- cutboxpt0 - pt0n.directionc
      } else if (ptsn.directionc < ptsn.directionl) {
        lbord     <- lbord - 1
        cutboxpt0 <- cutboxpt0 - pt0n.directionl
      }
    }
  }
  return(c(cbord, lbord))
} # end function cutbox

# computation of joint distributions up to 200,
# from script for https://doi.org/10.1371/journal.pone.0223233
# with time information from the actual run
starttime200symm.p240 <- Sys.time()
cr200symm.p240 <- crossrunsymm(nmax=200, printn=TRUE, prec=240)$pt
endtime200symm.p240 <- Sys.time()
endtime200symm.p240 - starttime200symm.p240 # about 21 minutes
# probability 0.8, precision 240:
starttime200.8.p240 <- Sys.time()
cr200.8.p240 <- crossrunbin(nmax=200, prob=.8, printn=TRUE, prec=240)$pt
endtime200.8.p240 <- Sys.time()
endtime200.8.p240 - starttime200.8.p240 # almost 1 hour (other computer)

# bounds for Anhøj rules
bounds <- data.frame(n=10:200, row.names=10:200)
bounds$ca <- qbinom(p=.05, size=bounds$n-1, prob=.5)
bounds$la <- round(log2(bounds$n)+3)

# bounds for bestbox rules (code adapted from the "crossrunbox" article)
bounds$cb <- NA
bounds$lb <- NA

for (nn in 10:200) {
  print(nn)
  bounds[bounds$n == nn, c('cb', 'lb')] <- 
  bestbox(pt0=cr200symm.p240, pts=cr200.8.p240, 
          n1 = nn, target = .925, prec=240)
  }

# find cut  boxes
bounds$cbord <- NA
bounds$lbord <- NA

for (nn in 10:200) {
  print(paste('cutbox', nn))
  bounds[bounds$n == nn, c('cbord', 'lbord')] <-
    cutbox(pt0=cr200symm.p240, pts=cr200.8.p240,
           n1 = nn,
           target = .925,
           prec=240,
           c1 = bounds$cb[bounds$n == nn],
           l1 = bounds$lb[bounds$n == nn])
}

# joint distributions with default precision, up to n=100
starttime100symm.p120 <- Sys.time()
cr100symm.p120 <- crossrunsymm(nmax=100, printn=TRUE, prec=120)$pt
endtime100symm.p120 <- Sys.time()
endtime100symm.p120 - starttime100symm.p120
# probability 0.8, precision 120:
starttime100.8.p120 <- Sys.time()
cr100.8.p120 <- crossrunbin(nmax=100, prob=.8, printn=TRUE, prec=120)$pt
endtime100.8.p120 <- Sys.time()
endtime100.8.p120 - starttime100.8.p120 

# c omputations as in https://doi.org/10.1371/journal.pone.0233920 ,
# with precision 120
bounds100 <- data.frame(n=10:100, row.names=10:100)
bounds100$ca <- qbinom(p=.05, size=bounds100$n-1, prob=.5)
bounds100$la <- round(log2(bounds100$n)+3)

# bounds for bestbox rules (code adapted from the "crossrunbox" article),
# up to n=100, default precision
bounds100$cb <- NA
bounds100$lb <- NA

for (nn in 10:100) {
  print(nn)
  bounds100[bounds100$n == nn, c('cb', 'lb')] <- 
    bestbox(pt0=cr100symm.p120, pts=cr100.8.p120, 
            n1 = nn, target = .925)
}

# find cut  boxes, n=100, default precision
bounds100$cbord <- NA
bounds100$lbord <- NA

for (nn in 10:100) {
  print(paste('cutbox', nn))
  bounds100[bounds100$n == nn, c('cbord', 'lbord')] <-
    cutbox(pt0=cr100symm.p120, pts=cr100.8.p120,
           n1 = nn,
           target = .925,
           c1 = bounds100$cb[bounds100$n == nn],
           l1 = bounds100$lb[bounds100$n == nn])
}

summary(bounds[bounds$n<=100,1:7] - bounds100)
# exactly the same

# closer investigations for the discrepancies with
# https://doi.org/10.1371/journal.pone.0233920

# n=35
# bestbox probabilities, previous bestbox
sum(cr200symm.p240[[35]][(12+1):35,1:8])/sum(cr200symm.p240[[35]]) # symmetry
sum(cr200.8.p240[[35]][(12+1):35,1:8])/sum(cr200.8.p240[[35]]) # shift 0.8
# bestbox probabilities, bestbox computed here
sum(cr200symm.p240[[35]][(13+1):35,1:10])/sum(cr200symm.p240[[35]]) # symmetry
sum(cr200.8.p240[[35]][(13+1):35,1:10])/sum(cr200.8.p240[[35]]) # shift 0.8
# in both cases above target with symmetry, while a bit lower box probability
# as computed here for shift 0.8

# n=43
# bestbox probabilities, previous bestbox
sum(cr200symm.p240[[43]][(14+1):43,1:8])/sum(cr200symm.p240[[43]]) # symmetry
sum(cr200.8.p240[[43]][(14+1):43,1:8])/sum(cr200.8.p240[[43]]) # shift 0.8
# bestbox probabilities, bestbox computed here
sum(cr200symm.p240[[43]][(16+1):43,1:9])/sum(cr200symm.p240[[43]]) # symmetry
sum(cr200.8.p240[[43]][(16+1):43,1:9])/sum(cr200.8.p240[[43]]) # shift 0.8
# in both cases above target with symmetry, while a bit lower box probability
# as computed here for shift 0.8

# n=87
# bestbox probabilities, previous bestbox
sum(cr200symm.p240[[87]][(35+1):87,1:10])/sum(cr200symm.p240[[87]]) # symmetry
sum(cr200.8.p240[[87]][(35+1):87,1:10])/sum(cr200.8.p240[[87]]) # shift 0.8
# bestbox probabilities, bestbox computed here
sum(cr200symm.p240[[87]][(36+1):87,1:11])/sum(cr200symm.p240[[87]]) # symmetry
sum(cr200.8.p240[[87]][(36+1):87,1:11])/sum(cr200.8.p240[[87]]) # shift 0.8
# in both cases above target with symmetry, while a little bit lower box 
# probability as computed here for shift 0.8

# n=91
# bestbox probabilities, previous bestbox
sum(cr200symm.p240[[91]][(37+1):91,1:10])/sum(cr200symm.p240[[91]]) # symmetry
sum(cr200.8.p240[[91]][(37+1):91,1:10])/sum(cr200.8.p240[[91]]) # shift 0.8
# bestbox probabilities, bestbox computed here
sum(cr200symm.p240[[91]][(38+1):91,1:11])/sum(cr200symm.p240[[91]]) # symmetry
sum(cr200.8.p240[[91]][(38+1):91,1:11])/sum(cr200.8.p240[[91]]) # shift 0.8
# in both cases above target with symmetry, while a little bit lower box 
# probability as computed here for shift 0.8

# other targets up to n=200, precision 240

# target 0.95
# bounds for bestbox rules (code adapted from the "crossrunbox" article)
bounds$cb.95 <- NA
bounds$lb.95 <- NA

for (nn in 10:200) {
  print(nn)
  print(Sys.time())
  bounds[bounds$n == nn, c('cb.95', 'lb.95')] <- 
    bestbox(pt0=cr200symm.p240, pts=cr200.8.p240, 
            n1 = nn, target = .95, prec=240)
}

# find cut  boxes

bounds$cbord.95 <- NA
bounds$lbord.95 <- NA

for (nn in 10:200) {
  print(paste('cutbox', nn))
  print(Sys.time())
  bounds[bounds$n == nn, c('cbord.95', 'lbord.95')] <-
    cutbox(pt0=cr200symm.p240, pts=cr200.8.p240,
           n1 = nn,
           target = .95,
           prec=240,
           c1 = bounds$cb.95[bounds$n == nn],
           l1 = bounds$lb.95[bounds$n == nn])
}

# target 0.92
# bounds for bestbox rules (code adapted from the "crossrunbox" article)
bounds$cb.92 <- NA
bounds$lb.92 <- NA

for (nn in 10:200) {
  print(nn)
  print(Sys.time())
  bounds[bounds$n == nn, c('cb.92', 'lb.92')] <- 
    bestbox(pt0=cr200symm.p240, pts=cr200.8.p240, 
            n1 = nn, target = .92, prec=240)
}

# find cut  boxes

bounds$cbord.92 <- NA
bounds$lbord.92 <- NA

for (nn in 10:200) {
  print(paste('cutbox', nn))
  print(Sys.time())
  bounds[bounds$n == nn, c('cbord.92', 'lbord.92')] <-
    cutbox(pt0=cr200symm.p240, pts=cr200.8.p240,
           n1 = nn,
           target = .92,
           prec=240,
           c1 = bounds$cb.92[bounds$n == nn],
           l1 = bounds$lb.92[bounds$n == nn])
}

# target 0.91
# bounds for bestbox rules (code adapted from the "crossrunbox" article)
bounds$cb.91 <- NA
bounds$lb.91 <- NA

for (nn in 10:200) {
  print(nn)
  print(Sys.time())
  bounds[bounds$n == nn, c('cb.91', 'lb.91')] <- 
    bestbox(pt0=cr200symm.p240, pts=cr200.8.p240, 
            n1 = nn, target = .91, prec=240)
}

# find cut  boxes

bounds$cbord.91 <- NA
bounds$lbord.91 <- NA

for (nn in 10:200) {
  print(paste('cutbox', nn))
  print(Sys.time())
  bounds[bounds$n == nn, c('cbord.91', 'lbord.91')] <-
    cutbox(pt0=cr200symm.p240, pts=cr200.8.p240,
           n1 = nn,
           target = .91,
           prec=240,
           c1 = bounds$cb.91[bounds$n == nn],
           l1 = bounds$lb.91[bounds$n == nn])
}

# target 0.93
# bounds for bestbox rules (code adapted from the "crossrunbox" article)
bounds$cb.93 <- NA
bounds$lb.93 <- NA

for (nn in 10:200) {
  print(nn)
  print(Sys.time())
  bounds[bounds$n == nn, c('cb.93', 'lb.93')] <- 
    bestbox(pt0=cr200symm.p240, pts=cr200.8.p240, 
            n1 = nn, target = .93, prec=240)
}

# find cut  boxes

bounds$cbord.93 <- NA
bounds$lbord.93 <- NA

for (nn in 10:200) {
  print(paste('cutbox', nn))
  print(Sys.time())
  bounds[bounds$n == nn, c('cbord.93', 'lbord.93')] <-
    cutbox(pt0=cr200symm.p240, pts=cr200.8.p240,
           n1 = nn,
           target = .93,
           prec=240,
           c1 = bounds$cb.93[bounds$n == nn],
           l1 = bounds$lb.93[bounds$n == nn])
}


# target 0.94
# bounds for bestbox rules (code adapted from the "crossrunbox" article)
bounds$cb.94 <- NA
bounds$lb.94 <- NA

for (nn in 10:200) {
  print(nn)
  print(Sys.time())
  bounds[bounds$n == nn, c('cb.94', 'lb.94')] <- 
    bestbox(pt0=cr200symm.p240, pts=cr200.8.p240, 
            n1 = nn, target = .94, prec=240)
}

# find cut  boxes

bounds$cbord.94 <- NA
bounds$lbord.94 <- NA

for (nn in 10:200) {
  print(paste('cutbox', nn))
  print(Sys.time())
  bounds[bounds$n == nn, c('cbord.94', 'lbord.94')] <-
    cutbox(pt0=cr200symm.p240, pts=cr200.8.p240,
           n1 = nn,
           target = .94,
           prec=240,
           c1 = bounds$cb.94[bounds$n == nn],
           l1 = bounds$lb.94[bounds$n == nn])
}


# target 0.96
# bounds for bestbox rules (code adapted from the "crossrunbox" article)
bounds$cb.96 <- NA
bounds$lb.96 <- NA

for (nn in 10:200) {
  print(nn)
  print(Sys.time())
  bounds[bounds$n == nn, c('cb.96', 'lb.96')] <- 
    bestbox(pt0=cr200symm.p240, pts=cr200.8.p240, 
            n1 = nn, target = .96, prec=240)
}

# find cut  boxes

bounds$cbord.96 <- NA
bounds$lbord.96 <- NA

for (nn in 10:200) {
  print(paste('cutbox', nn))
  print(Sys.time())
  bounds[bounds$n == nn, c('cbord.96', 'lbord.96')] <-
    cutbox(pt0=cr200symm.p240, pts=cr200.8.p240,
           n1 = nn,
           target = .96,
           prec=240,
           c1 = bounds$cb.96[bounds$n == nn],
           l1 = bounds$lb.96[bounds$n == nn])
}


# target 0.915
# bounds for bestbox rules (code adapted from the "crossrunbox" article)
bounds$cb.915 <- NA
bounds$lb.915 <- NA

for (nn in 10:200) {
  print(nn)
  print(Sys.time())
  bounds[bounds$n == nn, c('cb.915', 'lb.915')] <- 
    bestbox(pt0=cr200symm.p240, pts=cr200.8.p240, 
            n1 = nn, target = .915, prec=240)
}

# find cut  boxes

bounds$cbord.915 <- NA
bounds$lbord.915 <- NA

for (nn in 10:200) {
  print(paste('cutbox', nn))
  print(Sys.time())
  bounds[bounds$n == nn, c('cbord.915', 'lbord.915')] <-
    cutbox(pt0=cr200symm.p240, pts=cr200.8.p240,
           n1 = nn,
           target = .915,
           prec=240,
           c1 = bounds$cb.915[bounds$n == nn],
           l1 = bounds$lb.915[bounds$n == nn])
}


# target 0.935
# bounds for bestbox rules (code adapted from the "crossrunbox" article)
bounds$cb.935 <- NA
bounds$lb.935 <- NA

for (nn in 10:200) {
  print(nn)
  print(Sys.time())
  bounds[bounds$n == nn, c('cb.935', 'lb.935')] <- 
    bestbox(pt0=cr200symm.p240, pts=cr200.8.p240, 
            n1 = nn, target = .935, prec=240)
}

# find cut  boxes

bounds$cbord.935 <- NA
bounds$lbord.935 <- NA

for (nn in 10:200) {
  print(paste('cutbox', nn))
  print(Sys.time())
  bounds[bounds$n == nn, c('cbord.935', 'lbord.935')] <-
    cutbox(pt0=cr200symm.p240, pts=cr200.8.p240,
           n1 = nn,
           target = .935,
           prec=240,
           c1 = bounds$cb.935[bounds$n == nn],
           l1 = bounds$lb.935[bounds$n == nn])
}


# target 0.945
# bounds for bestbox rules (code adapted from the "crossrunbox" article)
bounds$cb.945 <- NA
bounds$lb.945 <- NA

for (nn in 10:200) {
  print(nn)
  print(Sys.time())
  bounds[bounds$n == nn, c('cb.945', 'lb.945')] <- 
    bestbox(pt0=cr200symm.p240, pts=cr200.8.p240, 
            n1 = nn, target = .945, prec=240)
}

# find cut  boxes

bounds$cbord.945 <- NA
bounds$lbord.945 <- NA

for (nn in 10:200) {
  print(paste('cutbox', nn))
  print(Sys.time())
  bounds[bounds$n == nn, c('cbord.945', 'lbord.945')] <-
    cutbox(pt0=cr200symm.p240, pts=cr200.8.p240,
           n1 = nn,
           target = .945,
           prec=240,
           c1 = bounds$cb.945[bounds$n == nn],
           l1 = bounds$lb.945[bounds$n == nn])
}


# target 0.955
# bounds for bestbox rules (code adapted from the "crossrunbox" article)
bounds$cb.955 <- NA
bounds$lb.955 <- NA

for (nn in 10:200) {
  print(nn)
  print(Sys.time())
  bounds[bounds$n == nn, c('cb.955', 'lb.955')] <- 
    bestbox(pt0=cr200symm.p240, pts=cr200.8.p240, 
            n1 = nn, target = .955, prec=240)
}

# find cut  boxes

bounds$cbord.955 <- NA
bounds$lbord.955 <- NA

for (nn in 10:200) {
  print(paste('cutbox', nn))
  print(Sys.time())
  bounds[bounds$n == nn, c('cbord.955', 'lbord.955')] <-
    cutbox(pt0=cr200symm.p240, pts=cr200.8.p240,
           n1 = nn,
           target = .955,
           prec=240,
           c1 = bounds$cb.955[bounds$n == nn],
           l1 = bounds$lb.955[bounds$n == nn])
}

# adjusting variable names and order
names(bounds)[4:7]
names(bounds)[4:7] <- paste0(names(bounds)[4:7], ".925")
names(bounds)[c(1:3,16:19,32:35,12:15,4:7,20:23,36:39,
                24:27,40:43,8:11,44:47,28:31)]
bounds <-bounds[,c(1:3,16:19,32:35,12:15,4:7,20:23,36:39,
                   24:27,40:43,8:11,44:47,28:31)]
# saved as boundsall.Rdata

# separately save parts
save(x=bounds, file="bounds.Rdata")
save(x=bounds100, file="bounds100.Rdata")
save(x=cr200symm.p240, file="cr200symm.p240.Rdata")
save(x=cr200.8.p240, file="cr200.8.p240.Rdata")
save(x=cr100symm.p120, file="cr100symm.p120.Rdata")
save(x=cr100.8.p120, file="cr200.8.p120.Rdata")
