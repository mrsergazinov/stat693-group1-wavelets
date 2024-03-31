### STAT693 Group 1 Project: Wavelets
### Trisha Dawn, Dan Drennan, Renat Sergazinov*, WeiWei Wang
### Spring 2024
require(wavelets)

Wavmat <- function(h,N,k0=log(N,2),shift=2){
  # WavMat -- Transformation Matrix of FWT_PO
  #  Usage
  #    W <- WavMat(h, N, k0, shift)
  #  Inputs
  #    h      low-pass filter corresponding to orthogonal WT
  #    N      size of matrix/length of data. Should be power of 2.
  #
  #    k0     depth of transformation. Ranges from 1 to J=log2(N).
  #           Default is J.
  #    shift  the matrix is not unique an any integer shift gives
  #           a valid transformation. Default is 2.
  #  Outputs
  #    W      N x N transformation matrix
  #
  #  Description
  #    For a quadrature mirror filter h (low pass) the wavelet
  #    matrix is formed. The algorithm is described in
  #    [BV99] Vidakovic, B. (1999). Statistical Modeling By Wavelets, Wiley,
  #    on pages 115-116.
  #    Any shift is valid.  Size N=2048 is still managable on a standard PC.
  #
  #  Usage
  #    We will mimic the example 4.3.1 from [BV99] page 112.
  #   > dat <- c(1 0 -3 2 1 0 1 2);
  #   > W <- WavMat(MakeONFilter('Haar',99),2^3,3,2);
  #   > wt <- W %*% t(dat) #should be [sqrt(2)  |  -sqrt(2) |   1 -1  | ...
  #              #  1/sqrt(2) -5/sqrt(2) 1/sqrt(2) - 1/sqrt(2) ]
  #   > data <- t(W) %*% wt # should return you to the 'dat'
  #
  #


  J <- log(N,2);
  g <- rev(h*(-1)^(1:length(h)));

  if(J != floor(J)){
    stop("N has to be a power of 2.");
  }

  h <- c(h, rep(0,N));
  g <- c(g, rep(0,N));

  oldmat <- diag(2^(J-k0));

  for(k in seq(k0,1,by=-1)){

    ubJk <- 2^(J-k);
    ubJk1 <- 2^(J-k+1);

    gmat <- matrix(data=0,nrow=ubJk1,ncol=ubJk);
    hmat <- matrix(data=0,nrow=ubJk1,ncol=ubJk);

    for(jj in 1:ubJk){
      for(ii in 1:ubJk1){
        modulus <- (N+ii-2*jj+shift) %% ubJk1;
        modulus <- modulus + (modulus ==0)*ubJk1;
        hmat[ii,jj] <- h[modulus];
        gmat[ii,jj] <- g[modulus];
      }
    }

    W <- rbind(oldmat %*% t(hmat),t(gmat));
    oldmat <- W;
  }

  return(W)
  
  #
  # Copyright (c) 2004. Brani Vidakovic
  #
  # ver 1.0 Built 8/24/04; ver 1.2 Built 12/1/2004
  # This is Copyrighted Material
  # Comments? e-mail brani@isye.gatech.edu
  #
  # updated: R code translation (Seong-joon Kim, 2015)
}

# Step 1: load data

# danjdrennan: We should not need to set the working directory for this problem.
# Commenting out the `setwd` command so we don't use it.
# setwd('D:/brani')
Xcal <- t(as.matrix(read.table("Xcal.txt")))
Xpre <- t(as.matrix(read.table("Xpre.txt")))
Ycal <- as.matrix(read.table("Ycal.txt"))
Ypre <- as.matrix(read.table("Ypre.txt"))


# Step 2: Define the orthonormal wavelet filter (Coiflet with 6 taps)
# More taps means more complicated wavelet filter, which can approximate higher order polynomial functions.
filtercoe =  wt.filter('c6')@h

# Step 3: Form the wavelet matrix of size 256 x 256 with 7 detail levels
# More detail levels mean it will capture finer details of signals.
W = Wavmat(filtercoe,2^8,7)

# Step 4: Transform the spectra to the wavelet domain, for both calibration and prediction
# sets. Transpose the results to get the observations as rows and wavelet coefs as columns
D1= W %*% Xcal
D2 = W %*% Xpre
D1 = t(D1)
D2 = t(D2)

# Step 5: Predict sugar (second column in matrix Y)
# First get true sugar levels from lab measurements for calibration and prediction sets
Y1 <- Ycal[, 2]
Y2 <- Ypre[, 2]

# Step 6: Find all correlations between Y1 and the columns in D1 corresponding
# the wavelet coefficients.
corrs <- apply(D1, 2, function(x) cor(Y1, x), simplify=TRUE)

# Step 7: Find wavelet coefficient in calibration data maximizing the correlation
# with Y1.
largest_correlation <- max(abs(corrs))
largest_idx <- which(abs(corrs) == largest_correlation)

correlation_analysis_results <- list(
  "max_corr"=largest_correlation,
  "max_corr_where"=largest_idx
)

# Step 8: Fit the regression between Y1 and the "best" wavelet coefficient from the calibration set
# Choose the "best" wavelet coefficient
best_coef = D1[, correlation_analysis_results$max_corr_where]

# Run linear regression 
model_cal = lm(Y1 ~ best_coef)

# Estimated regression coefficients for the calibration set
coeff_cal = model_cal$coefficients

# Step 9: Using the regression equation from calibration set we plug the "best" wavelet coefficient 
# but from the prediction set and predict Y2.
Y2_est = coeff_cal[1] + coeff_cal[2]*D2[, correlation_analysis_results$max_corr_where]

