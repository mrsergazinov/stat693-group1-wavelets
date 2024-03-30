# Step 1: load data
setwd('D:/brani')
Xcal_lines <- readLines("Xcal.txt")
Xpre_lines <- readLines("Xpre.txt")
Ycal_lines <- readLines("Ycal.txt")
Ypre_lines <- readLines("Ypre.txt")

# Split each line into elements (assuming space-separated values) and then put it in a matrix
Xcal_elements <- lapply(Xcal_lines, function(line) strsplit(line, " ")[[1]])

Xcal_matrix <- matrix(0, nrow = 256, ncol = 39)

for (i in 1:39) {
  Xcal_matrix[, i] = as.numeric( Xcal_elements[[i]] )
}



Xpre_elements <- lapply(Xpre_lines, function(line) strsplit(line, " ")[[1]])

Xpre_matrix <- matrix(0, nrow = 256, ncol = 32)

for (i in 1:32) {
  Xpre_matrix[, i] = as.numeric( Xpre_elements[[i]] )
}


Ycal_elements <- lapply(Ycal_lines, function(line) strsplit(line, " ")[[1]])

Ycal_matrix <- matrix(0, nrow = 4, ncol = 39)

for (i in 1:39) {
  Ycal_matrix[, i] = as.numeric( Ycal_elements[[i]] )
}



Ypre_elements <- lapply(Ypre_lines, function(line) strsplit(line, " ")[[1]])

Ypre_matrix <- matrix(0, nrow = 4, ncol = 32)

for (i in 1:32) {
  Ypre_matrix[, i] = as.numeric( Ypre_elements[[i]] )
}




# Step 2: Define the orthonormal wavelet filter (Coiflet with 6 taps)
# More taps means more complicated wavelet filter, which can approximate higher order polynomial functions.
filtercoe =  wt.filter('c6')@h

# Step 3: Form the wavelet matrix of size 256 x 256 with 7 detail levels
# More detail levels mean it will capture finer details of signals.

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
  # 
  # Copyright (c) 2004. Brani Vidakovic
  #        
  # ver 1.0 Built 8/24/04; ver 1.2 Built 12/1/2004
  # This is Copyrighted Material
  # Comments? e-mail brani@isye.gatech.edu
  #   
  # updated: R code translation (Seong-joon Kim, 2015)
}  

W = Wavmat(filtercoe,2^8,7)


# Step 4: Transform the spectra to the wavelet domain, for both calibration and prediction
# sets. Transpose the results to get the observations as rows and wavelet coefs as columns
D1= W %*% Xcal_matrix
D2 = W %*% Xpre_matrix
D1 = t(D1)
D2 = t(D2)