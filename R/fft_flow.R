


#' Computes flow in the vicinity of the object.
#'
#' Takes in object info (major axis and center) and two images to estimate ambient flow using FFT phase correlation.
#' Margin is the additional region arround the object used to comput the flow vectors.
get_objAmbientFlow <- function(obj_extent, img1, img2, margin) {
    #coordinates of the flowregion
    r1 <- obj_extent$obj_center[1] - obj_extent$width - margin
    r2 <- obj_extent$obj_center[1] + obj_extent$width + margin
    c1 <- obj_extent$obj_center[2] - obj_extent$width - margin
    c2 <- obj_extent$obj_center[2] + obj_extent$width + margin

    dims <- dim(img1)
    if(r1<=0 || c1 <=0 || r2>dims[1] || c2 > dims[2]){
        return(c(0, 0))                         #if echo is at the image boundary
    }

    flow_region1 <- img1[r1:r2, c1:c2]
    flow_region2 <- img2[r1:r2, c1:c2]

    return(fft_flowVectors(flow_region1, flow_region2))
}


#' Alternative to get_objAmbientFlow.
#'
#' Flow vectors magnitude is clipped to given magnitude to controll erratic output from FFT flow.
get_std_flowVector<-function(obj_extent, img1, img2, margin, magnitude){
    shift <- get_objAmbientFlow(obj_extent, img1, img2, margin)
    shift <- replace(shift, shift > magnitude, magnitude)
    magnitude_negative <- magnitude * -1
    shift <- replace(shift, shift < magnitude_negative, magnitude_negative)
    return(shift)
}


#' Estimates flow vectors in two images using cross covariance of the images.
#'
#' Leese, John A., Charles S. Novak, and Bruce B. Clark.
#'           "An automated technique for obtaining cloud motion from geosynchronous
#'           satellite data using cross correlation."
#'           Journal of applied meteorology 10.1 (1971): 118-132.
fft_flowVectors <- function (im1, im2) {
    if(max(im1)==0 || max(im2)==0){
        return(c(0, 0))                         #if no object found
    }

    #im1 <- replace(im1, im1==0, runif(1))       # add noise to image background
    #im2 <- replace(im2, im1==0, runif(1))       # when objects are big and smooth

    crossCov <- fft_crossCov(im1, im2)
    cov_smooth <- spatstat::blur(spatstat::as.im(crossCov))

    dims<-dim(im1)

    pshift <- which(cov_smooth$v==max(cov_smooth$v),arr.ind=TRUE)
    pshift <- pshift-(dims[1]/2)

    return(c(pshift[1], pshift[2]))
}


#' Computes cross-covariance using FFT method, returns shifted covariance image
fft_crossCov <- function (img1, img2) {
    fft1_conj <- Conj(fft(img1)) #complex conjugate
    fft2 <- fft(img2)

    C <- (fft2*fft1_conj)/abs(fft2*fft1_conj) #crossCov in Freq domain

    crossCov <- fft(C, inv=TRUE)/length(C)
    crossCov <- Re(crossCov)
    return(fft_shift(crossCov))
}


#' Rearranges the crossCov matrix
#'
#' Rearranges the crossCov matrix so that 'zero' frequency or DC component
#'  is in the middle of the matrix. Taken from stackoverflow Que. 30630632
fft_shift <- function(fft_mat) {
    if(class(fft_mat)=='matrix') {
        rd2 <- floor(nrow(fft_mat)/2)
        cd2 <- floor(ncol(fft_mat)/2)

        # Identify the first, second, third, and fourth quadrants
        q1 <- fft_mat[1:rd2,1:cd2]
        q2 <- fft_mat[1:rd2,(cd2+1):ncol(fft_mat)]
        q3 <- fft_mat[(rd2+1):nrow(fft_mat),(cd2+1):ncol(fft_mat)]
        q4 <- fft_mat[(rd2+1):nrow(fft_mat),1:cd2]

        # rearrange the quadrants
        centered.t <- rbind(q4,q1)
        centered.b <- rbind(q3,q2)
        centered <- cbind(centered.b,centered.t)

        invisible(Re(centered))
    } else {
        stop("input to fft_shift() should be a matrix")
    }
}

