#' Plots image with objects labels. 
#' 
#' This is used in development stage to test images when processed 1-by-1.
plot_objects_label <- function(labeled_image, xvalues, yvalues){
    image2D(replace(labeled_image, labeled_image==0, NA), x=xvalues, y=yvalues)
    grid()

    for(object in 1:max(labeled_image)) {
        #center indices of the object assuming it is a rectangle
        obj_id1 <- object

        obj_index <- which(labeled_image==obj_id1, arr.ind = TRUE)

        r1 <- min(obj_index[, 1], na.rm = TRUE)
        r2 <- max(obj_index[, 1], na.rm = TRUE)
        c1 <- min(obj_index[, 2], na.rm = TRUE)
        c2 <- max(obj_index[, 2], na.rm = TRUE)

        obj_centerIndex <- c((r1+r2)/2, (c1+c2)/2)
        text(x=xvalues[obj_centerIndex[1]], y=yvalues[obj_centerIndex[2]],
             toString(object), cex=0.7)
    }
}



#' standard Euclidean distance.
#' 
#' Returns  Euclidean distance between two vectors or matrices.
euclidean_dist <- function(vec1, vec2){
    sqrt(sum((vec1-vec2)^2))
}
