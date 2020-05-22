
#' FFT shifts are corrected using last headings.
#'
#' takes in flow vector based shift and current_object dataframe which has last
#' headings, and check if they are resonably close if not rejects or modify shift and return.
#' Note:  frame2 of last timestep is now frame1, but current_objects still has it as frame2.
#' So id2 in the last frame2 are actually ids related to frame1 now.
correct_shift<-function(this_shift, current_objects, object_id1){
    last_heads <- c(current_objects$xhead[current_objects$id2==object_id1],
                    current_objects$yhead[current_objects$id2==object_id1])

    #for small shifts and empty last shifts
    if((length(last_heads)<2) || any(is.na(last_heads)))
        return(this_shift)
    else if(any(abs(this_shift-last_heads)>4))  #if they are too different
        return(last_heads)                      #then trust last_heads
    else return((this_shift+last_heads)/2) #else retun the average of both
}


#' Returns object extent properties.
#'
#' Takes in a labeled image and finds the radius, area and the center of the given object.
get_objExtent <- function(labeled_image, obj_label) {
    #center indices of the object assuming it is a rectangle
    obj_index <- which(labeled_image==obj_label, arr.ind = TRUE)
    x_len <- (max(obj_index[, 1]) - min(obj_index[, 1]) + 1)
    y_len <- (max(obj_index[, 2]) - min(obj_index[, 2]) + 1)

    obj_width <- max(c(x_len, y_len)) #maximum possible object radius

    #definition of object center based on median, This is working better.
    obj_center <- get_objectCenter(obj_label, labeled_image)
    obj_area <- length(obj_index[, 1])  #size in pixels


    obj_extent<-list(obj_center=obj_center, width=obj_width,
                     obj_area=obj_area, obj_index=obj_index)

    return(obj_extent)
}



#' predict search region.
#'
#' Predicts search extent/region for the object in image2 given the image shift.
predict_searchExtent <- function(obj1_extent, shift){
    shifted_center <- obj1_extent$obj_center + shift
    search_radius <- ceiling(sqrt(obj1_extent$obj_area))

    if(search_radius< min_size) search_radius <- min_size

    x1 <- shifted_center[1] - search_radius
    x2 <- shifted_center[1] + search_radius
    y1 <- shifted_center[2] - search_radius
    y2 <- shifted_center[2] + search_radius

    return(list(x1=x1, x2=x2, y1=y1, y2=y2, center_pred=shifted_center))
}


#' checks that the search box is in the domain.
#'
#' Returns NA if search box  outside the image or search box is very small.
check_searchBox <- function(search_box, img_dims){

    if(search_box$x1 <= 0){
        search_box$x1 <- 1
    }
    if(search_box$y1 <= 0){
        search_box$y1 <- 1
    }
    if(search_box$x2 > img_dims[1]){
        search_box$x2 <- img_dims[1]
    }
    if(search_box$y2 > img_dims[2]){
        search_box$y2 <- img_dims[2]
    }

    #search box should be large enough
    if(search_box$x2-search_box$x1 < 4 || search_box$y2-search_box$y1 < 4 ){
        return(NA)
    } else {
        return(search_box)
    }
}


#' Returns vector of objects ids in the given reion.
#'
#' Given the search box and image2, returns objects in the region.
find_objects <- function(search_box, image2) {
    #if search box is NA then object left the image
    if(is.na(search_box[1])){
        obj_found <- NA
    } else {
        search_area <- image2[search_box$x1:search_box$x2, search_box$y1:search_box$y2]
        obj_found <- unique(as.vector(search_area))
    }
    return(obj_found)
}




#' return center indices of the object.
#'
#' Returns indices of center pixel of the given object id from a labeled image.
#' This may be done in better way for non-oval objects.
get_objectCenter<-function(obj_id, labeled_image){
    obj_index <- which(labeled_image==obj_id, arr.ind = TRUE)
    center_x <- median(obj_index[, 1]) #center column
    center_y <- median(obj_index[, 2]) #center row

    pix_value <- labeled_image[center_x, center_y]

    if(pix_value!=obj_id){
        dist_pix <- NULL
        for(apix in 1:length(obj_index[, 1])){
            dist_pix <-append(dist_pix, euclidean_dist(obj_index[apix,], c(center_x, center_y)))
        }
        closest_obj_pix <- which.min(dist_pix)
        center_x <- obj_index[closest_obj_pix, 1]
        center_y <- obj_index[closest_obj_pix, 2]
    }

    return(c(center_x, center_y))
}



#' Return all the object's size, location and classification info.
#'
#' xyDist should be a list of x_dist and y_dist in km.
#'
#'@export
get_objectProp <- function(image1, xyDist){
    objprop <- c(NULL)
    nobj <- max(image1)

    for(obj in seq(nobj)){
        obj_index <- which(image1==obj, arr.ind = TRUE)

        x_len <- (max(obj_index[, 1]) - min(obj_index[, 1]) + 1)
        y_len <- (max(obj_index[, 2]) - min(obj_index[, 2]) + 1)

        obj_width <- max(c(x_len, y_len)) #maximum possible object radius
        obj_breadth <- min(c(x_len, y_len)) #maximum possible object radius

        ellipse_par <- fitEllipse(obj_index)
        circularity <- ellipse_par$minor/ellipse_par$major #minor_axis/major_axis

        objprop$id1 <- append (objprop$id1, obj)  #id in frame1
        objprop$x <- append(objprop$x, floor(median(obj_index[, 1]))) #center column
        objprop$y <- append(objprop$y, floor(median(obj_index[, 2]))) #center row
        objprop$area <- append(objprop$area, length(obj_index[, 1]))
        objprop$width <- append(objprop$width, obj_width)
        objprop$breadth <- append(objprop$breadth, obj_breadth)
        objprop$circularity<- append(objprop$circularity, circularity)
        objprop$orientation<- append(objprop$orientation, ellipse_par$angle)
    }
    objprop <- attach_xyDist(objprop, xyDist$x, xyDist$y)
    invisible(objprop)
}




#'
#'
#' Attaches y and x distance from radar in km to object location indices
attach_xyDist<-function(obj_props, xdist, ydist){
    obj_props$xdist <- xdist[obj_props$x]
    obj_props$ydist <- ydist[obj_props$y]
    invisible(obj_props)
}


#'
#'
#'Try to fit optimum ellipse for the object circularity, given the object index.
#'The ellipse fitted with this method does not enclosed the object but it provide optimum ellipse parameters.
#'The ratio of major and minor axis and eccentricity are correctly estimated for tracking purpose.
fitEllipse <- function(object_index){
    x_len <- (max(object_index[, 1]) - min(object_index[, 1]) + 1)
    y_len <- (max(object_index[, 2]) - min(object_index[, 2]) + 1)

    ellipseGPar <- NULL

    if(x_len < 3 | y_len < 3){ # can not fit ellipse
        ellipseGPar <- list(center = c(NA, NA), major = NA,
                            minor = NA, angle = NA)
    }else{
        tryCatch({
        ellipDirect <- conicfit::EllipseDirectFit(object_index)
        ellipDirectG <- conicfit::AtoG(ellipDirect)$ParG
        major_ax <- max(ellipDirectG[3, ], ellipDirectG[4, ])
        minor_ax <- min(ellipDirectG[3, ], ellipDirectG[4, ])

        ellipseGPar <- list(center = c(ellipDirectG[1:2, ]), major = major_ax,
                            minor = minor_ax, angle = ellipDirectG[5, ])
        },error=function(err){
            ellipseGPar <<- list(center = c(NA, NA), major = NA, minor = NA, angle = NA)
            return(NULL)
        }

        )}
    return(ellipseGPar)
}

