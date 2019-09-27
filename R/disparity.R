
#' returns size change between the frames.
#'
#' Returns change in size of the eacho as ratio of bigger number to smaller
#' number when given two number, minus 1.
get_sizeChange<-function(x, y){
    if(x < 5 && y <5 ) # if too small, return zero
        return(0)
    else if(x>=y)
        return(x/y - 1)
    else
        return(y/x - 1)
}



#' Returns Disparity of all the objects in the region.
#'
#' Returns disparities of all the objects found within the search box or NA if
#' no object is present.
get_disparity_all <- function(obj_found, image2, search_box, obj1_extent) {

    if(is.na(obj_found[1]) || max(obj_found)==0) {
        obj_id2 <- 0
        dist_pred <- NA
        dist_actual <- NA
        disparity <- NA
    } else {
        obj_found <- obj_found[obj_found>0] #remove 0

        if(length(obj_found)==1){ # if this is the only object
            disparity <- get_disparity(obj_found, image2, search_box, obj1_extent)
            if(disparity <= 5) disparity <- 0 #lower the disparity if not too large

        } else { # when more than one objects to match
            disparity <- get_disparity(obj_found, image2, search_box, obj1_extent)
        }
    }
    return(disparity)
}




#' Actually computes desparity for a single object.
#'
#' Check how it is computed for detail.
#' This parameter has most effect on the acccuracy of the tracks.
get_disparity <- function(obj_found, image2, search_box, obj1_extent) {
    dist_pred <- c(NULL)
    dist_initial <- c(NULL)
    change <- c(NULL)
    overlap <- c(NULL)
    for(target_obj in obj_found){
        target_extent <- get_objExtent(image2, target_obj)
        overlap_area <- check_bigOverlap(obj1_extent, target_extent)
        overlap <- append(overlap, overlap_area)

        euc_dist<- euclidean_dist(target_extent$obj_center, search_box$center_pred)
        dist_pred <- append(dist_pred, euc_dist)

        euc_dist<- euclidean_dist(target_extent$obj_center, obj1_extent$obj_center)
        dist_initial <- append(dist_initial, euc_dist)
        size_changed <- get_sizeChange(target_extent$obj_area, obj1_extent$obj_area) #change in size
        change <- append(change, size_changed)

    }

    #This is crucial parameter that affect the results
    disparity <- dist_pred + change + dist_initial - overlap


    return(disparity)
}




#' checks overlapping regoin for big objects in both the frames.
#'
#' Checks overlapping area in pixels, size of the object and return if overlapping is considerable.
check_bigOverlap <- function(obj_extend, target_extend){
    duplicates <- duplicated(rbind(obj_extend$obj_index, target_extend$obj_index))
    overlap_area <- length(duplicates[duplicates==TRUE])
    if(obj_extend$obj_area > big_obj_size & overlap_area >= obj_extend$obj_area/2){
        return(overlap_area)
    } else
        return(0)
}



