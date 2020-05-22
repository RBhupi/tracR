
#' returns size change between the frames.
#'
#' Returns change in size of the echo.
get_sizeChange<-function(x, y){
    return(abs(x-y))
}



#' Returns Disparity of all the objects in the region.
#'
#' Returns disparities of all the objects found within the search box or NA if
#' no object is present.
get_disparity_all <- function(obj_found, image2, search_box, obj1_extent) {

    if(is.na(obj_found[1]) || max(obj_found)==0) {
        obj_id2 <- 0
        disparity <- NA
    } else {
        obj_found <- obj_found[obj_found>0] #remove 0
        disparity <- get_disparity(obj_found, image2, search_box, obj1_extent)
    }

    disparity <- replace(disparity, disparity<0, 0)
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
    #sink(file = "./output/disparity_param.txt", append = TRUE)
    for(target_obj in obj_found){
        target_extent <- get_objExtent(image2, target_obj)
        overlap_area <- check_bigOverlap(obj1_extent, target_extent)
        overlap <- append(overlap, overlap_area)
        euc_dist_pred<- euclidean_dist(target_extent$obj_center, search_box$center_pred)
        dist_pred <- append(dist_pred, euc_dist_pred)

        euc_dist_cent<- euclidean_dist(target_extent$obj_center, obj1_extent$obj_center)
        dist_initial <- append(dist_initial, euc_dist_cent)
        size_changed <- get_sizeChange(target_extent$obj_area, obj1_extent$obj_area) #change in size
        change <- append(change, size_changed)
        #cat(euc_dist_cent, euc_dist_pred, size_changed, overlap_area, "\n", sep = " ")
    }

    #This is crucial parameter that affect the results. this is found to be stable.
    #the idea is to combined change in location (distances) and change in area in disparity.
    disparity <- dist_pred + dist_initial + sqrt(change) - sqrt(overlap)

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



