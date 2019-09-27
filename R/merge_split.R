

#' Checks possible merging of the dead objects.
#'
#'This function takes in two R-lists containing information about current objects
#' in the frame1 and their properties, such as center location and area. If the
#' I am using an arbitrary crieterion for merging. If euclidean distance between centers
#' of the two objects c_dist < or = to r=sqrt(area), then merging is considered.
#' Here, if we assume square objects, then the r is length of a sides of the square.
check_merging<-function(dead_obj_id1, current_objects, object_props){
    nobj_frame1 <- length(current_objects$id1)
    c_dist_all <- NULL
    checked_id1 <- NULL

    #' If all objects are dead in frame2 then no merging happened.
    if(all(current_objects$id2==0)) return(0)

    for(check_obj in seq(nobj_frame1)){
        # skip checking dead objects of frame2
        if(current_objects$id2[check_obj]!=0){
            dead_xy <- c(obj_props$x[dead_obj_id1], obj_props$y[dead_obj_id1])
            merge_xy <- c(obj_props$x[check_obj], obj_props$y[check_obj])
            c_dist <- euclidean_dist(merge_xy, dead_xy)
            if(c_dist < sqrt(object_props$area[check_obj])){
                c_dist_all <-append(c_dist_all, c_dist)
                checked_id1 <- append(checked_id1, check_obj)
            }

        }
    }
    #if this is null, then no merging, else nearest object is the product of merging.
    if(is.null(c_dist_all)) return(0)
    else {
        product_id1 <- checked_id1[which(c_dist_all==min(c_dist_all))]
        return(current_objects$uid[product_id1])
    }


}

#' Find id of the parent of the new born object.
#'
#' returns unique id of the origin (or zero) for given object in frame1.
#' Also remember that old object id2 is actual id1 in frame1, as we still have
#' to update the object_ids.
get_origin_uid<-function(obj, frame1, old_objects){
    origin_id <- find_origin(obj, frame1)
    if (origin_id==0) return(0)

    origin_index <- which(old_objects$id2==origin_id)

    # If it is first observation of the object then it will not be recorded in
    # old_objects$id2, ans it will not be suitable as the origin.
    if(!(origin_id %in% old_objects$id2)) return(0)

    origin_uid <- old_objects$uid[origin_index]
    return(origin_uid)
}




#' Checks for parent in the vicinity.
#'
#' This function checks near by objects in the frame for the given new-born object.
#' origin is an object which existed before the new born objects,
#' has comparable or larger size and is close enough to the offspring.
find_origin <- function(id1_newObj, frame1){
    if(max(frame1)==1) return(0) # If there is only one object, then dont look for origin





    frame1_edges <- get_object_edges(frame1)

    #get length and indices of the given object's boundary pixels
    object_ind <- which(frame1_edges==id1_newObj, arr.ind = TRUE)
    #Do the same for all other objects in the frame
    neighbour_ind <- which(frame1_edges>0 & frame1_edges!=id1_newObj, arr.ind = T)


    nearest_object_id <- find_nearby_objects(object_ind, neighbour_ind)

    if (spatstat::is.empty(nearest_object_id)) #if no close neighbour return 0
        return(0)

    return(big_unique_obj(id1_newObj, nearest_object_id, frame1))

    # NOTE: 1. At this time we are calling big_diff_obj as origin in all the situations.
    # This looks like a good first guess. But if needed we can make it more
    # complex and use ratio and size_diff as cost function.
    # 2. We are not considering the possibility of multiple potential origins
    # beyond this point.
}






#This function returns an image/matrix with only edges of the objects.
get_object_edges <- function(frame) {
    frame_distmap <- EBImage::distmap(frame)
    frame_edges <- ifelse(frame_distmap==1, frame, 0)
    return(frame_edges)
}



find_nearby_objects <- function(object_ind, neighbour_ind) {

    object_size <- length(object_ind[,1])
    neighbour_size <- length(neighbour_ind[, 1])

    # make empty vectors
    neighbour_dist <- NULL
    neighbour_id <- NULL

    # We are chekcing for all object pixels and finding the nearest pixel.
    for(pix in seq(object_size)){
        for(neighbour in seq(neighbour_size)){
            euc_dist <- euclidean_dist(as.vector(object_ind[pix, ]),
                                       as.vector(neighbour_ind[neighbour, ]))
            neighbour_dist <- append(neighbour_dist, euc_dist)

            pix_id <- as.vector(neighbour_ind[neighbour, ])
            neighbour_id <- append(neighbour_id, frame1[pix_id[1], pix_id[2]])
        }
    }

    nearest_object_id <- neighbour_id[which(neighbour_dist < split_distance)]
    #the_nearest_object <- neighbour_id[which(neighbour_dist==min(neighbour_dist))]
    return(nearest_object_id)
}


# select the unique object that is large in size from the nearby objects
big_unique_obj <- function(new_obj_id, nearest_object_id, frame) {
    # This is to take care of multiple objects in the neighbouring region.
    neigh_objects <- unique(nearest_object_id)

    object_size <- length(frame1[frame1==new_obj_id])
    size_ratio <- NULL
    size_diff <- NULL
    for(object in neigh_objects){
        nearest_object_size <- length(frame1[frame==object])
        size_ratio <- append(size_ratio, nearest_object_size/object_size)
        size_diff <- append(size_diff, nearest_object_size - object_size)
    }

    # id of the object which has max size_ratio
    big_ratio_obj <- neigh_objects[which(size_ratio==max(size_ratio))]
    big_diff_obj  <- neigh_objects[which(size_diff==max(size_diff))]

    #if both are same call it the origin
    if(big_ratio_obj==big_diff_obj)
        return(big_diff_obj[1])
    else
        return(big_diff_obj[1])

}



#' standard Euclidean distance.
#'
#' Returns  Euclidean distance between two vectors or matrices.
euclidean_dist <- function(vec1, vec2){
    sqrt(sum((vec1-vec2)^2))
}





