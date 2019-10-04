

#' Checks possible merging of the dead objects.
#'
#'This function takes in two R-lists containing information about current objects
#' in the frame1 and their properties, such as center location and area. If the
#' I am using an arbitrary crieterion for merging. If euclidean distance between centers
#' of the two objects c_dist < or = to r=sqrt(area), then merging is considered.
#' Here, if we assume square objects, then the r is length of a sides of the square.
check_merging<-function(dead_obj_id1, current_objects, frame1, frame2){
    nobj_frame1 <- length(current_objects$id1)

    #' If all objects are dead in frame2 then no merging happened.
    if(all(current_objects$id2==0)) return(0)

    dead_obj_ind1 <- which(frame1==dead_obj_id1, arr.ind = TRUE)
    overlap_ind2 <- frame2[dead_obj_ind1]
    merge_id2 <- which.max(table(overlap_ind2[overlap_ind2>0]))
    merge_id2 <- names(merge_id2)
    merge_uid <- current_objects$uid[current_objects$id2==merge_id2]

    if(length(merge_uid)==0) return(0)
    return(merge_uid)
}

#' Find id of the parent of the new born object.
#'
#' returns unique id of the origin (or zero) for given object in frame1.
#' Also remember that old object id2 is actual id1 in frame1, as we still have
#' to update the object_ids.
get_origin_uid<-function(obj, frame1, old_objects, old_frame1){
    origin_id <- find_origin(obj, frame1, old_frame1)
    if (origin_id==0) return(0)

    origin_index <- which(old_objects$id1==origin_id)
    origin_uid <- old_objects$uid[origin_index]
    return(origin_uid)
}




#' Checks for parent in the vicinity.
#'
#' This function checks near by objects in the frame for the given new-born object.
#' origin is an object which existed before the new born objects,
#' has comparable or larger size and is close enough to the offspring.
find_origin <- function(id1_newObj, frame1, old_frame1){
    if(max(frame1)==1 || max(old_frame1)==0) return(0) # If there is only one object, then dont look for origin

    overlap_old_id1 <- old_frame1[which(frame1==id1_newObj, arr.ind = TRUE)]
    origin_old_id1 <- which.max(table(overlap_old_id1[overlap_old_id1>0]))
    if(length(origin_old_id1)==0) return(0)
    #OR else
    return(names(origin_old_id1))
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





