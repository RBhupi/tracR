
#'returns pairs of object ids from both frames.
#'
#' Given two images, the function identifies the matching
#' objects and pair them appropriatly using Hungarian method and desparity.
#'
#' @export
#' @seealso  code{disparity function}.
get_matchPairs <- function(image1, image2) {
    nObjects1 <- max(image1) #objects in first image
    nObjects2 <- max(image2) #objects in second image

    if(nObjects1==0){ #error if first image is empty
        stop("No echoes found in the first scan.")
    } else if (nObjects2==0){ #all objects will be zero if second image is empty
        zero_pairs <- rep(0, nObjects1)
        return(zero_pairs)
    }

    obj_match <- locate_allObjects(image1, image2)
    pairs <- match_pairs(obj_match) #1-to-1
    return(as.vector(pairs))
}

#' Matches objects into pairs and removes bad matching.
#'
#' The bad matching is when disparity is more than the set value.
match_pairs <- function(obj_match) {
    obj_match[obj_match<0] <- 0
    pairs <- clue::solve_LSAP(obj_match)
    pairs <- as.vector(pairs)
    # remove bad matching
    for(pair in 1:length(pairs)){
        if(obj_match[pair, pairs[pair]] > max_desparity){
            pairs[pair] <- 0
        }
    }
    return(pairs)
}


#' Matches all the obejects in image1 to the objects in image2.
#'
#' This is the main function to be called on two sets of radar images, for tracking.
locate_allObjects <- function(image1, image2) {
    nObjects1 <- max(image1) #objects in first image
    nObjects2 <- max(image2) #objects in second image


    if(nObjects2==0 || nObjects1==0){
        stop("No echoes to track!!!")
    }

    obj_match <- matrix(large_num, nrow = nObjects1, ncol = max(nObjects1, nObjects2))

    # here we match each object in image1 to all the near-by objects in image2.
    for(obj_id1 in 1:nObjects1) {
        obj1_extent <- get_objExtent(image1, obj_id1) #location, ind and radius
        shift <- get_std_flowVector(obj1_extent, image1, image2, flow_margin, maxFlow_mag)

        if(exists("current_objects"))
            shift <- correct_shift(shift, current_objects, obj_id1)

        search_box <- predict_searchExtent(obj1_extent, shift)
        search_box <- check_searchBox(search_box, dim(image2)) #search within the image
        obj_found <- find_objects(search_box, image2)  # gives possible candidates
        disparity <- get_disparity_all(obj_found, image2, search_box, obj1_extent)

        obj_match <- save_objMatch(obj_id1, obj_found, disparity, obj_match)
    }

    invisible(obj_match)
}




#' Corrects and saves disparity for the object matching stage.
#'
#' If disparity is large then it saves a large number for the value to reduce
#' the chances of this pairing to zero, else it save the value in the obj_match array.
save_objMatch <- function(obj_id1, obj_found, disparity, obj_match) {
    if(disparity > max_desparity || is.na(disparity)){
        obj_match[obj_id1, obj_found] <- large_num
    } else {
        obj_match[obj_id1, obj_found] <- disparity
    }
    return(obj_match)
}


