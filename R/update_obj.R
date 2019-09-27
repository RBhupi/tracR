
#' saves last x y movements of the objects.
#'
#' Attaches last xyheads to current objects for future use.
attach_xyheads <- function(frame1, frame2, current_objects) {

    nobj <- length(current_objects$uid)
    xhead <- yhead <- NULL

    for (obj in seq(nobj)) {
        if(current_objects$id1[obj]>0 && current_objects$id2[obj]>0){
            center1<- get_objectCenter(current_objects$id1[obj], frame1)
            center2<- get_objectCenter(current_objects$id2[obj], frame2)
            xhead <- append(xhead, center1[1]-center2[1])
            yhead <- append(yhead, center1[2]-center2[2])
        }else{ #if object is dead write NA
            xhead<-append(xhead, NA)
            yhead<-append(yhead, NA)
        }
    }
    return(cbind(current_objects, xhead, yhead)) #attach values
}


#' Removes dead objects, updates living objects and assign new uids to new born objects.
#'
#' Also, updates number of valid observations for each echo.
#' This function is called when rain continues from the last frame.
#' This is a complicated function to understand.
#'
#' @details See how the pairs vector looks like for a real case. The pairs
#' shows mapping of the current frame1 and frame2. This shows that frame2 has 4 objects.
#' The objects [1, 2, 3, 4] in current frame2 are mapped with objects [0, 1, 2, 3]
#' in current frame1. Thus, object 1 in frame2 is new born. Others can be traced back to
#' frame1.
#'
#' pairs>>
#'
#' 0, 1, 2, 3
#'
#' Now check old_objects and remember that at this instant, id2 (in the old_objects)
#' correspond to the objects in current frame1 which was frame2 in the earlier
#' time-step, and that they are the same frame.
#'
#' old_objects>>
#'
#' id1, uid, id2, obs_num, xhead, yhead
#'
#'  1, 1, 1,   2,   1,  0
#'
#'  2,  12,  3, 2, 0,  -1
#'
#' So the object 1 and 3 in current frame1 (earlier it was frame2 with id2) existed
#' before and has "uid" (11 and 12). We will copy their "uid" to our object_matrix
#' and increament the observation number (obs_num).
#' For object 2 and 4 in current frame2 which do not exist in frame1,
#' we will ask for new uids. This information will be written  in
#' current_objects and return for writting in to the output file.
#' @export
update_current_objects <- function(frame1, frame2, pairs, old_objects){
    nobj <- max(frame1)
    objects_mat <- matrix(data = NA, ncol = 5, nrow = nobj)

    objects_mat[, 1] <- seq(nobj) # this is id1 at current step

    for (obj in seq(nobj)){
        if(obj %in% old_objects$id2){ # but same was id2 in the last step
            # so they should get same uid as last time
            objects_mat[obj, 2] <- old_objects$uid[old_objects$id2==obj]
            objects_mat[obj, 4] <- old_objects$obs_num[old_objects$id2==obj] + 1
            objects_mat[obj, 5] <- old_objects$origin[old_objects$id2==obj]
        } else {
            objects_mat[obj, 2] <- next_uid()
            objects_mat[obj, 4] <- 1 #first observation of the echo
            objects_mat[obj, 5] <- get_origin_uid(obj, frame1, old_objects)
        }
    }

    objects_mat[, 3] <- as.vector(pairs) #match as they are in frame2

    current_objects <- data.frame(objects_mat, row.names = NULL)
    colnames(current_objects) <- c("id1", "uid", "id2", "obs_num", "origin")
    current_objects <- attach_xyheads(frame1, frame2, current_objects)
    invisible(current_objects)
}





#' Returns a dataframe for objects with unique ids and their corresponding ids in frame1 and frame2.
#'
#' This function is called when new rainy scan is seen after the period of no rain or the first time.
#' @param first_frame First image for tracking.
#' @param second_frame Second image for tracking.
#' @param pairs output of \code{get_match_pairs()}
#' @export
init_uids <- function(first_frame, second_frame, pairs){
    nobj <- max(first_frame) #number of objects in frame1
    objects_mat <- matrix(data = NA, ncol = 5, nrow = nobj)

    objects_mat[, 1] <- seq(nobj)               #id1
    objects_mat[, 2] <- next_uid(count = nobj) #unique ids
    objects_mat[, 3] <- as.vector(pairs) #as they are in frame2
    objects_mat[, 4] <-rep(1, nobj)     #observation number for the echo
    objects_mat[, 5] <-rep(0, nobj)

    current_objects <- data.frame(objects_mat, row.names = NULL)
    colnames(current_objects) <- c("id1", "uid", "id2", "obs_num", "origin")
    current_objects <- attach_xyheads(first_frame, second_frame, current_objects)
    return(current_objects)
}



#'
#'
#' Returns sequence of next unique ids and increament the uid_counter.
next_uid<-function(count=1){
    this_uid <- uid_counter + 1:count
    uid_counter <<- uid_counter + count
    return(this_uid)
}


#' Returns a list with number of objects lived, died and born between the current and the previousstep.
#'
#' @export
survival_stats <- function(pairs, num_obj2) {
    pairs_vec <- as.vector(pairs)
    obj_lived <- length(pairs_vec[pairs_vec>0])
    obj_died <- length(pairs_vec)-obj_lived
    obj_born <- num_obj2 - obj_lived
    return(list(lived=obj_lived, died=obj_died, born=obj_born, total=num_obj2))
}





