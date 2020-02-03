

#' Returns a single radar scan from the input netcdf file.
#'
#' Smaller objects are removed and the rest are lebeled.
#' @param ncfile input netcdf file object from ncdf4 package.
#' @param var_id string name of the varibale in the file.
#' @param scan_num index of frame to be read from the file.
#' @param min_size in pixels objects smaller than this will be removed. Default 2.
#' @export
get_filteredFrame <- function(ncfile, var_id, scan_num, min_size=2) {
    frame <- ncvar_get(nc=ncfile, varid = var_id,
                        start = c(1, 1, scan_num), count = c(-1, -1, 1))
    frame <- filterFrame(frame, min_size)
    invisible(frame)
}

#' filters frame for small objects when frame is already available.
#'
#' @export
filterFrame <- function(frame, min_size = 2)
{
    labeled_echo <- EBImage::bwlabel(frame)
    frame <- clear_smallEchoes(labeled_echo, min_size)
    invisible(frame)
}



#' Removed objects smaller than the given size.
#'
#' Takes in labeled image removes objects smaller than min_size and returns re-labeled image.
#' Pixels attached diagonally are not considered as continuouse part of the object.
clear_smallEchoes <- function(label_image, min_size) {
    size_table <- table(label_image[label_image>0]) # remove zero values
    onePix_objects <- as.numeric(names(which(size_table < min_size)))

    for(obj in onePix_objects){
        label_image <- replace(label_image, label_image == obj, 0.0)
    }

    label_image <- EBImage::bwlabel(label_image)
    invisible(label_image)
}



