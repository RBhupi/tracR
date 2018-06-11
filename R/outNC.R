

#' Creates output netcdf file for radar echo trajectories. 
#'
#' This is the longest function.
#' @param ofile name and path of the output file.
#' @param max_obs longest possible track record. Deafult maximum 65 observation per track.
#' @export
create_outNC_track <- function(ofile, max_obs) {
    if(file.exists(ofile)){
        print(paste("removing existing file", basename(ofile)))
        file.remove(ofile)
    }
    deflat <- 9

    dim_echo <- ncdim_def("echo_id", vals=1, units = "", unlim = TRUE,
                          longname = "unique id of convection echo", create_dimvar = TRUE)

    dim_obs <- ncdim_def("records", vals = seq(max_obs), units="",
                         longname = "observation records")

    dim_time <- ncdim_def("time", vals=1, units = "seconds since 1970-01-01 00:00:00 UTC",
                          longname = "time of the scan", unlim = TRUE, create_dimvar = TRUE)

    dim_stat <- ncdim_def("stat", vals = seq(4), units="", longname = "object survival vector; lived, died, born, total")

    # Define Variables
    var_survival <- ncvar_def("survival", units = "", longname = "survival from the last scan",
                              dim=list(dim_stat, dim_time), missval = -999, prec="integer",
                              compression = deflat, shuffle = TRUE)

    var_dur <- ncvar_def("duration", units = "", longname = "duration of echo in time-steps",
                         dim=dim_echo, missval = -999, prec="integer",
                         compression = deflat, shuffle = TRUE)

    var_origin <- ncvar_def("origin", units="", longname = "id from which the echo split up.",
                            dim=dim_echo, missval = 0, prec = "integer")

    var_merged <- ncvar_def("merged", units="", longname = "id in which the echo merged.",
                            dim=dim_echo, missval = 0, prec = "integer")


    var_time <- ncvar_def("record_time", units = "seconds since 1970-01-01 00:00:00 UTC",
                          longname = "time of the scan for each record",
                          dim = list(dim_obs, dim_echo), missval = -999, prec = "integer",
                          compression = deflat, shuffle = TRUE)

    var_xdist <- ncvar_def("x_dist", units = "Km", longname = "distance from Radar",
                           dim = list(dim_obs, dim_echo), missval = -999.0, prec = "float",
                           compression = deflat)

    var_ydist <- ncvar_def("y_dist", units = "Km", longname = "distance from Radar",
                           dim = list(dim_obs, dim_echo), missval = -999.0, prec = "float",
                           compression = deflat)


    var_x <- ncvar_def("x", units = "", longname = "index along x-coordinate",
                       dim = list(dim_obs, dim_echo), missval = -999, prec = "integer",
                       compression = deflat, shuffle = TRUE)

    var_y <- ncvar_def("y", units = "", longname = "index along y-coordinate",
                       dim = list(dim_obs, dim_echo), missval = -999, prec = "integer",
                       compression = deflat, shuffle = TRUE)

    var_npix <- ncvar_def("area", units = "pixels", longname = "area of the echo in pixels",
                          dim = list(dim_obs, dim_echo), missval = -999, prec = "integer",
                          compression = deflat, shuffle = TRUE)

    #var_ncg <- ncvar_def("Cu_cong", units = "pixels", longname = "num of Cu Congestus pixels",
    #                     dim = list(dim_obs, dim_echo), missval = -999, prec = "integer",
    #                     compression = deflat, shuffle = TRUE)

    #var_ncb <- ncvar_def("Cu_deep", units = "pixels", longname = "num of deep convection pixels",
    #                     dim = list(dim_obs, dim_echo), missval = -999, prec = "integer",
    #                     compression = deflat, shuffle = TRUE)

    #var_nco <- ncvar_def("Cu_over", units = "pixels", longname = "num of overshooting convection pixels",
    #                     dim = list(dim_obs, dim_echo), missval = -999, prec = "integer",
    #                     compression = deflat, shuffle = TRUE)

    var_list <- list(var_time, var_survival, var_dur, var_origin, var_merged, var_xdist, var_ydist,
                     var_x, var_y, var_npix) #, var_ncg, var_ncb, var_nco)


    outNC <- nc_create(filename = ofile, vars = var_list)

    write_settingParms_toNC(outNC)

    #for CF standards
    ncatt_put(outNC, varid = "echo_id", attname = "cf_role", attval = "trajectory_id")
    ncatt_put(outNC, varid = 0, attname = "featureType", attval = "trajectory")

    description <- paste("The CPOL (Darwin) radar echoes of convective types were separated using Steiner classification scheme and tracked.",
                         "Merging and splitting is added with echo ids.")

    ncatt_put(outNC, varid = 0, attname = "_description",
              attval = description, prec = "text")
    ncatt_put(outNC, varid = 0, attname = "_creator",
              attval = "Bhupendra Raut", prec = "text")
    ncatt_put(outNC, varid = 0, attname = "_url",
              attval = "www.baraut.info", prec = "text")
    ncatt_put(outNC, varid = 0, attname = "_date_created",
              attval = date(), prec = "text")

    invisible(outNC)
}


#' Writes all the setting parameters (as attributes) for the tracking. 
#' 
#' These parameters affect the sensitivity of the tracks, mergers and split definitions etc.
write_settingParms_toNC <- function(outNC){
    ncatt_put(outNC, varid = 0, attname = "search_margin", attval =search_margin, prec = "short")
    ncatt_put(outNC, varid = 0, attname = "flow_margin", attval =flow_margin, prec = "short")
    ncatt_put(outNC, varid = 0, attname = "maxFlow_magnitude", attval =maxFlow_mag, prec = "short")
    ncatt_put(outNC, varid = 0, attname = "min_echoSize_toTrack", attval =min_size, prec = "short")
    ncatt_put(outNC, varid = 0, attname = "max_desparity", attval =max_desparity, prec = "short")
}



#' Writes properties and uids for all objects into output netcdf file.
#' 
#' @param outNC output netcdf file object from function \code{\link{create_outNC_track}}
#' @param current_objects output of \code{update_current_objects}
#' @param obj_props output of \code{get_object_prop()}
#' @param obs_time time of first scan in POSIX format. units="seconds since 1970-01-01".
#' @export
write_update<-function(outNC, current_objects, obj_props, obs_time){
    nobj <- length(current_objects$id1) #num of objects in frame1

    for(object in seq(nobj)){
        nc_start <- c(current_objects$obs_num[object], current_objects$uid[object])
        nc_count <- c(1, 1)

        ncvar_put(outNC, varid = "record_time", obs_time, start = nc_start, count = nc_count)
        ncvar_put(outNC, varid = "x", obj_props$x[object], start = nc_start, count = nc_count)
        ncvar_put(outNC, varid = "y", obj_props$y[object], start = nc_start, count = nc_count)
        ncvar_put(outNC, varid = "x_dist", obj_props$xdist[object], start = nc_start, count = nc_count)
        ncvar_put(outNC, varid = "y_dist", obj_props$ydist[object], start = nc_start, count = nc_count)

        ncvar_put(outNC, varid = "area", obj_props$area[object],  start = nc_start, count = nc_count)
    }

    write_duration(outNC, current_objects)
}


#' Write duration of dead objects in output NC file.
#' 
#' Writes number of observations for dead objects. Duration is in time-steps.
write_duration <- function(outNC, current_objects){
    nobj <- length(current_objects$id1)

    for (obj in seq(nobj)){
        if(current_objects$id2[obj]==0){
            ncvar_put(outNC, varid = "duration", current_objects$obs_num[obj],
                      start=current_objects$uid[obj], count=1)
            ncvar_put(outNC, varid = "origin", current_objects$origin[obj],
                      start=current_objects$uid[obj], count=1)

            #check for merging
            merged_in <- check_merging(obj, current_objects, obj_props)
            ncvar_put(outNC, varid="merged", merged_in,
                      start=current_objects$uid[obj], count=1)
        }
    }
}


#' Write survival stats 
#' 
#' write number of lived, dead and born objects to the file for each scan.
#' 
#' @export
write_survival <- function(outNC, survival_stat, time, scan){
    if(!is.atomic(survival_stat)){
        survival_stat <- unlist(survival_stat, use.names = FALSE)
    }

    ncvar_put(outNC, varid = "survival", vals = survival_stat, start = c(1, scan), count = c(4, 1))
    ncvar_put(outNC, varid = "time", vals = time, start = scan, count=1)
}

