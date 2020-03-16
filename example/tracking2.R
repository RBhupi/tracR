library(trackR)
library(ncdf4)

start_time <- proc.time()

#+ echo=TRUE
#'----------------------- Settings for tracking method ------------------------#
#search_margin <- 4          #pixels
flow_margin <- 5            #pixels
maxFlow_mag <- 5            #fft_flow will not be faster than this
large_num <- 1000           #a large number for Hungarian method
max_obs<- 65                #longest recoreded track (else show error).
min_size <- 4               #objects smaller than this will be filtered out
max_desparity <- 18         # two objects with more desparity than this value, are not the same.
big_obj_size <- 100          # Use overlap method
#split_distance <- 3         # new object's neighboure closer than this distance is its origin
#==============================================================================#


setwd("/home/bhupendra/projects/Tracking/")

#get file names and keep only last few files.
#file_list_dbz <- Sys.glob(paths = "./data/reflectivity/*.nc")
file_class <- "./data/Darwin_REClass_2017-2017.nc"



uid_counter <- 0            #(a global variable) to start uid with 1, count from 0.
var_name<-"ATWT_ECHO_CLASSIFICATION"

ncfile <- nc_open(file_class)
x <- ncvar_get(ncfile, varid = "x")
y <- ncvar_get(ncfile, varid = "y")
nc_close(ncfile)

start_scan <- 1 #19300
newRain <- TRUE         #is this new rainy scan after dry period?

#read x, y and time from the file
ncfile <- nc_open(file_class)
time <- ncvar_get(ncfile, varid="time")
#time <- change_baseEpoch(time, From_epoch = as.Date("2004-01-01"))


outfile_name <- "~/Desktop/delete.nc" #stringr::str_replace(file_list_steiner, ".nc", "_test_delete.nc")
#outfile_name <- "~/Desktop/test_tracks.nc"
print(paste("Opening output file", outfile_name))

outNC <- create_outNC_track(outfile_name, max_obs)


print(paste("start scan = ", start_scan))
end_scan <- length(time)
print(paste("Total scans in this file", end_scan-start_scan+1))
pb = txtProgressBar(min =start_scan, max = end_scan, initial = 1, style = 3) #progress bar

#' To continuousely track the echoes in the images,
#' we read the first scan and save it in to array named frame2 then in the loop,
#' this frame will be copied to the array named frame1 and next scan will be frame2.
#'
#' frame1 <-- frame2 <-- next scan  (repeat)
#+ echo=TRUE, eval=FALSE, warning=FALSE, error=FALSE, message=FALSE
#frame2 <- get_filteredFrame(ncfile, start_scan, min_size)
#class2 <- get_classFrame(ncfile, start_scan) #classifictaion

frame2 <- ncvar_get(nc = ncfile, varid = "ATWT_ECHO_CLASSIFICATION", start = c(1, 1, start_scan),
                    count = c(-1, -1, 1))
frame2 <- replace(frame2, is.na(frame2), 0)
frame2 <- replace(frame2, frame2<2, 0)
frame2 <- filterFrame(frame2, min_size)



for(scan_ind in (start_scan+1):end_scan){
    setTxtProgressBar(pb, scan_ind) #progress bar

    #this frame0 is for checking overlap for origin of the echo.
    if(exists("frame1"))
        frame0 <- frame1

    frame1 <- frame2
    frame2 <- ncvar_get(nc = ncfile, varid = "ATWT_ECHO_CLASSIFICATION", start = c(1, 1, scan_ind),
                        count = c(-1, -1, 1))

    frame2 <- replace(frame2, is.na(frame2), 0)
    frame2 <- replace(frame2, frame2<2, 0)
    frame2 <- filterFrame(frame2, min_size)

    #if this is the last scan make it zero. This kills all the objects.
    if(scan_ind==end_scan)
        frame2 <- replace(frame2, frame2>0, 0)



    #skip if no echoes in frame 1
    if(max(frame1, na.rm = TRUE)==0){       #if no echoes in frame1
        newRain = TRUE                      #next rain will be newRain
        if(exists("current_objects")){
            rm(current_objects)
        }

        #write zeros for empty frames
        write_survival(outNC, survival_stat = rep(0, 4),
                       time = time[scan_ind], scan = scan_ind)
        next
    }

    #if(scan_ind==19399) browser()
    # track when echoes are present in frame1
    pairs <- get_matchPairs(frame1, frame2)
    obj_props <- get_objectProp(frame1, list(x=x, y=y)) #of frame1

    #when echoes are found in frame1 and it is newRain, init uids
    if(newRain){
        current_objects <- init_uids(frame1, frame2, pairs) #init ids and return list of objects

        #for newRain, all the objects are born in this frame1
        num_obj1 <- max(frame1)
        survival <- c(rep(0, 2), num_obj1, num_obj1)

        #except for the first scan where values should be missing
        if(scan_ind==start_scan+1)
            survival <- c(rep(-999, 3), num_obj1)

        write_survival(outNC, survival_stat = survival,
                       time = time[scan_ind-1], scan = scan_ind-1)

        newRain <- FALSE

    } else {            #else update old ids
        current_objects <- update_current_objects(frame0, frame1, frame2, pairs, current_objects)
    }

    #print(scan_ind)
    write_update(outNC, current_objects, obj_props, time[scan_ind-1], frame1, frame2) #for frame1

    #Survival for frame2
    num_obj2 <- max(frame2)
    obj_survival <- survival_stats(pairs, num_obj2)
    write_survival(outNC, survival_stat = obj_survival, time = time[scan_ind], scan = scan_ind)
}

#dev.off()

cat("\n") #new line required for progress bar

print("closing files")
nc_close(ncfile)
#write unlimited dim and close
ncvar_put(outNC, varid = "echo_id", vals = seq(uid_counter), start = 1, count = uid_counter)


# Stop the clock and print the time elapsed
time_elapsed <- (proc.time() - start_time)
print(paste("time elapsed", round(time_elapsed[3]/60), "minutes"))

nc_close(outNC)

