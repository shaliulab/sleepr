#' Motion detector for Ethocope data
#'
#' Defines whether a *single animal* is moving according to:
#'
#' * Validated and corrected subpixel velocity ([max_velocity_detector]), the most rigorous
#' * Uncorrected subpixel velocity ([max_velocity_detector_legacy])
#' * Crossing a virtual beam in the middle of the region of interest ([virtual_beam_cross_detector])
#'
#' [max_velocity_detector] is the default movement classification for real-time ethoscope experiments.
#' It is benchmarked against human-generated ground truth.
#' @name motion_detectors
#' @param data [data.table::data.table] containing behavioural variables of *a single animal* (no id).
#' It must have the columns `xy_dist_log10x1000`(for computing subpixel velocity),
#' `x`(beam cross), `t` and `has_interacted` (whether a stimulus was delivered).
#' @param velocity_correction_coef an empirical coefficient to correct velocity with respect
#'  to variable framerate.
#' @inheritParams sleep_annotation
#' @inheritParams prepare_data_for_motion_detector
#' @param masking_duration number of seconds during which any movement is ignored (velocity is set to 0) after
#' a stimulus is delivered (a.k.a. interaction).
#' @param velocity_threshold uncorrected velocity above which an animal is classified as `moving' (for the legacy version).
#' @return an object of the same type as `data` (i.e. [data.table::data.table] or [behavr::behavr])  with additional columns:
#' * `moving` Logical, TRUE iff. motion was detected.
#' * `beam_crosses` The number of beam crosses
#' (when the animal crosses x = 0.5 -- that is the midpoint of the region of interest) within the time window
#' * `max_velocity` The maximal velocity within the time window.
#' The resulting data is sampled at a period equals to `time_window_length`.
#' @details
#'  These functions are *rarely called directly*, but typically used is in the context of [sleep_annotation].
#' @seealso
#' * [sleep_annotation] -- which requieres a motion detector
#' @export
#' @importFrom dplyr right_join summarise
#' @importFrom magrittr `%>%`
  max_velocity_detector  <- function(data,
                                   time_window_length,
                                   velocity_correction_coef = 3e-3,
                                   masking_duration = 6,
                                   curate=TRUE
                                   ){
  dt = x = .N = . = velocity = moving = dist = beam_cross = has_interacted = NULL
  dt = beam_crossed =  interaction_id = masked = interactions =  NULL
  xy_dist_log10x1000 = max_velocity = velocity_corrected = NULL

  message(sprintf("Velocity correction coefficient: %#.4f", velocity_correction_coef))


  ## Preprocessing
  ## ----
  # get only the columns that the motion detector needs
  # bin the data by windows of time_window_length seconds
  # remove sparse data i.e.
  # data in windows with very few data already
  # very few means < 20 datapoints in 1 minute
  d <- prepare_data_for_motion_detector(data,
                                        c("t", "xy_dist_log10x1000", "x"),
                                        time_window_length,
                                        "has_interacted", curate=curate)

  ## Define velocity as the distance traversed
  ## between two consecutive frames and the time
  ## between them
  ## Unit: fraction_of_roi_width / s
  ## TODO Change the code so it becomes mm / s
  ## or something with a direct physical interpretation
  ## ----
  # compute the inter frame time
  d$dt <-  c(NA,diff(d$t))

  #d[,surface_change := xor_dist * 1e-3]

  # restore the distance from the log-transformed variable
  d[,dist := 10^(xy_dist_log10x1000/1000) ]
  d$velocity <- 10^(d$xy_dist_log10x1000/1000)/d$dt

  a = velocity_correction_coef
  ## ----

  ## Compute beam crosses for DAM-like data generation
  ## ----
  # compute the side the fly is in on each frame
  # -1 for right and 1 for left (does not matter which is which)
  # changes in side are thus encoded as a change in sign
  # making a diff of the side thus returns 0 when the side does not change
  # and 1/-1 when it does
  # we dont care about the direction i.e. whether it is 1 or -1
  # so take the absolute value and take that as a beam cross
  d$beam_cross <- as.logical(abs(c(0,diff(sign(.5 -d$x)))))
  ## ----


  ## Masking
  ## Intended to ignore activity within masking_duration (default 6) seconds
  ## after an interaction happens
  ## ----
  if(!"has_interacted" %in% colnames(d)){
    if(masking_duration > 0)
      warning("Data does not contain an `has_interacted` column.
              Cannot apply masking!.
              Set `masking_duration = 0` to ignore masking")
    d$has_interacted <- 0
  }

  # create a unique identifier for every interaction
  # this is done so we can split d
  # so rows within the same block share last interaction
  d$interaction_id <- cumsum(d$has_interacted)

  # masked becomes TRUE if t is within masking_duration
  # after t[1]
  # since the first row of the block
  # represents the last interaction
  # the time of the interaction is t[1]
  # d[,
  #   masked := t < (t[1] + masking_duration),
  #   by=interaction_id
  #   ]

  t0 <- dplyr::group_by(d, interaction_id) %>% dplyr::summarise(t0 = min(t))
  d <- right_join(d, t0[, c("t0", "interaction_id")])
  # browser()

  d$masked <- ifelse(d$t < (d$t0 + masking_duration), TRUE, FALSE)
  d$t0 <- NULL

  # velocity is 0 if the mask is TRUE
  # all the data up to the first interaction in the whole experiment time course
  # is masked with the above protocol (under interaction_id 0)
  # so ignore the mask if interaction_id is 0 because there was no interaction
  d$velocity <- ifelse(d$masked & d$interaction_id != 0, 0, d$velocity)

  # in the same way, beam cross can only be TRUE if masked is FALSE
  d$beam_cross <- !d$masked & d$beam_cross

  # remove the interaction_id and masked columns
  # to preserve state
  d$interaction_id <- NULL
  d$masked <- NULL
  ## ----

  # velocity correction to handle
  # velocity being dependent on FPS
  # See quentin's Thesis
  # PDF -> https://spiral.imperial.ac.uk:8443/handle/10044/1/69514
  # First paragraph in https://github.com/rethomics/sleepr/issues/7#issuecomment-579297206
  d$velocity_corrected <-  d$velocity  * d$dt  / a



  # Get a central summary value for variables of interest
  # for each window given by t_round
  # See prepare_data_for_motion_detector to learn
  # how is t_round computed
  # velocity_corrected -> max
  # has_interacted -> sum
  # beam_cross -> sum
  stats <- dplyr::group_by(d, t_round) %>% summarise(
      max_velocity = max(velocity[2:length(velocity)]),
      interactions = as.integer(sum(has_interacted)),
      beam_crosses = as.integer(sum(beam_cross))
    )

  d_small <- dplyr::right_join(d, stats[, c("max_velocity", "interactions", "beam_crosses", "t_round")])

  # Gist of the program!!
  # Score movement as TRUE/FALSE value for every window
  # Score is TRUE if max_velocity of the window is > 1
  # Score FALSE otherwise
  d_small$moving <-  ifelse(d_small$max_velocity > 1, TRUE,FALSE)

  # Set t_round as the representative time of the window
  # i.e. t becomes the begining of the window and not the t
  # of the first frame in the window
  d_small <- data.table::as.data.table(d_small)
  data.table::setnames(d_small, "t_round", "t")
  setkey(d_small, t)
  return(d_small)
}


attr(max_velocity_detector, "needed_columns") <- function(...){
  c("xy_dist_log10x1000", "x", "has_interacted")
}




#' @export
#' @rdname motion_detectors
virtual_beam_cross_detector <- function(data, time_window_length){
  beam_crosses = x=.= beam_cross= NULL
  d <- prepare_data_for_motion_detector(data,
                                        c("t", "x"),
                                        time_window_length)
  d[,beam_cross := abs(c(0,diff(sign(.5 - x))))]
  d[,beam_cross := as.logical(beam_cross)]

  d_small <- d[,
               .(moving = any(beam_cross)),
               by = "t_round"]
  data.table::setnames(d_small, "t_round", "t")
  d_small
}

attr(virtual_beam_cross_detector, "needed_columns") <- function(...){
  c("x")
}

#' Copy needed columns, and add t_round var for downsampling
#'
#' @inheritParams sleep_annotation
#' @param needed_columns Columns that must be present in the data
#' @param optional_columns Optional, columns that can be used if available
#' @param curate If true, remove sparse ROI datasets (probably spurious fly detections)
#' @importFrom data.table data.table setkeyv
#' @export
prepare_data_for_motion_detector <- function(data,
                                             needed_columns,
                                             time_window_length,
                                             optional_columns = NULL,
                                             curate=TRUE){
  # todo assert no key/unique
  t_round <- NULL
  if(! all(needed_columns %in% names(data)))
    stop(sprintf(
      "data from ethoscope should have columns named %s!",
      paste(needed_columns, collapse = ", ")
    ))

  # include the optional columns if available
  # and make sure the resulting colums are unique
  # optional_columns is tipically columns that the motion detector can use
  # but they are not required i.e. has_interacted for instance
  needed_columns <- unique(c(needed_columns, intersect(names(data), optional_columns)))

  # retrieve data containing only the needed columns
  d <- data.table::copy(data[, needed_columns, with = FALSE])

  # compute t_round
  # t_round becomes the same for all datapoints in the same window
  # even if their exact t is not
  # each window is time_window_length seconds long
  # the first window starts at t = 0 s
  # the windows do not overlap
  # for t in [300, 309.999] t_round = 300
  # for t 310, it is 310 i.e. data form 310 and thereafter is analyzed as another window
  d$t_round <-  time_window_length * floor(d$t / time_window_length)

  # remove datapoints belonging to windows (of size 60 seconds by default)
  # where the number of datapoints is less than a default of 20.
  if (curate) {
    before_n <- nrow(d)
    d <- curate_sparse_roi_data(d)
    after_n <- nrow(d)
    if (before_n > 0 & after_n == 0) {
      message("Data is too sparse (datapoints per bin are too low)")
    }
  }

  # TODO in order to ...
  data.table::setkeyv(d, "t_round")
}
