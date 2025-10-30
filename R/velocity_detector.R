

log10x1000_inv <- function(x) {
  return(10 ^ (x / 1000))
}

compute_velocity_old <- function(data) {
  x <- data$feature / data$dt
  if (length(x) > 1) x <- mean(x)
  return(x)
}

compute_velocity <- function(data) {
  distance <- sum(data$feature)
  deltaT <- tail(data[, t], 1) - head(data[, t], 1)
  velocity <- distance / deltaT
  return(velocity)
}

compute_distance <- function(data) {
  return(sum(data$feature))
}

process_timeseries <- function(data, time_window_length=10, threshold=1) {

  dt <- . <- NULL

  message(paste0("Movement detector - running.\ntime_window_length = ", time_window_length))

  feature <- "xy_dist_log10x1000"
  preproc_FUN <- log10x1000_inv
  func <- compute_velocity
  statistic <- "mean_velocity"

  # data$body_movement <- data$xy_dist_log10x1000
  d <- sleepr::prepare_data_for_motion_detector(data,
                                        c("t", feature, "x"),
                                        time_window_length,
                                        "has_interacted")

  d[, dt := c(NA, diff(t))]

  data.table::setnames(d, feature, "feature")
  # restore the distance from the log-transformed variable
  if (! is.null(preproc_FUN)) d[, feature := preproc_FUN(feature)]

  # Get a central summary value for variables of interest
  # for each window given by t_round
  # See prepare_data_for_motion_detector to learn
  # how is t_round computed
  # velocity_corrected -> max
  # has_interacted -> sum
  # beam_cross -> sum

  key <- data.table::key(d)
  d_small <- d[, .(
    statistic = func(.SD)
  ), by = "t_round"]
  data.table::setkeyv(d_small, key)

  # Gist of the program!!
  # Score movement as TRUE/FALSE value for every window
  # Score is TRUE if max_velocity of the window is > 1
  # Score FALSE otherwise
  data.table::setnames(d_small, "statistic", statistic)

  # Set t_round as the representative time of the window
  # i.e. t becomes the begining of the window and not the t
  # of the first frame in the window
  return(d_small)
}
attr(process_timeseries, "needed_columns") <- function() {
  c("t", "mean_velocity")
}
attr(process_timeseries, "parameters") <- function() {
  return(names(formals(compute_velocity)))
}

attr(process_timeseries, "variables") <- function() {
  "mean_velocity"
}


mean_velocity_detector_single_fly <- function(d, ...) {

  columns_to_keep <- c("t", attr(process_timeseries, 'needed_columns')())

  is_interpolated <- NULL
  if(nrow(d) < 100)
    return(NULL)
  # todo if t not unique, stop
  d_small <- process_timeseries(d, time_window_length, ...)
  data.table::setnames(d_small, "t_round", "t")

  if(data.table::key(d_small) != "t")
    stop("Key in output of motion_classifier_FUN MUST be `t'")

  if(nrow(d_small) < 1)
    return(NULL)
  # the times to  be queried
  time_map <- data.table::data.table(
      t = seq(from = d_small[1,t], to = d_small[.N, t], by = time_window_length),
      key = "t"
    )
  missing_val <- time_map[!d_small]

  d_small <- d_small[time_map, roll = TRUE]
  d_small[, is_interpolated := FALSE]
  d_small[missing_val, is_interpolated := TRUE]
  d_small <- stats::na.omit(d[d_small, on = c("t"), roll = TRUE])
  d_small <- d_small[, intersect(columns_to_keep, colnames(d_small)), with = FALSE]
  return(d_small)
}

#' @export
mean_velocity_detector <- function(
  data,
  time_window_length = 10, #s
  ...
){
  .SD <- NULL
  # all columns likely to be needed.
  if(is.null(data.table::key(data)))
    return(mean_velocity_detector(data))

  stopifnot(!is.null(data.table::key(data)))
  key <- data.table::key(data)
  return(
    data[,
      mean_velocity_detector_single_fly(.SD),
      by = eval(key)
    ]
  )
}


attr(mean_velocity_detector, "needed_columns") <- function() {attr(process_timeseries, 'needed_columns')()}
attr(mean_velocity_detector, "parameters") <- function() {
  args <- names(formals(mean_velocity_detector))
  args <- c(args, attr(process_timeseries, "parameters")())
  args <- unique(args)
  args <- args[args != "..."]
  args <- args[args != "data"]
  return(args)
}

attr(mean_velocity_detector, "variables") <- function() {
  attr(mean_velocity_detector, "variables")()
}
