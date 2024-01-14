
## PROCESS SERIES AND ASSOCIATED FUNCTIONS

rotate_quads <- function(x, quadrant){
  if (quadrant == "LR"){
    return(x)
  } else if (quadrant == "LL") {
    x[,2] <- -x[,2]
    return(x)
  } else if (quadrant == "UR") {
    x[,3] <- -x[,3]
    return(x)
  } else if (quadrant == "UL") {
    x[,2:3] <- -x[,2:3]
    return(x)
  }
}

process_series <- function(path){

  series <- read.csv(path, sep = "\t")

  quad_id <- strsplit(path, ".txt")[[1]][1]
  quad_id <- str_sub(quad_id, start = -2)

  series <- rotate_quads(series, quad_id)
  series <- cbind.data.frame(series,
                             integer(length = length(rownames(series))),
                             integer(length = length(rownames(series))),
                             integer(length = length(rownames(series))),
                             integer(length = length(rownames(series))),
                             integer(length = length(rownames(series))),
                             integer(length = length(rownames(series))))
  colnames(series) <- c("Time", "XCoordinate", "YCoordinate",
                        "Quadrant", "Distance", "Freezing",
                        "CenterDistance", "PortDistance", "OpenField")

  # needs to be for loop because it requires values ahead and behind it
  for (i in 1:length(rownames(series))){
    # get rid of signal dropouts
    if (i > 1){
      if (abs(series[i,2]) > 30){
        series[i,2:3] <- series[i-1,2:3]
      } else if (abs(series[i,3]) > 30){
        series[i,2:3] <- c(mean(c(series[i - 1,2], series[i + 1,2])),
                           mean(c(series[i - 1,3], series[i + 1,3])))
      }
    } else if (i == 1){
      if (abs(series[i,2]) > 30){
        series[i,2:3] <- 0
      } else if (abs(series[i,3]) > 30){
        series[i,2:3] <- 0
      }
    }

    # rotate the data if quadrant is not LR to align properly

    if (series[i,2] > 0){
      if (series[i,3] > 0){
        series[i,4] <- "UR"
      } else {
        series[i,4] <- "LR"
      }
    } else {
      if (series[i,3] > 0){
        series[i,4] <- "UL"
      } else {
        series[i,4] <- "LL"
      }
    }

    # give instantaneous velocity at each time point
    if (i > 1){
      series[i,5] <- sqrt((series[i,2] - series[(i-1),2])^2 + (series[i,3] - series[(i-1),3])^2)
    } else if (i == 0){
      series[i,5] <- 0
    }

    # assign freezing
    if (series[i,5] < 1) {
      # need to assign 1 to this to ensure a long string of unique but repeated values
      series[i,6] <- 1
    } else {
      series[i,6] <- runif(1, max = .99)
    }

    # find center distance
    series[i,7] <- sqrt(series[i,2]^2 + series[i,3]^2)

    # find distance from the port
    series[i,8] <- sqrt(abs(series[i,2] - max(series$XCoordinate))^2 +
                          abs(series[i,3] - min(series$YCoordinate))^2)

    if (series[i,2] > -15) {
      if (series[i,2] < 15) {
        if (series[i,3] > -15) {
          if (series[i,3] < 15) {
            series[i,9] <- "Center"
          }
          else {
            series[i,9] <- "Surround"
          }
        }
        else {
          series[i,9] <- "Surround"
        }
      }
      else {
        series[i,9] <- "Surround"
      }
    }
    else {
      series[i,9] <- "Surround"
    }
  }

  # add in freezing lasting for at least one second
  mintime <- rle(series[, 6]) # only count freezing that lasts one second or more

  newreps <- NULL
  for (i in 1:length(mintime$lengths)){
    repeat_block <- rep(mintime$lengths[i], times = mintime$lengths[i])
    newreps <- c(newreps, repeat_block)
  }
  series[, 6] <- ifelse(newreps >= 4, "Immobile", "Mobile")

  # break into baseline and treatment at 10 minute mark
  breakpoint <- findInterval(600, series$Time)
  startpoint <- findInterval(start_min*60, series$Time)
  endpoint <- findInterval(end_min*60, series$Time)

  baseline_series <- series[1:breakpoint,]
  treatment_series <- series[startpoint:endpoint,]

  list_output <- list(series, baseline_series, treatment_series)
  names(list_output) <- c("full", "baseline", "treatment")
  return(list_output)
}

process_series_epm <- function(path, mode){

  series <- read.csv(path, sep = "\t")

  series <- cbind.data.frame(series,
                             integer(length = length(rownames(series))),
                             integer(length = length(rownames(series))),
                             integer(length = length(rownames(series))),
                             integer(length = length(rownames(series))))
  colnames(series) <- c("Time", "XCoordinate", "YCoordinate",
                        "Region", "Distance", "OpenEntry", "ClosedEntry")

  # needs to be for loop because it requires values ahead and behind it
  for (i in 1:length(rownames(series))){

    # get rid of signal dropouts
    if (i > 1){
      if (sqrt((series[i,2] - series[(i-1),2])^2 + (series[i,3] - series[(i-1),3])^2) > 10){
        series[i,2:3] <- series[i-1,2:3]
      }
    }

    if (i > 1){
      if (abs(series[i,2]) > 2.5) {
        if (abs(series[i,3]) > 2.5) {
          series[i,2:3] <- series[i-1,2:3]
        }
      }
    }

    if (abs(series[i,2]) > 2.5){
      series[i,4] <- "Open"
    } else if (abs(series[i,3]) > 2.5){
      series[i,4] <- "Closed"
    }
    else {
      series[i,4] <- "Center"
    }

    # give instantaneous velocity at each time point
    if (i > 1){
      series[i,5] <- sqrt((series[i,2] - series[(i-1),2])^2 + (series[i,3] - series[(i-1),3])^2)
    } else if (i == 1){
      series[i,5] <- 0
    }

    # mark open arm entries
    if (i > 1){
      if (series[i,4] == "Open") {
        if (series[i-1,4] != "Open") {
          series[i,6] <- 1
        } else {
          series[i,6] <- 0
        }
      } else {
        series[i,6] <- 0
      }
    } else if (i == 1){
      series[i,6] <- 0
    }

    # mark closed arm entries
    if (i > 1){
      if (series[i,4] == "Closed") {
        if (series[i-1,4] != "Closed") {
          series[i,7] <- 1
        } else {
          series[i,7] <- 0
        }
      } else {
        series[i,7] <- 0
      }
    } else if (i == 1){
      series[i,7] <- 0
    }
  }
  # break into baseline, treatment, and baseline again at 5 and 10 minute mark
  if (mode == "stim"){

    baseline1_series <- series[1:1200,]
    treatment_series <- series[1201:2400,]
    baseline2_series <- series[2401:3600,]

    baseline_all_series <- rbind.data.frame(baseline1_series, baseline2_series)

    minbymin <- vector("list", length = 15)
    for (i in 1:15){
      minbymin[[i]] <- series[(((i-1)*240)+1):(i*240),]
      names(minbymin)[i] <- i
    }

    list_output <- list(series, baseline1_series, treatment_series, baseline2_series, baseline_all_series, minbymin)
    names(list_output) <- c("full", "baseline_pre", "treatment", "baseline_post", "baseline_all", "each_min")

  } else if (mode == "silence"){

    series <- series[1:2400,]
    list_output <- list(series)
    names(list_output) <- "full"

  }

  return(list_output)
}

process_series_openfield <- function(path, mode){

  series <- read.csv(path, sep = "\t")

  series <- cbind.data.frame(series,
                             integer(length = length(rownames(series))),
                             integer(length = length(rownames(series))),
                             integer(length = length(rownames(series))),
                             integer(length = length(rownames(series))))
  colnames(series) <- c("Time", "XCoordinate", "YCoordinate",
                        "Region", "Distance", "Freezing", "Corner")

  # needs to be for loop because it requires values ahead and behind it
  for (i in 1:length(rownames(series))){

    # fix other coordinates from labview code, bring to zero than transform to centimeters
    series[i,2] <- series[i,2] * 25.2 # old conversion factor
    series[i,2] <- series[i,2] + 531 # old center
    series[i,2] <- series[i,2] - 585 # new center, midpoint varies between 330 and 340
    series[i,2] <- series[i,2] * (27.3/500) # conversion factor, is 500 pixels wide and is 27.3 cm wide

    series[i,3] <- series[i,3] * -25.2
    series[i,3] <- series[i,3] + 520
    series[i,3] <- series[i,3] - 550
    series[i,3] <- series[i,3] * (27.3/500)
  }

  for (i in 1:length(rownames(series))){
    # get rid of signal dropouts
    if (i > 1){
      if (sqrt((series[i,2] - series[(i-1),2])^2 + (series[i,3] - series[(i-1),3])^2) > 10){
        series[i,2:3] <- series[i-1,2:3]
      }
    }

    # sometimes counts on walls, need to adjust to control amount of motion
    if (i > 1){
      if (abs(series[i,2]) > 13.65){
        series[i,2] <- series[i-1,2]
      }
      if (abs(series[i,3]) > 13.65) {
        series[i,3] <- series[i-1,3]
      }
    } else if (i == 1){
      if (abs(series[i,2]) > 13.65){
        series[i,2] <- 13.64
      }
      if (abs(series[i,3]) > 13.65) {
        series[i,3] <- 13.64
      }
    }

    # center is 50% of area, or 19.3 cm square within 27.3 cm square
    if (abs(series[i,2]) < 6.83){
      if (abs(series[i,3]) < 6.83){
        series[i,4] <- "Center"
      } else {
        series[i,4] <- "Surround"
      }
    } else {
      series[i,4] <- "Surround"
    }

    # give instantaneous velocity at each time point
    if (i > 1){
      series[i,5] <- sqrt((series[i,2] - series[(i-1),2])^2 + (series[i,3] - series[(i-1),3])^2)
    } else if (i == 1){
      series[i,5] <- 0
    }

    # assign freezing
    if (series[i,5] < 0.5) {
      # need to assign 1 to this to ensure a long string of unique but repeated values
      series[i,6] <- 1
    } else {
      series[i,6] <- runif(1, max = .99)
    }

    # add corners for additional information
    if (abs(series[i,2]) > 6.83){
      if (abs(series[i,3]) > 6.83){
        series[i,7] <- "Corner"
      } else {
        series[i,7] <- "Non-corner"
      }
    } else {
      series[i,7] <- "Non-corner"
    }

    # sometimes yields 0, need to fix
    if (series[i,7] == 0){
      series[i,7] <- series[i-1,7]
    }
  }

  # add in freezing lasting for at least one second
  mintime <- rle(series[, 6]) # only count freezing that lasts one second or more

  newreps <- NULL
  for (i in 1:length(mintime$lengths)){
    repeat_block <- rep(mintime$lengths[i], times = mintime$lengths[i])
    newreps <- c(newreps, repeat_block)
  }
  series[, 6] <- ifelse(newreps >= 4, "Immobile", "Mobile")

  # break into baseline, treatment, and baseline again at 5 and 10 minute mark
  if (mode == 'stim'){

    baseline1_series <- series[1:1200,]
    treatment1_series <- series[1201:2400,]
    baseline2_series <- series[2401:3600,]
    treatment2_series <- series[3601:4800,]
    baseline3_series <- series[4801:6000,]

    baseline_all_series <- rbind.data.frame(baseline1_series, baseline2_series, baseline3_series)
    treatment_all_series <- rbind.data.frame(treatment1_series, treatment2_series)

    minbymin <- vector("list", length = 25)
    for (i in 1:25){
      minbymin[[i]] <- series[(((i-1)*240)+1):(i*240),]
      names(minbymin)[i] <- i
    }

    list_output <- list(series,
                        baseline1_series, treatment1_series, baseline2_series, treatment2_series, baseline3_series,
                        baseline_all_series, treatment_all_series, minbymin)
    names(list_output) <- c("full", "baseline_1", "treatment_1", "baseline_2", "treatment_2", "baseline_3",
                            "baseline_all", "treatment_all", "each_min")

  } else if (mode == 'silence'){

    series <- series[1:2400,]
    list_output <- list(series)
    names(list_output) <- "full"
  }

  return(list_output)
}

add_custom_names <- function(array, path, group_id, group_shorthand){

  # get names of datasets within groups
  name_vector <- c(Sys.glob(paste0(path, group_id, "/*")))
  name_vector <- unlist(strsplit(name_vector, paste0(path, group_id, "/")))[c(FALSE, TRUE)]

  # add group identifier for downstream transformation
  name_vector <- paste0(group_shorthand, "|", name_vector)

  names(array[[group_id]]) <- name_vector

  return(array)
}

merge_all <- function(x, y) {
  merge(x, y, all = TRUE, by= "Time")
}

## READ DATA

read_stim_data <- function(assay, groups, path, labels){

  data_structure <- vector(mode = "list", length = 4)

  if (assay == "4quad"){
    data_structure <- lapply(groups, function(x) lapply(Sys.glob(paste0(path, x, "/*")), process_series))
  } else if (assay == "epm"){
    data_structure <- lapply(groups, function(x) lapply(Sys.glob(paste0(path, x, "/*")),
                                                        function(y) process_series_epm(y, "stim")))
  } else if (assay == "oft"){
    data_structure <- lapply(groups, function(x) lapply(Sys.glob(paste0(path, x, "/*")),
                                                        function(y) process_series_openfield(y, "stim")))
  }

  names(data_structure) <- Sys.glob(paste0(path, "*"))
  names(data_structure) <- unlist(strsplit(names(data_structure), path))[c(FALSE, TRUE)]

  for (i in 1:length(groups)){
    data_structure <- add_custom_names(data_structure, path, groups[i], labels[i])
  }

  return(data_structure)
}

## STATS

display_anova <- function(anova_dataframe, metric){

  group_anova <- get_anova_table(anova_test(data = anova_dataframe,
                                            dv = .data[[metric]],
                                            wid = id,
                                            between = group,
                                            within = time),
                                 correction = "auto")

  return(group_anova)
}

pairwise_comp_anova <- function(anova_dataframe, metric){

  pairwise_formula <- reformulate("time", metric)

  pairwise_group_anova <- anova_dataframe %>%
    group_by(group) %>%
    tukey_hsd(pairwise_formula, paired = TRUE,
              p.adjust.method = "bonferroni")

  return(pairwise_group_anova)

}

get_pi_summary_stats <- function(array){

  # isolate base metric of interest
  baseline_array <- lapply(vg_data, function(x) lapply(x, function(y) y[["baseline"]][,c("Quadrant")]))
  treatment_array <- lapply(vg_data, function(x) lapply(x, function(y) y[["treatment"]][,c("Quadrant")]))

  baseline_array <- lapply(baseline_array, function(x) lapply(x, function(y) sum(y == "LR") / sum(y != 0)))
  treatment_array <- lapply(treatment_array, function(x) lapply(x, function(y) sum(y == "LR") / sum(y != 0)))

  baseline_array <- lapply(baseline_array, function(x) lapply(x, function(y) (y - 0.25) / 0.25))
  treatment_array <- lapply(treatment_array, function(x) lapply(x, function(y) (y - 0.25) / 0.25))

  pi_dataframe_anova <- cbind.data.frame(as.factor(unlist(lapply(vg_data, names))),
                                         unlist(lapply(vg_data,
                                                       function(x) strsplit(names(x), "[|]")))[c(TRUE, FALSE)],
                                         unlist(baseline_array), unlist(treatment_array))
  colnames(pi_dataframe_anova) <- c("id", "group", "baseline_pi", "treatment_pi")

  pi_dataframe_anova <- pi_dataframe_anova %>%
    gather(key = "time", value = "PI", baseline_pi, treatment_pi) %>%
    convert_as_factor(group, time)

  pi_dataframe_anova %>%
    group_by(group, time) %>%
    get_summary_stats(PI, type = "mean_se")

  return(pi_dataframe_anova)
}

## TRANSFORM FUNCTIONS

transform_pi <- function(array_list, bin, metric, start, end){

  # isolate base metric of interest
  array_list <- lapply(array_list, function(x) lapply(x, function(y) y[[bin]][start:end, metric]))
  # figure out proportion of time spent in the quadrant of interest
  array_list <- lapply(array_list, function(x) lapply(x, function(y) sum(y == "LR") / sum(y != 0)))
  # transform into PI
  array_list <- lapply(array_list, function(x) lapply(x, function(y) (y - 0.25) / 0.25))

  return(array_list)
}

transform_of <- function(array_list, bin, metric, start, end){

  # isolate base metric of interest
  array_list <- lapply(array_list, function(x) lapply(x, function(y) y[[bin]][start:end, metric]))
  # figure out proportion of time spent in the quadrant of interest
  array_list <- lapply(array_list, function(x) lapply(x, function(y) sum(y == "Center") / sum(y != 0)))

  return(array_list)
}

transform_pd <- function(array_list, bin, metric, start, end){

  # isolate base metric of interest
  array_list <- lapply(array_list, function(x) lapply(x, function(y) y[[bin]][start:end, metric]))
  # figure out proportion of time spent in the quadrant of interest
  array_list <- lapply(array_list, function(x) lapply(x, mean))

  return(array_list)
}

transform_epm_opentime <- function(array_list, bin, metric, start, end){

  # isolate base metric of interest
  array_list <- lapply(array_list, function(x) lapply(x, function(y) y[[bin]][start:end, metric]))
  # figure out proportion of time spent in the open arm
  array_list <- lapply(array_list, function(x) lapply(x, function(y) sum(y == "Open") / sum(y != 0)))

  return(array_list)
}

transform_epm_opentime_bymin <- function(array_list, metric){

  # isolate base metric of interest
  array_list <- lapply(array_list, function(x) lapply(x, function(y) lapply(y$each_min, function(z) z[[metric]])))
  # figure out proportion of time spent in the open arm
  array_list <- lapply(array_list, function(x) lapply(x, function(y) lapply(y, function(z) sum(z == "Open") / sum(z != 0))))
  # change to vectors
  array_list <- lapply(array_list, function(x) lapply(x, unlist))
  # change to arrays
  array_list <- lapply(array_list, function(x) do.call(cbind.data.frame, x))
  mins <- c(1:15)
  array_list <- lapply(array_list, function(x) cbind.data.frame(mins, x))
  # change to ggplot-ready arrays
  array_list <- lapply(array_list, function(x) x %>% pivot_longer(cols = !mins, names_to = "id", values_to = "OpenTime"))

  return(array_list)
}

transform_epm_openentry <- function(array_list, bin, metric, start, end){

  # isolate base metric of interest
  array_list <- lapply(array_list, function(x) lapply(x, function(y) y[[bin]][start:end, metric]))
  # figure out proportion of time spent in the open arm
  array_list <- lapply(array_list, function(x) lapply(x, function(y) sum(y)))

  return(array_list)
}

transform_epm_openentry_bymin <- function(array_list, metric){

  # isolate base metric of interest
  array_list <- lapply(array_list, function(x) lapply(x, function(y) lapply(y$each_min, function(z) z[[metric]])))
  # figure out proportion of time spent in the open arm
  array_list <- lapply(array_list, function(x) lapply(x, function(y) lapply(y, function(z) sum(z))))
  # change to vectors
  array_list <- lapply(array_list, function(x) lapply(x, unlist))
  # change to arrays
  array_list <- lapply(array_list, function(x) do.call(cbind.data.frame, x))
  mins <- c(1:15)
  array_list <- lapply(array_list, function(x) cbind.data.frame(mins, x))
  # change to ggplot-ready arrays
  array_list <- lapply(array_list, function(x) x %>% pivot_longer(cols = !mins, names_to = "id", values_to = "OpenEntry"))

  return(array_list)
}

transform_openfield_centertime <- function(array_list, bin, metric, start, end){

  # isolate base metric of interest
  array_list <- lapply(array_list, function(x) lapply(x, function(y) y[[bin]][start:end, metric]))
  # figure out proportion of time spent in the open arm
  array_list <- lapply(array_list, function(x) lapply(x, function(y) sum(y == "Center") / sum(y != 0)))

  return(array_list)
}

transform_openfield_centertime_bymin <- function(array_list, metric){

  # isolate base metric of interest
  array_list <- lapply(array_list, function(x) lapply(x, function(y) lapply(y$each_min, function(z) z[[metric]])))
  # figure out proportion of time spent in the center region
  array_list <- lapply(array_list, function(x) lapply(x, function(y) lapply(y, function(z) sum(z == "Center") / sum(z != 0))))
  # change to vectors
  array_list <- lapply(array_list, function(x) lapply(x, unlist))
  # change to arrays
  array_list <- lapply(array_list, function(x) do.call(cbind.data.frame, x))
  mins <- c(1:25)
  array_list <- lapply(array_list, function(x) cbind.data.frame(mins, x))
  # change to ggplot-ready arrays
  array_list <- lapply(array_list, function(x) x %>% pivot_longer(cols = !mins, names_to = "id", values_to = "CenterTime"))

  return(array_list)
}

transform_openfield_cornertime <- function(array_list, bin, metric, start, end){

  # isolate base metric of interest
  array_list <- lapply(array_list, function(x) lapply(x, function(y) y[[bin]][start:end, metric]))
  # figure out proportion of time spent in the open arm
  array_list <- lapply(array_list, function(x) lapply(x, function(y) sum(y == "Corner") / sum(y != 0)))

  return(array_list)
}

transform_openfield_cornertime_bymin <- function(array_list, metric){

  # isolate base metric of interest
  array_list <- lapply(array_list, function(x) lapply(x, function(y) lapply(y$each_min, function(z) z[[metric]])))
  # figure out proportion of time spent in the center region
  array_list <- lapply(array_list, function(x) lapply(x, function(y) lapply(y, function(z) sum(z == "Corner") / sum(z != 0))))
  # change to vectors
  array_list <- lapply(array_list, function(x) lapply(x, unlist))
  # change to arrays
  array_list <- lapply(array_list, function(x) do.call(cbind.data.frame, x))
  mins <- c(1:25)
  array_list <- lapply(array_list, function(x) cbind.data.frame(mins, x))
  # change to ggplot-ready arrays
  array_list <- lapply(array_list, function(x) x %>% pivot_longer(cols = !mins, names_to = "id", values_to = "CornerTime"))

  return(array_list)
}

transform_openfield_mobility <- function(array_list, bin, metric, start, end){

  # isolate base metric of interest
  array_list <- lapply(array_list, function(x) lapply(x, function(y) y[[bin]][start:end, metric]))
  # figure out proportion of time spent in the open arm
  array_list <- lapply(array_list, function(x) lapply(x, function(y) sum(y == "Immobile") / sum(y != 0)))

  return(array_list)
}

transform_openfield_mobility_bymin <- function(array_list, metric){

  # isolate base metric of interest
  array_list <- lapply(array_list, function(x) lapply(x, function(y) lapply(y$each_min, function(z) z[[metric]])))
  # figure out proportion of time spent in the center region
  array_list <- lapply(array_list, function(x) lapply(x, function(y) lapply(y, function(z) sum(z == "Immobile") / sum(z != 0))))
  # change to vectors
  array_list <- lapply(array_list, function(x) lapply(x, unlist))
  # change to arrays
  array_list <- lapply(array_list, function(x) do.call(cbind.data.frame, x))
  mins <- c(1:25)
  array_list <- lapply(array_list, function(x) cbind.data.frame(mins, x))
  # change to ggplot-ready arrays
  array_list <- lapply(array_list, function(x) x %>% pivot_longer(cols = !mins, names_to = "id", values_to = "Freezing"))

  return(array_list)
}

transform_openfield_distance_bymin <- function(array_list, metric){

  # isolate base metric of interest
  array_list <- lapply(array_list, function(x) lapply(x, function(y) lapply(y$each_min, function(z) z[[metric]])))
  # figure out proportion of time spent in the open arm
  array_list <- lapply(array_list, function(x) lapply(x, function(y) lapply(y, function(z) sum(z))))
  # change to vectors
  array_list <- lapply(array_list, function(x) lapply(x, unlist))
  # change to arrays
  array_list <- lapply(array_list, function(x) do.call(cbind.data.frame, x))
  mins <- c(1:25)
  array_list <- lapply(array_list, function(x) cbind.data.frame(mins, x))
  # change to ggplot-ready arrays
  array_list <- lapply(array_list, function(x) x %>% pivot_longer(cols = !mins, names_to = "id", values_to = "Distance"))

  return(array_list)
}

## GROUP FUNCTIONS

linear_array <- function(array_list, baseline_array, treatment_array, label){

  line_array <- cbind.data.frame(unlist(lapply(array_list,
                                               function(x) strsplit(names(x), "[|]")))[c(TRUE, FALSE)],
                                 as.numeric(paste0("-", str_sub(unlist(lapply(data_topography,
                                                                              function(x) strsplit(names(x), "[|]")))[c(FALSE, TRUE)], start = 1, end = 4))),
                                 unlist(baseline_array), unlist(treatment_array),
                                 unlist(treatment_array) - unlist(baseline_array))
  colnames(line_array) <- c("group", "ap_coords",
                            paste0("baseline_", label), paste0("treatment_", label), paste0("difference_", label))

  line_array <- cbind.data.frame(line_array, unlist(strsplit(line_array$group, " "))[c(FALSE, TRUE)])
  colnames(line_array)[6] <- "type"
  line_array$type <- as.factor(line_array$type)

  return(line_array)
}

group_array <- function(array_list, baseline_array, treatment_array, label){

  grouped_array <- cbind.data.frame(as.factor(unlist(lapply(array_list, names))),
                                    unlist(lapply(array_list,
                                                  function(x) strsplit(names(x), "[|]")))[c(TRUE, FALSE)],
                                    unlist(baseline_array), unlist(treatment_array))
  colnames(grouped_array) <- c("id", "group", paste0("baseline_", label), paste0("treatment_", label))

  return(grouped_array)
}

group_array_epm <- function(array_list, baseline_pre_array, treatment_array, baseline_post_array, label){

  grouped_array <- cbind.data.frame(as.factor(unlist(lapply(array_list, names))),
                                    unlist(lapply(array_list,
                                                  function(x) strsplit(names(x), "[|]")))[c(TRUE, FALSE)],
                                    unlist(baseline_pre_array), unlist(treatment_array), unlist(baseline_post_array))
  colnames(grouped_array) <- c("id", "group",
                               paste0("baseline_pre_", label), paste0("treatment_", label), paste0("baseline_post_", label))

  return(grouped_array)
}

group_array_openfield <- function(array_list, baseline_1_array, treatment_1_array,
                                  baseline_2_array, treatment_2_array, baseline_3_array, label){

  grouped_array <- cbind.data.frame(as.factor(unlist(lapply(array_list, names))),
                                    unlist(lapply(array_list,
                                                  function(x) strsplit(names(x), "[|]")))[c(TRUE, FALSE)],
                                    unlist(baseline_1_array), unlist(treatment_1_array), unlist(baseline_2_array),
                                    unlist(treatment_2_array), unlist(baseline_3_array))
  colnames(grouped_array) <- c("id", "group",
                               paste0("baseline_1_", label), paste0("treatment_1_", label), paste0("baseline_2_", label),
                               paste0("treatment_2_", label), paste0("baseline_3_", label))

  return(grouped_array)
}
