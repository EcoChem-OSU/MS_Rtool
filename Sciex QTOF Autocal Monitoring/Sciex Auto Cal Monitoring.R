###################################################
## Visualizing the Sciex QTOF instrument (7600 or else) auto-calibration results during a specified running period
#PPM, Resolution, intensity, etc for different calibrants 
####################################################
## Sciex Autocal Monitoring CODE
## version:1.0
## Date: 2025-08-21
## Author: Jiaju Fu @ Oregon State University
###############################################################################
## Description (points to be clear before using this script):

  ### point1- ## This script is using the auto-cal results files (.csv docs) generated from SCIEX QTOF during batches and batches acquisition as input files
  #(For 7600, the auto-cal results files are all already stored in the folder of D:\SCIEX OS Data\7600\Data\Cal
  ### point2- ## If you want to check auto-cal directly in the instrument OS, please download and store this script in the same folder of auto-cal results files. 
  #Or pick up some results files to you own PC
  ### point3- ## All the output plots would be stored in a new sub-folder named "Output_Polarity_Username_Checkdate"


############# (Need to be specified by user)#############################################################################################
 # Set1: Working directory (i.e., autocal file and script storage directory)
workdir <- "/Volumes/research/EMT/ALS1113/share/Jiaju Fu/Field LAB/7600 autocal monitoring/Sciex Auto Cal Monitoring"
setwd(workdir)
data_dir <- "."
 # Set2: Choose which polarity to analyze: "positive" or "negative" ---
polarity <- "negative"

 # Set3: User name for output folder (leave "" for no name) ----
Username <- "Jiaju"   # e.g., "Jiaju"; if "", username part is omitted

################### File folder naming format: Not suggest to change ##########################
# ---- Output folder (all exports go here) with polarity + username + today's date ----
checkdate_str <- format(Sys.Date(), "%Y%m%d")

# sanitize username to be filesystem-friendly
username_clean <- gsub("[^A-Za-z0-9]+", "_", Username)
username_clean <- gsub("^_|_$", "", username_clean)
username_part  <- if (nzchar(username_clean)) paste0("_", username_clean) else ""

# Use the existing `polarity` variable ("positive" / "negative") and title-case it
polarity_title <- if (tolower(polarity) == "positive") "Positive" else "Negative"

# Build output directory name: e.g., Output_Negative_Jiaju_Checkdate_20250821
out_dir_name <- paste0("Output_", polarity_title, username_part, "_Checkdate_", checkdate_str)
out_dir <- file.path(data_dir, out_dir_name)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
###############################################################################################

 # Set4: Date range you want to check for detail
start_date <- as.Date("2024-02-01")  # inclusive
end_date   <- as.Date("2100-12-31")  # inclusive

 # Set5: Explicitly declare whether to draw normal and cumulative plots for each variable 
##(JUST need to specify the TRUE or FALSE options below)

# Available column names (case-sensitive):
# "m/z Before Calibration (Da)",
# "Delta m/z Before Calibration (Da)",
# "PPM Error Before Calibration",
# "m/z After Calibration (Da)",
# "Delta m/z After Calibration (Da)",
# "PPM Error After Calibration",
# "Intensity (cps)",
# "Resolution"

plot_variables <- tibble::tibble(
  variable = c(
    "m/z Before Calibration (Da)",
    "Delta m/z Before Calibration (Da)",
    "PPM Error Before Calibration",
    "m/z After Calibration (Da)",
    "Delta m/z After Calibration (Da)",
    "PPM Error After Calibration",
    "Intensity (cps)",
    "Resolution"
  ),
  plot        = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE,  TRUE,  TRUE),#####Specify here for the basic plot
  plot_cumsum = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE,  FALSE,  FALSE)#####Specify here for the cumulative plot
)

 # Set6: x-axis tick label font size (user configurable) ======
x_axis_text_size <- 9   # e.g., 9/10/11/12...

#######################################################################################################################################


## SCRIPT STARTS HERE## To run anyway
#########################################################################################################################################################

# ====== Libraries ======
library(dplyr)
library(stringr)
library(tidyr)
library(lubridate)
library(purrr)
library(ggplot2)
library(writexl)
library(patchwork)

# ====== I/O paths ======
data_dir <- "."   # current working directory; change to your CSV folder if needed
pattern_str <- "^Cal-\\d{4}-\\d{2}-\\d{2}-\\d{2}-\\d{2}-\\d{2}\\.csv$"


# ====== List files & date filter ======
files <- list.files(path = data_dir, pattern = pattern_str, full.names = TRUE)

get_file_date <- function(filename) {
  parts <- stringr::str_match(basename(filename), "Cal-(\\d{4})-(\\d{2})-(\\d{2})")
  if (any(is.na(parts))) return(NA)
  as.Date(sprintf("%s-%s-%s", parts[2], parts[3], parts[4]))
}
file_dates <- sapply(files, get_file_date)
files <- files[!is.na(file_dates) & file_dates >= start_date & file_dates <= end_date]

# print summary of files to read
cat(sprintf("Found %d file(s) between %s and %s. Target polarity: %s\n",
            length(files), format(start_date), format(end_date), tolower(polarity)))

# ====== Hold extracted tables ======
all_data <- list(TOF_MS = list(), TOF_MSMS = list(), TOF_Zeno = list())

# ====== Extract one block between two anchors ======
extract_block <- function(lines_split, start_line, end_line, file_time) {
  block <- lines_split[(start_line + 1):(end_line - 1)]
  # header row starts with "Used"
  header_idx <- which(sapply(block, function(x) length(x) > 0 && x[1] == "Used"))
  if (length(header_idx) == 0) return(NULL)
  
  header <- block[[header_idx[1]]]
  data_lines <- block[(header_idx[1] + 1):length(block)]
  # keep rows whose first column is "Yes"
  data_lines <- Filter(function(x) length(x) > 0 && x[1] == "Yes", data_lines)
  if (length(data_lines) == 0) return(NULL)
  
  df <- as.data.frame(do.call(rbind, data_lines), stringsAsFactors = FALSE)
  colnames(df) <- header
  df$Time <- file_time
  df
}

# helper to get first match index or NA
first_idx <- function(v) {
  w <- which(v)
  if (length(w)) w[1] else NA_integer_
}

# ====== Progress bar over files (with per-file message) ======
n_files <- length(files)
if (n_files == 0) {
  warning("No CSV matched the pattern/date range in: ", normalizePath(data_dir))
} else {
  pb <- utils::txtProgressBar(min = 0, max = n_files, style = 3)
  on.exit(close(pb), add = TRUE)
  
  # choose parent m/z anchors by desired polarity
  parent_mz <- if (tolower(polarity) == "positive") "609.2807" else "928.8344"
  parent_key <- paste0("TOF MS/MS Parent = ", parent_mz, " Da")
  zeno_key   <- paste0("TOF MS/MS Zeno Parent = ", parent_mz, " Da")
  
  for (i in seq_along(files)) {
    file <- files[i]
    cat(sprintf("Reading (%d/%d): %s\n", i, n_files, basename(file)))
    
    # timestamp for this file
    time_str <- stringr::str_match(
      basename(file),
      "Cal-(\\d{4})-(\\d{2})-(\\d{2})-(\\d{2})-(\\d{2})-(\\d{2})"
    )[, 2:7]
    file_time <- lubridate::ymd_hms(paste(time_str, collapse = "-"))
    
    # read lines and split by comma
    lines <- readLines(file, warn = FALSE)
    lines_trimmed <- trimws(lines)
    lines_split <- strsplit(lines_trimmed, ",", fixed = TRUE)
    
    # first column strings of each line
    h1 <- sapply(lines_split, function(x) x[1])
    
    # detect this file's polarity from content
    file_pol <- if ("Polarity: Positive" %in% h1) {
      "positive"
    } else if ("Polarity: Negative" %in% h1) {
      "negative"
    } else {
      NA_character_
    }
    
    if (is.na(file_pol)) {
      warning("Polarity not found in file, skipped: ", basename(file))
      utils::setTxtProgressBar(pb, i); next
    }
    if (tolower(file_pol) != tolower(polarity)) {
      cat(sprintf("  Skipped (polarity %s, want %s)\n", file_pol, tolower(polarity)))
      utils::setTxtProgressBar(pb, i); next
    }
    
    # anchor lines (first occurrence), using m/z that matches chosen polarity
    idxs <- list(
      tof  = first_idx(h1 == "TOF MS"),
      msms = first_idx(h1 == parent_key),
      zeno = first_idx(h1 == zeno_key),
      end  = first_idx(h1 == "Calibration Parameters Summary")
    )
    
    if (any(is.na(unlist(idxs)))) {
      warning("Missing anchors for chosen polarity, skipped: ", basename(file))
      utils::setTxtProgressBar(pb, i)
      next
    }
    
    # extract each section
    tof   <- extract_block(lines_split, idxs$tof,  idxs$msms, file_time)
    msms  <- extract_block(lines_split, idxs$msms, idxs$zeno, file_time)
    zeno  <- extract_block(lines_split, idxs$zeno, idxs$end,  file_time)
    
    if (!is.null(tof))  all_data$TOF_MS[[length(all_data$TOF_MS) + 1]]     <- tof
    if (!is.null(msms)) all_data$TOF_MSMS[[length(all_data$TOF_MSMS) + 1]] <- msms
    if (!is.null(zeno)) all_data$TOF_Zeno[[length(all_data$TOF_Zeno) + 1]] <- zeno
    
    utils::setTxtProgressBar(pb, i)
  }
}

# ====== Export detail workbook: one sheet per m/z ======
export_split_by_mz <- function(data_list, out_file) {
  if (length(data_list) == 0) return(invisible(NULL))
  all_df <- dplyr::bind_rows(data_list)
  split_list <- split(all_df, all_df[["m/z Expected (Da)"]])
  names(split_list) <- make.names(names(split_list))
  writexl::write_xlsx(split_list, out_file)
}

export_split_by_mz(all_data$TOF_MS,   file.path(out_dir, "TOF_MS_Detail.xlsx"))
export_split_by_mz(all_data$TOF_MSMS, file.path(out_dir, "TOF_MSMS_Detail.xlsx"))
export_split_by_mz(all_data$TOF_Zeno, file.path(out_dir, "TOF_Zeno_Detail.xlsx"))

# ====== Filename cleaner (remove Da/cps and non-alnum) ======
clean_var_name <- function(var) {
  var %>%
    stringr::str_replace_all(stringr::regex("\\((da|cps)\\)", ignore_case = TRUE), "") %>%
    stringr::str_replace_all(stringr::regex("\\b(da|cps)\\b", ignore_case = TRUE), "") %>%
    stringr::str_replace_all("[^A-Za-z0-9]", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "") %>%
    stringr::str_trim()
}


# ---------- helper: build top year header (full coverage + small gaps only at year boundaries) ----------
year_header_layers <- function(t_min, t_max, y_line, y_text, gap_days = 6) {
  yrs <- seq(lubridate::floor_date(t_min, "year"),
             lubridate::floor_date(t_max, "year"),
             by = "1 year")
  if (length(yrs) == 0) return(list())
  
  gap <- lubridate::days(gap_days)
  layers <- list()
  
  for (k in seq_along(yrs)) {
    x0 <- max(yrs[k], t_min)
    x1 <- min(if (k < length(yrs)) yrs[k + 1] else lubridate::ceiling_date(t_max, "year"), t_max)
    
    # leave gap ONLY at internal year boundaries
    seg_start <- if (k == 1) x0 else x0 + gap
    seg_end   <- if (k == length(yrs)) x1 else x1 - gap
    
    if (seg_end > seg_start) {
      layers <- append(layers, list(
        ggplot2::annotate("segment", x = seg_start, xend = seg_end,
                          y = y_line, yend = y_line,
                          colour = "#888888", linewidth = 0.4)
      ))
    }
    
    lab_x <- as.POSIXct((as.numeric(x0) + as.numeric(x1)) / 2,
                        origin = "1970-01-01", tz = lubridate::tz(t_min))
    layers <- append(layers, list(
      ggplot2::annotate("text", x = lab_x, y = y_text,
                        label = format(x0, "%Y"),
                        size = 3.2, colour = "#666666", vjust = 0.5)
    ))
  }
  layers
}



# ====== Plot: value vs time (bottom month labels + extra ticks + TOP year line/label) ======
plot_variable_by_mz <- function(data_list, variable, title_prefix, out_pdf) {
  if (length(data_list) == 0) return(invisible(NULL))
  
  df <- dplyr::bind_rows(data_list)
  mz_var <- "m/z Expected (Da)"
  df$Time <- lubridate::as_datetime(df$Time)
  
  # "Not Found" -> NA -> impute by mean for plotting (but color as Not Found)
  df[[variable]][df[[variable]] == "Not Found"] <- NA
  df[[variable]] <- as.numeric(df[[variable]])
  
  split_df <- split(df, df[[mz_var]])
  
  plots <- purrr::imap(split_df, function(subdf, mz_val) {
    subdf <- dplyr::arrange(subdf, Time)
    mean_val <- mean(subdf[[variable]], na.rm = TRUE)
    subdf$status <- factor(ifelse(is.na(subdf[[variable]]), "Not Found", "Normal"),
                           levels = c("Normal", "Not Found"))
    subdf$Value <- ifelse(is.na(subdf[[variable]]), mean_val, subdf[[variable]])
    
    # dynamic y-limits + headroom for year header
    yr <- range(subdf$Value, na.rm = TRUE); if (!all(is.finite(yr))) yr <- c(0, 1)
    span <- diff(yr)
    base_buf   <- if (span > 0) 0.05 * span else max(0.05 * abs(yr[1]), 0.1)
    header_pad <- if (span > 0) 0.10 * span else 0.2
    y_limits   <- c(yr[1] - base_buf, yr[2] + base_buf + header_pad)
    
    # month ticks along bottom
    mb_seq <- seq(lubridate::floor_date(min(subdf$Time, na.rm = TRUE), "month"),
                  lubridate::floor_date(max(subdf$Time, na.rm = TRUE), "month"),
                  by = "1 month")
    
    # year header placements
    t_min <- min(subdf$Time, na.rm = TRUE)
    t_max <- max(subdf$Time, na.rm = TRUE)
    y_line <- yr[2] + base_buf + 0.15 * header_pad
    y_text <- yr[2] + base_buf + 0.55 * header_pad
    year_layers <- year_header_layers(t_min, t_max, y_line, y_text, gap_days = 6)
    
    p <- ggplot(subdf, aes(x = Time, y = Value)) +
      geom_line(aes(group = 1), color = "#2E86AB", show.legend = FALSE) +
      geom_point(aes(color = status), shape = 19, size = 1.8, show.legend = TRUE) +
      scale_color_manual(
        name   = NULL,
        breaks = "Not Found",
        labels = c("Not Found" = "Not Found (replaced by aver.)"),
        values = c("Normal" = "#2E86AB", "Not Found" = "#E74C3C"),
        drop   = FALSE
      ) +
      guides(color = guide_legend(override.aes = list(linetype = NA, shape = 19, size = 3))) +
      scale_x_datetime(
        breaks = scales::date_breaks("1 month"),
        labels = scales::label_date("%b")
      ) +
      geom_rug(
        data = data.frame(Time = mb_seq),
        mapping = aes(x = Time),
        inherit.aes = FALSE,
        sides = "b",
        length = grid::unit(5, "pt"),
        color = "#666666",
        show.legend = FALSE
      ) +
      scale_y_continuous(limits = y_limits) +
      coord_cartesian(clip = "off") +
      year_layers +
      labs(title = paste0(title_prefix, " - m/z ", mz_val), x = "Time", y = variable) +
      theme_minimal(base_size = 12) +
      theme(
        plot.margin        = margin(t = 18, r = 10, b = 10, l = 10),
        axis.text.x        = element_text(size = x_axis_text_size, margin = margin(t = 2)),
        axis.title.x       = element_text(margin = margin(t = 4)),
        legend.position    = "bottom",
        legend.title       = element_blank(),
        legend.text        = element_text(size = 9),
        legend.key         = element_blank()
      )
    
    p
  })
  
  combined <- patchwork::wrap_plots(plots, ncol = 2, guides = "collect")
  ggsave(out_pdf, plot = combined,
         width = 12, height = ceiling(length(plots) / 2) * 3)
  invisible(NULL)
}

# ====== Plot: cumulative sum vs time (bottom month labels + extra ticks + TOP year line/label) ======
plot_variable_cumsum_by_mz <- function(data_list, variable, title_prefix, out_pdf) {
  if (length(data_list) == 0) return(invisible(NULL))
  
  df <- dplyr::bind_rows(data_list)
  mz_var <- "m/z Expected (Da)"
  df$Time <- lubridate::as_datetime(df$Time)
  
  df[[variable]][df[[variable]] == "Not Found"] <- NA
  df[[variable]] <- as.numeric(df[[variable]])
  
  split_df <- split(df, df[[mz_var]])
  
  plots <- purrr::imap(split_df, function(subdf, mz_val) {
    subdf <- dplyr::arrange(subdf, Time)
    mean_val <- mean(subdf[[variable]], na.rm = TRUE)
    subdf$status <- factor(ifelse(is.na(subdf[[variable]]), "Not Found", "Normal"),
                           levels = c("Normal", "Not Found"))
    subdf$Value  <- ifelse(is.na(subdf[[variable]]), mean_val, subdf[[variable]])
    subdf$Cumsum <- cumsum(subdf$Value)
    
    yr <- range(subdf$Cumsum, na.rm = TRUE); if (!all(is.finite(yr))) yr <- c(0, 1)
    span <- diff(yr)
    base_buf   <- if (span > 0) 0.05 * span else max(0.05 * abs(yr[1]), 0.1)
    header_pad <- if (span > 0) 0.10 * span else 0.2
    y_limits   <- c(yr[1] - base_buf, yr[2] + base_buf + header_pad)
    
    mb_seq <- seq(lubridate::floor_date(min(subdf$Time, na.rm = TRUE), "month"),
                  lubridate::floor_date(max(subdf$Time, na.rm = TRUE), "month"),
                  by = "1 month")
    
    t_min <- min(subdf$Time, na.rm = TRUE)
    t_max <- max(subdf$Time, na.rm = TRUE)
    y_line <- yr[2] + base_buf + 0.15 * header_pad
    y_text <- yr[2] + base_buf + 0.55 * header_pad
    year_layers <- year_header_layers(t_min, t_max, y_line, y_text, gap_days = 6)
    
    p <- ggplot(subdf, aes(x = Time, y = Cumsum)) +
      geom_line(aes(group = 1), color = "#2E86AB", show.legend = FALSE) +
      geom_point(aes(color = status), shape = 19, size = 1.8, show.legend = TRUE) +
      scale_color_manual(
        name   = NULL,
        breaks = "Not Found",
        labels = c("Not Found" = "Not Found (replaced by aver.)"),
        values = c("Normal" = "#2E86AB", "Not Found" = "#E74C3C"),
        drop   = FALSE
      ) +
      guides(color = guide_legend(override.aes = list(linetype = NA, shape = 19, size = 3))) +
      scale_x_datetime(
        breaks = scales::date_breaks("1 month"),
        labels = scales::label_date("%b")
      ) +
      geom_rug(
        data = data.frame(Time = mb_seq),
        mapping = aes(x = Time),
        inherit.aes = FALSE,
        sides = "b",
        length = grid::unit(5, "pt"),
        color = "#666666",
        show.legend = FALSE
      ) +
      scale_y_continuous(limits = y_limits) +
      coord_cartesian(clip = "off") +
      year_layers +
      labs(title = paste0(title_prefix, " - m/z ", mz_val),
           x = "Time", y = paste("Cumulative", variable)) +
      theme_minimal(base_size = 12) +
      theme(
        axis.title.y     = element_text(size = 10),
        plot.margin      = margin(t = 18, r = 10, b = 10, l = 10),
        axis.text.x      = element_text(size = x_axis_text_size, margin = margin(t = 2)),
        axis.title.x     = element_text(margin = margin(t = 4)),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.text      = element_text(size = 9),
        legend.key       = element_blank()
      )
    
    p
  })
  
  combined <- patchwork::wrap_plots(plots, ncol = 2, guides = "collect")
  ggsave(out_pdf, plot = combined,
         width = 12, height = ceiling(length(plots) / 2) * 3)
  invisible(NULL)
}

# ====== Iterate variables and export plots ======
for (i in seq_len(nrow(plot_variables))) {
  var <- plot_variables$variable[i]
  var_file <- clean_var_name(var)
  
  # Regular (non-cumulative) plots — only when plot == TRUE
  if (isTRUE(plot_variables$plot[i])) {
    plot_variable_by_mz(all_data$TOF_MS,   var, "TOF MS",
                        file.path(out_dir, paste0("TOF_MS_",   var_file, ".pdf")))
    plot_variable_by_mz(all_data$TOF_MSMS, var, "TOF MS/MS",
                        file.path(out_dir, paste0("TOF_MSMS_", var_file, ".pdf")))
    plot_variable_by_mz(all_data$TOF_Zeno, var, "TOF MS/MS Zeno",
                        file.path(out_dir, paste0("TOF_Zeno_", var_file, ".pdf")))
  }
  
  # Cumulative plots — evaluated independently;
  # even if plot == FALSE, draw cumulative plots when plot_cumsum == TRUE
  if (isTRUE(plot_variables$plot_cumsum[i])) {
    plot_variable_cumsum_by_mz(all_data$TOF_MS,   var, "TOF MS",
                               file.path(out_dir, paste0("TOF_MS_",   var_file, "_Cumsum.pdf")))
    plot_variable_cumsum_by_mz(all_data$TOF_MSMS, var, "TOF MS/MS",
                               file.path(out_dir, paste0("TOF_MSMS_", var_file, "_Cumsum.pdf")))
    plot_variable_cumsum_by_mz(all_data$TOF_Zeno, var, "TOF MS/MS Zeno",
                               file.path(out_dir, paste0("TOF_Zeno_", var_file, "_Cumsum.pdf")))
  }
}

#########################################################################################################################################################