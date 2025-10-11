###################################################
## Monitoring the QTOF instrument sensitivity across months or years (absolute response for target/internal standard) 
## Using the Area txt files generated from every obtained batch quantification results in SCIEX OS Analytics molecule
####################################################
## Sciex 7600 QTOF Sensitivity Monitoring CODE
## version:1.0
## Date: 2025-08-19
## Author: Jiaju Fu @ Oregon State University
###############################################################################
## Description (points to be clear before using this script):

### point1- In the workdir folder, require one subfolder named "Area logs" with .txt files inside 
##### for a Area file name example:20241003_Bag1_d1000000_DJM_AFFFs_A_Area.txt, alway making sure using the year/month/date (XXXX/XX/XX) as the name beginning (first 8 digits) and _Area as the name ending, other parts are totally customized by the user

### point2- making sure your QC samples info in the txt files is described as below
##### Sample Name should contain CCV for the high concentration continuing calibration verification standards sample and LLCCV for the low level continuing calibration verification standards sample
######## Sample Name suffix is acceptable such as LLCCV_01 CCV_01 LLCCV_02 CCV_02
##### Sample Type only for the CCV and LLCCV should be specified as Quality Control

#############Set working directory (need to be specified by user)#####################################################
workdir <- "/Volumes/research/EMT/ALS1113/share/Jiaju Fu/Field LAB/7600 sensitivity log/test"
###############################################################################

## START SCRIPT ## To run anyway
#################
setwd(workdir)

# ===== Load required packages =====
library(readr)
library(dplyr)
library(stringr)
library(readxl)
library(ggplot2)
library(tidyr)
library(writexl)

# ===== Define Output Folder =====
output_dir <- "Output"
dir.create(output_dir, showWarnings = FALSE)

# ===== List and sort input files in "Area logs" folder =====
all_files <- list.files("Area logs", pattern = "_Area\\.txt$", full.names = TRUE)

# Extract date from filename (first 8 digits as YYYYMMDD)
get_date_from_filename <- function(fname) {
  basename(fname) |> str_extract("^\\d{8}") |> as.character()
}

# Sort files by date
sorted_files <- all_files[order(sapply(all_files, get_date_from_filename))]

# ===== Initialize container and sample counter =====
combined_df <- data.frame()
sample_counter <- 0

# ===== Loop through each file to extract LLCCV and CCV samples =====
for (file in sorted_files) {
  df <- read_tsv(file, na = c("N/A", "NA", ""), show_col_types = FALSE)
  pfas_cols <- setdiff(colnames(df), c("Sample Name", "Sample Type"))
  df[pfas_cols] <- lapply(df[pfas_cols], as.numeric)
  df_qc <- df %>% filter(`Sample Type` == "Quality Control")
  
  df_llccv <- df_qc %>% filter(str_detect(`Sample Name`, "LLCCV"))
  df_ccv <- df_qc %>% filter(str_detect(`Sample Name`, "CCV"), !str_detect(`Sample Name`, "LLCCV"))
  
  if (nrow(df_llccv) == 0 & nrow(df_ccv) == 0) next
  
  file_date <- get_date_from_filename(file)
  file_mmdd <- format(as.Date(file_date, "%Y%m%d"), "%m-%d")
  file_year <- substr(file_date, 1, 4)
  
  max_n <- max(nrow(df_llccv), nrow(df_ccv))
  for (i in seq_len(max_n)) {
    sample_counter <- sample_counter + 1
    if (i <= nrow(df_llccv)) {
      df_llccv[i, "SampleIndex"] <- sample_counter
      df_llccv[i, "Group"] <- "LLCCV"
      df_llccv[i, "Date_MMDD"] <- file_mmdd
      df_llccv[i, "Year"] <- file_year
      df_llccv[i, "FileName"] <- basename(file)
    }
    if (i <= nrow(df_ccv)) {
      df_ccv[i, "SampleIndex"] <- sample_counter
      df_ccv[i, "Group"] <- "CCV"
      df_ccv[i, "Date_MMDD"] <- file_mmdd
      df_ccv[i, "Year"] <- file_year
      df_ccv[i, "FileName"] <- basename(file)
    }
  }
  
  combined_df <- bind_rows(combined_df, df_llccv, df_ccv)
}

# ===== Get PFAS compound column names =====
compound_cols <- setdiff(colnames(combined_df), c("Sample Name", "Sample Type", "Group", "SampleIndex", "Date_MMDD", "Year", "FileName"))

# ===== Plotting function for each compound =====
plot_combined_llccv_ccv <- function(df, compound, out_dir) {
  dir.create(out_dir, showWarnings = FALSE)
  
  plot_df <- df %>%
    select(SampleIndex, Group, Date_MMDD, Year, FileName, all_of(compound)) %>%
    rename(Area = all_of(compound))
  
  # Grouping info for axis
  group_info <- plot_df %>%
    group_by(FileName, Date_MMDD, Year) %>%
    summarise(
      xstart = min(SampleIndex),
      xend = max(SampleIndex),
      xmid = mean(SampleIndex),
      .groups = "drop"
    )
  
  year_lines <- group_info %>%
    group_by(Year) %>%
    summarise(
      xstart = min(xstart),
      xend = max(xend),
      xmid = mean(xmid),
      .groups = "drop"
    )
  
  # Y axis structure
  max_y <- max(plot_df$Area, na.rm = TRUE)
  y_line_group <- 0
  y_text_mmdd <- -0.05 * max_y
  y_year_line <- max_y * 1.05
  y_year_text <- max_y * 1.08
  
  # Draw
  p <- ggplot(plot_df, aes(x = SampleIndex, y = Area, color = Group, group = Group)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.4, na.rm = TRUE) +
    geom_segment(data = group_info, aes(x = xstart, xend = xend, y = y_line_group, yend = y_line_group),
                 inherit.aes = FALSE, color = "grey30", linewidth = 0.6) +
    geom_text(data = group_info, aes(x = xmid, y = y_text_mmdd, label = Date_MMDD),
              inherit.aes = FALSE, size = 2.2, angle = 45, hjust = 1, vjust = 1.2, color = "black") +
    geom_segment(data = year_lines, aes(x = xstart, xend = xend, y = y_year_line, yend = y_year_line),
                 inherit.aes = FALSE, color = "black", linewidth = 0.4) +
    geom_text(data = year_lines, aes(x = xmid, y = y_year_text, label = Year),
              inherit.aes = FALSE, size = 3.2, fontface = "bold", color = "black") +
    scale_x_continuous(
      breaks = df$SampleIndex,
      labels = rep("", nrow(df)),
      expand = expansion(mult = c(0.01, 0.07))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = paste0("LLCCV + CCV: ", compound),
      x = "Data Acquisition Date",
      y = "Area"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_text(size = 10, face = "bold", margin = margin(t = 20)),
      panel.grid.minor.x = element_blank(),
      plot.margin = margin(t = 20, r = 10, b = 20, l = 10, unit = "pt")
    )
  
  ggsave(file.path(out_dir, paste0("LLCCV_CCV_", compound, ".pdf")), p, width = 7, height = 4)
}

# ===== Generate plots for all compounds =====
plot_output_dir <- file.path(output_dir, "LLCCV_CCV_MergedPlot")
for (compound in compound_cols) {
  plot_combined_llccv_ccv(combined_df, compound, out_dir = plot_output_dir)
}

# ===== Export Excel: LLCCV and CCV in separate sheets =====
llccv_export <- combined_df %>%
  filter(Group == "LLCCV") %>%
  select(`Sample Name`, Group, FileName, all_of(compound_cols))

ccv_export <- combined_df %>%
  filter(Group == "CCV") %>%
  select(`Sample Name`, Group, FileName, all_of(compound_cols))

write_xlsx(list(
  LLCCV = llccv_export,
  CCV   = ccv_export
), file.path(output_dir, "QC_LLCCV_CCV_AreaMatrix.xlsx"))