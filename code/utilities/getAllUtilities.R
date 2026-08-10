paths <- yaml::yaml.load_file(here::here("config/paths.yaml"))
setwd(paths$root)

.util_dir       <- file.path(paths$root, paths$code, "utilities")
.file_list      <- list.files(.util_dir, full.names = TRUE)
.files          <- .file_list[endsWith(.file_list, ".R")]
.R_script_files <- .files[!sapply(.files, grepl, pattern = "getAllUtilities")]

for (.script_file in .R_script_files) {
  source(.script_file)
}
