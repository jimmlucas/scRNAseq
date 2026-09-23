# Export the exact installed environment, including packages installed outside Conda.
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1]] else "/opt/environment-manifest.yml"

command_lines <- function(command, arguments) {
  result <- suppressWarnings(system2(command, arguments, stdout = TRUE, stderr = FALSE))
  if (!is.null(attr(result, "status")) && attr(result, "status") != 0L)
    stop(command, " failed with status ", attr(result, "status"))
  result
}

conda <- jsonlite::fromJSON(
  paste(command_lines("conda", c("list", "--json", "-n", "scrnaseq")), collapse = "\n"),
  simplifyVector = FALSE
)
conda <- lapply(conda, function(p) {
  p[c("name", "version", "build_string", "channel")]
})
explicit <- command_lines("conda", c("list", "--explicit", "-n", "scrnaseq"))
explicit <- explicit[grepl("^https?://", explicit)]

fields <- c("Repository", "RemoteSha", "RemoteRepo", "RemoteUsername", "RemoteRef")
r <- installed.packages(fields = fields, noCache = TRUE)
r_packages <- lapply(seq_len(nrow(r)), function(i) {
  row <- r[i, ]
  info <- list(name = unname(row[["Package"]]), version = unname(row[["Version"]]),
               built = unname(row[["Built"]]))
  for (field in fields) {
    value <- unname(row[[field]])
    if (!is.na(value) && nzchar(value)) info[[tolower(field)]] <- value
  }
  info
})

dpkg <- command_lines("dpkg-query", "-W")
system_packages <- lapply(dpkg, function(line) {
  parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
  list(name = parts[[1]], version = parts[[2]])
})

manifest <- list(
  schema = "scrnaseq-environment-v1",
  generated_at_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  platform = R.version$platform,
  r_version = as.character(getRversion()),
  bioconductor_version = as.character(BiocManager::version()),
  conda = list(packages = conda, explicit_package_urls = as.list(explicit)),
  r_packages = r_packages,
  system_packages = system_packages
)
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
yaml::write_yaml(manifest, out)
cat("Wrote ", out, "\n", sep = "")
