# Helper for the record_*/dataset-building scripts in this directory.
#
# Those scripts read from lab project directories outside the package. Writing
# the full path out records the project's directory name and the raw file
# names, both of which routinely contain a collaborator's name, and this
# directory is public. `source_path()` globs for the file instead, so a script
# names only the project number and the part of the filename that describes
# what the file is.
#
# The glob must match exactly one path; anything else is an error rather than a
# silent choice of the first match.
#
# Usage:
#   source("data-raw/source_paths.R")
#   psm_inf <- source_path("2026/2026_11_*/raw/*_TMT_PSMs.txt.gz")

projects_root <- "~/git_repos/projects"

source_path <- function(pattern, root = projects_root) {
  hits <- Sys.glob(path.expand(file.path(root, pattern)))

  if (length(hits) != 1) {
    stop(sprintf(
      "Expected exactly 1 path matching '%s', found %i.",
      pattern, length(hits)
    ))
  }

  hits
}
