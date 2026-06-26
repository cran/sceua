# Note: variables prefixed with `.` are used for text replacement in
# Makevars.in and Makevars.win.in.

# Adapted from
# https://github.com/belian-earth/a5R/blob/main/tools/config.R

source("tools/msrv.R")

env_debug <- Sys.getenv("DEBUG")
env_not_cran <- Sys.getenv("NOT_CRAN")

vendor_exists <- file.exists("src/rust/vendor.tar.xz")
is_not_cran <- nzchar(env_not_cran)
is_debug <- nzchar(env_debug)

if (is_debug) {
  is_not_cran <- TRUE
  message("Creating DEBUG build.")
}

if (!is_not_cran) {
  message("Building for CRAN.")
}

# Use --offline whenever vendored crates are present. This prevents Cargo
# network access during R CMD check and CRAN installation.
.cran_flags <- if (vendor_exists) {
  "-j 2 --offline"
} else {
  ""
}

.profile <- if (is_debug) {
  ""
} else {
  "--release"
}

.clean_targets <- if (is_debug) {
  ""
} else {
  "$(TARGET_DIR)"
}

webr_target <- "wasm32-unknown-emscripten"
is_wasm <- identical(R.version$platform, webr_target)

if (is_wasm) {
  message("Building for WebR")
}

target_libpath <- if (is_wasm) {
  webr_target
} else {
  NULL
}

cfg <- if (is_debug) {
  "debug"
} else {
  "release"
}

.libdir <- paste(c(target_libpath, cfg), collapse = "/")

.target <- if (is_wasm) {
  paste0("--target=", webr_target)
} else {
  ""
}

.panic_exports <- if (is_wasm) {
  "CARGO_PROFILE_DEV_PANIC=\"abort\" CARGO_PROFILE_RELEASE_PANIC=\"abort\" "
} else {
  ""
}

is_windows <- .Platform[["OS.type"]] == "windows"

mv_fp <- if (is_windows) {
  "src/Makevars.win.in"
} else {
  "src/Makevars.in"
}

mv_ofp <- if (is_windows) {
  "src/Makevars.win"
} else {
  "src/Makevars"
}

if (file.exists(mv_ofp)) {
  message("Cleaning previous `", mv_ofp, "`.")
  invisible(file.remove(mv_ofp))
}

mv_txt <- readLines(mv_fp)
new_txt <- gsub("@CRAN_FLAGS@", .cran_flags, mv_txt, fixed = TRUE)
new_txt <- gsub("@PROFILE@", .profile, new_txt, fixed = TRUE)
new_txt <- gsub("@CLEAN_TARGET@", .clean_targets, new_txt, fixed = TRUE)
new_txt <- gsub("@LIBDIR@", .libdir, new_txt, fixed = TRUE)
new_txt <- gsub("@TARGET@", .target, new_txt, fixed = TRUE)
new_txt <- gsub("@PANIC_EXPORTS@", .panic_exports, new_txt, fixed = TRUE)

message("Writing `", mv_ofp, "`.")
con <- file(mv_ofp, open = "wb")
writeLines(new_txt, con, sep = "\n")
close(con)

message("`tools/config.R` has finished.")
