#!/usr/bin/env Rscript

# Render the GitBook, then inline its search dependency for offline use.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Run with Rscript")
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."))
setwd(root)

bookdown::render_book(
  input = ".", output_format = "bookdown::gitbook", output_dir = "docs"
)

# Bookdown 0.46 references Fuse.js from a CDN even with self_contained = TRUE.
# Reuse the exact 6.4.6 copy embedded in the previous edition of this book.
dependency <- paste0(
  '<script src="https://cdn.jsdelivr.net/npm/',
  'fuse.js@6.4.6/dist/fuse.min.js"></script>'
)
fuse <- paste(readLines("assets/fuse-6.4.6.min.js", warn = FALSE),
              collapse = "\n")
embedded <- paste0("<script>\n", fuse, "\n</script>")
pages <- list.files("docs", pattern = "[.]html$", full.names = TRUE)
replaced <- 0L
for (page in pages) {
  html <- paste(readLines(page, warn = FALSE), collapse = "\n")
  if (!grepl(dependency, html, fixed = TRUE)) next
  occurrences <- lengths(regmatches(html, gregexpr(dependency, html,
                                                    fixed = TRUE)))
  html <- paste(strsplit(html, dependency, fixed = TRUE)[[1L]],
                collapse = embedded)
  writeLines(html, page, useBytes = TRUE)
  replaced <- replaced + occurrences
}
cat("Embedded Fuse.js in ", replaced, " script tags across ",
    length(pages), " HTML pages.\n", sep = "")
