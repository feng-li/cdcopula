rtags(getwd(),
      pattern = "[.]*\\.[Rr]$",
      keep.re = "/R/",
      verbose = TRUE,
      ofile = "TAGS",
      append = FALSE,
      recursive = TRUE)
