#' Make the introduction slides
#' Uses text content from input scenario
library(magick)
library(purrr)
library(dplyr)

# source("input/kermadec_94.in.R")
source("input/solomon_90.in.R")

fps <- 24
n_slides <- c(9 * fps, 10 * fps, 12 * fps)

logo <- image_read("input/intro_images/geoscience_inline.png")
logo <- image_scale(logo, "x100")
logo <- image_extent(logo, "1260x100", gravity = "east", color = "white")
right_logo_pad <- image_blank(20, 100, "white")
logo <- image_append(c(logo, right_logo_pad))

intro_dir <- file.path(scenario$dir, "intro")

# first slide
wrapped_text <- paste(strwrap(scenario$info_text, width = 42), collapse = "\n")
intro1 <- image_blank(1280, 700, "white") %>%
  image_annotate(
    scenario$intro_title,
    font = "C059",
    size = 64,
    color = "#00718B",
    gravity = "north",
    location = "+0+10"
  ) %>%
  image_annotate(
    wrapped_text,
    font = "C059",
    size = 42,
    color = "black",
    gravity = "west",
    location = "+200+0"
  )
intro1 <- image_append(c(intro1, logo), stack = TRUE)

# second slide
wrapped_text <- paste(strwrap(scenario$info_text_2, width = 44), collapse = "\n")
intro2 <- image_blank(1280, 700, "white") %>%
  image_annotate(
    scenario$title,
    font = "C059",
    size = 64,
    color = "#00718B",
    gravity = "north"
  ) %>%
  image_annotate(
    wrapped_text,
    font = "C059",
    size = 42,
    color = "black",
    gravity = "west",
    location = "+200+0"
  )
intro2 <- image_append(c(intro2, logo), stack = TRUE)

# third slide
wrapped_text <- paste(strwrap(scenario$info_text_3, width = 35), collapse = "\n")
goff_image <- image_read("input/intro_images/goff2022_kermadec.jpg") %>%
  image_scale("500") %>%
  image_extent("500x600", gravity = "center", color = "white")
header <- image_blank(1260, 100, color = "white") %>%
  image_annotate(
    scenario$title,
    font = "C059",
    size = 64,
    color = "#00718B",
    gravity = "north",
    location = "+0+10"
  )
intro3 <- image_blank(780, 600, "white") %>%
  image_annotate(
    wrapped_text,
    font = "C059",
    size = 42,
    color = "black",
    gravity = "west",
    location = "+50+0"
  )
citation <- paste(strwrap("Goff, J. (2022). In search of Holocene trans-Pacific palaeotsunamis. Earth-Science Reviews. https://doi.org/10.1016/j.earscirev.2022.104194", width = 70), collapse = "\n")
logo_with_citation <- logo %>%
  image_annotate(
    citation,
    font = "C059",
    size = 20,
    color = "black",
    gravity = "northwest"
  )

goff_image
main <- image_append(c(intro3, goff_image))
intro3 <- image_append(c(header, main, logo_with_citation), stack = TRUE)

dir.create(intro_dir, showWarnings = FALSE)

intro_id <- seq_len(n_slides[1])
map(intro_id, ~ {
  filename <- file.path(intro_dir, sprintf("intro_%03d.png", .x))
  image_write(intro1, path = filename)
})

intro_id <- seq_len(n_slides[2])
map(intro_id, ~ {
  filename <- file.path(intro_dir, sprintf("intro_%03d.png", .x + n_slides[1]))
  image_write(intro2, path = filename)
})

intro_id <- seq_len(n_slides[3])
map(intro_id, ~ {
  filename <- file.path(intro_dir, sprintf("intro_%03d.png", .x + n_slides[1] + n_slides[2]))
  image_write(intro3, path = filename)
})
