#' Make the introduction slides
#' Uses text content from input scenario
library(magick)
library(purrr)
library(dplyr)
library(qrcode)

fps <- 24
n_slides <- 8 * fps

outro_text <-
  "Main text here."
acknowledgement <-
  "Thanks to the funders."

# Make qr code
# if (!file.exists("input/intro_images/report_qr.png")) {
png(filename = "input/intro_images/report_qr.png", width = 450, height = 500)
qr <- qr_code("https://dx.doi.org/xxx")
plot(qr)
dev.off()
# }

# header
header <- image_blank(1260, 100, color = "white") %>%
  image_annotate(
    "Further Information",
    font = "C059",
    size = 64,
    color = "#00718B",
    gravity = "north",
    location = "+0+10"
  )

# main content
wrapped_text <- paste(strwrap(outro_text, width = 35), collapse = "\n")
content <- image_blank(780, 600, "white") %>%
  image_annotate(
    wrapped_text,
    font = "C059",
    size = 40,
    color = "black",
    gravity = "west",
    location = "+50+0"
  )
qr_code <- image_read("input/intro_images/report_qr.png") %>%
  image_extent("x550", gravity = "south", color = "white") %>%
  image_annotate(
    "Scan to read.",
    font = "C059",
    size = 42,
    gravity = "north",
    location = "+0+50"
  ) %>%
  image_annotate(
    "https://dx.doi.org/xxxx",
    font = "C059",
    size = 24,
    gravity = "south"
    )

# add qr code and footer
main <- image_append(c(content, qr_code)) %>%
  image_annotate(
    acknowledgement,
    font = "C059",
    size = 24,
    color = "black",
    gravity = "south",
    location = "-50+0"
  )

# Process logos
logo <- c(
  image_read("input/intro_images/logo_1.png"),
  image_read("input/intro_images/logo_2.jpg"),
  image_read("input/intro_images/logo_3.png"),
  image_read("input/intro_images/logo_4.png")
)
# convert the qfd logo from CMYK to RGB
logo[2] <- image_convert(logo[2], colorspace = "sRGB")
# make them all 100 pixels high
logo <- image_scale(logo, "x100")
# pad them all on the right side
right_logo_pad <- image_blank(40, 100, "white")
logo[1] <- image_append(c(logo[1], right_logo_pad))
logo[2] <- image_append(c(logo[2], right_logo_pad))
logo[3] <- image_append(c(logo[3], right_logo_pad))

# combine all logos into a footer
footer <- image_append(logo)
footer <- image_extent(footer, "1260x100", gravity = "center", color = "white")

# add footer
full_slide <- image_append(c(header, main, footer), stack = TRUE)

outro_dir <- file.path("outro")
dir.create(outro_dir, showWarnings = FALSE)

intro_id <- seq_len(n_slides[1])
map(intro_id, ~ {
  filename <- file.path(outro_dir, sprintf("outro_%03d.png", .x))
  image_write(full_slide, path = filename)
})
