library(magick)


combine_frame <- function(i, scenario, scene) {
  # Make header
  # check available fonts with magick_fonts()
  header <- image_blank(1920, 130, "#00718B") %>%
    image_annotate(
      scenario$title,
      font = "C059",
      size = 64,
      color = "white",
      gravity = "center"
    )
  image_write(header, file.path(scene$dir_write, "header.png"))


  # get time for side panel
  # round input in case floating point errors
  time_s <- round(scenario$time$seconds[[i]], 3)
  hours <- floor(time_s / 3600)
  time_s <- time_s - hours * 3600
  minutes <- floor(time_s / 60)
  time_s <- time_s - minutes * 60
  seconds <- round(time_s)

  # make side panels
  space <- ifelse(minutes < 10, "\u2007", "")
  side_panel <- image_blank(620, 960, color = "white") %>%
    image_annotate(
      sprintf("%s%2d h %s%2.1d m", "Time post-earthquake\n", hours, space, minutes),
      font = "C059",
      size = 42,
      color = "black",
      gravity = "center"
    ) %>%
    image_annotate(
      sprintf("\n%s", scene$title),
      font = "C059",
      size = 54,
      color = "black",
      gravity = "north"
    )

  # place logo at bottom of side panel
  logo <- image_read(
    "input/logos/geoscience_inline.png"
  ) %>%
  image_resize("620x100")
  side_panel <- image_composite(side_panel, logo, gravity = "south", offset = "+0+20")

  # merge header, side panel and image
  merge_frames(i, scene, side_panel)
}

merge_frames <- function(i, scene, side) {
  header <- image_read(file.path(scene$dir_write, "header.png"))
  main <- image_read(file.path(scene$dir_write, sprintf("stage%03d.png", i)))
  # resize to ensure 1300x950
  main <- image_resize(main, "1300x950")

  body <- image_append(c(main, side))
  full <- image_append(c(header, body), stack = TRUE)

  image_write(full, file.path(scene$dir_write, sprintf("frame_%03d.png", i)))
}
