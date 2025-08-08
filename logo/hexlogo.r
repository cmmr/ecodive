library(ggplot2)
library(ggpattern)
library(magick)

ggplot() + 
  geom_polygon_pattern(
    mapping = aes(
      x = 1200 * sqrt(3)/2 * c(0, 1, 1, 0, -1, -1), 
      y = 1200 * .5 * c(2, 1, -1, -2, -1, 1) ),
    pattern          = 'image',
    pattern_type     = 'expand',
    pattern_filename = 'logo/cartoon.png',
    color            = '#0078a5',
    linewidth        = NA ) + 
  coord_fixed(ratio = 1) +
  theme_void() +
  theme(rect = element_rect(fill = 'transparent')) +
  annotate(
    geom   = 'text', 
    label  = 'ecodive',
    family = 'Super Ocean',
    alpha  = 0.8,
    angle  = 30,
    size   = 20,
    x      = 290, 
    y      = -500 )

ggsave(
  path     = 'logo',
  filename = 'ecodive.png', 
  device   = 'png',
  width    = 2000, 
  height   = 2000, 
  dpi      = 380, 
  units    = 'px',
  bg       = 'transparent' )


image_read('logo/ecodive.png') |>
  image_trim() |>
  image_resize('x200') |>
  image_write('man/figures/logo.png')
