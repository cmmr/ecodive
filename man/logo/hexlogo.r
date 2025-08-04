library(ggplot2)
library(ggpattern)

ggplot() + 
  geom_polygon_pattern(
    mapping = aes(
      x = 100 * sqrt(3)/2 * c(0, 1, 1, 0, -1, -1), 
      y = 100 * .5 * c(2, 1, -1, -2, -1, 1) ),
    pattern          = 'image',
    pattern_type     = 'expand',
    pattern_filename = 'logo/ecodiver.png',
    color         = '#0C356A',
    linewidth     = 0.5 ) + 
  coord_fixed(ratio = 1) +
  theme_void() +
  theme(rect = element_rect(fill = 'transparent')) +
  annotate(
    geom   = 'text', 
    label  = 'ecodive',
    family = 'Super Ocean',
    alpha  = 0.8,
    angle  = 30,
    size   = 2,
    x      = 25, 
    y      = -40 )

ggsave(
  path     = 'logo',
  filename = 'ecodive.png', 
  device   = 'png',
  width    = 200, 
  height   = 200, 
  units    = 'px',
  bg       = 'transparent' )
