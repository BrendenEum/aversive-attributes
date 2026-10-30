
.color_e = 'lightcyan1'
.color_c = 'mistyrose1'
.color_j = 'lightcyan2'

# Figure size. Every figure has two columns (Study 1 | Study 2) and a column can
# hold up to three facets (N, B, P), so .figw is the width of ONE column and
# .figh is the height of ONE row. A figure with r rows is saved at
#   width = 2 * .figw, height = r * .figh
# so every figure has the same width and every row has the same height.
.figw = 9
.figh = 6.5

# Text size. The figures are drawn on a big canvas (2 * .figw inches wide) and
# then shrunk to .print_width inches in the paper (include them at \textwidth;
# letter paper with 1 in margins gives 6.5 in). To keep the text the same size
# on the page in every figure, ggplot has to draw it .fontscale times bigger
# than it will print. Set the sizes below in PRINTED points and leave the rest.
.print_width = 6.5
.fontscale   = (2 * .figw) / .print_width
.tick_pt     = 8     # axis tick labels and legend labels
.title_pt    = 10    # axis titles, legend title, Study titles and Study strips
.strip_pt    = 7     # attribute-type facet labels (Aversive (N), ...), not bold

.myPlot = list(
  theme_bw(),
  coord_cartesian(expand=T, default=T),
  theme(
    legend.position="None",
    legend.background=element_blank(),
    legend.key = element_rect(fill = NA),
    legend.key.size = unit(1.5 * .tick_pt * .fontscale, "pt"),
    #legend.spacing.x = unit(0.1, 'cm'),
    #legend.spacing.y = unit(0.1, 'cm'),
    plot.margin = unit(c(.5,.5,.5,.5), "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = .title_pt * .fontscale),
    axis.title = element_text(size = .title_pt * .fontscale),
    axis.text = element_text(size = .tick_pt * .fontscale),
    legend.title = element_text(size = .title_pt * .fontscale),
    legend.text = element_text(size = .tick_pt * .fontscale)
  )
)


.linewidth = 1.5
.markersize = .01
.errsize = 1.5
.ribbonalpha = 0.4

# Named scales shared across notebooks (Rotman palette).
# Pass directly to scale_*_manual(values = ...).
.type_colors  = c("P" = "#004990", "N" = "#E20778", "B" = "#FCAF17")   # attribute type
.wn_colors    = c("Narrow" = "#09AD87", "Wide" = "#3E77B7")             # range condition
.wn_linetypes = c("Narrow" = "solid",   "Wide" = "dashed")
.wn_shapes    = c("Narrow" = 16,        "Wide" = 17)                    # filled circle / triangle
