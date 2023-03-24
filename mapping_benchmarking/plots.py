from snakehelp.plotting import PlotType, Plot
from .parameters import parameters, result

# Plot types

f1_score_vs_read_length = PlotType(
    type="line",
    y: f1_score
    x: read_length
    color: method
    facet_col: variant_filter
    facet_row: read_type
