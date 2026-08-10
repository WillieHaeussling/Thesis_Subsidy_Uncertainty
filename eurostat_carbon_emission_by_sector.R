# 1. Install required packages if missing
library(eurostat)
library(dplyr)

a4_width <- 11.69   # Landscape orientation
a4_height <- 8.27
resolution <- 300

# 2. Download the Air Emissions dataset
cat("Downloading data from Eurostat... This might take a moment.\n")
df <- get_eurostat("env_ac_ainah_r2", time_format = "num")

# 3. Filter for CO2, the EU-27, and Tonnes
df_filtered <- df %>%
  filter(
    geo == "EU27_2020",      # European Union - 27 countries (from 2020)
    airpol == "CO2",         # Carbon dioxide 
    unit == "T"              # Measured in Tonnes
  )

# 4. Filter for top-level NACE sectors to avoid double counting
df_sectors <- df_filtered %>%
  filter(nchar(nace_r2) == 1)

# 5. Convert Eurostat codes to readable labels
df_labelled <- label_eurostat(df_sectors)

# 6. Group sectors into Top 5 and 'Other'
latest_year <- max(df_labelled$TIME_PERIOD, na.rm = TRUE)
top_sectors <- df_labelled %>%
  filter(TIME_PERIOD == latest_year) %>%
  slice_max(values, n = 5) %>%
  pull(nace_r2)

df_grouped <- df_labelled %>%
  mutate(
    sector_group = ifelse(nace_r2 %in% top_sectors, as.character(nace_r2), "Other")
  ) %>%
  group_by(TIME_PERIOD, sector_group) %>%
  summarise(values = sum(values, na.rm = TRUE), .groups = "drop")

# Calculate the percentage breakdown per year
df_pct <- df_grouped %>%
  group_by(TIME_PERIOD) %>%
  mutate(percentage = values / sum(values) * 100) %>%
  ungroup()

# 7. Plot the trends
png("eu_carbon_emissions_by_sector.png", width = a4_width, height = a4_height, units = "in", res = resolution)

# Set up a layout: 2 columns in the top row (plots), 1 merged column in the bottom row (legend)
# heights = c(5, 1.5) means the plot row gets about 77% of the height, legend gets 23%
layout_matrix <- matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE)
layout(layout_matrix, heights = c(5, 1.5))

# Outer margins to leave room for the source text at the very bottom edge
par(oma = c(2, 1, 2, 1))

plot_sectors <- c(top_sectors, "Other")
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7", "#999999")

# ---------------------------------------------------------
# PLOT 1: Absolute Emissions (Million Tonnes)
# ---------------------------------------------------------
par(mar = c(4, 4, 3, 2)) # Standard margins for the plots
x_limits <- range(df_grouped$TIME_PERIOD, na.rm = TRUE)
y_limits <- range(df_grouped$values / 1e6, na.rm = TRUE)

plot(
  x_limits, 
  y_limits, 
  type = "n",
  xlab = "Year", 
  ylab = "CO2 Emissions (Million Tonnes)",
  main = "EU CO2 Emitting Sectors (Top 5 + Other)",
  las = 1,                 
  bty = "l"                
)

grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dashed")

for (i in seq_along(plot_sectors)) {
  sector_data <- subset(df_grouped, sector_group == plot_sectors[i])
  sector_data <- sector_data[order(sector_data$TIME_PERIOD), ]
  
  lines(
    sector_data$TIME_PERIOD, 
    sector_data$values / 1e6,
    col = cb_palette[i], 
    lwd = 2
  )
}

# ---------------------------------------------------------
# PLOT 2: Percentage Breakdown (%)
# ---------------------------------------------------------
par(mar = c(4, 4, 3, 2))
y_limits_pct <- c(0, max(df_pct$percentage, na.rm = TRUE)) 

plot(
  x_limits, 
  y_limits_pct, 
  type = "n",
  xlab = "Year", 
  ylab = "Share of Total CO2 Emissions (%)",
  main = "Sector Breakdown across EU (%)",
  las = 1,                 
  bty = "l"                
)

grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dashed")

for (i in seq_along(plot_sectors)) {
  sector_data_pct <- subset(df_pct, sector_group == plot_sectors[i])
  sector_data_pct <- sector_data_pct[order(sector_data_pct$TIME_PERIOD), ]
  
  lines(
    sector_data_pct$TIME_PERIOD, 
    sector_data_pct$percentage,
    col = cb_palette[i], 
    lwd = 2
  )
}

# ---------------------------------------------------------
# PANEL 3: Merged Legend at the Bottom
# ---------------------------------------------------------
par(mar = c(0, 0, 0, 0)) # Remove margins completely so legend is centered in the space
plot.new()               # Create an empty plot to hold the legend

legend(
  "center",              # Center it in the new bottom panel
  legend = plot_sectors, 
  col = cb_palette, 
  lwd = 2, 
  ncol = 2,              # Split into 2 columns so long labels fit
  bty = "n",             # Remove the box around the legend
  cex = 0.9              # Slightly smaller text
)

# ---------------------------------------------------------
# ADD SOURCE TEXT
# ---------------------------------------------------------
# outer = TRUE places it relative to the whole canvas (oma), below the legend panel
mtext("Source: Eurostat (env_ac_ainah_r2)", side = 1, line = 0, outer = TRUE, adj = 0.98, cex = 0.8)

dev.off()