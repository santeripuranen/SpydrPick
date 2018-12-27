# ----------------------------------------------------------------------
# Genome-wide epistasis studies (GWES) Manhattan plot script
# Copyright (c) 2018 Juri Kuronen and Santeri Puranen.
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# 1) Settings - Fill in missing values according to SpydrPick's output.
# ----------------------------------------------------------------------

# Number of edges to read from the files.
n_edges_read <- 
# Number of edges to draw.
n_edges <- 
# Linkage disequilibrium distance (0 - not drawn).
ld_dist <- 0
# Outlier thresholds (0 - not drawn).
outlier_threshold <- 0.0
extreme_outlier_threshold <- 0.0
# File type (0 - binary file, 1 - csv file).
file_type <- 0
# File paths.
edges_file_path <- ""
mi_file_path <- ""
aracne_file_path <- ""
# Colors.
color_direct <- rgb(0, 115, 190, maxColorValue = 255)
color_indirect <- rgb(192, 192, 192, maxColorValue = 255)
# Compression (FALSE - draw all values, TRUE - compresses input, if necessary, to speed-up drawing).
compression <- FALSE
# Plot sizes
width <- 1200
height <- 800

# ----------------------------------------------------------------------
# Initialization - Create the necessary functions
# ----------------------------------------------------------------------

#
read_file <- function(file_path, n_edges_read) {
    file <- file(file_path, "rb"); ret <- readBin(file, "double", n_edges_read); close(file); return(ret);
}

#
calculate_distances <- function(edges, n_edges) {
    distances <- numeric(n_edges)
    for (i in 1:n_edges) {distances[i] <- abs(edges[i * 2 - 1] - edges[i * 2])}
    return(distances)
}

# 
calculate_plot_indices <- function(edges, mi, n_edges, width, height) {
    indices = logical(n_edges)
    # Implement compression code here.
    return(indices)
}


# ----------------------------------------------------------------------
# File reading
# ----------------------------------------------------------------------

start_reading <- proc.time()
if (file_type == 0) { # Binary file.
    # Read edges and then replace them with distances.
    edges <- calculate_distances(read_file(edges_file_path, 2 * n_edges_read), n_edges_read)
    aracne <- read_file(aracne_file_path, n_edges_read) != 0
    mi <- read_file(mi_file_path, n_edges_read)
} else { # Csv file.
    # Implement csv reading logic here.
}
end_reading <- proc.time()
time_reading <- end_reading - start_reading

start_compression <- proc.time()
# Get edge indices to plot.
if (compression == TRUE) {
    indices <- calculate_plot_indices(edges, mi, n_edges, width, height);
} else {
    indices <- numeric(n_edges) == 0
}
end_compression <- proc.time()
time_compression <- end_compression - start_compression

options(scipen=999)

#1/1000 ~1000x1000 images

max_dist <- max(edges[1:n_edges])
min_mi <- min(mi[1:n_edges])
max_mi <- max(mi[1:n_edges])
multiplier <- ceiling(log10(max_dist))

start_plotting <- proc.time()
indirect_edges <- edges[1:n_edges][indices & !aracne[1:n_edges]]
direct_edges <- edges[1:n_edges][indices & aracne[1:n_edges]]
indirect_mi <- mi[1:n_edges][indices & !aracne[1:n_edges]]
direct_mi <- mi[1:n_edges][indices & aracne[1:n_edges]]
png(paste("plots/plot_", width, "x", height, "_", n_edges/ 1e6, "M_comp", compression, ".png", sep=""), width = width, height = height, pointsize = 12)
plot(indirect_edges, indirect_mi, col = color_indirect, type = "p", cex = 0.1, 
     xlim = c(0, max_dist), ylim = c(min_mi, 1),
     xlab = "Distance between positions (bp)", ylab = "Mutual information", 
     xaxt = "n", yaxt = "n")
axis(1, at = seq(0, ceiling(max_dist), 10^(multiplier - 1)), labels = seq(0, max_dist / 10^(multiplier - 1)))
axis(2, at = seq(0.05, 1, 0.05), labels = FALSE)
axis(2, at = seq(0.1, 1, 0.1), las = 1)
lines(direct_edges, direct_mi, col = color_direct, type = "p", cex = 0.1, 
      xlim = c(0, max_dist), ylim = c(min_mi, 1))
segments(0, extreme_outlier_threshold, max_dist, extreme_outlier_threshold, col = "red", lty = 2)
segments(0, outlier_threshold, max_dist, outlier_threshold, col = "red", lty = 2)
segments(ld_dist, min_mi, ld_dist, 1, col = "red", lty = 2)
end_plotting <- proc.time()
time_plotting <- end_plotting - start_plotting

# (DEBUG) Draw run times in the plot.
debug <- 1
if (debug > 0) {
  text(max_dist, max_mi, adj = 1, paste(n_edges / 1e6, "(x 1e6) edges, ", y_compression, "y compression", x_compression, "x compression"), cex=1)
  x <- time_reading[3]
  time_str <- "seconds"
  if (x > 60) {time_str <- "minutes"; x <- x / 60} else if (x > 3600) {time_str <- "hours"; x <- x / 3600}
  text(max_dist, max_mi - 0.05, adj = 1, paste(round(x, 2), time_str, "for file reading"), cex=1)
  x <- time_compression[3]
  time_str <- "seconds"
  if (x > 60) {time_str <- "minutes"; x <- x / 60} else if (x > 3600) {time_str <- "hours"; x <- x / 3600}
  text(max_dist, max_mi - 0.10, adj = 1, paste(round(x, 2), time_str, "for compressing"), cex=1)
  x <- time_plotting[3]
  time_str <- "seconds"
  if (x > 60) {time_str <- "minutes"; x <- x / 60} else if (x > 3600) {time_str <- "hours"; x <- x / 3600}
  text(max_dist, max_mi - 0.15, adj = 1, paste(round(x, 2), time_str, "for plotting"), cex=1)
  x <- time_plotting[3]
  text(max_dist, max_mi - 0.20, adj = 1, paste("Compression", 100 * round(1 - sum(indices) / n_edges, 5), "%"), cex=1)
  dev.off()
}
