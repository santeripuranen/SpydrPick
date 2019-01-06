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
ld_dist <- 
# Outlier thresholds (0 - not drawn).
outlier_threshold <- 
extreme_outlier_threshold <- 
# If chromosome type is circular, set genome_length, otherwise leave as 0.
genome_length <- 
# File type ("binary" or "csv").
file_type <- "binary"
# File paths.
edges_file_path <- ""
mi_file_path <- ""
aracne_file_path <- ""
output_directory_path <- ""
# Colors.
color_direct <- rgb(0, 115, 190, maxColorValue = 255)
color_indirect <- rgb(192, 192, 192, maxColorValue = 255)
# Plot sizes
width <- 1200
height <- 450
pointsize <- 16

# ----------------------------------------------------------------------
# Initialization - Create the necessary functions
# ----------------------------------------------------------------------

# Read <n_edges_read> instances of <data_type> from a file.
read_binary_file <- function(file_path, n_edges_read, data_type) {
    file <- file(file_path, "rb"); res <- readBin(file, data_type, n_edges_read); close(file); return(res)
}

# Read edges from file and calculate edge distances. 
read_binary_edges_and_calculate_distances <- function(file_path, n_edges_read, data_type, genome_length) {
    edges <- read_binary_file(file_path, n_edges_read * 2, data_type)
    edge_distances <- abs(edges[seq(1, n_edges_read, 2)] - edges[seq(2, n_edges_read, 2)])
    edge_distances <- pmin(genome_length - edge_distances, edge_distances)
    return(edge_distances)
}

# ----------------------------------------------------------------------
# File reading
# ----------------------------------------------------------------------

start_reading <- proc.time()
if (file_type == "binary") { # Binary file.
    edge_distances <- read_binary_edges_and_calculate_distances(edges_file_path, n_edges_read, "integer", genome_length)
    aracne <- read_binary_file(aracne_file_path, n_edges_read, "integer") != 0
    mi <- read_binary_file(mi_file_path, n_edges_read, "double")
} else { # Csv file.
    # Implement csv reading logic here.
}
end_reading <- proc.time()
time_reading <- end_reading - start_reading

# Disable scientific notation
options(scipen=999)

start_plotting <- proc.time()
max_distance <- max(edge_distances[1:n_edges])
min_mi <- min(mi[1:n_edges])
max_mi <- max(mi[1:n_edges])
exponent <- round(log10(max_distance)) - 1

indirect_edge_distances <- edge_distances[1:n_edges][!aracne[1:n_edges]]
direct_edge_distances <- edge_distances[1:n_edges][aracne[1:n_edges]]
indirect_mi <- mi[1:n_edges][!aracne[1:n_edges]]
direct_mi <- mi[1:n_edges][aracne[1:n_edges]]

png(paste(output_directory_path, "/plot_", width, "x", height, "_", n_edges, "edges.png", sep=""), 
    width = width, height = height, pointsize = pointsize)
plot(indirect_edge_distances, indirect_mi, col = color_indirect, type = "p", cex = 0.1, 
     xlim = c(0, max_distance), ylim = c(min_mi, 1), xaxs = "i", yaxs = "i",
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
lines(direct_edge_distances, direct_mi, col = color_direct, type = "p", cex = 0.1)
axis(1, at = seq(0, max_distance, 10^exponent), tick = FALSE, labels = seq(0, max_distance / 10^exponent), line = -0.8)
title(xlab = "Distance between positions (bp)", line = 1.2)
title(xlab = substitute(x10^exp, list(exp = exponent)), line = 1.4, adj = 1)
axis(2, at = seq(0.05, 1, 0.05), labels = FALSE, tcl = -0.5)
axis(2, at = seq(0.1, 1, 0.1), labels = seq(0.1, 1, 0.1), las = 1, tcl = -0.5)
title(ylab = "Mutual information", line = 2.5)
segments(0, extreme_outlier_threshold, max_distance, extreme_outlier_threshold, col = "red", lty = 2)
segments(0, outlier_threshold, max_distance, outlier_threshold, col = "red", lty = 2)
segments(ld_dist, min_mi, ld_dist, 1, col = "red", lty = 2)
text(0, extreme_outlier_threshold, "**", col = "red", pos = 2, offset = 0.2, cex = 1, xpd = NA)
text(0, outlier_threshold, "*", col = "red", pos = 2, offset = 0.2, cex = 1, xpd = NA)
legend(x = max_distance, y = max_mi, cex = 1.2, pch = 19, bty = "n", xjust = 1.2, yjust = 1.2,
       c("Indirect", "Direct"), col = c(color_indirect, color_direct))
dev.off()
end_plotting <- proc.time()
time_plotting <- end_plotting - start_plotting

