#!/usr/bin/Rscript

#
# LOADING REQUIRED PACKAGES
#

# Install and load mixtools
if (!requireNamespace("mixtools", quietly = TRUE)) {
  install.packages("mixtools")
}
library(mixtools)

# Install and load miscTools
if (!requireNamespace("miscTools", quietly = TRUE)) {
  install.packages("miscTools")
}
library(miscTools)

# Install and load miscTools
if (!requireNamespace("bestNormalize", quietly = TRUE)) {
  install.packages("bestNormalize")
}
library(bestNormalize)

# Install and load miscTools
if (!requireNamespace("stringr", quietly = TRUE)) {
  install.packages("stringr")
}
library(stringr)


# Install and load Wakefield
if (!requireNamespace("MCMCpack", quietly = TRUE)) {
  install.packages("MCMCpack")
}
library(MCMCpack)

# Install and load miscTools
if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}
library(optparse)


#
# SOURCING USER DEFINED FUNCTIONS
#

source("scripts/utils.R")

#
# CONFIGURE AND READING COMMAND LINE ARGUMENTS
#

# Defining command line argument list
option_list <- list(
  make_option(c("-i", "--input_file"), type="character", help="Path to the input file whose phenotyping is intended to be refined"),
  make_option(c("-o", "--output_path"), type="character", default=".", help="Path to generate output file"),
  make_option(c("-d", "--dapi_threshold"), type="double", default=1.9, help="Threshold of DAPI intensity to be used to identify folds in the image"),
  make_option(c("-s", "--size_threshold"), type="double", default=600.0, help="Threshold of nuxlear size to remove segmentation errors"),
  make_option(c("-t", "--tumour_marker"), type="character", default="PAN.CK", help="Name of the tumour marker"),
  make_option(c("-p", "--tumour_phenotype"), type="character", default="PAN-CK", help="Name of the tumour phenotype generated in the threshold based clustering"),
  make_option(c("-e", "--expected_phenotypes_and_markers"), type="character", default = "/Users/isasiain/PhD/Projects/immune_spatial/PhenoImager_automatic_phenotyping/data/expected_phenotypes_and_markers.tsv", help="Data frame with the combination of markers in the expected phenotypes. Phenotypes as rows and markers as columns.")
)


# Parse command line arguments
args <- parse_args(OptionParser(option_list=option_list))

# Reading input files
expected_markers <- read.table(args$"expected_phenotypes_and_markers")
intensities_file <- read.csv(args$"input_file")

#
# CELL FILTER
#

# Filtering cells based on dapi
dapi_filter <- scale(intensities_file$"DAPI") < args$"dapi_threshold"
# Filtering cells based on nuclear size
size_filter <- intensities_file$"Nuclear_area" < args$"size_threshold"


# Removing cells from data frame
intensities_file <- intensities_file[dapi_filter & size_filter, ]


#
# REFINING PHENOTYPING: TUMOUR VS STROMA
#

# Getting_data
to_cluster <- intensities_file[, c("Nuclear_area", args$"tumour_marker", "Phenotype")]

# Perform winsorization and archsinh(x/5) transformation of pan-ck intensities
to_cluster[, args$"tumour_marker"] <- arcsinh_x(winsorize(to_cluster[, args$"tumour_marker"])/5)


# Run Gaussian mixture modelling for two groups setting initial centroid parameters based on thresholding based phenotyping
my_model <- mvnormalmixEM(x = to_cluster,
                          lambda = c(nrow(to_cluster[to_cluster$Phenotype == "PAN-CK",]) / nrow(to_cluster),
                                     nrow(to_cluster[to_cluster$Phenotype != "PAN-CK",]) / nrow(to_cluster)),
                          sigma = list(cov(to_cluster[to_cluster$Phenotype == "PAN-CK", c("Nuclear_area", args$"tumour_marker")]),
                                       cov(to_cluster[to_cluster$Phenotype != "PAN-CK", c("Nuclear_area", args$"tumour_marker")])),
                          mu = list(colMeans(to_cluster[to_cluster$Phenotype == "PAN-CK", c("Nuclear_area", args$"tumour_marker")]),
                                    colMeans(to_cluster[to_cluster$Phenotype != "PAN-CK", c("Nuclear_area", args$"tumour_marker")]))
)

# Getting clusters based on posterior probabilities
cluster_assignments <- apply(my_model$posterior, 1, which.max)

# Assign "Tumour" or "Stroma" categories to clusters based on mean panCK expression
cl_1_panck <- mean(intensities_file[cluster_assignments == 1, args$"tumour_marker"])
cl_2_panck <- mean(intensities_file[cluster_assignments == 2, args$"tumour_marker"])

if (cl_1_panck > cl_2_panck) {
  cluster_assignments[cluster_assignments == 1] <- "Tumour"
  cluster_assignments[cluster_assignments == 2] <- "Stroma"
} else {
  cluster_assignments[cluster_assignments == 1] <- "Stroma"
  cluster_assignments[cluster_assignments == 2] <- "Tumour"
}

# Reassigning possibly misassigned stromal cells
sum_of_stroma = intensities_file[cluster_assignments == "Tumour", "CD8"] +
  intensities_file[cluster_assignments == "Tumour", "CD4"] +
  intensities_file[cluster_assignments == "Tumour", "CD20"] +
  intensities_file[cluster_assignments == "Tumour", "CD68"] +
  intensities_file[cluster_assignments == "Tumour", "FOXP3"]

# Determining a maximum level of strimal marker expression
IQR <- IQR(sum_of_stroma)
Q3 <- quantile(sum_of_stroma, 0.75)
outlier_limit <- Q3 + 1.5 * IQR


# Reassign the stromal cells classified as tumour
cluster_assignments[cluster_assignments == "Tumour"][sum_of_stroma > outlier_limit] <- "Stroma"


#
# REFINING PHENOTYPING: STROMAL SUBTYPES
#

# Determine stromal markers
markers <- colnames(expected_markers)
markers <- markers[markers != args$tumour_marker]

# Getting relevant data
to_cluster <- intensities_file[cluster_assignments == "Stroma", c(markers, "Phenotype")]


# Determining mean residual intensity for each markers based on the thresholding-based phenotyping
res_intensities <- sapply(markers,
                          function(marker) {

                            # Determine the phenotypes to take into account to calculate resiudal intensity
                            phenotypes_to_dismiss <- rownames(expected_markers)[expected_markers[,marker]]

                            # Calculating residual intensity
                            residual_intensity <- mean(to_cluster[!(to_cluster$"Phenotype" %in% phenotypes_to_dismiss), marker])
                            return(residual_intensity)

                          })

# Substract the residual intensity to the intensity of each cell
for (marker in markers) {

  to_cluster[, marker] = to_cluster[, marker] - res_intensities[marker]

}

# Arcsinh(x/5) transform the data
to_cluster[, markers] <- sapply(to_cluster[, markers], function(x) arcsinh_x(as.numeric(x/5))$x.t)


# Generating centroid for each identified phenotype (except PAN-CK)
# Get unique phenotypes
unique_phenotypes <- unique(to_cluster$Phenotype)

# THIS IS HARDCODED. CHANGE IT!!!
unique_phenotypes <- unique_phenotypes[unique_phenotypes != "PAN-CK"]

unique_phenotypes <- unique_phenotypes[which(unique_phenotypes != args$tumour_phenotype)]

# Initialize an empty matrix to store centroids
num_markers <- ncol(to_cluster) - 1  # Number of markers (excluding the Phenotype column)
centroids_matrix <- matrix(NA, nrow = length(unique_phenotypes), ncol = num_markers)
rownames(centroids_matrix) <- unique_phenotypes
colnames(centroids_matrix) <- colnames(to_cluster)[which(colnames(to_cluster) != "Phenotype")]  # Marker names

# Loop through each phenotype and calculate the centroid
for (phenotype in unique_phenotypes) {
  # Subset the data for the current phenotype
  subset_data <- to_cluster[to_cluster$Phenotype == phenotype, ]

  # Calculate the median of each column (excluding the Phenotype column)
  # Convert the subset_data to a matrix for consistency
  subset_matrix <- as.matrix(subset_data[, -ncol(subset_data)])

  # Calculate medians for each marker (column)
  centroids <- apply(subset_matrix, 2, median)

  # Store the result in the matrix
  centroids_matrix[phenotype, ] <- centroids
}

# Initialize empty vector to store best fitness value
fitness_eval <- seq(1:20)

# GENETIC ALGORITHM I: GENERATE INITIAL CHROMOSOMES

# Define the weights assignment function
assign_weights <- function(row) {
  num_true <- sum(row)

  if (num_true == 1) {
    weights_set1 <- ifelse(row, 1.4, 0.9)
    weights_set2 <- ifelse(row, 3, 0.5)
  } else if (num_true == 2) {
    weights_set1 <- ifelse(row, 1.225, 0.85)
    weights_set2 <- ifelse(row, 1.825, 0.45)
  } else {
    weights_set1 <- rep(1, length(row))
    weights_set2 <- rep(1, length(row))
  }

  return(list(weights_set1, weights_set2))
}

# Define the phenotypes
phenotypes <- rownames(expected_markers)[rownames(expected_markers) != "Tumour"]

# Generate all possible combinations (64 combinations for 7 phenotypes)
total_combinations <- 2^length(phenotypes)
combinations <- expand.grid(rep(list(1:2), length(phenotypes)))

# Filter to ensure we only have valid combinations
valid_combinations <- combinations[1:64, ]

# Initialize matrix to store initial values
matrix_of_weights <- matrix(nrow = nrow(valid_combinations), ncol = length(markers) * length(phenotypes))
colnames(matrix_of_weights) <- unlist(sapply(phenotypes, function(p) paste(p, markers, sep = ";")))

# Iterate through all valid combinations
for (i in 1:nrow(valid_combinations)) {
  # Initialize empty vector of weights for the chromosome
  vector_of_weights <- c()

  # Iterate through each phenotype in the combination
  for (j in 1:length(phenotypes)) {
    phenotype <- phenotypes[j]
    weight_set <- valid_combinations[i, j]

    # Get the weights based on the chosen set
    weights_list <- assign_weights(expected_markers[phenotype, markers])
    chosen_weights <- weights_list[[weight_set]]

    # Append the weights to the vector
    vector_of_weights <- c(vector_of_weights, chosen_weights)
  }

  # Add the vector of weights to the initial_values matrix
  matrix_of_weights[i, ] <- vector_of_weights
}

# Defining number of chromosomes
number_of_chromosomes <- nrow(matrix_of_weights)

# Determining the number of parent chromosomes to keep (15% of total). Even number
number_of_parents <- round(number_of_chromosomes * 0.15)

if (number_of_parents %% 2 != 0) {

  number_of_parents = number_of_parents + 1

}

# Determining the number of child chromosomes needed
number_of_childs <- number_of_chromosomes - number_of_parents


# GENETIC ALGORITHM II: STARTING THE OPTIMIZATION LOOP

stromal_markers <- markers[markers != args$tumour_marker]

# Using boolean to indicate when to activate or stop the genetic algorithm
activated_ga <- TRUE
counter = 1

while(isTRUE(activated_ga)) {

  cat("\n\nGenetic algorithm iteration: ", counter)
  counter <- counter + 1

  # GENETIC ALGORITHM III: DETERMINING FITNESS
  fitness <- apply(X = matrix_of_weights,
                   MARGIN = 1,
                   FUN = function (row) evaluation_function (
                     X = to_cluster[, markers],
                     centroids_matrix = centroids_matrix,
                     alpha=2,
                     beta=0.5,
                     weights = row
                   ))

  cat("\nMean fitness: ", mean(fitness))
  cat("\nBest fitness: ", min(fitness))
  cat("\n===================================\n")


  # GENETIC ALGORITHM VI: EVALUATION

  # Ordering chromosomes based on fitness
  fitness_order <- order(fitness)

  fitness <- fitness[fitness_order]
  matrix_of_weights <- matrix_of_weights[fitness_order, ]

  # Use custom function to add an element to a vector and shift the other elements to the left
  fitness_eval <- shift_and_add(fitness_eval, min(fitness))

  # Breaking the loop if there is no improvement in the best solution
  if (length(unique(round(fitness_eval, digits = 2))) == 1) {

    activated_ga = FALSE
    break

  }


  # GENETIC ALGORITHM V: PERFORMING SELECTION

  # Applying selection to the chromosomes with the worst fitness
  fitness[(number_of_chromosomes - number_of_childs + 1):number_of_chromosomes] <- NaN
  matrix_of_weights[(number_of_chromosomes - number_of_childs + 1):number_of_chromosomes,] <- NaN

  # GENETIC ALGORITHM VI: PERFORMING CROSSOVER

  # Perform crossover using a loop
  for (element in 1:(number_of_childs/2)) {

    # Select randomly parent chromosomes
    selected_chromosomes <- sample(1:number_of_parents, 2)

    # Generating childs
    my_child_chrs <- crossover_function(matrix_of_weights[selected_chromosomes[1],], matrix_of_weights[selected_chromosomes[2],])

    # Append child chromosomes to matrix_of_weights
    matrix_of_weights[number_of_parents + element*2-1,] <- my_child_chrs[[1]]
    matrix_of_weights[number_of_parents + element*2,] <- my_child_chrs[[2]]

  }


  # GENETIC ALGORITHM VII: MUTATION

  for (chr in 1:number_of_chromosomes) {

    matrix_of_weights[chr, ] <- mutation_function(matrix_of_weights[chr, ])

  }

}

# GENETIC ALGORITHM VIII: GETTING OPTIMAL SOLUTION

# Getting optimal solution from matrix
optimal_weights <- matrix_of_weights[1, ]

# Transform vector into weights matrix
weights_matrix <- matrix(NA, ncol=length(markers), nrow=length(phenotypes), dimnames = list(phenotypes, markers))

for (element in names(optimal_weights)) {

  # Getting row and colnames
  row <- gsub(";.*", "", element)
  col <-gsub(".*;", "", element)

  # Adding value to matrix
  weights_matrix[row, col] <- optimal_weights[element]

}

#
# GENERATING OUTPUT
#

# Generating phenotypes using the optimized weights
final_cell_assignment <- sapply(1:nrow(to_cluster), function(row) {
                                nearest_centroid_assigner(
                                    data_point = to_cluster[row, stromal_markers],
                                    centroid_matrix = centroids_matrix,
                                    weights_matrix = weights_matrix)
  }
)

# Updating refined phenotype of stromal cells
cluster_assignments[cluster_assignments == "Stroma"] <- final_cell_assignment

# Adding refined clusters to intensity dataframe
intensities_file$"Refined_phenotype" <- cluster_assignments

match <- str_extract(args$"input_file", "Core\\[\\d+,\\d+,[A-Za-z]+\\]")

# Save updated intensity dataframe as output
write.csv(intensities_file, paste0(args$"output_path", "/", match, "_refined.csv"))


# ggplot(data=intensities_file, aes(x=X_coor, y=Y_coor)) +
#   geom_point(aes(colour = as.factor(Refined_phenotype)), cex=0.7)  +
#   labs(x="X", y="Y", colour="Phenotype") +
#   theme_classic()
#
# ggplot(data=intensities_file, aes(x=X_coor, y=Y_coor)) +
#   geom_point(aes(colour = as.factor(Phenotype)), cex=0.7)  +
#   labs(x="X", y="Y", colour="Phenotype") +
#   theme_classic()
#
# ggplot(data=intensities_file[intensities_file$Phenotype == "CD4",], aes(x=X_coor, y=Y_coor)) +
#   geom_point(aes(colour = as.factor(Refined_phenotype)), cex=0.7)
#
# ggplot(data=intensities_file[intensities_file$Refined_phenotype == "CD8_FOXP3",], aes(x=X_coor, y=Y_coor)) +
#   geom_point(aes(colour = as.factor(Refined_phenotype)), cex=0.7)


