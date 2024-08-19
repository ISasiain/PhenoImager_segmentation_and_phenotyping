#!/usr/bin/Rscript

#
# DEFINING FUNCTIONS FOR PHENOTYPE REFINEMENT
#

# Define function to winsorize data
winsorize <- function(x, lower_percentile = 0, upper_percentile = 0.96) {

  # Ensure valid percentiles
  if (lower_percentile < 0 || lower_percentile >= 1 || upper_percentile <= lower_percentile || upper_percentile > 1) {
    stop("Percentiles must be between 0 and 1, with lower_percentile < upper_percentile")
  }

  # Compute the lower and upper bounds for winsorization
  lower_bound <- quantile(x, lower_percentile, na.rm = TRUE)
  upper_bound <- quantile(x, upper_percentile, na.rm = TRUE)

  # Apply winsorization
  x_winsorized <- pmin(pmax(x, lower_bound), upper_bound)

  # Return winsorized values
  return(x_winsorized)
}


# Define function to assign cells to phenotype centroids based on weighted distance
nearest_centroid_assigner <- function(data_point, centroid_matrix, weights_matrix) {

  # Calculate distance between cell and centroids
  distances <- sapply(rownames(centroid_matrix), function(phenotype) {

    sqrt(sum(weights_matrix[phenotype, markers] * ((as.numeric(data_point) - as.numeric(centroid_matrix[phenotype,]))^2)))

  })

  # Determine and return nearest's centroid ID (phenotype name). If the distance is exactly the same return the first one
  nearest <- names(which.min(distances))[1]
  return(nearest)
}



# Defining fitness evaluation function for the genetic algorithm
evaluation_function <- function(X, centroids_matrix, weights, alpha=2, beta=0.5) {

  # DEFINING INTERNAL FUNCTIONS

  # Function to calculate inter-cluster distance. Distance between cluster means
  inter_cluster_distance <- function(X, cell_phenotypes) {
    # Getting all the identified cell phenotypes
    total_phenotypes <- as.character(unique(cell_phenotypes))
    cell_phenotypes <- as.character(cell_phenotypes)

    # Initialize variable to store inter-cluster distance
    distance <- 0

    for (i in 1:(length(total_phenotypes) - 1)) {
      for (j in (i + 1):length(total_phenotypes)) {
        distance <- distance + sqrt(sum((colMeans(X[cell_phenotypes == total_phenotypes[i],]) - colMeans(X[cell_phenotypes == total_phenotypes[j],]))^2))
      }
    }
    return(distance)
  }

  # Function to calculate intra-cluster variability. Sum of squares
  intra_cluster_variability <- function(X, cell_phenotypes) {
    unique_phenotypes <- unique(cell_phenotypes)
    total_wcss <- 0

    for (phenotype in unique_phenotypes) {
      cluster_data <- X[cell_phenotypes == phenotype, ]
      cluster_center <- colMeans(cluster_data)
      wcss <- sum(rowSums((cluster_data - cluster_center)^2))
      total_wcss <- total_wcss + wcss
    }

    return(total_wcss)
  }

  # TRANSFORMING WEIGHTS INTO A MATRIX

  # Building the matrix
  weights_matrix <- matrix(NA, ncol=length(markers), nrow=length(phenotypes), dimnames = list(phenotypes, markers))

  for (element in names(weights)) {

    # Getting row and colnames
    row <- gsub(";.*", "", element)
    col <-gsub(".*;", "", element)

    # Adding value to matrix
    weights_matrix[row, col] <- weights[element]

  }

  # DEFINE CELL PHENOTYPES FOR SELECTED WEIGHTS
  cell_phenotypes <- sapply(1:nrow(X), function(row) {
    nearest_centroid_assigner(
      data_point = to_cluster[row, markers],
      centroid_matrix = centroids_matrix,
      weights_matrix = weights_matrix
    )
  })

  # DEFINE CONSTRAINTS
  unique_phenotypes <- unique(cell_phenotypes)
  sum_of_weights_per_phenotype <- sapply(unique_phenotypes, function(phenotype) {
    sum(weights_matrix[phenotype, ])
  })

  # Calculate the total number of phenotypes
  num_phenotypes <- length(unique_phenotypes)


  # Determining the two elements of the fitness function
  inter_cluster_distance_value <- inter_cluster_distance(X = X, cell_phenotypes = cell_phenotypes)
  intra_cluster_variability_value <- intra_cluster_variability(X = X, cell_phenotypes = cell_phenotypes)

  # Determining fitness value
  fitness_value <- (alpha * intra_cluster_variability_value) / (beta * inter_cluster_distance_value)

  # Add penalisation
  fitness_value <- fitness_value + abs((sum(sum_of_weights_per_phenotype) - 5 * num_phenotypes)) * 10000000


  # Returning fitness value
  return(fitness_value)

}


# Define function to perform crossover in the genetic algorithm
crossover_function <- function(parent_chr_1, parent_chr_2, blocks = 7, block_size = 5, crossover_rate = 0.8) {

  # Generating boolean to set if the crossover is performed or not
  crossover <- sample(c(TRUE, FALSE), size = 1, prob = c(crossover_rate, 1-crossover_rate))

  # Performing chromosome crossover if crossover == TRUE
  if (isTRUE(crossover)) {

    # Determine cross_over points
    crossover_site <- sample(1:(blocks-1), 1)

    # Generate childs
    child_chr_1 <- c(parent_chr_1[1:(crossover_site * block_size)], parent_chr_2[(crossover_site * block_size + 1):(blocks * block_size)])
    child_chr_2 <- c(parent_chr_2[1:(crossover_site * block_size)], parent_chr_1[(crossover_site * block_size + 1):(blocks * block_size)])

    # Return output
    return(list(child_chr_1, child_chr_2))

  } else {

    # Return output
    return(list(parent_chr_1, parent_chr_2))

  }
}


# Generate function to perform mutation in the genetic algorithm
mutation_function <- function(chromosome,
                              blocks = 7,
                              block_size = 5,
                              mutation_rate = 0.15,
                              mutation_size = 0.075,
                              min = 0,
                              max = 5) {

  for (block in 1:blocks) {

    # Generating boolean to set if the mutation happens or not
    mutation <- sample(c(TRUE, FALSE), size = 1, prob = c(mutation_rate, 1-mutation_rate))

    if (isTRUE(mutation)) {

      # Determining the position of the main mutation in the block
      position <- sample(1:5, size = 1)

      # Determine the sign of the mutation
      sign <- sample(c("-", "+"), size = 1)

      # Mutating the gene. Compensate the weight change with the remaining ones
      if (sign == "-") {

        # Substract 0,05 to the position
        chromosome[(block-1)*block_size + position] <- chromosome[(block-1)*block_size + position] - mutation_size

        # Compensate the values of the remaining positions
        random_numbers <- runif(block_size -1)
        weight_change <- (random_numbers / sum(random_numbers)) * mutation_size

        # Determining indices of positions to update
        indices <- seq((block-1)*block_size + 1, block*block_size)
        indices <- indices[which(indices != (block-1)*block_size + position)]

        # Updating values
        chromosome[indices] <- chromosome[indices] + weight_change

      } else {

        # Substract 0,05 to the position
        chromosome[(block-1)*block_size + position] <- chromosome[(block-1)*block_size + position] + mutation_size

        # Compensate the values of the remaining positions
        random_numbers <- runif(block_size -1)
        weight_change <- - (random_numbers / sum(random_numbers)) * mutation_size

        # Determining indices of positions to update
        indices <- seq((block-1)*block_size + 1, block*block_size)
        indices <- indices[which(indices != (block-1)*block_size + position)]

        # Updating values
        chromosome[indices] <- chromosome[indices] + weight_change

      }

    }

  }

  # Return chromosome
  return(chromosome)

}


# Define the function to add an element and shift the vector
shift_and_add <- function(vec, new_element) {
  # Concatenate the new element to the end of the vector
  vec <- c(vec, new_element)
  # Remove the first element by taking a subset from the second element to the end
  vec <- vec[-1]
  return(vec)
}
