################################################################
#### Supplements for Simultaneous Detection of Cheaters and ####
########## Cheating Items Using a Biclustering Approach ########
############# Hyeryung Lee & Walter P. Vispoel (2025) ##########
################################################################

#### Library
library(LNIRT)
library(QUBIC)
library(DescTools)


#### Data: Cizek & Wollack (2016) - Form 1
data <- CredentialForm1
K = 170  # number of items
N = nrow(data)  # number of examinees
true_cheater <- which(data$Flagged==1) # true cheaters
true_ci <- c(1,3,4,7,8,9,10,14,15,21,
             22,25,26,29,30,35,38,40,41,43,
             50,53,54,58,59,61,62,67,68,72,
             74,75,81,84,87,89,96,97,98,106,
             110,116,121,127,128,129,134,135,136,140,
             144,145,151,152,154,159,160,161,162,163,
             164,165,169,170) # true cheating items


#### Input matrix generation 
mat <- matrix(0, nrow = N, ncol = K)
time <- data[, 411:580]
time[time == 0] <- NA

# Standardize time for each item to create a time matrix
item_mean_time <- colMeans(time, na.rm = TRUE) # Compute mean response time for each item (column)
time_mat <- matrix(NA, nrow = N, ncol = K)
for (i in 1:K) {
  time_mat[, i] <- (time[, i] - item_mean_time[i]) / sd(time[, i], na.rm = T)
}

# Compute each examinee's average standardized response time
person_mean_time <- apply(time_mat, 1, mean, na.rm = TRUE)

# Compute time difference for each examinee and item
time_dif_mat <- matrix(NA, nrow = N, ncol = K)
for (p in 1:N) {
  time_dif_mat[p, ] <- (time_mat[p, ] - person_mean_time[p]) / sd(time_mat[p, ], na.rm = T)
}

# Adjust responses in 'mat' based on response time
for (p in 1:N) {
  for (i in 1:K) {
    iraw  <- paste0("iraw.", i)   # Response accuracy
    idur  <- paste0("idur.", i)   # Response time
    iresp <- paste0("iresp.", i)  # Response choice
    
    # Create a temporary data frame for the current item
    comp <- data.frame(
      accuracy = as.numeric(data[[iraw]]),
      time     = as.numeric(data[[idur]]),
      choice   = as.factor(data[[iresp]])
    )
    # Set a threshold: 50% of the median response time for this item
    median2 <- median(comp$time, na.rm = T) * 0.5
    
    if (is.na(time_dif_mat[p, i])) {
      mat[p, i] <- 9 # If time difference is NA, mark as 9 (flag value)
    } else if (comp$time[p] < median2 && time_dif_mat[p, i] < 0 && comp$accuracy[p] == 1) {
      mat[p, i] <- 10 # If response time is fast (below threshold), time is below average, and answer is correct, code 10
    } else if (comp$time[p] < median2 && time_dif_mat[p, i] < 0 && comp$accuracy[p] == 0) {
      mat[p, i] <- as.numeric(as.character(comp$choice[p])) # If response time is fast, time is below average, and answer is incorrect, use the original response choice
    } else if (is.na(comp$time[p]) || is.na(comp$accuracy[p])) {
      mat[p, i] <- 9
    }
  }
}
# Define a customized QUBIC function (excluding discretization)
qubic_no_disc <- function(x, r = 1L, q = 0.06, c = 0.95, o = 100, f = 1,
                          k = max(ncol(x) %/% 20, 2), type = "default", P = F,
                          C = F, verbose = T, weight = NULL, seedbicluster = NULL) {
  x_d <- x
  return(qubiclust_d(x_d, c, o, f, k, type, P, C, verbose, weight, seedbicluster))
}



#### Cheating detection using QUBIC
bicluster_n <- 10 # adjust the desired number of outputted biclusters
res <- qubic_no_disc(mat, c = 1, o= bicluster_n, k = 2, verbose = F)

# Collect biclusters from the QUBIC result
biclusters <- list()
for (h in 1:res@Number) {
  bicluster <- list(
    rows = which(res@RowxNumber[, h] == T),  # Examinees 
    cols = which(res@NumberxCol[h, ] == T)   # Items 
  )
# Filter biclusters based on the proportion(50%) of correct, fast responses
  if (length(bicluster$rows) >= 1 &&
      length(which(colMeans(mat[bicluster$rows, bicluster$cols]) == 10)) > length(bicluster$cols) * 0.5) {
    biclusters[[h]] <- bicluster
  }
}

# Calculate p-values for each bicluster and filter by p < 0.05
p_values <- sapply(biclusters, function(bic) {
  if (is.null(bic) || length(bic$rows) == 0 || length(bic$cols) == 0) return(NA)
  ori_pattern <- apply(mat[bic$rows, bic$cols], 2, Mode)  
  pattern <- rep(NA, ncol(mat))
  pattern[bic$cols] <- ori_pattern
  column_probs <- sapply(bic$cols, function(col) mean(mat[, col] == pattern[col], na.rm = T))
  pattern_prob <- prod(column_probs, na.rm = T)
  expected <- N * pattern_prob
  observed <- length(bic$rows)
  p_value <- ppois(observed - 1, lambda = expected, lower.tail = F)
  return(p_value)
})
filtered_indices <- which(p_values < 0.05)
biclusters_filtered <- biclusters[filtered_indices]


#### Flag possible cheaters and cheating items
person_det <- integer()
item_det <- integer()
for (bic in biclusters_filtered) {
  person_det <- union(person_det, bic$rows)
  item_det <- union(item_det, bic$cols)
}
flag_person <- unique(person_det)

thre_t <- median(unlist(time), na.rm=T)/4
flag_items <- integer()
for (bic in biclusters_filtered){
  time_bic <- time_data[bic$rows, bic$cols] 
  flag_items <- c(flag_items ,bic$cols[which(colMeans(time_bic) < thre_t)]) 
}
flag_items <- unique(flag_items)


#### Cheating detection results
res_tab <- matrix(nrow = 1, ncol = 4)
colnames(res_tab) <- c("examinee_sensitivity", "examinee_specificity", "item_sensitivity", "item_specificity")
res_tab[1] <- length(intersect(flag_person, true_cheater)) / length(true_cheater)
res_tab[2] <- length(intersect(setdiff(1:N, flag_person), setdiff(1:N, true_cheater))) / 
  length(setdiff(1:N, true_cheater))
res_tab[3] <- length(intersect(flag_items, true_ci)) / length(true_ci)
res_tab[4] <- length(intersect(setdiff(1:K, flag_items), setdiff(1:K, true_ci))) / 
  length(setdiff(1:K, true_ci))

# Output the final detection results
round(res_tab, 3)
