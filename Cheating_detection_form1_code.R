################################################################
#### Supplements for Simultaneous Detection of Cheaters and ####
########## Cheating Items Using Biclustering Approach ##########
############# Hyeryung Lee & Walter P. Vispoel (2024) ##########
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


#### Input Matrix Generation 
mat <- matrix(0, nrow=nrow(data), ncol=170)
for (i in 1:170){
  iraw <- paste0("iraw.", i) # response accuracy
  idur <- paste0("idur.", i) # response time
  iresp <- paste0("iresp.", i) # response choice
  
  comp <- data.frame(
    accuracy = as.numeric(data[[iraw]]), 
    time = as.numeric(data[[idur]]),
    choice = as.factor(data[[iresp]]))
  
  # Code 10 for correct and fast responses
  correct <- which(comp[,1]==1 & comp[,2]< median(comp[,2])/2)
  mat[correct, i] <- 10  
  # Use response choice for incorrect and fast responses
  incorrect <- which(comp[,1]==0 & comp[,2]<median(comp[,2])/2)
  mat[incorrect, i] <- comp[incorrect, 3] 
}

# Customized QUBIC function : exclude discretization 
qubic_no_disc <- function (x, r = 1L, q = 0.06, c = 0.95, o = 100, f = 1, 
                           k = max(ncol(x)%/%20, 2), type = "default", P = FALSE, 
                           C = FALSE, verbose = TRUE, 
                           weight = NULL, seedbicluster = NULL) 
{x_d <- x
return(qubiclust_d(x_d, c, o, f, k, type, P, C, verbose, 
                   weight, seedbicluster))
}

# Result table format
res_tab <- matrix(nrow = 1, ncol = 4)
colnames(res_tab) <- c("examinee_sensitivity", "examinee_specificity",  
                       "item_sensitivity","item_specificity")


#### Cheating detection using QUBIC
range_k <- 3  # Initial value for minimum column width of a bicluster
detect <- matrix(0, nrow= N, ncol= K) # Detection result matrix 
biclusters_all <- list()
while (TRUE) {  # Loop that will break when no bicluster is found
  res <- tryCatch(
    # Bicluster generation with the desired number of output biclusters = 100
    qubic_no_disc(mat, P = T, c = 0.8, o = 100, k = range_k, verbose = F),
    error = function(e) { NULL } 
    )
  
  if (is.null(res)) {
    print(paste("No bicluster found at k =", range_k))
    break 
  } else {
    biclusters <- list() 
    for (h in 1:res@Number) { # Collect biclusters
      bicluster <- list(
        rows = which(res@RowxNumber[, h] == TRUE), # row = examinee
        cols = which(res@NumberxCol[h, ] == TRUE)  # col = item
      )
      # Filter biclusters with > 2/3 correct responses
      if (length(bicluster$rows) < 1 ||
          length(which(colMeans(mat[bicluster$rows, bicluster$cols]) == 10)) < 
          length(bicluster$cols)/1.5) { 
        biclusters[[h]] <- NULL
      } else {
        biclusters[[h]] <- bicluster  # Collect filtered biclusters
      }
    }
    # Calculate p-values for biclusters and filter by p < 0.05
    if (length(biclusters) > 0) {
      p_mat <- matrix(NA, nrow = length(biclusters), ncol = 1) 
      for (z in 1:length(biclusters)) {
        if (is.null(biclusters[[z]]) || 
            length(biclusters[[z]]$rows) == 0 || 
            length(biclusters[[z]]$cols) == 0) {
          p_mat[z] <- NA
        } else {
          rowss <- biclusters[[z]]$rows
          colss <- biclusters[[z]]$cols
          
          
          ori_pattern <- apply(mat[rowss, colss], 2, Mode) # Pattern within a bicluster
          pattern <- rep(NA, ncol(mat))  
          pattern[colss] <- ori_pattern 
          
          column_probs <- sapply(colss, function(col) {
            mean(mat[, col] == pattern[col], na.rm = T) 
          })
          pattern_prob <- prod(column_probs, na.rm = T) # Probability of the same pattern
          expected <- nrow(mat) * pattern_prob # Expected number of rows with the same pattern
          
          # Calculate p-value using Poisson distribution
          observed <- length(rowss)
          p_value <- ppois(observed - 1, lambda = expected, lower.tail = F)
          p_mat[z] <- p_value
          low_p <- which(p_mat < 0.05) # Filter biclusters by p-value < 0.05
          
          # Collect the filtered biclusters
          person_det <- integer()
          item_det <- integer()
          for (c in low_p) {
            per_current <- biclusters[[c]]$rows
            person_det <- union(person_det, per_current)
            item_current <- biclusters[[c]]$cols
            item_det <- union(item_det, item_current)
          }
        }
      }
      detect_current <- matrix(0, nrow = N, ncol = K) # Cheating detection matrix
      detect_current[person_det, item_det] <- 1 # Cells appearing in the filtered biclusters
      detect <- detect + detect_current # Collect the cells
      
      # Collect biclusters which has p-value < 0.05
      biclusters_all <- c(biclusters[low_p], biclusters_all) 
    }
  }
  range_k <- range_k + 1 # Increment range_k for the next iteration
}

#### Flag possbile cheaters
max_detect <- max(as.numeric(names(table(detect))))
flag_person <- unique(which(detect > round(max_detect/5), arr.ind=T)[,"row"])

#### Flag possible cheating items
# Cheating item time threshold: median of all response time/4
time_start <- which(colnames(data)=="idur.1")
time <- unlist(data[,time_start:(time_start+169)])
time <- time[-which(time==0)] # Exclude the cases with time=0
median_t <- median(time)
thre_t <- median_t/4

# Collect all filtered biclusters
bic_unique <- unique(biclusters_all)
bic_unique <- Filter(Negate(is.null), bic_unique)

bic_col_time <- list()
for (i in 1:length(bic_unique)){
  # biclusters comprising only flagged examinees
  if(length(intersect(bic_unique[[i]]$rows, flag_person))==length(bic_unique[[i]]$rows)){
    row_time <- data[bic_unique[[i]]$rows, time_start:(time_start+169)] # time variables
    bic_time <- row_time[,bic_unique[[i]]$cols] # response time of a bicluster
    bic_col_time[[i]] <- (colMeans(bic_time)) # average response time for each column within a bicluster
  }else{
    bic_col_time[[i]] <- NA
  }
}
# Filter: column which has average response time within a bicluster less than median of all response time/4
filter_item <- unique(names(which(unlist(bic_col_time) < thre_t))) 
flag_item <- gsub("idur\\.","", filter_item)


#### Cheating detection results
# Examinee sensitivity
res_tab[1] <- length(intersect(flag_person, true_cheater)) / length(true_cheater)

# Examinee specificity
res_tab[2] <- length(intersect(setdiff(1:N, flag_person), setdiff(1:N, true_cheater))) / length(setdiff(1:N, true_cheater))

# Item sensitivity
res_tab[3] <- length(intersect(flag_item, true_ci)) / length(true_ci)

# Item specificity
res_tab[4] <- length(intersect(setdiff(1:K, flag_item), setdiff(1:K, true_ci))) / length(setdiff(1:K, true_ci))

# Final cheating detection results
round(res_tab, 3) 

