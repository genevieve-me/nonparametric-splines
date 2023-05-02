library(KernSmooth)
library(logspline)
library(splines)
library(pspline)

# SSE
SSE <- function(y_actual, y_fitted) {
  SSE <- sum((y_actual - y_fitted)^2, na.rm = TRUE)
}

# approximate MISE using Riemann sum
IMSE <- function(y_actual, y_fitted_mat, x_spacing) {
  # different trails are arranged as different column in y_fitted_mat
  IMSE <- x_spacing * sum(rowMeans((y_fitted_mat - y_actual)^2))
  return(IMSE)
}

# split data
k_fold_data_split <- function(data, k) {
  ran_ind <- sample(seq(1, nrow(data)), replace = FALSE)
  grouping_ind <- cut(seq(1, nrow(data)), breaks = k, labels = FALSE)
  output_list <- list(ran_ind, grouping_ind)
  names(output_list) <- c("ran_ind", "grouping_ind")
  return(output_list)
}

# KDE evaluate by the definition
density_predict <- function(olddatax, newdatax, h) {
  return(mean(dnorm(newdatax, mean = olddatax, sd = h)))
}

# # KDE alternative; slow
# density_predict = function(olddatax, newdatax, h) {
#   kernelValues <- rep(0, length(olddatax))
#     for(i in 1:length(olddatax)){
#         transformed = (newdatax - olddatax[i]) / h
#         kernelValues[i] <- dnorm(transformed, mean = 0, sd = 1) / h
#     }
#     return(sum(kernelValues) / length(olddatax))
# }

### Create simulated data for density estimation

# set RNG seed
set.seed(345)

# analytical expression for the mixture
density_expression <- function(x, p1, p2, mu1, mu2, mu3, sigma1, sigma2, sigma3) {
  p_of_x <- p1 * dnorm(x, mu1, sigma1) + p2 * dnorm(x, mu2, sigma2) +
    (1 - p1 - p2) * dnorm(v, mu3, sigma3)

  return(p_of_x)
}


# define the Gaussian mixture
density_simulator <- function(n, p1, p2, mu1, mu2, mu3, sigma1, sigma2, sigma3) {
  x <- matrix(0, n, 1)

  for (i in 1:n) {
    u <- runif(1)
    if (u < p1) {
      x[i] <- rnorm(1, mu1, sigma1)
    } else if (u < (p1 + p2)) {
      x[i] <- rnorm(1, mu2, sigma2)
    } else {
      x[i] <- rnorm(1, mu3, sigma3)
    }
  }

  return(x)
}

# generate simulated training and test data set
n_test <- 1000
p1 <- 1 / 4
p2 <- 1 / 4
mu1 <- -5
mu2 <- -2
mu3 <- 4
sigma1 <- 0.6
sigma2 <- 0.6
sigma3 <- 1.5

x_test <- density_simulator(n_test, p1, p2, mu1, mu2, mu3, sigma1, sigma2, sigma3)

# evaluation grid
lbound <- -10
ubound <- 10.46
spacing <- 0.02
num_of_points <- (ubound - lbound) / spacing + 1

### Kernel density estimation
bcv_bw <- bw.bcv(x_test)
kde_sim_display <- density(x_test,
  bw = bcv_bw,
  from = lbound, to = ubound, n = num_of_points
)
v <- kde_sim_display$x

### Spline density estimation
spline_sim_obj <- logspline(x_test, lbound, ubound)
spline_sim_display <- dlogspline(v, spline_sim_obj)

# compute approximated IMSE for kde
p_actual <- density_expression(v, p1, p2, mu1, mu2, mu3, sigma1, sigma2, sigma3)
Monte_Carlo_repetition <- 1000
p_fitted_mat <- matrix(0, nrow = length(v), ncol = Monte_Carlo_repetition)
for (i in 1:Monte_Carlo_repetition) {
  x_IMSE <- density_simulator(n_test, p1, p2, mu1, mu2, mu3, sigma1, sigma2, sigma3)

  bcv_bw <- bw.bcv(x_IMSE)
  kde_sim <- density(x_IMSE, bw = bcv_bw, from = lbound, to = ubound, n = num_of_points)
  p_fitted_mat[, i] <- kde_sim$y
}

kde_IMSE <- IMSE(p_actual, p_fitted_mat, spacing)

# compute approximated IMSE for spline
Monte_Carlo_repetition <- 1000
p_fitted_mat <- matrix(0, nrow = length(v), ncol = Monte_Carlo_repetition)
for (i in 1:Monte_Carlo_repetition) {
  x_IMSE <- density_simulator(n_test, p1, p2, mu1, mu2, mu3, sigma1, sigma2, sigma3)

  spline_sim_obj <- logspline(x_test, lbound, ubound)
  p_fitted_mat[, i] <- dlogspline(v, spline_sim_obj)
}

spline_IMSE <- IMSE(p_actual, p_fitted_mat, spacing)

#### visualization
hist(x_test,
  freq = FALSE, ylim = c(0, 0.25),
  main = "Simulated data for the density estimation problem", xlab = "x"
)
rug(x_test, ticksize = 0.05)
lines(v, density_expression(v, p1, p2, mu1, mu2, mu3, sigma1, sigma2, sigma3),
  lty = 1, lwd = 1, col = "black"
)
lines(v, kde_sim_display$y, lty = 1, lwd = 2, col = "purple")
lines(v, spline_sim_display, lty = 1, lwd = 2, col = "red")

legend("topright", c(
  "Groundtruth", paste("KDE IMSE =", round(kde_IMSE, 7)),
  paste("logspline IMSE =", round(spline_IMSE, 7))
),
col = c("black", "purple", "red"), lty = c(1, 1, 1), lwd = c(1, 2, 2)
)

rm(p_fitted_mat)


# create simulated data for regression

# set RNG seed
set.seed(345)

# define the regression function
r1 <- function(x) {
  y <- as.numeric(length(x)) # y=matrix(0,1,length(x))
  for (i in 1:length(x)) {
    if (x[i] < 0) {
      y[i] <- 0
    } else if (x[i] < 1) {
      y[i] <- 1
    } else if (x[i] < 2.5) {
      y[i] <- 3
    } else if (x[i] < 3) {
      y[i] <- 0.5
    } else if (x[i] < 5) {
      y[i] <- 2
    } else if (x[i] < 10) {
      y[i] <- 6
    } else {
      y[i] <- 0
    }
  }
  return(y)
}

# Generate simulated data
spacing <- .01
noise_mean <- 0
noise_sigma <- 0.5
x2 <- seq(from = -2, to = 12, by = spacing)
y2 <- r1(x2) + rnorm(length(x2), noise_mean, noise_sigma)
sim_staircase_data <- data.frame(x2, y2)

# x grid for evaluation
grid_x_2 <- sim_staircase_data$x2

# split data into training and test data sets
train_indices <- sample(1:nrow(sim_staircase_data), round(nrow(sim_staircase_data) * 0.75),
  replace = FALSE
)
sim_staircase_data_train_data <- sim_staircase_data[train_indices, ]
sim_staircase_data_test_data <- sim_staircase_data[-train_indices, ]

# split training data
k <- 3
training_folds <- k_fold_data_split(sim_staircase_data_train_data, k)

##### 3rd degree local polynomial model (bandwidth selection)
# y2_locpoly = locpoly(x2, y2, degree = 3, bandwidth = 0.25)
bw_search <- seq(0.01, 0.05, 0.0025)

L2_error_mat <- matrix(0, length(bw_search), k)
for (i in 1:length(bw_search)) {
  for (j in 1:k) {
    sim_staircase_data_train_data_not_j <-
      sim_staircase_data_train_data[training_folds$ran_ind[training_folds$grouping_ind != j], ]
    y2_loess <-
      loess(y2 ~ x2, span = bw_search[i], degree = 2, data = sim_staircase_data_train_data_not_j)
    # this loess predict will produce NA sometimes
    y2_loess_predict <- predict(
      y2_loess,
      sim_staircase_data_train_data$x2[training_folds$ran_ind[training_folds$grouping_ind == j]]
    )
    L2_error_mat[i, j] <-
      SSE(sim_staircase_data_train_data$y2[training_folds$ran_ind[training_folds$grouping_ind == j]], y2_loess_predict)
  }
}

L2_error <- rowMeans(L2_error_mat)
cv_bw_loess <- bw_search[which.min(L2_error)]
plot(bw_search, L2_error, main = "local poly bandwidth selection")
abline(v = cv_bw_loess, col = "red")
legend("topright", c("L2 error", "Min error bandwidth"),
  pch = c(1, NA), lty = c(NA, 1), col = c("black", "red")
)

# 3rd degree local polynomial model
y2_loess <- loess(y2 ~ x2,
  span = cv_bw_loess, degree = 2,
  data = sim_staircase_data_train_data
)
y2_loess_predict <- predict(y2_loess, grid_x_2)

##### Splines model (df selection)
bw_search <- seq(60, 80, 1)

L2_error_mat <- matrix(0, length(bw_search), k)
for (i in 1:length(bw_search)) {
  for (j in 1:k) {
    sim_staircase_data_train_data_not_j <-
      sim_staircase_data_train_data[training_folds$ran_ind[training_folds$grouping_ind != j], ]

    y2_spline <- lm(y2 ~ bs(x2, df = bw_search[i], degree = 3),
      data = sim_staircase_data_train_data_not_j
    )
    y2_spline_predict <- predict(
      y2_spline,
      data.frame(x2 = sim_staircase_data_train_data$x2[training_folds$ran_ind[training_folds$grouping_ind == j]])
    )

    L2_error_mat[i, j] <-
      SSE(
        sim_staircase_data_train_data$y2[training_folds$ran_ind[training_folds$grouping_ind == j]],
        y2_spline_predict
      )
  }
}

L2_error <- rowMeans(L2_error_mat)
cv_bw_spline <- bw_search[which.min(L2_error)]
plot(bw_search, L2_error, main = "Spline bandwidth selection")
abline(v = cv_bw_spline, col = "red")
legend("topright", c("L2 error", "Min error bandwidth"),
  pch = c(1, NA), lty = c(NA, 1), col = c("black", "red")
)

# Splines model
y2_spline <- lm(y2 ~ bs(x2, df = cv_bw_spline, degree = 3), data = sim_staircase_data_train_data)
y2_spline_predict <- predict(y2_spline, data.frame(x2 = grid_x_2))

# compute approximated IMSE for locpoly
Monte_Carlo_repetition <- 1000
y_fitted_mat <- matrix(0, nrow = length(grid_x_2), ncol = Monte_Carlo_repetition)
for (i in 1:Monte_Carlo_repetition) {
  spacing <- .01
  noise_mean <- 0
  noise_sigma <- 0.5
  x2 <- seq(from = -2, to = 12, by = spacing)
  y2 <- r1(x2) + rnorm(length(x2), noise_mean, noise_sigma)
  sim_staircase_data_IMSE <- data.frame(x2, y2)

  y2_loess <- loess(y2 ~ x2, span = cv_bw_loess, degree = 2, data = sim_staircase_data_IMSE)
  y_fitted_mat[, i] <- predict(y2_loess, grid_x_2)
}

loess_IMSE <- IMSE(r1(x2), y_fitted_mat, spacing)

# compute approximated IMSE for spline
Monte_Carlo_repetition <- 1000
y_fitted_mat <- matrix(0, nrow = length(grid_x_2), ncol = Monte_Carlo_repetition)
for (i in 1:Monte_Carlo_repetition) {
  spacing <- .01
  noise_mean <- 0
  noise_sigma <- 0.5
  x2 <- seq(from = -2, to = 12, by = spacing)
  y2 <- r1(x2) + rnorm(length(x2), noise_mean, noise_sigma)
  sim_staircase_data_IMSE <- data.frame(x2, y2)

  y2_spline <- lm(y2 ~ bs(x2, df = cv_bw_spline, degree = 3), data = sim_staircase_data_IMSE)
  y_fitted_mat[, i] <- predict(y2_spline, data.frame(x2 = grid_x_2))
}

spline_IMSE <- IMSE(r1(x2), y_fitted_mat, spacing)

# Visualization
plot(sim_staircase_data_test_data$x2, sim_staircase_data_test_data$y2,
  xlab = "Predictor x", ylab = "Response y",
  main = "Simulated data for the regression problem", cex = 0.4
)
lines(grid_x_2, r1(grid_x_2), type = "l", lty = 1, lwd = 1, col = "black")
# lines(y2_locpoly$x, y2_locpoly$y, type="l", lty=1, lwd = 2,  col="green")
lines(grid_x_2, y2_loess_predict, type = "l", lty = 1, lwd = 2, col = "purple")
lines(grid_x_2, y2_spline_predict, lty = 1, lwd = 2, col = "red")
# lines(grid_x_2, y2_penalized_spline_predict, lty=1, lwd = 2, col="yellow")
legend("topleft", c(
  "Testset data", "Groundtruth",
  paste("Local poly; MSE =", round(loess_IMSE, 4)),
  paste("Spline; MSE =", round(spline_IMSE, 4))
),
lty = c(NA, 1, 1, 1), pch = c(1, NA, NA, NA),
col = c("black", "black", "purple", "red"), lwd = c(NA, 1, 2, 2)
)

##### penalized spline

# penalized splines model (bandwidth selection)
bw_search <- seq(0, 0.04, 0.001)

L2_error_mat <- matrix(0, length(bw_search), k)
for (i in 1:length(bw_search)) {
  for (j in 1:k) {
    sim_staircase_data_train_data_not_j <- sim_staircase_data_train_data[training_folds$ran_ind[training_folds$grouping_ind != j], ]

    sim_staircase_data_train_data_not_j <- sim_staircase_data_train_data_not_j[order(sim_staircase_data_train_data_not_j$x2), ]

    y2_penalized_spline <- smooth.Pspline(sim_staircase_data_train_data_not_j$x2, sim_staircase_data_train_data_not_j$y2, norder = 2, spar = bw_search[i])
    y2_penalized_spline_predict <- predict(y2_penalized_spline, sim_staircase_data_train_data$x2[training_folds$ran_ind[training_folds$grouping_ind == j]])

    L2_error_mat[i, j] <- SSE(sim_staircase_data_train_data$y2[training_folds$ran_ind[training_folds$grouping_ind == j]], y2_penalized_spline_predict)
  }
}

L2_error <- rowMeans(L2_error_mat)
cv_bw_pspline <- bw_search[which.min(L2_error)]
plot(bw_search, L2_error, main = "Penalized spline lambda selection")
abline(v = cv_bw_pspline, col = "red")
legend("topright", c("L2 error", "Min error lambda"), pch = c(1, NA), lty = c(NA, 1), col = c("black", "red"))

# penalized splines model
sim_staircase_data_train_data <- sim_staircase_data_train_data[order(sim_staircase_data_train_data$x2), ]

y2_penalized_spline <- smooth.Pspline(sim_staircase_data_train_data$x2, sim_staircase_data_train_data$y2, norder = 2, spar = cv_bw_pspline)
y2_penalized_spline_predict <- predict(y2_penalized_spline, grid_x_2)

# compute IMSE for pspline
Monte_Carlo_repetition <- 1000
y_fitted_mat <- matrix(0, nrow = length(grid_x_2), ncol = Monte_Carlo_repetition)
for (i in 1:Monte_Carlo_repetition) {
  spacing <- .01
  noise_mean <- 0
  noise_sigma <- 0.5
  x2 <- seq(from = -2, to = 12, by = spacing)
  y2 <- r1(x2) + rnorm(length(x2), noise_mean, noise_sigma)
  sim_staircase_data_IMSE <- data.frame(x2, y2)

  sim_staircase_data_IMSE <- sim_staircase_data_IMSE[order(sim_staircase_data_IMSE$x2), ]

  y2_penalized_spline <- smooth.Pspline(sim_staircase_data_IMSE$x2, sim_staircase_data_IMSE$y2, norder = 2, spar = cv_bw_pspline)
  y_fitted_mat[, i] <- predict(y2_penalized_spline, grid_x_2)
}

Pspline_IMSE <- IMSE(r1(x2), y_fitted_mat, spacing)

# Visualization
plot(sim_staircase_data_test_data$x2, sim_staircase_data_test_data$y2,
  xlab = "Predictor x", ylab = "Response y",
  main = "Simulated data for the regression problem", cex = 0.4
)
lines(grid_x_2, r1(grid_x_2), type = "l", lty = 1, lwd = 1, col = "black")
lines(grid_x_2, y2_spline_predict, lty = 1, lwd = 2, col = "red")
lines(grid_x_2, y2_penalized_spline_predict, lty = 1, lwd = 2, col = "orange")
legend("topleft", c(
  "Testset data", "Groundtruth", paste("Spline; MSE =", round(spline_IMSE, 4)),
  paste("Pspline; MSE =", round(Pspline_IMSE, 4))
),
lty = c(NA, 1, 1, 1),
pch = c(1, NA, NA, NA),
col = c("black", "black", "red", "orange"), lwd = c(NA, 1, 2, 2)
)

# pspline lambda effect

# need variables from previous code chunk
lambda_vec <- c(0, 0.001, 0.01, 0.1, 1)
y2_penalized_spline_predict_mat <- matrix(0, nrow = length(grid_x_2), ncol = length(lambda_vec))

for (i in 1:length(lambda_vec)) {
  y2_penalized_spline <- smooth.Pspline(sim_staircase_data_train_data$x2,
    sim_staircase_data_train_data$y2,
    norder = 2, spar = lambda_vec[i]
  )
  y2_penalized_spline_predict_mat[, i] <- predict(y2_penalized_spline, grid_x_2)
}

# Visualization
plot(sim_staircase_data_test_data$x2, sim_staircase_data_test_data$y2,
  xlab = "Predictor x", ylab = "Response y",
  main = "Simulated data for the regression problem", cex = 0.4
)
lines(grid_x_2, r1(grid_x_2), type = "l", lty = 1, lwd = 1, col = "black")
lines(grid_x_2, y2_spline_predict, lty = 1, lwd = 2, col = "red")
lines(grid_x_2, y2_penalized_spline_predict_mat[, 1], lty = 1, lwd = 2, col = "pink")
lines(grid_x_2, y2_penalized_spline_predict_mat[, 2], lty = 1, lwd = 2, col = "yellow")
lines(grid_x_2, y2_penalized_spline_predict_mat[, 3], lty = 1, lwd = 2, col = "green")
lines(grid_x_2, y2_penalized_spline_predict_mat[, 4], lty = 1, lwd = 2, col = "cyan")
lines(grid_x_2, y2_penalized_spline_predict_mat[, 5], lty = 1, lwd = 2, col = "purple")

legend("topleft", c(
  "Testset data", "Groundtruth", "Spline", "Pspline lambda = 0",
  "Pspline lambda = 0.001", "Pspline lambda = 0.01", "Pspline lambda = 0.1", "Pspline lambda = 1"
),
lty = c(NA, 1, 1, 1, 1, 1, 1, 1), pch = c(1, NA, NA, NA, NA, NA, NA, NA),
col = c("black", "black", "red", "pink", "yellow", "green", "cyan", "purple"),
lwd = c(NA, 1, 2, 2, 2, 2, 2, 2), cex = 0.6
)


# real data for density estimation

# set RNG seed
set.seed(345)

# load data
Fiji_data <- read.table("/Users/Samson/Documents/Github/496-splines/fijiquakes.dat.txt", header = TRUE)

# split data into training and test data sets
train_indices <- sample(1:nrow(Fiji_data), round(nrow(Fiji_data) * 0.25),
  replace = FALSE
)

Fiji_train_data <- Fiji_data[train_indices, ]
Fiji_test_data <- Fiji_data[-train_indices, ]

# evaluation grid
lbound <- -30
ubound <- 750
spacing <- 0.75
num_of_points <- 1024

# kernel density estimation
# ucv_bw = bw.ucv(Fiji_train_data$depth) # failed
bcv_bw <- bw.bcv(Fiji_train_data$depth)
kde_fiji <- density(Fiji_train_data$depth, bw = bcv_bw, from = lbound, to = ubound, n = num_of_points)
fiji_grid <- kde_fiji$x

# spline model
spline_fiji_obj <- logspline(Fiji_train_data$depth, lbound, ubound)
spline_sim <- dlogspline(fiji_grid, spline_fiji_obj)

# find empirical cdf for kde, spline and the test data
test_ecdf <- ecdf(Fiji_test_data$depth)
kde_cdf <- cumsum(kde_fiji$y) * (fiji_grid[2] - fiji_grid[1])
spline_cdf <- cumsum(spline_sim) * (fiji_grid[2] - fiji_grid[1])

kde_l2 <- sqrt(sum((kde_cdf - test_ecdf(fiji_grid))^2))
spline_l2 <- sqrt(sum((spline_cdf - test_ecdf(fiji_grid))^2))

plot(fiji_grid, test_ecdf(fiji_grid),
  type = "l",
  main = "Test set empirical cdf and estimated cdf", lwd = 2
)
lines(fiji_grid, kde_cdf, col = "purple", lwd = 2)
lines(fiji_grid, spline_cdf, col = "red", lwd = 2)
legend("topleft", c("Testset ecdf", paste(
  "KDE cdf; L2 error =",
  round(kde_l2, 4)
), paste("logspline cdf; L2 error =", round(spline_l2, 4))),
lty = c(1, 1, 1), col = c("black", "purple", "red"), lwd = c(2, 2, 2)
)

# compute the log likelihood of each density given the test data
kde_pdf_val <- rep(0, length(Fiji_test_data$depth))
for (i in 1:length(Fiji_test_data$depth)) {
  kde_pdf_val[i] <- density_predict(Fiji_train_data$depth, Fiji_test_data$depth[i], bcv_bw)
}
LL_kde <- sum(log(kde_pdf_val))

spline_pdf_val <- dlogspline(Fiji_test_data$depth, spline_fiji_obj)
LL_spline <- sum(log(spline_pdf_val))

LLR <- LL_spline / LL_kde

# visualization
hist(Fiji_test_data$depth,
  freq = FALSE,
  main = "Fiji earthquake depth data (test set)",
  xlab = "Depth of the Fiji earthquake", ylim = c(0, 0.0085)
)
rug(Fiji_test_data$depth, ticksize = 0.05)
lines(kde_fiji, lty = 1, lwd = 2, col = "purple")
lines(fiji_grid, spline_sim, lty = 1, lwd = 2, col = "red")
legend("top", c(
  paste("KDE; log-likelihood =", round(LL_kde, 3)),
  paste("logspline; log-likelihood =", round(LL_spline, 3))
),
lty = c(1, 1), col = c("purple", "red"), lwd = c(2, 2)
)

# real data for regression

# set RNG seed
set.seed(345)

# load data
lidar_data <- read.table("/Users/Samson/Documents/Github/496-splines/lidar.dat.txt", header = TRUE)

# evaluation grid
grid_lidar <- seq(from = 395, to = 715, by = 1)

# split data into training and testing data
train_indices <- sample(1:nrow(lidar_data), round(nrow(lidar_data) * 0.7),
  replace = FALSE
)

lidar_train_data <- lidar_data[train_indices, ]
lidar_test_data <- lidar_data[-train_indices, ]

# split training data
k <- 3
training_folds <- k_fold_data_split(lidar_train_data, k)

# local polynomial model (bandwidth selection)
bw_search <- seq(0.02, 1, 0.02)

L2_error_mat <- matrix(0, length(bw_search), k)
for (i in 1:length(bw_search)) {
  for (j in 1:k) {
    lidar_train_data_not_j <- lidar_train_data[training_folds$ran_ind[training_folds$grouping_ind != j], ]
    y2_loess <- loess(logratio ~ range, span = bw_search[i], degree = 2, data = lidar_train_data_not_j)
    y2_loess_predict <- predict(y2_loess, lidar_train_data$range[training_folds$ran_ind[training_folds$grouping_ind == j]])
    L2_error_mat[i, j] <- SSE(lidar_train_data$logratio[training_folds$ran_ind[training_folds$grouping_ind == j]], y2_loess_predict)
  }
}

L2_error <- rowMeans(L2_error_mat)
cv_bw_loess <- bw_search[which.min(L2_error)]
plot(bw_search, L2_error, main = "local poly bandwidth selection")
abline(v = cv_bw_loess, col = "red")
legend("topright", c("L2 error", "Min error bandwidth"), pch = c(1, NA), lty = c(NA, 1), col = c("black", "red"))


# local polynomial model
y2_loess <- loess(logratio ~ range, span = cv_bw_loess, degree = 2, data = lidar_train_data)
y2_lidar_loess_predict <- predict(y2_loess, grid_lidar)

# spline model (bandwidth selection)
bw_search <- seq(10, 30, 1)

L2_error_mat <- matrix(0, length(bw_search), k)
for (i in 1:length(bw_search)) {
  for (j in 1:k) {
    lidar_train_data_not_j <- lidar_train_data[training_folds$ran_ind[training_folds$grouping_ind != j], ]
    y2_spline <- lm(logratio ~ bs(range, df = bw_search[i], degree = 3), data = lidar_train_data_not_j)
    y2_spline_predict <- predict(y2_spline, data.frame(range = lidar_train_data$range[training_folds$ran_ind[training_folds$grouping_ind == j]]))
    L2_error_mat[i, j] <- SSE(lidar_train_data$logratio[training_folds$ran_ind[training_folds$grouping_ind == j]], y2_spline_predict)
  }
}

L2_error <- rowMeans(L2_error_mat)
cv_bw_spline <- bw_search[which.min(L2_error)]
plot(bw_search, L2_error, main = "Spline bandwidth selection")
abline(v = cv_bw_spline, col = "red")
legend("topright", c("L2 error", "Min error bandwidth"), pch = c(1, NA), lty = c(NA, 1), col = c("black", "red"))


# spline model
lidar_spline <- lm(logratio ~ bs(range, df = cv_bw_spline, degree = 3), data = lidar_train_data)
y2_lidar_spline_predict <- predict(lidar_spline, data.frame(range = grid_lidar))

# find SSE for kde
SSE_loess <- SSE(lidar_test_data$logratio, y2_lidar_loess_predict)

# find SSE for spline
SSE_spline <- SSE(lidar_test_data$logratio, y2_lidar_spline_predict)

# Visualization
plot(lidar_test_data$range, lidar_test_data$logratio,
  main = "LIDAR data", ylab = "logratio", xlab = "range", cex = 0.4
)
lines(grid_lidar, y2_lidar_loess_predict, lty = 1, lwd = 2, col = "purple")
lines(grid_lidar, y2_lidar_spline_predict, lty = 1, lwd = 2, col = "red")
legend("bottomleft", c(
  "Test data", paste("Local poly; SSE =", round(SSE_loess, 5)),
  paste("Spline; SSE =", round(SSE_spline, 5))
),
lty = c(NA, 1, 1),
pch = c(1, NA, NA), lwd = c(NA, 2, 2),
col = c("black", "purple", "red")
)
