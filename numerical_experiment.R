# Load required library for sparse matrix operations
if (!require("Matrix")) install.packages("Matrix")
library(Matrix)

## ==============================================================================
# PHASE 1: INITIALIZATION & PARAMETERS ==========================================
## ==============================================================================

T_end <- 10
dt    <- 1/252
Nt    <- floor(T_end / dt)

X_min <- -1000;  X_max <- 1000; Nx <- 8001

dx <- (X_max - X_min) / (Nx - 1)

X <- seq(X_min, X_max, length.out = Nx)

# Parameters required for the simulation (using the values from Colaneri, Frey, and Koeck (2024) Section 6.3)
params <- list(
  q_max     = 10,
  p_g       = 0.2,
  x_bar     = 20,
  theta_bar = 3,
  theta_min = 0,
  theta_max = 5.1 ,
  nu        = 1,
  py        = 2.1,
  r         =  0.04,
  sigma     = 0.2,
  delta     = 0.02,
  kappa     = 0.5,
  terminal_discount = 0.7
)

tol      <- 1e-6
max_iter <- 50

# Define A4 dimensions in inches (Landscape orientation)
a4_width  <- 11.69
a4_height <- 8.27
resolution <- 300 

## ==============================================================================
# PHASE 1.1: Functions ==========================================================
## ==============================================================================

# Function to calculate q_tilde given x and the parameter list
calculate_q_tilde <- function(x, prm, P_gx, e_x) {
  with(prm, {
    
    # Case 1: py >= theta_bar (or e_x is 0 to avoid division by zero)
    if ((py >= theta_bar) || e_x == 0) {
      q_tilde <- (4/9) * py^2 + P_gx
      return(min(q_tilde, q_max))
    } 
    
    # Case 2: py < theta_bar
    # Calculate q_22 for condition B1. 
    q_22 <- (4/9) * ((max(py, theta_min) - py) * e_x + py)^2 + P_gx
    
    if (e_x * q_22 > 2 * nu * (theta_bar - py)) {
      # Condition B1 is met
      q_tilde <- q_22
      
    } else {
      # Calculate q_212 for condition B2
      q_212 <- (4/9) * (theta_min * e_x + (1 - e_x) * py)^2 + P_gx
      
      if ((2 * nu * (theta_bar - theta_min) < e_x * q_212) && (e_x * q_212 <= 2 * nu * (theta_bar - py))) {
        # Condition B2 is met
        q_tilde <- q_212
        
      } else {
        # Condition B3 (Fallback to the quadratic solution)
        a <- (e_x^4) / (4 * nu^2)
        b <- -(9/4  + (e_x^2 / nu) * ((theta_bar - py) * e_x + py))
        c_val <- ((theta_bar - py) * e_x + py)^2 + (9/4) * P_gx
        
        discriminant <- b^2 - 4 * a * c_val
        
        # Protect against floating point inaccuracies resulting in tiny negatives
        if (discriminant < 0) {
          discriminant <- 0
        }
        
        # Calculate the strict minus branch of the quadratic
        q_cand <- (-b - sqrt(discriminant)) / (2 * a)
        
        # Apply the P_g(x) \vee ... logic from the math derivation
        q_tilde <- max(P_gx, q_cand)
      }
    }
    
    # Apply the final upper bound q_max for all Case 2 branches
    return(min(q_tilde, q_max))
  })
}

calculate_theta_tilde <- function(x, q, prm, Pg, e_x){
  if (prm$py > prm$theta_bar){
    return(prm$theta_bar)
  } else{
    if ( e_x*q <= 2*prm$nu*(prm$theta_bar - prm$py)){
      return(prm$theta_bar - e_x*q/(2*prm$nu))
    } else {
      return(prm$py)
    }
  }
}

saddle_value_i <- function(x, prm){
  # Evaluate the function g(x, y, q, theta)
  P_gx <- prm$p_g * max(x - prm$x_bar, 0)
  e_x <- min(P_gx / prm$q_max, 1)
  
  q <- calculate_q_tilde(x, prm, P_gx, e_x)
  theta <- calculate_theta_tilde(x, q, prm , P_gx, e_x)
  
  term1 <- prm$py * q                         # Revenue
  term2 <- -max(q - P_gx, 0)^1.5              # Production Cost C(q, x, y)
  term3 <- max(theta - prm$py, 0) * e_x * q   # Subsidy paid
  term4 <- prm$nu * (theta - prm$theta_bar)^2 # Penalty
  
  g_val <- term1 + term2 + term3 + term4
  
  return(c(q,theta,g_val))
}

simulate_fixed_p <- function(theta_bar_input, nu_input, params){
  
  ## ==============================================================================
  # PHASE 2: GRID SETUP & PRE-COMPUTATION =========================================
  ## ==============================================================================
  
  # Adjust the parameters for uncertainty nu and price floor benchmark theta_bar
  params$nu <- nu_input
  params$theta_bar <- theta_bar_input
  
  # Pre-compute heavy constants outside the time loop
  c_sig2_dx2  <- dt * (params$sigma^2) / (dx^2)
  c_sig2_2dx2 <- c_sig2_dx2 / 2
  c_dt_dx <- dt/dx
  
  # set interior X
  X_int <- X[c(-1,-Nx)]
  
  # Save game values
  U_res <- as.list(1:(Nt+1))
  Gamma_res <- as.list(1:Nt)
  Howard_error<-as.numeric(Nt)
  Howard_iter<-as.numeric(Nt)
  
  # Terminal condition
  U_res[[Nt +1 ]] <- pmax(params$terminal_discount * X, 0)    
  
  # Saddle Points
  cat("Pre-computing Saddle point and value Matrices...\n")
  q_matrix <- matrix(0, nrow = Nx, ncol = 1)
  theta_matrix <- matrix(0, nrow = Nx, ncol = 1)
  G_matrix <- matrix(0, nrow = Nx, ncol = 1)
  for (i in 1:Nx) {
    res <- saddle_value_i(X[i], params)
    theta_matrix[i] <- res[2]
    q_matrix[i] <- res[1]
    G_matrix[i] <- res[3]
  }
  cat("... Saddle point and value computation complete.\n")
  
  
  
  ## ==============================================================================
  # PHASE 3: BACKWARD TIME-STEPPING ===============================================
  ## ==============================================================================
  
  cat("Starting backward time integration...\n")
  start_time <- Sys.time()
  
  # First is terminal condition 
  U_prev <- U_res[[Nt+1]]
  for (i in Nt:1){
    U_k <- U_prev
    
    for (k in 1:max_iter){
      
      ### Step 3.1: Calculate optimal gamma for this iteration ====
      
      u_x_fwd <- (U_prev[-(1:2)]-U_prev[-c(1,Nx)])/dx
      u_x_bwd <- (U_prev[-c(1,Nx)]-U_prev[-c(Nx-1,Nx)])/dx
      
      gamma_fwd <- pmax(u_x_fwd - 1,0)/(2*params$kappa)
      gamma_bwd <- pmax(u_x_bwd - 1,0)/(2*params$kappa)
      
      sup_g_fwd <- pmax(gamma_fwd - params$delta * X_int,0)*u_x_fwd + pmin(gamma_fwd - params$delta * X_int,0)*u_x_bwd - gamma_fwd - params$kappa * gamma_fwd^2
      sup_g_bwd <- pmax(gamma_bwd - params$delta * X_int,0)*u_x_fwd + pmin(gamma_bwd - params$delta * X_int,0)*u_x_bwd - gamma_bwd - params$kappa * gamma_bwd^2
      
      gamma_opt <-ifelse(sup_g_fwd>sup_g_bwd, gamma_fwd, gamma_bwd)
      
      
      ### Step 3.2: Vectorized M-Matrix Assembly (Interior + Boundaries) ====
      
      # Calculate interior diagonals
      lower_diag <- -c_sig2_2dx2 + c_dt_dx*pmin(gamma_opt - params$delta*X[-c(1,Nx)],0)
      lower_coord_i <- 2:(Nx-1)
      lower_coord_j <- 1:(Nx-2)
      
      upper_diag <- -c_sig2_2dx2 - c_dt_dx*pmax(gamma_opt - params$delta*X[-c(1,Nx)],0)
      upper_coord_i <- 2:(Nx-1)
      upper_coord_j <- 3:Nx
      
      center_diag <- rep(1+dt*params$r,Nx-2) - lower_diag - upper_diag
      center_coord_i <- 2:(Nx-1)
      center_coord_j <- 2:(Nx-1)
      
      #### LOWER BOUNDARY: LINEAR EXTRAPOLATION (u_xx = 0) ====
      # Enforces: 1*U[1] - 2*U[2] + 1*U[3] = 0 
      left_i <- c(1, 1, 1)
      left_j <- c(1, 2, 3)
      left_v <- c(1, -2, 1)
      
      #### UPPER BOUNDARY: LINEAR EXTRAPOLATION (u_xx = 0) ====
      # Enforces: 1*U[Nx] - 2*U[Nx-1] + 1*U[Nx-2] = 0
      right_i <- c(Nx, Nx, Nx)
      right_j <- c(Nx, Nx-1, Nx-2)
      right_v <- c(1, -2, 1)
      
      #### COMBINE ALL TRIPLETS ====
      coord_values <- c(left_v, right_v, lower_diag, upper_diag, center_diag)
      coord_full_i <- c(left_i, right_i, lower_coord_i, upper_coord_i, center_coord_i)
      coord_full_j <- c(left_j, right_j, lower_coord_j, upper_coord_j, center_coord_j)
      
      A <- sparseMatrix(i = coord_full_i, 
                        j = coord_full_j,
                        x = coord_values)
      
      #### ASSEMBLE RIGHT-HAND SIDE VECTOR B ====
      # Both boundary rows equal 0 on the right-hand side of their algebraic equations
      B <- c(0, 
             U_res[[i+1]][-c(1,Nx)] - (gamma_opt + params$kappa * gamma_opt^2 - G_matrix[-c(1,Nx)])*dt, 
             0)
      
      U_new <- solve(A, B)
      
      error <- max(abs(U_new - U_prev))
      U_prev <- U_new
      if (error < tol) {break}
    }
    
    Howard_error[i] <- error
    Howard_iter[i] <- k
    Gamma_res[[i]] <- gamma_opt
    U_res[[i]] <- U_prev
    
    if (i %% max(1, floor(Nt / 10)) == 0) {
      cat(sprintf("... Time Step %d / %d completed. (t = %.2f)\n", i, Nt, i*dt))
    }
  }
  return(list("G"=G_matrix,
              "q"=q_matrix,
              "theta"=theta_matrix,
              "gamma"=Gamma_res,
              "value"=U_res,
              "h_error"=Howard_error,
              "h_iter"=Howard_iter))
  
}

## ==============================================================================
# EVALUATE FINITE DIFFERENCE SCHEME =============================================
## ==============================================================================

th_vals <- c(5, 3, 2.3, 2.1)
uncertainty_nu <- c(0.1, 1, 100)

cat("Calculating Base Scenarios (Varying Nu) sequentially...\n")
res_high_unc <- lapply(th_vals, function(th) simulate_fixed_p(th, uncertainty_nu[1], params))
res_mid_unc  <- lapply(th_vals, function(th) simulate_fixed_p(th, uncertainty_nu[2], params))
res_low_unc  <- lapply(th_vals, function(th) simulate_fixed_p(th, uncertainty_nu[3], params))

## ==============================================================================
# PHASE 4.0: Saddle point processes =============================================
## ==============================================================================

### Plot for Optimal price floor (theta) ====
png("optimal_subsidy_theta.png",  width = a4_width, height = a4_height, units = "in", res = resolution)
#width = 800, height = 600)
par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
for (k in 1:4) {
  idx_low <- round(Nx / 2,0)
  idx_high <- idx_low + round(100/dx,0) + 1
  x_sub <- X[idx_low:idx_high]
  t_low <- res_low_unc[[k]]$theta[idx_low:idx_high]
  t_mid <- res_mid_unc[[k]]$theta[idx_low:idx_high]
  t_high <- res_high_unc[[k]]$theta[idx_low:idx_high]
  
  plot(x_sub, t_low, type = "l", lwd = 2, col = "darkgreen",
       main = paste("theta_bar =", th_vals[k]),
       xlab = "Investment (x)", ylab = "Price Floor (theta)",
       ylim = range(c(t_low, t_mid, t_high, params$py, th_vals[k])))
  
  lines(x_sub, t_mid, col = "black", lwd = 2)
  lines(x_sub, t_high, col = "orange", lwd = 2)
  
  abline(h = params$py, lty = 2, col = "darkgrey")
  abline(h = th_vals[k], lty = 2, col = "darkgrey")
  
  if (k == 1) {
    legend("bottomright", legend = paste("nu =", uncertainty_nu), 
           col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n")
  }
}
mtext("Optimal Price Floor (theta) across scenarios", outer = TRUE, cex = 1.2, font = 2)
dev.off()

### Plot for Optimal Quantity (q) ====
png("optimal_production_q.png",  width = a4_width, height = a4_height, units = "in", res = resolution)
##width = 800, height = 600)
par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
for (k in 1:4) {
  idx_low <- round(Nx / 2,0)
  idx_high <- idx_low + round(100/dx,0) + 1
  x_sub <- X[idx_low:idx_high]
  q_low <- res_low_unc[[k]]$q[idx_low:idx_high]
  q_mid <- res_mid_unc[[k]]$q[idx_low:idx_high]
  q_high <- res_high_unc[[k]]$q[idx_low:idx_high]
  
  plot(x_sub, q_low, type = "l", lwd = 2, col = "darkgreen",
       main = paste("theta_bar =", th_vals[k]),
       xlab = "Investment (x)", ylab="Quantity (q)",
       ylim = range(c(q_low, q_mid, q_high, params$q_max)))
  
  lines(x_sub, q_mid, col = "black", lwd = 2)
  lines(x_sub, q_high, col = "orange", lwd = 2)
  
  abline(h = params$q_max, lty = 2, col = "darkgrey")
  abline(h = 1.96, lty = 2, col = "darkgrey")
  
  if (k == 1) {
    legend("bottomright", legend = paste("nu =", uncertainty_nu), 
           col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n")
  }
}
mtext("Optimal Quantity (q) across scenarios", outer = TRUE, cex = 1.2, font = 2)
dev.off()

## ==============================================================================
# PHASE 4.0.1: Phase Plots (theta vs q and q vs theta) ==========================
## ==============================================================================

### Optimal Price Floor (theta) given Quantity (q) ====
png("optimal_subsidy_given_q.png",  width = a4_width, height = a4_height, units = "in", res = resolution)
    #width = 800, height = 600)
par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
for (k in 1:4) {
  q_low  <- res_low_unc[[k]]$q
  t_low  <- res_low_unc[[k]]$theta
  
  q_mid  <- res_mid_unc[[k]]$q
  t_mid  <- res_mid_unc[[k]]$theta
  
  q_high <- res_high_unc[[k]]$q
  t_high <- res_high_unc[[k]]$theta
  
  plot(q_low, t_low, type = "l", lwd = 2, col = "darkgreen",
       main = paste("theta_bar =", th_vals[k]),
       xlab = "Quantity (q)", ylab = "Optimal Price Floor (theta)",
       xlim = range(c(q_low, q_mid, q_high)),
       ylim = range(c(t_low, t_mid, t_high, params$py, th_vals[k])))
  
  lines(q_mid, t_mid, col = "black", lwd = 2)
  lines(q_high, t_high, col = "orange", lwd = 2)
  
  abline(h = params$py, lty = 2, col = "darkgrey")
  abline(h = th_vals[k], lty = 2, col = "darkgrey")
  
  if (k == 1) {
    legend("topright", legend = paste("nu =", uncertainty_nu), 
           col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n")
  }
}
mtext("Optimal Price Floor (theta) given Quantity (q)", outer = TRUE, cex = 1.2, font = 2)
dev.off()

### Optimal Quantity (q) given Price Floor (theta) ====
png("optimal_production_given_theta.png",  width = a4_width, height = a4_height, units = "in", res = resolution)
    #width = 800, height = 600)
par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
for (k in 1:4) {
  q_low  <- res_low_unc[[k]]$q
  t_low  <- res_low_unc[[k]]$theta
  
  q_mid  <- res_mid_unc[[k]]$q
  t_mid  <- res_mid_unc[[k]]$theta
  
  q_high <- res_high_unc[[k]]$q
  t_high <- res_high_unc[[k]]$theta
  
  plot(t_low, q_low, type = "l", lwd = 2, col = "darkgreen",
       main = paste("theta_bar =", th_vals[k]),
       xlab = "Price Floor (theta)", ylab = "Optimal Quantity (q)",
       xlim = range(c(t_low, t_mid, t_high, params$py, th_vals[k])),
       ylim = range(c(q_low, q_mid, q_high)))
  
  lines(t_mid, q_mid, col = "black", lwd = 2)
  lines(t_high, q_high, col = "orange", lwd = 2)
  
  abline(v = params$py, lty = 2, col = "darkgrey")
  abline(v = th_vals[k], lty = 2, col = "darkgrey")
  
  if (k == 1) {
    legend("bottomright", legend = paste("nu =", uncertainty_nu), 
           col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n")
  }
}
mtext("Optimal Quantity (q) given Price Floor (theta)", outer = TRUE, cex = 1.2, font = 2)
dev.off()

# CONFIDENCE ====================================================================

## ==============================================================================
# PHASE 4.1: MC of INVESTMENT PROCESS (Memory Optimized & Sequential) ===========
## ==============================================================================

Npaths <- 10000
set.seed(42) 
Z_std <- matrix(rnorm(Npaths * Nt), nrow = Npaths) 
W <- params$sigma * sqrt(dt) * Z_std

simulate_paths <- function(res_obj, params, Npaths, Nt, dx, dt, W, X0) {
  S_current <- rep(X0, Npaths)
  
  S_avg     <- numeric(Nt + 1)
  S_lower   <- numeric(Nt + 1)
  S_upper   <- numeric(Nt + 1)
  S_avg[1]   <- X0
  S_lower[1] <- X0
  S_upper[1] <- X0
  
  gamma_avg   <- numeric(Nt)
  gamma_lower <- numeric(Nt)
  gamma_upper <- numeric(Nt)
  
  q_avg     <- numeric(Nt)
  q_lower   <- numeric(Nt)
  q_upper   <- numeric(Nt)
  
  theta_avg   <- numeric(Nt)
  theta_lower <- numeric(Nt)
  theta_upper <- numeric(Nt)
  
  for (t_step in 1:Nt) {
    idx <- pmax(1, pmin(Nx - 2, floor((S_current - X_min) / dx) + 1))
    
    gamma_step <- res_obj$gamma[[t_step]][idx]
    q_step     <- res_obj$q[idx]
    theta_step <- res_obj$theta[idx]
    
    gamma_avg[t_step]   <- mean(gamma_step)
    gamma_lower[t_step] <- quantile(gamma_step, 0.05, names = FALSE)
    gamma_upper[t_step] <- quantile(gamma_step, 0.95, names = FALSE)
    
    q_avg[t_step]       <- mean(q_step)
    q_lower[t_step]     <- quantile(q_step, 0.05, names = FALSE)
    q_upper[t_step]     <- quantile(q_step, 0.95, names = FALSE)
    
    theta_avg[t_step]   <- mean(theta_step)
    theta_lower[t_step] <- quantile(theta_step, 0.05, names = FALSE)
    theta_upper[t_step] <- quantile(theta_step, 0.95, names = FALSE)
    
    S_current <- S_current + (gamma_step - params$delta * S_current) * dt + W[, t_step]
    
    S_avg[t_step+1]   <- mean(S_current)
    S_lower[t_step+1] <- quantile(S_current, 0.05, names = FALSE)
    S_upper[t_step+1] <- quantile(S_current, 0.95, names = FALSE)
  }
  
  list(
    S_avg = S_avg, S_lower = S_lower, S_upper = S_upper,
    gamma_avg = gamma_avg, gamma_lower = gamma_lower, gamma_upper = gamma_upper,
    q_avg = q_avg, q_lower = q_lower, q_upper = q_upper,
    theta_avg = theta_avg, theta_lower = theta_lower, theta_upper = theta_upper
  )
}

X0_vals <- c(params$x_bar, params$x_bar - 1 , params$x_bar - 2)
X0_titles <- c("x_bar", "x_bar - 1", "x_bar - 2")
X0_tags <- c("x_bar", "x_bar_minus_1", "x_bar_minus_2")

cat("Running Base Monte Carlo paths sequentially...\n")
mc_results <- list()
for (v in 1:length(X0_vals)) {
  X0_current <- X0_vals[v]
  mc_results[[v]] <- list(
    low  = lapply(1:4, function(k) simulate_paths(res_low_unc[[k]], params, Npaths, Nt, dx, dt, W, X0_current)),
    mid  = lapply(1:4, function(k) simulate_paths(res_mid_unc[[k]], params, Npaths, Nt, dx, dt, W, X0_current)),
    high = lapply(1:4, function(k) simulate_paths(res_high_unc[[k]], params, Npaths, Nt, dx, dt, W, X0_current))
  )
}

## ==============================================================================
# PHASE 4.2.1: Process plots combined (Rows = X0, Columns = theta_bar) ==========
## ==============================================================================
time_seq <- (0:Nt) * dt
time_seq_Nt <- (1:Nt) * dt

# Rows are X0, Columns are theta_bar
n_rows <- length(X0_vals)
n_cols <- length(th_vals)

### Panel Plot for Avg. Investment Paths ====
png("avg_investment_paths_combined.png", width = a4_width, height = a4_height, units = "in", res = resolution)
par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))

col_ylim <- list()
for (k in 1:n_cols) {
  row_vals <- c()
  for(v in 1:n_rows) {
    row_vals <- c(row_vals, 
                  mc_results[[v]]$low[[k]]$S_lower, mc_results[[v]]$low[[k]]$S_upper,
                  mc_results[[v]]$mid[[k]]$S_lower, mc_results[[v]]$mid[[k]]$S_upper,
                  mc_results[[v]]$high[[k]]$S_lower, mc_results[[v]]$high[[k]]$S_upper)
  }
  col_ylim[[k]] <- range(row_vals, na.rm = TRUE)
}

for (v in 1:n_rows) {
  for (k in 1:n_cols) {
    mc_low  <- mc_results[[v]]$low
    mc_mid  <- mc_results[[v]]$mid
    mc_high <- mc_results[[v]]$high
    
    plot(time_seq, mc_low[[k]]$S_avg, type = "n", 
         xlab = "Time (Years)", ylab = "Investment Value (X)",
         ylim = col_ylim[[k]])
    
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    
    polygon(c(time_seq, rev(time_seq)), c(mc_high[[k]]$S_lower, rev(mc_high[[k]]$S_upper)), col = adjustcolor("orange", alpha.f = 0.2), border = NA)
    polygon(c(time_seq, rev(time_seq)), c(mc_mid[[k]]$S_lower, rev(mc_mid[[k]]$S_upper)), col = adjustcolor("black", alpha.f = 0.2), border = NA)
    polygon(c(time_seq, rev(time_seq)), c(mc_low[[k]]$S_lower, rev(mc_low[[k]]$S_upper)), col = adjustcolor("darkgreen", alpha.f = 0.2), border = NA)
    
    lines(time_seq, mc_high[[k]]$S_avg, col = "orange", lwd = 2)
    lines(time_seq, mc_mid[[k]]$S_avg, col = "black", lwd = 2)
    lines(time_seq, mc_low[[k]]$S_avg, col = "darkgreen", lwd = 2)
    
    if (v == 1) { 
      mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    if (k == 1) { 
      mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}
mtext("Avg. of Investment Value (X) Paths over Time", outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                            paste0("Mid Unc (", uncertainty_nu[2], ")"),
                            paste0("Low Unc (", uncertainty_nu[3], ")")), 
       col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
dev.off()


### Panel Plot for Avg. Price Floor (theta) ====
png("avg_subsidies_paths_combined.png", width = a4_width, height = a4_height, units = "in", res = resolution)
par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))

col_ylim <- list()
for (k in 1:n_cols) {
  row_vals <- c()
  for(v in 1:n_rows) {
    row_vals <- c(row_vals, 
                  mc_results[[v]]$low[[k]]$theta_lower, mc_results[[v]]$low[[k]]$theta_upper,
                  mc_results[[v]]$mid[[k]]$theta_lower, mc_results[[v]]$mid[[k]]$theta_upper,
                  mc_results[[v]]$high[[k]]$theta_lower, mc_results[[v]]$high[[k]]$theta_upper)
  }
  col_ylim[[k]] <- range(row_vals, na.rm = TRUE)
}

for (v in 1:n_rows) {
  for (k in 1:n_cols) {
    mc_low  <- mc_results[[v]]$low
    mc_mid  <- mc_results[[v]]$mid
    mc_high <- mc_results[[v]]$high
    
    plot(time_seq_Nt, mc_low[[k]]$theta_avg, type = "n", 
         xlab = "Time (Years)", ylab = "Price Floor (theta)",
         ylim = col_ylim[[k]])
    
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    abline(h = th_vals[k], lty = 2, col = "darkgrey")
    abline(h = params$py, lty = 2, col = "darkgrey")
    
    polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_high[[k]]$theta_lower, rev(mc_high[[k]]$theta_upper)), col = adjustcolor("orange", alpha.f = 0.2), border = NA)
    polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_mid[[k]]$theta_lower, rev(mc_mid[[k]]$theta_upper)), col = adjustcolor("black", alpha.f = 0.2), border = NA)
    polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_low[[k]]$theta_lower, rev(mc_low[[k]]$theta_upper)), col = adjustcolor("darkgreen", alpha.f = 0.2), border = NA)
    
    lines(time_seq_Nt, mc_high[[k]]$theta_avg, col = "orange", lwd = 2)
    lines(time_seq_Nt, mc_mid[[k]]$theta_avg, col = "black", lwd = 2)
    lines(time_seq_Nt, mc_low[[k]]$theta_avg, col = "darkgreen", lwd = 2)
    
    if (v == 1) {
      mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    if (k == 1) {
      mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}
mtext("Avg. of Price Floor (theta) Paths over Time", outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                            paste0("Mid Unc (", uncertainty_nu[2], ")"),
                            paste0("Low Unc (", uncertainty_nu[3], ")")), 
       col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
dev.off()


### Panel Plot for Avg. Quantities (q) ====
png("avg_quantities_paths_combined.png", width = a4_width, height = a4_height, units = "in", res = resolution)
par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))

col_ylim <- list()
for (k in 1:n_cols) {
  row_vals <- c()
  for(v in 1:n_rows) {
    row_vals <- c(row_vals, 
                  mc_results[[v]]$low[[k]]$q_lower, mc_results[[v]]$low[[k]]$q_upper,
                  mc_results[[v]]$mid[[k]]$q_lower, mc_results[[v]]$mid[[k]]$q_upper,
                  mc_results[[v]]$high[[k]]$q_lower, mc_results[[v]]$high[[k]]$q_upper)
  }
  col_ylim[[k]] <- range(row_vals, na.rm = TRUE)
}

for (v in 1:n_rows) {
  for (k in 1:n_cols) {
    mc_low  <- mc_results[[v]]$low
    mc_mid  <- mc_results[[v]]$mid
    mc_high <- mc_results[[v]]$high
    
    plot(time_seq_Nt, mc_low[[k]]$q_avg, type = "n", 
         xlab = "Time (Years)", ylab = "Quantity (q)",
         ylim = col_ylim[[k]])
    
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    abline(h = params$q_max, lty = 2, col = "darkgrey")
    abline(h = 1.96, lty = 2, col = "darkgrey")
    
    polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_high[[k]]$q_lower, rev(mc_high[[k]]$q_upper)), col = adjustcolor("orange", alpha.f = 0.2), border = NA)
    polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_mid[[k]]$q_lower, rev(mc_mid[[k]]$q_upper)), col = adjustcolor("black", alpha.f = 0.2), border = NA)
    polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_low[[k]]$q_lower, rev(mc_low[[k]]$q_upper)), col = adjustcolor("darkgreen", alpha.f = 0.2), border = NA)
    
    lines(time_seq_Nt, mc_high[[k]]$q_avg, col = "orange", lwd = 2)
    lines(time_seq_Nt, mc_mid[[k]]$q_avg, col = "black", lwd = 2)
    lines(time_seq_Nt, mc_low[[k]]$q_avg, col = "darkgreen", lwd = 2)
    
    if (v == 1) {
      mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    if (k == 1) {
      mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}
mtext("Avg. of Quantity (q) Paths over Time", outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                            paste0("Mid Unc (", uncertainty_nu[2], ")"),
                            paste0("Low Unc (", uncertainty_nu[3], ")")), 
       col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
dev.off()


### Panel Plot for Avg. Investment Rate (gamma) ====
png("avg_investment_rates_paths_combined.png", width = a4_width, height = a4_height, units = "in", res = resolution)
par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))

col_ylim <- list()
for (k in 1:n_cols) {
  row_vals <- c()
  for(v in 1:n_rows) {
    row_vals <- c(row_vals, 
                  mc_results[[v]]$low[[k]]$gamma_lower, mc_results[[v]]$low[[k]]$gamma_upper,
                  mc_results[[v]]$mid[[k]]$gamma_lower, mc_results[[v]]$mid[[k]]$gamma_upper,
                  mc_results[[v]]$high[[k]]$gamma_lower, mc_results[[v]]$high[[k]]$gamma_upper)
  }
  col_ylim[[k]] <- range(row_vals, na.rm = TRUE)
}

for (v in 1:n_rows) {
  for (k in 1:n_cols) {
    mc_low  <- mc_results[[v]]$low
    mc_mid  <- mc_results[[v]]$mid
    mc_high <- mc_results[[v]]$high
    
    plot(time_seq_Nt, mc_low[[k]]$gamma_avg, type = "n", 
         xlab = "Time (Years)", ylab = "Investment Rate (gamma)",
         ylim = col_ylim[[k]])
    
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    abline(h = 0, lty = 2, col = "darkgrey")
    
    polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_high[[k]]$gamma_lower, rev(mc_high[[k]]$gamma_upper)), col = adjustcolor("orange", alpha.f = 0.2), border = NA)
    polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_mid[[k]]$gamma_lower, rev(mc_mid[[k]]$gamma_upper)), col = adjustcolor("black", alpha.f = 0.2), border = NA)
    polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_low[[k]]$gamma_lower, rev(mc_low[[k]]$gamma_upper)), col = adjustcolor("darkgreen", alpha.f = 0.2), border = NA)
    
    lines(time_seq_Nt, mc_high[[k]]$gamma_avg, col = "orange", lwd = 2)
    lines(time_seq_Nt, mc_mid[[k]]$gamma_avg, col = "black", lwd = 2)
    lines(time_seq_Nt, mc_low[[k]]$gamma_avg, col = "darkgreen", lwd = 2)
    
    if (v == 1) {
      mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    if (k == 1) {
      mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}
mtext("Avg. of Investment Rate (gamma) Paths over Time", outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                            paste0("Mid Unc (", uncertainty_nu[2], ")"),
                            paste0("Low Unc (", uncertainty_nu[3], ")")), 
       col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
dev.off()


## ==============================================================================
# PHASE 4.2.2: Process plots for Avg. Quantities (q) Breakdown ==================
## ==============================================================================

### Panel Plot for Avg. Quantities (q) Breakdown ====
png("avg_quantities_areas_breakdown_combined.png", width = a4_width, height = a4_height, units = "in", res = resolution)

n_rows_nu <- 3 
par(mfrow = c(n_rows_nu, n_cols), oma = c(4, 6, 4, 1), mar = c(4, 4, 3, 1))

X0_state <- 1 # We plot for the first initial state X0

# Pre-calculate global ylim incorporating the upper and lower confidence intervals
col_ylim <- list()
for (k in 1:n_cols) {
  col_vals <- c(mc_results[[X0_state]]$low[[k]]$q_lower, mc_results[[X0_state]]$low[[k]]$q_upper,
                mc_results[[X0_state]]$mid[[k]]$q_lower, mc_results[[X0_state]]$mid[[k]]$q_upper,
                mc_results[[X0_state]]$high[[k]]$q_lower, mc_results[[X0_state]]$high[[k]]$q_upper)
  
  lims <- range(col_vals, na.rm = TRUE)
  col_ylim[[k]] <- c(0, max(lims[2], params$q_max)) 
}

for (r in 1:n_rows_nu) {
  for (k in 1:n_cols) {
    
    if (r == 1) curr_data <- mc_results[[X0_state]]$high[[k]]
    if (r == 2) curr_data <- mc_results[[X0_state]]$mid[[k]]
    if (r == 3) curr_data <- mc_results[[X0_state]]$low[[k]]
    
    q_val <- curr_data$q_avg
    q_lower <- curr_data$q_lower
    q_upper <- curr_data$q_upper
    
    # CALCULATE Pg_x
    x_val <- curr_data$S_avg[-1] 
    Pg_x <- params$p_g * pmax(x_val - params$x_bar, 0)
    
    # Prevent overlap bugs
    q_val <- pmax(q_val, Pg_x)
    
    # Draw Plot
    plot(time_seq_Nt, q_val, type = "n", 
         xlab = "Time (Years)", ylab = "Quantity (q)",
         ylim = col_ylim[[k]])
    
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    abline(h = params$q_max, lty = 2, col = "darkgrey")
    abline(h = 1.96, lty = 2, col = "darkgrey")
    
    # Draw Green Area
    polygon(x = c(time_seq_Nt, rev(time_seq_Nt)), 
            y = c(rep(0, length(time_seq_Nt)), rev(Pg_x)), 
            col = "lightgreen", border = NA)
    
    # Draw Brown Area
    polygon(x = c(time_seq_Nt, rev(time_seq_Nt)), 
            y = c(Pg_x, rev(q_val)), 
            col = "tan", border = NA)
    
    # Add Borders
    lines(time_seq_Nt, Pg_x, col = "darkgreen", lwd = 2)
    lines(time_seq_Nt, q_val, col = "saddlebrown", lwd = 2)
    
    # Add dashed line for brown overproduction
    lines(time_seq_Nt, Pg_x + 1.96, col = "saddlebrown", lty = 2, lwd = 1)
    
    # Margins and Titles
    if (r == 1) {
      mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    if (k == 1) {
      nu_label <- switch(r, 
                         paste0("High Unc (", uncertainty_nu[1], ")"),
                         paste0("Mid Unc (", uncertainty_nu[2], ")"),
                         paste0("Low Unc (", uncertainty_nu[3], ")"))
      mtext(nu_label, side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}

# Main Outer Title
mtext(paste0("Avg. of Quantity (q) Paths over Time with Initial State X0 = ", X0_titles[X0_state]), 
      outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)

# Outer Legend
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", 
       legend = c(expression("P"["g"]*"(x) (Green Area)"), 
                  expression("q - P"["g"]*"(x) (Brown Area)"), 
                  expression("P"["g"]*"(x) + 1.96 (Brown Production Lower Bound)")), 
       fill = c("lightgreen", "tan", NA), 
       border = c("darkgreen", "saddlebrown", NA),
       col = c(NA, NA, "saddlebrown"), 
       lty = c(0, 0, 2), # 0 = no line for boxes, 2 = dashed line for the 3rd item
       lwd = c(0, 0, 2), 
       bty = "n", horiz = TRUE, cex = 1.2)
dev.off()

## ==============================================================================
# PHASE 5: Panel Value plots ==================================================
## ==============================================================================

color_vec <- c("darkblue", "blue", "lightblue", "lightgrey", "grey")
step_values <- c(0.25, 0.5, 0.75, 1) * Nt

# Create a 3x4 grid: 3 rows (nu/uncertainty), 4 columns (theta_bar)
png("game_value_combined.png", width = a4_width, height = a4_height, units = "in", res = resolution)
par(mfrow = c(3, 4), oma = c(5, 6, 4, 1), mar = c(4, 4, 3, 1))

# Group lists and labels to easily loop through them
res_lists <- list(res_high_unc, res_mid_unc, res_low_unc)
nu_labels <- c(paste0("High Unc (", uncertainty_nu[1], ")"),
               paste0("Mid Unc (", uncertainty_nu[2], ")"),
               paste0("Low Unc (", uncertainty_nu[3], ")"))

# Pre-calculate global ylim per column so the Y-axis is scaled equally across rows
col_ylim <- list()
lims <- c()
for (k in 1:4) {
  for (r in 1:3) {
    U_val <- matrix(unlist(res_lists[[r]][[k]]$value), nrow = Nt + 1, ncol = Nx, byrow = TRUE)
    lims <- c(lims, range(U_val, na.rm = TRUE))
  }
}
col_ylim[[1]] <- range(lims)

# Generate the 3x4 plots
for (r in 1:3) {
  for (k in 1:4) {
    
    # Extract the current matrix
    U_val <- matrix(unlist(res_lists[[r]][[k]]$value), nrow = Nt + 1, ncol = Nx, byrow = TRUE)
    
    # Base Plot
    plot(X, U_val[1,], type = "l", col = color_vec[1], lwd = 2,
         xlab = "Investment X", ylab = "Game Value",
         ylim = col_ylim[[1]]) 
    
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    
    # Add step lines
    for (i in 1:4) {
      lines(X, U_val[step_values[i], ], col = color_vec[i+1], lwd = 2)
    }
    
    # Add Column Titles (theta_bar) on the top row
    if (r == 1) {
      mtext(paste("theta_bar =", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    
    # Add Row Titles (nu) on the first column
    if (k == 1) {
      mtext(nu_labels[r], side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}

# Main Outer Title
mtext("Game Value over X by Uncertainty and Theta Bar", outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)

# Outer Global Legend at the bottom
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", 
       legend = c("t=0", paste("t=", step_values*10/Nt)), 
       col = color_vec, lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)

dev.off()




## ==============================================================================
# PHASE 6: Benchmarking of Code =================================================
## ==============================================================================

# Group lists and labels to easily loop through them
res_lists <- list(res_high_unc, res_mid_unc, res_low_unc)
nu_labels <- c(paste0("High Unc (", uncertainty_nu[1], ")"),
               paste0("Mid Unc (", uncertainty_nu[2], ")"),
               paste0("Low Unc (", uncertainty_nu[3], ")"))

# Line colors for the different uncertainty levels (matching your original function calls)
row_colors <- c("orange", "black", "darkgreen")

## ==============================================================================
# PHASE 6.1: First Derivative ===================================================
## ==============================================================================

png("marginal_value_benchmark_combined.png", width = a4_width, height = a4_height, units = "in", res = resolution)
par(mfrow = c(3, 4), oma = c(5, 6, 4, 1), mar = c(4, 4, 3, 1))

# Pre-calculate global ylim so the Y-axis is scaled equally across all 12 panels
lims <- c()
for (r in 1:3) {
  for (k in 1:4) {
    Ux <- diff(res_lists[[r]][[k]]$value[[1]]) / dx
    lims <- c(lims, range(Ux, na.rm = TRUE))
  }
}
global_ylim <- range(lims)

# Generate the 3x4 plots
for (r in 1:3) {
  for (k in 1:4) {
    
    # Extract and calculate the derivative for the current matrix
    Ux <- diff(res_lists[[r]][[k]]$value[[1]]) / dx
    
    # Base Plot
    plot(X[-Nx], Ux, type = "l", col = row_colors[r], lwd = 2,
         xlab = "Investment Space (X)", ylab = "U_x",
         ylim = global_ylim) 
    
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    
    # Add benchmark vertical line
    abline(v = params$x_bar, lty = 3, col = "red")
    
    # Add Column Titles (theta_bar) on the top row
    if (r == 1) {
      mtext(paste("theta_bar =", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    
    # Add Row Titles (nu) on the first column
    if (k == 1) {
      mtext(nu_labels[r], side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}

# Main Outer Title
mtext("First Derivative (Marginal Value) Benchmark by Uncertainty and Theta Bar", outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)

# Outer Global Legend at the bottom
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", 
       legend = c("High Unc", "Mid Unc", "Low Unc", "x_bar Threshold"), 
       col = c(row_colors[1], row_colors[2], row_colors[3], "red"), 
       lwd = c(2, 2, 2, 1), lty = c(1, 1, 1, 3), 
       bty = "n", horiz = TRUE, cex = 1.2)

dev.off()

## ==============================================================================
# PHASE 6.2: Howard Policy Iteration Diagnostics ================================
## ==============================================================================

# Panel Plot for Howard Error (h_error)
png("howard_error_combined.png", width = a4_width, height = a4_height, units = "in", res = resolution)
par(mfrow = c(3, 4), oma = c(5, 6, 4, 1), mar = c(4, 4, 3, 1))

# Pre-calculate global ylim (safeguarded against flat lines by forcing a small max)
lims_err <- c()
for (r in 1:3) {
  for (k in 1:4) {
    lims_err <- c(lims_err, range(res_lists[[r]][[k]]$h_error, na.rm = TRUE))
  }
}
global_ylim_err <- c(min(0, min(lims_err, na.rm = TRUE)), 
                     max(1e-8, max(lims_err, na.rm = TRUE)))

# Generate the 3x4 plots
for (r in 1:3) {
  for (k in 1:4) {
    
    err_val <- res_lists[[r]][[k]]$h_error
    steps <- 1:length(err_val)
    
    # Base Plot
    plot(steps, err_val, type = "l", col = row_colors[r], lwd = 1,
         xlab = "Time Step Count in dt", ylab = "Howard Error",
         ylim = global_ylim_err) 
    
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    
    # Add Column Titles (theta_bar) on the top row
    if (r == 1) {
      mtext(paste("theta_bar =", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    
    # Add Row Titles (nu) on the first column
    if (k == 1) {
      mtext(nu_labels[r], side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}

# Main Outer Title
mtext("Howard Policy Iteration Error by Uncertainty and Theta Bar", outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)

# Outer Global Legend at the bottom
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", 
       legend = c("High Unc", "Mid Unc", "Low Unc"), 
       col = row_colors, 
       lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)

dev.off()


# Panel Plot for Howard Iterations (h_iter)
png("howard_iter_combined.png", width = a4_width, height = a4_height, units = "in", res = resolution)
par(mfrow = c(3, 4), oma = c(5, 6, 4, 1), mar = c(4, 4, 3, 1))

# Pre-calculate global ylim (safeguarded so perfectly flat iterations won't crash ylim)
lims_iter <- c()
for (r in 1:3) {
  for (k in 1:4) {
    lims_iter <- c(lims_iter, range(res_lists[[r]][[k]]$h_iter, na.rm = TRUE))
  }
}
global_ylim_iter <- c(min(0, min(lims_iter, na.rm = TRUE)), 
                      max(1, max(lims_iter, na.rm = TRUE)))

# Generate the 3x4 plots
for (r in 1:3) {
  for (k in 1:4) {
    
    iter_val <- res_lists[[r]][[k]]$h_iter
    steps <- 1:length(iter_val)
    
    # Base Plot
    plot(steps, iter_val, type = "p", col = row_colors[r], lwd = 0.5,
         xlab = "Time Step Count in dt", ylab = "Iterations",
         ylim = global_ylim_iter) 
    
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    
    # Add Column Titles (theta_bar) on the top row
    if (r == 1) {
      mtext(paste("theta_bar =", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    
    # Add Row Titles (nu) on the first column
    if (k == 1) {
      mtext(nu_labels[r], side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}

# Main Outer Title
mtext("Howard Policy Iterations by Uncertainty and Theta Bar", outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)

# Outer Global Legend at the bottom
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", 
       legend = c("High Unc", "Mid Unc", "Low Unc"), 
       col = row_colors, 
       lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)

dev.off()


## ==============================================================================
# PHASE 7: SENSITIVITY ANALYSIS (DELTA & SIGMA) (Sequential) ====================
## ==============================================================================

base_nu <- uncertainty_nu[2]
delta_vals <- c(0.002, params$delta, 0.2) 
sigma_vals <- c(0.002, params$sigma, 2)   
kappa_vals <- c(0.05, params$kappa, 5)   

cat("Calculating Delta Sensitivity (Depreciation) sequentially...\n")
res_delta_low <- lapply(th_vals, function(th) {
  p <- params; p$delta <- delta_vals[1]
  simulate_fixed_p(th, base_nu, p)
})
res_delta_high <- lapply(th_vals, function(th) {
  p <- params; p$delta <- delta_vals[3]
  simulate_fixed_p(th, base_nu, p)
})

cat("Calculating Sigma Sensitivity (Volatility) sequentially...\n")
res_sigma_low <- lapply(th_vals, function(th) {
  p <- params; p$sigma <- sigma_vals[1]
  simulate_fixed_p(th, base_nu, p)
})
res_sigma_high <- lapply(th_vals, function(th) {
  p <- params; p$sigma <- sigma_vals[3]
  simulate_fixed_p(th, base_nu, p)
})

cat("Calculating kappa Sensitivity (Transaction Costs) sequentially...\n")
res_kappa_low <- lapply(th_vals, function(th) {
  p <- params; p$kappa <- kappa_vals[1]
  simulate_fixed_p(th, base_nu, p)
})
res_kappa_high <- lapply(th_vals, function(th) {
  p <- params; p$kappa <- sigma_vals[3]
  simulate_fixed_p(th, base_nu, p)
})

cat("Running MC Simulations for Sensitivity sequentially...\n")
W_sigma_low  <- sigma_vals[1] * sqrt(dt) * Z_std
W_sigma_high <- sigma_vals[3] * sqrt(dt) * Z_std

mc_sens_results <- list()
for(v in 1:length(X0_vals)) {
  X0_current <- X0_vals[v]
  
  mc_sens_results[[v]] <- list(
    delta_low  = lapply(1:4, function(k) { p <- params; p$delta <- delta_vals[1]; simulate_paths(res_delta_low[[k]], p, Npaths, Nt, dx, dt, W, X0_current) }),
    delta_high = lapply(1:4, function(k) { p <- params; p$delta <- delta_vals[3]; simulate_paths(res_delta_high[[k]], p, Npaths, Nt, dx, dt, W, X0_current) }),
    sigma_low  = lapply(1:4, function(k) { p <- params; p$sigma <- sigma_vals[1]; simulate_paths(res_sigma_low[[k]], p, Npaths, Nt, dx, dt, W_sigma_low, X0_current) }),
    sigma_high = lapply(1:4, function(k) { p <- params; p$sigma <- sigma_vals[3]; simulate_paths(res_sigma_high[[k]], p, Npaths, Nt, dx, dt, W_sigma_high, X0_current) }),
    kappa_low  = lapply(1:4, function(k) { p <- params; p$kappa <- kappa_vals[1]; simulate_paths(res_kappa_low[[k]], p, Npaths, Nt, dx, dt, W, X0_current) }),
    kappa_high = lapply(1:4, function(k) { p <- params; p$kappa <- kappa_vals[3]; simulate_paths(res_kappa_high[[k]], p, Npaths, Nt, dx, dt, W, X0_current) })
  )
}


## ==============================================================================
# PHASE 7.1: SENSITIVITY ANALYSIS (DELTA & SIGMA) Combined Plot ================
## ==============================================================================

### COMBINED PLOTS: DELTA & SIGMA SENSITIVITY Combined Panel Plots ====
time_seq <- (0:Nt) * dt

n_rows <- length(X0_vals)
n_cols <- length(th_vals)

#### DELTA SENSITIVITY PLOT ====
png("sens_delta_avg_investment_combined_nu_1.png", width = a4_width, height = a4_height, units = "in", res = resolution)
par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))

# Calculate limits by row (for each X0 across all theta_bar columns) incorporating CIs
row_ylim <- list()
for (v in 1:n_rows) {
  vals <- c()
  for (k in 1:n_cols) {
    vals <- c(vals, 
              mc_sens_results[[v]]$delta_low[[k]]$S_lower, mc_sens_results[[v]]$delta_low[[k]]$S_upper,
              mc_results[[v]]$mid[[k]]$S_lower, mc_results[[v]]$mid[[k]]$S_upper,
              mc_sens_results[[v]]$delta_high[[k]]$S_lower, mc_sens_results[[v]]$delta_high[[k]]$S_upper)
  }
  row_ylim[[v]] <- range(vals, na.rm = TRUE)
}

for (v in 1:n_rows) {
  for (k in 1:n_cols) {
    mc_delta_low  <- mc_sens_results[[v]]$delta_low
    mc_delta_high <- mc_sens_results[[v]]$delta_high
    mc_mid        <- mc_results[[v]]$mid
    
    # Initialize empty plot
    plot(time_seq, mc_delta_high[[k]]$S_avg, type = "n", 
         xlab = "Time (Years)", ylab = "Investment Value (X)",
         ylim = row_ylim[[v]])
    
    # Add grid
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    
    # Draw Confidence Intervals (Shading)
    polygon(c(time_seq, rev(time_seq)), c(mc_delta_high[[k]]$S_lower, rev(mc_delta_high[[k]]$S_upper)), col = adjustcolor("red", alpha.f = 0.2), border = NA)
    polygon(c(time_seq, rev(time_seq)), c(mc_mid[[k]]$S_lower, rev(mc_mid[[k]]$S_upper)), col = adjustcolor("black", alpha.f = 0.2), border = NA)
    polygon(c(time_seq, rev(time_seq)), c(mc_delta_low[[k]]$S_lower, rev(mc_delta_low[[k]]$S_upper)), col = adjustcolor("blue", alpha.f = 0.2), border = NA)
    
    # Draw lines
    lines(time_seq, mc_delta_high[[k]]$S_avg, col = "red", lwd = 2)
    lines(time_seq, mc_mid[[k]]$S_avg, col = "black", lwd = 2) 
    lines(time_seq, mc_delta_low[[k]]$S_avg, col = "blue", lwd = 2)
    
    # Headers
    if (v == 1) {
      mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    if (k == 1) {
      mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}
mtext("Sensitivity to Depreciation (Delta): Avg. of Investment (X) Paths over Time at nu = 1", outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", legend = c(paste0("High Delta (", delta_vals[3], ")"), 
                            paste0("Base Delta (", delta_vals[2], ")"),
                            paste0("Low Delta (", delta_vals[1], ")")), 
       col = c("red", "black", "blue"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
dev.off()

#### SIGMA SENSITIVITY PLOT ====
png("sens_sigma_avg_investment_combined_nu_1.png", width = a4_width, height = a4_height, units = "in", res = resolution)
par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))

# Calculate limits by row (for each X0 across all theta_bar columns) incorporating CIs
row_ylim <- list()
for (v in 1:n_rows) {
  vals <- c()
  for (k in 1:n_cols) {
    vals <- c(vals, 
              mc_sens_results[[v]]$sigma_low[[k]]$S_lower, mc_sens_results[[v]]$sigma_low[[k]]$S_upper,
              mc_results[[v]]$mid[[k]]$S_lower, mc_results[[v]]$mid[[k]]$S_upper,
              mc_sens_results[[v]]$sigma_high[[k]]$S_lower, mc_sens_results[[v]]$sigma_high[[k]]$S_upper)
  }
  row_ylim[[v]] <- range(vals, na.rm = TRUE)
}

for (v in 1:n_rows) {
  for (k in 1:n_cols) {
    mc_sigma_low  <- mc_sens_results[[v]]$sigma_low
    mc_sigma_high <- mc_sens_results[[v]]$sigma_high
    mc_mid        <- mc_results[[v]]$mid
    
    # Initialize empty plot
    plot(time_seq, mc_sigma_high[[k]]$S_avg, type = "n", 
         xlab = "Time (Years)", ylab = "Investment Value (X)",
         ylim = row_ylim[[v]])
    
    # Add grid
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    
    # Draw Confidence Intervals (Shading)
    polygon(c(time_seq, rev(time_seq)), c(mc_sigma_high[[k]]$S_lower, rev(mc_sigma_high[[k]]$S_upper)), col = adjustcolor("red", alpha.f = 0.2), border = NA)
    polygon(c(time_seq, rev(time_seq)), c(mc_mid[[k]]$S_lower, rev(mc_mid[[k]]$S_upper)), col = adjustcolor("black", alpha.f = 0.2), border = NA)
    polygon(c(time_seq, rev(time_seq)), c(mc_sigma_low[[k]]$S_lower, rev(mc_sigma_low[[k]]$S_upper)), col = adjustcolor("blue", alpha.f = 0.2), border = NA)
    
    # Draw lines
    lines(time_seq, mc_sigma_high[[k]]$S_avg, col = "red", lwd = 2)
    lines(time_seq, mc_mid[[k]]$S_avg, col = "black", lwd = 2) 
    lines(time_seq, mc_sigma_low[[k]]$S_avg, col = "blue", lwd = 2)
    
    # Headers
    if (v == 1) {
      mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    if (k == 1) {
      mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}
mtext("Sensitivity to Volatility (Sigma):  Avg. of Investment (X) Paths over Time at nu = 1", outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", legend = c(paste0("High Sigma (", sigma_vals[3], ")"), 
                            paste0("Base Sigma (", sigma_vals[2], ")"),
                            paste0("Low Sigma (", sigma_vals[1], ")")), 
       col = c("red", "black", "blue"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
dev.off()

par(mfrow = c(1, 1))


#### kappa SENSITIVITY PLOT ====
png("sens_kappa_avg_investment_combined_nu_1.png", width = a4_width, height = a4_height, units = "in", res = resolution)
par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))

# Calculate limits by row (for each X0 across all theta_bar columns) incorporating CIs
row_ylim <- list()
for (v in 1:n_rows) {
  vals <- c()
  for (k in 1:n_cols) {
    vals <- c(vals, 
              mc_sens_results[[v]]$kappa_low[[k]]$S_lower, mc_sens_results[[v]]$kappa_low[[k]]$S_upper,
              mc_results[[v]]$mid[[k]]$S_lower, mc_results[[v]]$mid[[k]]$S_upper,
              mc_sens_results[[v]]$kappa_high[[k]]$S_lower, mc_sens_results[[v]]$kappa_high[[k]]$S_upper)
  }
  row_ylim[[v]] <- range(vals, na.rm = TRUE)
}

for (v in 1:n_rows) {
  for (k in 1:n_cols) {
    mc_kappa_low  <- mc_sens_results[[v]]$kappa_low
    mc_kappa_high <- mc_sens_results[[v]]$kappa_high
    mc_mid        <- mc_results[[v]]$mid
    
    # Initialize empty plot
    plot(time_seq, mc_kappa_high[[k]]$S_avg, type = "n", 
         xlab = "Time (Years)", ylab = "Investment Value (X)",
         ylim = row_ylim[[v]])
    
    # Add grid
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    
    # Draw Confidence Intervals (Shading)
    polygon(c(time_seq, rev(time_seq)), c(mc_kappa_high[[k]]$S_lower, rev(mc_kappa_high[[k]]$S_upper)), col = adjustcolor("red", alpha.f = 0.2), border = NA)
    polygon(c(time_seq, rev(time_seq)), c(mc_mid[[k]]$S_lower, rev(mc_mid[[k]]$S_upper)), col = adjustcolor("black", alpha.f = 0.2), border = NA)
    polygon(c(time_seq, rev(time_seq)), c(mc_kappa_low[[k]]$S_lower, rev(mc_kappa_low[[k]]$S_upper)), col = adjustcolor("blue", alpha.f = 0.2), border = NA)
    
    # Draw lines
    lines(time_seq, mc_kappa_high[[k]]$S_avg, col = "red", lwd = 2)
    lines(time_seq, mc_mid[[k]]$S_avg, col = "black", lwd = 2) 
    lines(time_seq, mc_kappa_low[[k]]$S_avg, col = "blue", lwd = 2)
    
    # Headers
    if (v == 1) {
      mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    if (k == 1) {
      mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}
mtext("Sensitivity to Transaction Costs (kappa):  Avg. of Investment (X) Paths over Time at nu = 1", outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", legend = c(paste0("High kappa (", kappa_vals[3], ")"), 
                            paste0("Base kappa (", kappa_vals[2], ")"),
                            paste0("Low kappa (", kappa_vals[1], ")")), 
       col = c("red", "black", "blue"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
dev.off()

par(mfrow = c(1, 1))



## =============================================================================
# PHASE 8. EVALUATE again for different terminal discount ======================
## =============================================================================

prev_disc <- params$terminal_discount
prev_r <- params$r

cat("Calculating Base Scenarios (Varying Nu) sequentially...\n")
params$terminal_discount <- 0.9
res_high_unc <- lapply(th_vals, function(th) simulate_fixed_p(th, uncertainty_nu[1], params))
res_mid_unc  <- lapply(th_vals, function(th) simulate_fixed_p(th, uncertainty_nu[2], params))
res_low_unc  <- lapply(th_vals, function(th) simulate_fixed_p(th, uncertainty_nu[3], params))


### MC of INVESTMENT PROCESS (Memory Optimized & Sequential) ====


cat("Running Base Monte Carlo paths sequentially...\n")
mc_results <- list()
for (v in 1:length(X0_vals)) {
  X0_current <- X0_vals[v]
  mc_results[[v]] <- list(
    low  = lapply(1:4, function(k) simulate_paths(res_low_unc[[k]], params, Npaths, Nt, dx, dt, W, X0_current)),
    mid  = lapply(1:4, function(k) simulate_paths(res_mid_unc[[k]], params, Npaths, Nt, dx, dt, W, X0_current)),
    high = lapply(1:4, function(k) simulate_paths(res_high_unc[[k]], params, Npaths, Nt, dx, dt, W, X0_current))
  )
}

### Process plots combined (Rows = X0, Columns = theta_bar) ====

time_seq <- (0:Nt) * dt
time_seq_Nt <- (1:Nt) * dt

# Rows are X0, Columns are theta_bar
n_rows <- length(X0_vals)
n_cols <- length(th_vals)

### Panel Plot for Avg. Investment Paths ====
png("avg_investment_paths_combined_new_terminal_discount.png", width = a4_width, height = a4_height, units = "in", res = resolution)

par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))

col_ylim <- list()
for (k in 1:n_cols) {
  row_vals <- c()
  for(v in 1:n_rows) {
    row_vals <- c(row_vals, 
                  mc_results[[v]]$low[[k]]$S_lower, mc_results[[v]]$low[[k]]$S_upper,
                  mc_results[[v]]$mid[[k]]$S_lower, mc_results[[v]]$mid[[k]]$S_upper,
                  mc_results[[v]]$high[[k]]$S_lower, mc_results[[v]]$high[[k]]$S_upper)
  }
  col_ylim[[k]] <- range(row_vals, na.rm = TRUE)
}

for (v in 1:n_rows) {
  for (k in 1:n_cols) {
    mc_low  <- mc_results[[v]]$low
    mc_mid  <- mc_results[[v]]$mid
    mc_high <- mc_results[[v]]$high
    
    # Initialize empty plot for background
    plot(time_seq, mc_low[[k]]$S_avg, type = "n", 
         xlab = "Time (Years)", ylab = "Investment Value (X)",
         ylim = col_ylim[[k]])
    
    # Add grid
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    
    # Draw Confidence Intervals (Shading)
    polygon(c(time_seq, rev(time_seq)), c(mc_high[[k]]$S_lower, rev(mc_high[[k]]$S_upper)), col = adjustcolor("orange", alpha.f = 0.2), border = NA)
    polygon(c(time_seq, rev(time_seq)), c(mc_mid[[k]]$S_lower, rev(mc_mid[[k]]$S_upper)), col = adjustcolor("black", alpha.f = 0.2), border = NA)
    polygon(c(time_seq, rev(time_seq)), c(mc_low[[k]]$S_lower, rev(mc_low[[k]]$S_upper)), col = adjustcolor("darkgreen", alpha.f = 0.2), border = NA)
    
    # Draw data lines on top of grid
    lines(time_seq, mc_high[[k]]$S_avg, col = "orange", lwd = 2)
    lines(time_seq, mc_mid[[k]]$S_avg, col = "black", lwd = 2)
    lines(time_seq, mc_low[[k]]$S_avg, col = "darkgreen", lwd = 2)
    
    # Headers
    if (v == 1) { 
      mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    if (k == 1) { 
      mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}
mtext(paste0("Avg. of Investment Value (X) Paths over Time at ", 100*(1- params$terminal_discount), "% Terminal Depreciation"),
      outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                            paste0("Mid Unc (", uncertainty_nu[2], ")"),
                            paste0("Low Unc (", uncertainty_nu[3], ")")), 
       col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
dev.off()


### Panel Plot for Avg. Price Floor (theta) ====
png("avg_subsidies_paths_combined_new_terminal_discount.png", width = a4_width, height = a4_height, units = "in", res = resolution)
par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))

col_ylim <- list()
for (k in 1:n_cols) {
  row_vals <- c()
  for(v in 1:n_rows) {
    row_vals <- c(row_vals, 
                  mc_results[[v]]$low[[k]]$theta_lower, mc_results[[v]]$low[[k]]$theta_upper,
                  mc_results[[v]]$mid[[k]]$theta_lower, mc_results[[v]]$mid[[k]]$theta_upper,
                  mc_results[[v]]$high[[k]]$theta_lower, mc_results[[v]]$high[[k]]$theta_upper)
  }
  col_ylim[[k]] <- range(row_vals, na.rm = TRUE) 
}

for (v in 1:n_rows) {
  for (k in 1:n_cols) {
    mc_low  <- mc_results[[v]]$low
    mc_mid  <- mc_results[[v]]$mid
    mc_high <- mc_results[[v]]$high
    
    # Initialize empty plot
    plot(time_seq_Nt, mc_low[[k]]$theta_avg, type = "n", 
         xlab = "Time (Years)", ylab = "Price Floor (theta)",
         ylim = col_ylim[[k]])
    
    # Add grid and background reference lines
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    abline(h = th_vals[k], lty = 2, col = "darkgrey")
    abline(h = params$py, lty = 2, col = "darkgrey")
    
    # Draw Confidence Intervals (Shading)
    polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_high[[k]]$theta_lower, rev(mc_high[[k]]$theta_upper)), col = adjustcolor("orange", alpha.f = 0.2), border = NA)
    polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_mid[[k]]$theta_lower, rev(mc_mid[[k]]$theta_upper)), col = adjustcolor("black", alpha.f = 0.2), border = NA)
    polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_low[[k]]$theta_lower, rev(mc_low[[k]]$theta_upper)), col = adjustcolor("darkgreen", alpha.f = 0.2), border = NA)
    
    # Draw data lines on top
    lines(time_seq_Nt, mc_high[[k]]$theta_avg, col = "orange", lwd = 2)
    lines(time_seq_Nt, mc_mid[[k]]$theta_avg, col = "black", lwd = 2)
    lines(time_seq_Nt, mc_low[[k]]$theta_avg, col = "darkgreen", lwd = 2)
    
    if (v == 1) {
      mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    if (k == 1) {
      mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}
mtext(paste0("Avg. of Price Floor (theta) Paths over Time at ", 100*(1- params$terminal_discount), "% Terminal Depreciation"), 
      outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                            paste0("Mid Unc (", uncertainty_nu[2], ")"),
                            paste0("Low Unc (", uncertainty_nu[3], ")")), 
       col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
dev.off()


### Panel Plot for Avg. Quantities (q) ====
png("avg_quantities_paths_combined_new_terminal_discount.png", width = a4_width, height = a4_height, units = "in", res = resolution)
par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))

col_ylim <- list()
for (k in 1:n_cols) {
  row_vals <- c()
  for(v in 1:n_rows) {
    row_vals <- c(row_vals, 
                  mc_results[[v]]$low[[k]]$q_lower, mc_results[[v]]$low[[k]]$q_upper,
                  mc_results[[v]]$mid[[k]]$q_lower, mc_results[[v]]$mid[[k]]$q_upper,
                  mc_results[[v]]$high[[k]]$q_lower, mc_results[[v]]$high[[k]]$q_upper)
  }
  col_ylim[[k]] <- range(row_vals, na.rm = TRUE)
}

for (v in 1:n_rows) {
  for (k in 1:n_cols) {
    mc_low  <- mc_results[[v]]$low
    mc_mid  <- mc_results[[v]]$mid
    mc_high <- mc_results[[v]]$high
    
    # Initialize empty plot
    plot(time_seq_Nt, mc_low[[k]]$q_avg, type = "n", 
         xlab = "Time (Years)", ylab = "Quantity (q)",
         ylim = col_ylim[[k]])
    
    # Add grid and background reference lines
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    abline(h = params$q_max, lty = 2, col = "darkgrey")
    abline(h = 1.96, lty = 2, col = "darkgrey")
    
    # Draw Confidence Intervals (Shading)
    polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_high[[k]]$q_lower, rev(mc_high[[k]]$q_upper)), col = adjustcolor("orange", alpha.f = 0.2), border = NA)
    polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_mid[[k]]$q_lower, rev(mc_mid[[k]]$q_upper)), col = adjustcolor("black", alpha.f = 0.2), border = NA)
    polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_low[[k]]$q_lower, rev(mc_low[[k]]$q_upper)), col = adjustcolor("darkgreen", alpha.f = 0.2), border = NA)
    
    # Draw data lines on top
    lines(time_seq_Nt, mc_high[[k]]$q_avg, col = "orange", lwd = 2)
    lines(time_seq_Nt, mc_mid[[k]]$q_avg, col = "black", lwd = 2)
    lines(time_seq_Nt, mc_low[[k]]$q_avg, col = "darkgreen", lwd = 2)
    
    if (v == 1) {
      mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    if (k == 1) {
      mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}
mtext(paste0("Avg. of Quantity (q) Paths over Time at ", 100*(1- params$terminal_discount), "% Terminal Depreciation"),
      outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                            paste0("Mid Unc (", uncertainty_nu[2], ")"),
                            paste0("Low Unc (", uncertainty_nu[3], ")")), 
       col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
dev.off()


### Panel Plot for Avg. Investment Rate (gamma) ====
png("avg_investment_rates_paths_combined_new_terminal_discount.png", width = a4_width, height = a4_height, units = "in", res = resolution)
par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))

col_ylim <- list()
for (k in 1:n_cols) {
  row_vals <- c()
  for(v in 1:n_rows) {
    row_vals <- c(row_vals, 
                  mc_results[[v]]$low[[k]]$gamma_lower, mc_results[[v]]$low[[k]]$gamma_upper,
                  mc_results[[v]]$mid[[k]]$gamma_lower, mc_results[[v]]$mid[[k]]$gamma_upper,
                  mc_results[[v]]$high[[k]]$gamma_lower, mc_results[[v]]$high[[k]]$gamma_upper)
  }
  col_ylim[[k]] <- range(row_vals, na.rm = TRUE)
}

for (v in 1:n_rows) {
  for (k in 1:n_cols) {
    mc_low  <- mc_results[[v]]$low
    mc_mid  <- mc_results[[v]]$mid
    mc_high <- mc_results[[v]]$high
    
    # Initialize empty plot
    plot(time_seq_Nt, mc_low[[k]]$gamma_avg, type = "n", 
         xlab = "Time (Years)", ylab = "Investment Rate (gamma)",
         ylim = col_ylim[[k]])
    
    # Add grid and background reference lines
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    abline(h = 0, lty = 2, col = "darkgrey")
    
    # Draw Confidence Intervals (Shading)
    polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_high[[k]]$gamma_lower, rev(mc_high[[k]]$gamma_upper)), col = adjustcolor("orange", alpha.f = 0.2), border = NA)
    polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_mid[[k]]$gamma_lower, rev(mc_mid[[k]]$gamma_upper)), col = adjustcolor("black", alpha.f = 0.2), border = NA)
    polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_low[[k]]$gamma_lower, rev(mc_low[[k]]$gamma_upper)), col = adjustcolor("darkgreen", alpha.f = 0.2), border = NA)
    
    # Draw data lines on top
    lines(time_seq_Nt, mc_high[[k]]$gamma_avg, col = "orange", lwd = 2)
    lines(time_seq_Nt, mc_mid[[k]]$gamma_avg, col = "black", lwd = 2)
    lines(time_seq_Nt, mc_low[[k]]$gamma_avg, col = "darkgreen", lwd = 2)
    
    if (v == 1) {
      mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    if (k == 1) {
      mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}
mtext(paste0("Avg. of Investment Rate (gamma) Paths over Time at ", 100*(1- params$terminal_discount), "% Terminal Depreciation"),
      outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                            paste0("Mid Unc (", uncertainty_nu[2], ")"),
                            paste0("Low Unc (", uncertainty_nu[3], ")")), 
       col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
dev.off()

### Set back to original ====
params$terminal_discount <- prev_disc



## ==============================================================================
# PHASE 9: EXTENDED SENSITIVITY ANALYSIS (DELTA, KAPPA, SIGMA) ==================
## ==============================================================================

# Store original base values explicitly so they never get corrupted by previous runs
base_delta <- 0.02
base_sigma <- 0.2
base_kappa <- 0.5

param_grid <- data.frame(
  delta = c(0.002, 0.2, base_delta, base_delta, base_delta, base_delta),
  sigma = c(base_sigma, base_sigma, 0.002, 2, base_sigma, base_sigma),
  kappa = c(base_kappa, base_kappa, base_kappa, base_kappa, 0.05, 5)
)


# Values for evaluation
th_vals <- c(5, 3, 2.3, 2.1)
uncertainty_nu <- c(0.1, 1, 100)
X0_vals <- c(params$x_bar, params$x_bar - 1 , params$x_bar - 2)
X0_titles <- c("x_bar", "x_bar - 1", "x_bar - 2")
X0_tags <- c("x_bar", "x_bar_minus_1", "x_bar_minus_2")



# MAIN LOOP FOR THE 9 SCENARIOS

for (iter in 1:nrow(param_grid)) {
  
  # 1. Update params for this specific iteration
  params$delta <- param_grid$delta[iter]
  params$sigma <- param_grid$sigma[iter]
  params$kappa <- param_grid$kappa[iter]
  
  # 2. Recalculate W inside the loop because sigma might have changed
  W <- params$sigma * sqrt(dt) * Z_std
  
  # 3. Create a unique suffix for saving the plots
  suffix <- sprintf("_case_%d_d_%.3f_s_%.3f_k_%.3f.png", iter, params$delta, params$sigma, params$kappa)
  title_suffix <- sprintf("(Case %d: d=%.3f, s=%.3f, k=%.3f)", iter, params$delta, params$sigma, params$kappa)
  
  cat(sprintf("\n=================================================================\n"))
  cat(sprintf("Running Case %d of 6: delta = %.3f, sigma = %.3f, kappa = %.3f\n", 
              iter, params$delta, params$sigma, params$kappa))
  cat(sprintf("=================================================================\n\n"))
  
  
  
  ## 9.0 EVALUATE FINITE DIFFERENCE SCHEME =========================================
  
  cat("Calculating Base Scenarios (Varying Nu) sequentially...\n")
  res_high_unc <- lapply(th_vals, function(th) simulate_fixed_p(th, uncertainty_nu[1], params))
  res_mid_unc  <- lapply(th_vals, function(th) simulate_fixed_p(th, uncertainty_nu[2], params))
  res_low_unc  <- lapply(th_vals, function(th) simulate_fixed_p(th, uncertainty_nu[3], params))
  
  
  
  ## 9.4.1: MC of INVESTMENT PROCESS (Memory Optimized & Sequential) =============
  
  cat("Running Base Monte Carlo paths sequentially...\n")
  mc_results <- list()
  for (v in 1:length(X0_vals)) {
    X0_current <- X0_vals[v]
    mc_results[[v]] <- list(
      low  = lapply(1:4, function(k) simulate_paths(res_low_unc[[k]], params, Npaths, Nt, dx, dt, W, X0_current)),
      mid  = lapply(1:4, function(k) simulate_paths(res_mid_unc[[k]], params, Npaths, Nt, dx, dt, W, X0_current)),
      high = lapply(1:4, function(k) simulate_paths(res_high_unc[[k]], params, Npaths, Nt, dx, dt, W, X0_current))
    )
  }
  
  
  
  ## 9.4.2.1: Process plots combined (Rows = X0, Columns = theta_bar) ============
  time_seq <- (0:Nt) * dt
  time_seq_Nt <- (1:Nt) * dt
  
  # Rows are X0, Columns are theta_bar
  n_rows <- length(X0_vals)
  n_cols <- length(th_vals)
  
  ### Panel Plot for Avg. Investment Paths ====
  png(paste0("avg_investment_paths_combined", suffix), width = a4_width, height = a4_height, units = "in", res = resolution)
  par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))
  
  col_ylim <- list()
  for (k in 1:n_cols) {
    row_vals <- c()
    for(v in 1:n_rows) {
      row_vals <- c(row_vals, 
                    mc_results[[v]]$low[[k]]$S_lower, mc_results[[v]]$low[[k]]$S_upper,
                    mc_results[[v]]$mid[[k]]$S_lower, mc_results[[v]]$mid[[k]]$S_upper,
                    mc_results[[v]]$high[[k]]$S_lower, mc_results[[v]]$high[[k]]$S_upper)
    }
    col_ylim[[k]] <- range(row_vals, na.rm = TRUE)
  }
  
  for (v in 1:n_rows) {
    for (k in 1:n_cols) {
      mc_low  <- mc_results[[v]]$low
      mc_mid  <- mc_results[[v]]$mid
      mc_high <- mc_results[[v]]$high
      
      plot(time_seq, mc_low[[k]]$S_avg, type = "n", 
           xlab = "Time (Years)", ylab = "Investment Value (X)",
           ylim = col_ylim[[k]])
      
      grid(col = "lightgray", lty = "dotted", lwd = 1)
      
      polygon(c(time_seq, rev(time_seq)), c(mc_high[[k]]$S_lower, rev(mc_high[[k]]$S_upper)), col = adjustcolor("orange", alpha.f = 0.2), border = NA)
      polygon(c(time_seq, rev(time_seq)), c(mc_mid[[k]]$S_lower, rev(mc_mid[[k]]$S_upper)), col = adjustcolor("black", alpha.f = 0.2), border = NA)
      polygon(c(time_seq, rev(time_seq)), c(mc_low[[k]]$S_lower, rev(mc_low[[k]]$S_upper)), col = adjustcolor("darkgreen", alpha.f = 0.2), border = NA)
      
      lines(time_seq, mc_high[[k]]$S_avg, col = "orange", lwd = 2)
      lines(time_seq, mc_mid[[k]]$S_avg, col = "black", lwd = 2)
      lines(time_seq, mc_low[[k]]$S_avg, col = "darkgreen", lwd = 2)
      
      if (v == 1) { 
        mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
      }
      if (k == 1) { 
        mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
      }
    }
  }
  mtext(paste("Avg. Investment Value (X) Paths", title_suffix), outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend("bottom", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                              paste0("Mid Unc (", uncertainty_nu[2], ")"),
                              paste0("Low Unc (", uncertainty_nu[3], ")")), 
         col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
  dev.off()
  
  
  ### Panel Plot for Avg. Price Floor (theta) ====
  png(paste0("avg_subsidies_paths_combined", suffix), width = a4_width, height = a4_height, units = "in", res = resolution)
  par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))
  
  col_ylim <- list()
  for (k in 1:n_cols) {
    row_vals <- c()
    for(v in 1:n_rows) {
      row_vals <- c(row_vals, 
                    mc_results[[v]]$low[[k]]$theta_lower, mc_results[[v]]$low[[k]]$theta_upper,
                    mc_results[[v]]$mid[[k]]$theta_lower, mc_results[[v]]$mid[[k]]$theta_upper,
                    mc_results[[v]]$high[[k]]$theta_lower, mc_results[[v]]$high[[k]]$theta_upper)
    }
    col_ylim[[k]] <- range(row_vals, na.rm = TRUE)
  }
  
  for (v in 1:n_rows) {
    for (k in 1:n_cols) {
      mc_low  <- mc_results[[v]]$low
      mc_mid  <- mc_results[[v]]$mid
      mc_high <- mc_results[[v]]$high
      
      plot(time_seq_Nt, mc_low[[k]]$theta_avg, type = "n", 
           xlab = "Time (Years)", ylab = "Price Floor (theta)",
           ylim = col_ylim[[k]])
      
      grid(col = "lightgray", lty = "dotted", lwd = 1)
      abline(h = th_vals[k], lty = 2, col = "darkgrey")
      abline(h = params$py, lty = 2, col = "darkgrey")
      
      polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_high[[k]]$theta_lower, rev(mc_high[[k]]$theta_upper)), col = adjustcolor("orange", alpha.f = 0.2), border = NA)
      polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_mid[[k]]$theta_lower, rev(mc_mid[[k]]$theta_upper)), col = adjustcolor("black", alpha.f = 0.2), border = NA)
      polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_low[[k]]$theta_lower, rev(mc_low[[k]]$theta_upper)), col = adjustcolor("darkgreen", alpha.f = 0.2), border = NA)
      
      lines(time_seq_Nt, mc_high[[k]]$theta_avg, col = "orange", lwd = 2)
      lines(time_seq_Nt, mc_mid[[k]]$theta_avg, col = "black", lwd = 2)
      lines(time_seq_Nt, mc_low[[k]]$theta_avg, col = "darkgreen", lwd = 2)
      
      if (v == 1) {
        mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
      }
      if (k == 1) {
        mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
      }
    }
  }
  mtext(paste("Avg. Price Floor (theta) Paths", title_suffix), outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend("bottom", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                              paste0("Mid Unc (", uncertainty_nu[2], ")"),
                              paste0("Low Unc (", uncertainty_nu[3], ")")), 
         col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
  dev.off()
  
  
  ### Panel Plot for Avg. Quantities (q) ====
  png(paste0("avg_quantities_paths_combined", suffix), width = a4_width, height = a4_height, units = "in", res = resolution)
  par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))
  
  col_ylim <- list()
  for (k in 1:n_cols) {
    row_vals <- c()
    for(v in 1:n_rows) {
      row_vals <- c(row_vals, 
                    mc_results[[v]]$low[[k]]$q_lower, mc_results[[v]]$low[[k]]$q_upper,
                    mc_results[[v]]$mid[[k]]$q_lower, mc_results[[v]]$mid[[k]]$q_upper,
                    mc_results[[v]]$high[[k]]$q_lower, mc_results[[v]]$high[[k]]$q_upper)
    }
    col_ylim[[k]] <- range(row_vals, na.rm = TRUE)
  }
  
  for (v in 1:n_rows) {
    for (k in 1:n_cols) {
      mc_low  <- mc_results[[v]]$low
      mc_mid  <- mc_results[[v]]$mid
      mc_high <- mc_results[[v]]$high
      
      plot(time_seq_Nt, mc_low[[k]]$q_avg, type = "n", 
           xlab = "Time (Years)", ylab = "Quantity (q)",
           ylim = col_ylim[[k]])
      
      grid(col = "lightgray", lty = "dotted", lwd = 1)
      abline(h = params$q_max, lty = 2, col = "darkgrey")
      abline(h = 1.96, lty = 2, col = "darkgrey")
      
      polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_high[[k]]$q_lower, rev(mc_high[[k]]$q_upper)), col = adjustcolor("orange", alpha.f = 0.2), border = NA)
      polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_mid[[k]]$q_lower, rev(mc_mid[[k]]$q_upper)), col = adjustcolor("black", alpha.f = 0.2), border = NA)
      polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_low[[k]]$q_lower, rev(mc_low[[k]]$q_upper)), col = adjustcolor("darkgreen", alpha.f = 0.2), border = NA)
      
      lines(time_seq_Nt, mc_high[[k]]$q_avg, col = "orange", lwd = 2)
      lines(time_seq_Nt, mc_mid[[k]]$q_avg, col = "black", lwd = 2)
      lines(time_seq_Nt, mc_low[[k]]$q_avg, col = "darkgreen", lwd = 2)
      
      if (v == 1) {
        mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
      }
      if (k == 1) {
        mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
      }
    }
  }
  mtext(paste("Avg. Quantity (q) Paths", title_suffix), outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend("bottom", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                              paste0("Mid Unc (", uncertainty_nu[2], ")"),
                              paste0("Low Unc (", uncertainty_nu[3], ")")), 
         col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
  dev.off()
  
  
  ### Panel Plot for Avg. Investment Rate (gamma) ====
  png(paste0("avg_investment_rates_paths_combined", suffix), width = a4_width, height = a4_height, units = "in", res = resolution)
  par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))
  
  col_ylim <- list()
  for (k in 1:n_cols) {
    row_vals <- c()
    for(v in 1:n_rows) {
      row_vals <- c(row_vals, 
                    mc_results[[v]]$low[[k]]$gamma_lower, mc_results[[v]]$low[[k]]$gamma_upper,
                    mc_results[[v]]$mid[[k]]$gamma_lower, mc_results[[v]]$mid[[k]]$gamma_upper,
                    mc_results[[v]]$high[[k]]$gamma_lower, mc_results[[v]]$high[[k]]$gamma_upper)
    }
    col_ylim[[k]] <- range(row_vals, na.rm = TRUE)
  }
  
  for (v in 1:n_rows) {
    for (k in 1:n_cols) {
      mc_low  <- mc_results[[v]]$low
      mc_mid  <- mc_results[[v]]$mid
      mc_high <- mc_results[[v]]$high
      
      plot(time_seq_Nt, mc_low[[k]]$gamma_avg, type = "n", 
           xlab = "Time (Years)", ylab = "Investment Rate (gamma)",
           ylim = col_ylim[[k]])
      
      grid(col = "lightgray", lty = "dotted", lwd = 1)
      abline(h = 0, lty = 2, col = "darkgrey")
      
      polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_high[[k]]$gamma_lower, rev(mc_high[[k]]$gamma_upper)), col = adjustcolor("orange", alpha.f = 0.2), border = NA)
      polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_mid[[k]]$gamma_lower, rev(mc_mid[[k]]$gamma_upper)), col = adjustcolor("black", alpha.f = 0.2), border = NA)
      polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_low[[k]]$gamma_lower, rev(mc_low[[k]]$gamma_upper)), col = adjustcolor("darkgreen", alpha.f = 0.2), border = NA)
      
      lines(time_seq_Nt, mc_high[[k]]$gamma_avg, col = "orange", lwd = 2)
      lines(time_seq_Nt, mc_mid[[k]]$gamma_avg, col = "black", lwd = 2)
      lines(time_seq_Nt, mc_low[[k]]$gamma_avg, col = "darkgreen", lwd = 2)
      
      if (v == 1) {
        mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
      }
      if (k == 1) {
        mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
      }
    }
  }
  mtext(paste("Avg. Investment Rate (gamma) Paths", title_suffix), outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend("bottom", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                              paste0("Mid Unc (", uncertainty_nu[2], ")"),
                              paste0("Low Unc (", uncertainty_nu[3], ")")), 
         col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
  dev.off()
  
  
  
  ## 9.4.2.2: Process plots for Avg. Quantities (q) Breakdown ====================
  
  ### Panel Plot for Avg. Quantities (q) Breakdown ====
  png(paste0("avg_quantities_areas_breakdown_combined", suffix), width = a4_width, height = a4_height, units = "in", res = resolution)
  
  n_rows_nu <- 3 
  par(mfrow = c(n_rows_nu, n_cols), oma = c(4, 6, 4, 1), mar = c(4, 4, 3, 1))
  
  X0_state <- 1 # We plot for the first initial state X0
  
  # Pre-calculate global ylim incorporating the upper and lower confidence intervals
  col_ylim <- list()
  for (k in 1:n_cols) {
    col_vals <- c(mc_results[[X0_state]]$low[[k]]$q_lower, mc_results[[X0_state]]$low[[k]]$q_upper,
                  mc_results[[X0_state]]$mid[[k]]$q_lower, mc_results[[X0_state]]$mid[[k]]$q_upper,
                  mc_results[[X0_state]]$high[[k]]$q_lower, mc_results[[X0_state]]$high[[k]]$q_upper)
    
    lims <- range(col_vals, na.rm = TRUE)
    col_ylim[[k]] <- c(0, max(lims[2], params$q_max)) 
  }
  
  for (r in 1:n_rows_nu) {
    for (k in 1:n_cols) {
      
      if (r == 1) curr_data <- mc_results[[X0_state]]$high[[k]]
      if (r == 2) curr_data <- mc_results[[X0_state]]$mid[[k]]
      if (r == 3) curr_data <- mc_results[[X0_state]]$low[[k]]
      
      q_val <- curr_data$q_avg
      q_lower <- curr_data$q_lower
      q_upper <- curr_data$q_upper
      
      # CALCULATE Pg_x
      x_val <- curr_data$S_avg[-1] 
      Pg_x <- params$p_g * pmax(x_val - params$x_bar, 0)
      
      # Prevent overlap bugs
      q_val <- pmax(q_val, Pg_x)
      
      # Draw Plot
      plot(time_seq_Nt, q_val, type = "n", 
           xlab = "Time (Years)", ylab = "Quantity (q)",
           ylim = col_ylim[[k]])
      
      grid(col = "lightgray", lty = "dotted", lwd = 1)
      abline(h = params$q_max, lty = 2, col = "darkgrey")
      abline(h = 1.96, lty = 2, col = "darkgrey")
      
      # Draw Green Area
      polygon(x = c(time_seq_Nt, rev(time_seq_Nt)), 
              y = c(rep(0, length(time_seq_Nt)), rev(Pg_x)), 
              col = "lightgreen", border = NA)
      
      # Draw Brown Area
      polygon(x = c(time_seq_Nt, rev(time_seq_Nt)), 
              y = c(Pg_x, rev(q_val)), 
              col = "tan", border = NA)
      
      # Add Borders
      lines(time_seq_Nt, Pg_x, col = "darkgreen", lwd = 2)
      lines(time_seq_Nt, q_val, col = "saddlebrown", lwd = 2)
      
      # Add dashed line for brown overproduction
      lines(time_seq_Nt, Pg_x + 1.96, col = "saddlebrown", lty = 2, lwd = 1)
      
      # Margins and Titles
      if (r == 1) {
        mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
      }
      if (k == 1) {
        nu_label <- switch(r, 
                           paste0("High Unc (", uncertainty_nu[1], ")"),
                           paste0("Mid Unc (", uncertainty_nu[2], ")"),
                           paste0("Low Unc (", uncertainty_nu[3], ")"))
        mtext(nu_label, side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
      }
    }
  }
  
  # Main Outer Title
  mtext(sprintf("Avg. Quantity (q) Paths for X0 = %s %s", X0_titles[X0_state], title_suffix), 
        outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)
  
  # Outer Legend
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend("bottom", 
         legend = c(expression("P"["g"]*"(x) (Green Area)"), 
                    expression("q - P"["g"]*"(x) (Brown Area)"), 
                    expression("P"["g"]*"(x) + 1.96 (Brown Prod. Lower Bound)")), 
         fill = c("lightgreen", "tan", NA), 
         border = c("darkgreen", "saddlebrown", NA),
         col = c(NA, NA, "saddlebrown"), 
         lty = c(0, 0, 2), # 0 = no line for boxes, 2 = dashed line for the 3rd item
         lwd = c(0, 0, 2), 
         bty = "n", horiz = TRUE, cex = 1.2)
  dev.off()
  
  
  
  ## 9.8. EVALUATE again for different terminal discount ========================
  
  prev_disc <- params$terminal_discount
  
  cat("Calculating New Terminal Scenarios sequentially...\n")
  params$terminal_discount <- 0.9
  res_high_unc <- lapply(th_vals, function(th) simulate_fixed_p(th, uncertainty_nu[1], params))
  res_mid_unc  <- lapply(th_vals, function(th) simulate_fixed_p(th, uncertainty_nu[2], params))
  res_low_unc  <- lapply(th_vals, function(th) simulate_fixed_p(th, uncertainty_nu[3], params))
  
  
  ### MC of INVESTMENT PROCESS (Memory Optimized & Sequential) ====
  
  cat("Running Base Monte Carlo paths sequentially for New Terminal Discount...\n")
  mc_results <- list()
  for (v in 1:length(X0_vals)) {
    X0_current <- X0_vals[v]
    mc_results[[v]] <- list(
      low  = lapply(1:4, function(k) simulate_paths(res_low_unc[[k]], params, Npaths, Nt, dx, dt, W, X0_current)),
      mid  = lapply(1:4, function(k) simulate_paths(res_mid_unc[[k]], params, Npaths, Nt, dx, dt, W, X0_current)),
      high = lapply(1:4, function(k) simulate_paths(res_high_unc[[k]], params, Npaths, Nt, dx, dt, W, X0_current))
    )
  }
  
  ### Process plots combined (Rows = X0, Columns = theta_bar) ====
  
  ### Panel Plot for Avg. Investment Paths ====
  png(paste0("avg_investment_paths_combined_new_terminal_discount", suffix), width = a4_width, height = a4_height, units = "in", res = resolution)
  
  par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))
  
  col_ylim <- list()
  for (k in 1:n_cols) {
    row_vals <- c()
    for(v in 1:n_rows) {
      row_vals <- c(row_vals, 
                    mc_results[[v]]$low[[k]]$S_lower, mc_results[[v]]$low[[k]]$S_upper,
                    mc_results[[v]]$mid[[k]]$S_lower, mc_results[[v]]$mid[[k]]$S_upper,
                    mc_results[[v]]$high[[k]]$S_lower, mc_results[[v]]$high[[k]]$S_upper)
    }
    col_ylim[[k]] <- range(row_vals, na.rm = TRUE)
  }
  
  for (v in 1:n_rows) {
    for (k in 1:n_cols) {
      mc_low  <- mc_results[[v]]$low
      mc_mid  <- mc_results[[v]]$mid
      mc_high <- mc_results[[v]]$high
      
      # Initialize empty plot for background
      plot(time_seq, mc_low[[k]]$S_avg, type = "n", 
           xlab = "Time (Years)", ylab = "Investment Value (X)",
           ylim = col_ylim[[k]])
      
      # Add grid
      grid(col = "lightgray", lty = "dotted", lwd = 1)
      
      # Draw Confidence Intervals (Shading)
      polygon(c(time_seq, rev(time_seq)), c(mc_high[[k]]$S_lower, rev(mc_high[[k]]$S_upper)), col = adjustcolor("orange", alpha.f = 0.2), border = NA)
      polygon(c(time_seq, rev(time_seq)), c(mc_mid[[k]]$S_lower, rev(mc_mid[[k]]$S_upper)), col = adjustcolor("black", alpha.f = 0.2), border = NA)
      polygon(c(time_seq, rev(time_seq)), c(mc_low[[k]]$S_lower, rev(mc_low[[k]]$S_upper)), col = adjustcolor("darkgreen", alpha.f = 0.2), border = NA)
      
      # Draw data lines on top of grid
      lines(time_seq, mc_high[[k]]$S_avg, col = "orange", lwd = 2)
      lines(time_seq, mc_mid[[k]]$S_avg, col = "black", lwd = 2)
      lines(time_seq, mc_low[[k]]$S_avg, col = "darkgreen", lwd = 2)
      
      # Headers
      if (v == 1) { 
        mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
      }
      if (k == 1) { 
        mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
      }
    }
  }
  mtext(sprintf("Avg. Investment (X) Paths at %s%% Terminal Depreciation %s", 100*(1- params$terminal_discount), title_suffix),
        outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend("bottom", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                              paste0("Mid Unc (", uncertainty_nu[2], ")"),
                              paste0("Low Unc (", uncertainty_nu[3], ")")), 
         col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
  dev.off()
  
  
  ### Panel Plot for Avg. Price Floor (theta) ====
  png(paste0("avg_subsidies_paths_combined_new_terminal_discount", suffix), width = a4_width, height = a4_height, units = "in", res = resolution)
  par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))
  
  col_ylim <- list()
  for (k in 1:n_cols) {
    row_vals <- c()
    for(v in 1:n_rows) {
      row_vals <- c(row_vals, 
                    mc_results[[v]]$low[[k]]$theta_lower, mc_results[[v]]$low[[k]]$theta_upper,
                    mc_results[[v]]$mid[[k]]$theta_lower, mc_results[[v]]$mid[[k]]$theta_upper,
                    mc_results[[v]]$high[[k]]$theta_lower, mc_results[[v]]$high[[k]]$theta_upper)
    }
    col_ylim[[k]] <- range(row_vals, na.rm = TRUE) 
  }
  
  for (v in 1:n_rows) {
    for (k in 1:n_cols) {
      mc_low  <- mc_results[[v]]$low
      mc_mid  <- mc_results[[v]]$mid
      mc_high <- mc_results[[v]]$high
      
      # Initialize empty plot
      plot(time_seq_Nt, mc_low[[k]]$theta_avg, type = "n", 
           xlab = "Time (Years)", ylab = "Price Floor (theta)",
           ylim = col_ylim[[k]])
      
      # Add grid and background reference lines
      grid(col = "lightgray", lty = "dotted", lwd = 1)
      abline(h = th_vals[k], lty = 2, col = "darkgrey")
      abline(h = params$py, lty = 2, col = "darkgrey")
      
      # Draw Confidence Intervals (Shading)
      polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_high[[k]]$theta_lower, rev(mc_high[[k]]$theta_upper)), col = adjustcolor("orange", alpha.f = 0.2), border = NA)
      polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_mid[[k]]$theta_lower, rev(mc_mid[[k]]$theta_upper)), col = adjustcolor("black", alpha.f = 0.2), border = NA)
      polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_low[[k]]$theta_lower, rev(mc_low[[k]]$theta_upper)), col = adjustcolor("darkgreen", alpha.f = 0.2), border = NA)
      
      # Draw data lines on top
      lines(time_seq_Nt, mc_high[[k]]$theta_avg, col = "orange", lwd = 2)
      lines(time_seq_Nt, mc_mid[[k]]$theta_avg, col = "black", lwd = 2)
      lines(time_seq_Nt, mc_low[[k]]$theta_avg, col = "darkgreen", lwd = 2)
      
      if (v == 1) {
        mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
      }
      if (k == 1) {
        mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
      }
    }
  }
  mtext(sprintf("Avg. Price Floor (theta) Paths at %s%% Terminal Depreciation %s", 100*(1- params$terminal_discount), title_suffix), 
        outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend("bottom", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                              paste0("Mid Unc (", uncertainty_nu[2], ")"),
                              paste0("Low Unc (", uncertainty_nu[3], ")")), 
         col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
  dev.off()
  
  
  ### Panel Plot for Avg. Quantities (q) ====
  png(paste0("avg_quantities_paths_combined_new_terminal_discount", suffix), width = a4_width, height = a4_height, units = "in", res = resolution)
  par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))
  
  col_ylim <- list()
  for (k in 1:n_cols) {
    row_vals <- c()
    for(v in 1:n_rows) {
      row_vals <- c(row_vals, 
                    mc_results[[v]]$low[[k]]$q_lower, mc_results[[v]]$low[[k]]$q_upper,
                    mc_results[[v]]$mid[[k]]$q_lower, mc_results[[v]]$mid[[k]]$q_upper,
                    mc_results[[v]]$high[[k]]$q_lower, mc_results[[v]]$high[[k]]$q_upper)
    }
    col_ylim[[k]] <- range(row_vals, na.rm = TRUE)
  }
  
  for (v in 1:n_rows) {
    for (k in 1:n_cols) {
      mc_low  <- mc_results[[v]]$low
      mc_mid  <- mc_results[[v]]$mid
      mc_high <- mc_results[[v]]$high
      
      # Initialize empty plot
      plot(time_seq_Nt, mc_low[[k]]$q_avg, type = "n", 
           xlab = "Time (Years)", ylab = "Quantity (q)",
           ylim = col_ylim[[k]])
      
      # Add grid and background reference lines
      grid(col = "lightgray", lty = "dotted", lwd = 1)
      abline(h = params$q_max, lty = 2, col = "darkgrey")
      abline(h = 1.96, lty = 2, col = "darkgrey")
      
      # Draw Confidence Intervals (Shading)
      polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_high[[k]]$q_lower, rev(mc_high[[k]]$q_upper)), col = adjustcolor("orange", alpha.f = 0.2), border = NA)
      polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_mid[[k]]$q_lower, rev(mc_mid[[k]]$q_upper)), col = adjustcolor("black", alpha.f = 0.2), border = NA)
      polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_low[[k]]$q_lower, rev(mc_low[[k]]$q_upper)), col = adjustcolor("darkgreen", alpha.f = 0.2), border = NA)
      
      # Draw data lines on top
      lines(time_seq_Nt, mc_high[[k]]$q_avg, col = "orange", lwd = 2)
      lines(time_seq_Nt, mc_mid[[k]]$q_avg, col = "black", lwd = 2)
      lines(time_seq_Nt, mc_low[[k]]$q_avg, col = "darkgreen", lwd = 2)
      
      if (v == 1) {
        mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
      }
      if (k == 1) {
        mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
      }
    }
  }
  mtext(sprintf("Avg. Quantity (q) Paths at %s%% Terminal Depreciation %s", 100*(1- params$terminal_discount), title_suffix),
        outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend("bottom", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                              paste0("Mid Unc (", uncertainty_nu[2], ")"),
                              paste0("Low Unc (", uncertainty_nu[3], ")")), 
         col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
  dev.off()
  
  
  ### Panel Plot for Avg. Investment Rate (gamma) ====
  png(paste0("avg_investment_rates_paths_combined_new_terminal_discount", suffix), width = a4_width, height = a4_height, units = "in", res = resolution)
  par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))
  
  col_ylim <- list()
  for (k in 1:n_cols) {
    row_vals <- c()
    for(v in 1:n_rows) {
      row_vals <- c(row_vals, 
                    mc_results[[v]]$low[[k]]$gamma_lower, mc_results[[v]]$low[[k]]$gamma_upper,
                    mc_results[[v]]$mid[[k]]$gamma_lower, mc_results[[v]]$mid[[k]]$gamma_upper,
                    mc_results[[v]]$high[[k]]$gamma_lower, mc_results[[v]]$high[[k]]$gamma_upper)
    }
    col_ylim[[k]] <- range(row_vals, na.rm = TRUE)
  }
  
  for (v in 1:n_rows) {
    for (k in 1:n_cols) {
      mc_low  <- mc_results[[v]]$low
      mc_mid  <- mc_results[[v]]$mid
      mc_high <- mc_results[[v]]$high
      
      # Initialize empty plot
      plot(time_seq_Nt, mc_low[[k]]$gamma_avg, type = "n", 
           xlab = "Time (Years)", ylab = "Investment Rate (gamma)",
           ylim = col_ylim[[k]])
      
      # Add grid and background reference lines
      grid(col = "lightgray", lty = "dotted", lwd = 1)
      abline(h = 0, lty = 2, col = "darkgrey")
      
      # Draw Confidence Intervals (Shading)
      polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_high[[k]]$gamma_lower, rev(mc_high[[k]]$gamma_upper)), col = adjustcolor("orange", alpha.f = 0.2), border = NA)
      polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_mid[[k]]$gamma_lower, rev(mc_mid[[k]]$gamma_upper)), col = adjustcolor("black", alpha.f = 0.2), border = NA)
      polygon(c(time_seq_Nt, rev(time_seq_Nt)), c(mc_low[[k]]$gamma_lower, rev(mc_low[[k]]$gamma_upper)), col = adjustcolor("darkgreen", alpha.f = 0.2), border = NA)
      
      # Draw data lines on top
      lines(time_seq_Nt, mc_high[[k]]$gamma_avg, col = "orange", lwd = 2)
      lines(time_seq_Nt, mc_mid[[k]]$gamma_avg, col = "black", lwd = 2)
      lines(time_seq_Nt, mc_low[[k]]$gamma_avg, col = "darkgreen", lwd = 2)
      
      if (v == 1) {
        mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
      }
      if (k == 1) {
        mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
      }
    }
  }
  mtext(sprintf("Avg. Investment Rate (gamma) Paths at %s%% Terminal Depreciation %s", 100*(1- params$terminal_discount), title_suffix),
        outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend("bottom", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                              paste0("Mid Unc (", uncertainty_nu[2], ")"),
                              paste0("Low Unc (", uncertainty_nu[3], ")")), 
         col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
  dev.off()
  
  ### Set back to original terminal discount before moving to the next case ====
  params$terminal_discount <- prev_disc
  
} # END OF 9-CASE ITERATION LOOP