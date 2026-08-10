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
  term3 <- max(theta - prm$py, 0) * e_x * q   # Subsidy Gain
  term4 <- prm$nu * (theta - prm$theta_bar)^2 # Penalty
  
  g_val <- term1 + term2 + term3 + term4
  
  return(c(q,theta,g_val))
}

simulate_fixed_p <- function(theta_bar_input, nu_input, params){
  
  ## ==============================================================================
  # PHASE 2: GRID SETUP & PRE-COMPUTATION =========================================
  ## ==============================================================================
  
  # Adjust the parameters for uncertainty nu and subsidy benchmark theta_bar
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

### Plot for Optimal Subsidy (theta) ====
png("optimal_subsidy_theta.png", width = 800, height = 600)
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
       xlab = "Investment (x)", ylab = "Subsidy (theta)",
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
mtext("Optimal Subsidy (theta) across scenarios", outer = TRUE, cex = 1.2, font = 2)
dev.off()

### Plot for Optimal Production (q) ====
png("optimal_production_q.png", width = 800, height = 600)
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
       xlab = "Investment (x)", ylab = "Production (q)",
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
mtext("Optimal Production (q) across scenarios", outer = TRUE, cex = 1.2, font = 2)
dev.off()

## ==============================================================================
# PHASE 4.0.1: Phase Plots (theta vs q and q vs theta) ==========================
## ==============================================================================

### Optimal Subsidy (theta) given Production (q) ====
png("optimal_subsidy_given_q.png", width = 800, height = 600)
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
       xlab = "Optimal Production (q)", ylab = "Optimal Subsidy (theta)",
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
mtext("Optimal Subsidy (theta) given Production (q)", outer = TRUE, cex = 1.2, font = 2)
dev.off()

### Optimal Production (q) given Subsidy (theta) ====
png("optimal_production_given_theta.png", width = 800, height = 600)
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
       xlab = "Optimal Subsidy (theta)", ylab = "Optimal Production (q)",
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
mtext("Optimal Production (q) given Subsidy (theta)", outer = TRUE, cex = 1.2, font = 2)
dev.off()

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
  S_avg[1]  <- X0
  
  gamma_avg <- numeric(Nt)
  q_avg     <- numeric(Nt)
  theta_avg <- numeric(Nt)
  
  for (t_step in 1:Nt) {
    idx <- pmax(1, pmin(Nx - 2, floor((S_current - X_min) / dx) + 1))
    
    gamma_step <- res_obj$gamma[[t_step]][idx]
    q_step     <- res_obj$q[idx]
    theta_step <- res_obj$theta[idx]
    
    gamma_avg[t_step] <- mean(gamma_step)
    q_avg[t_step]     <- mean(q_step)
    theta_avg[t_step] <- mean(theta_step)
    
    S_current <- S_current + (gamma_step - params$delta * S_current) * dt + W[, t_step]
    S_avg[t_step+1] <- mean(S_current)
  }
  
  list(
    S_avg     = S_avg,
    gamma_avg = gamma_avg,
    q_avg     = q_avg,
    theta_avg = theta_avg
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
# PHASE 4.2: Process plots ======================================================
## ==============================================================================
time_seq <- (0:Nt) * dt
time_seq_Nt <- (1:Nt) * dt

for (v in 1:length(X0_vals)) {
  mc_low  <- mc_results[[v]]$low
  mc_mid  <- mc_results[[v]]$mid
  mc_high <- mc_results[[v]]$high
  
  ### Plot for Avg. Investment Paths ====
  png(paste0("avg_investment_paths_", X0_tags[v], ".png"), width = 800, height = 600)
  par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
  for (k in 1:4) {
    plot(time_seq, mc_low[[k]]$S_avg, type = "l", col = "darkgreen", lwd = 2,
         main = paste("theta_bar =", th_vals[k]),
         xlab = "Time (Years)", ylab = "Investment",
         ylim = range(c(mc_low[[k]]$S_avg, mc_mid[[k]]$S_avg, mc_high[[k]]$S_avg)))
    
    lines(time_seq, mc_mid[[k]]$S_avg, col = "black", lwd = 2)
    lines(time_seq, mc_high[[k]]$S_avg, col = "orange", lwd = 2)
    
    if(k == 1) legend("topleft", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                                            paste0("Mid Unc (", uncertainty_nu[2], ")"),
                                            paste0("Low Unc (", uncertainty_nu[3], ")")), 
                      col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n")
  }
  mtext(paste0("Avg. Investment (X) over Time (X0 = ", X0_titles[v], ")"), outer = TRUE, cex = 1.2, font = 2)
  dev.off()
  
  ### Plot for Avg. Subsidies (theta) ====
  png(paste0("avg_subsidies_paths_", X0_tags[v], ".png"), width = 800, height = 600)
  par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
  for (k in 1:4) {
    plot(time_seq_Nt, mc_low[[k]]$theta_avg, type = "l", col = "darkgreen", lwd = 2,
         main = paste("theta_bar =", th_vals[k]),
         xlab = "Time (Years)", ylab = "Subsidies",
         ylim = range(c(mc_low[[k]]$theta_avg, mc_mid[[k]]$theta_avg, mc_high[[k]]$theta_avg, params$py, th_vals[k])))
    
    lines(time_seq_Nt, mc_mid[[k]]$theta_avg, col = "black", lwd = 2)
    lines(time_seq_Nt, mc_high[[k]]$theta_avg, col = "orange", lwd = 2)
    
    abline(h = th_vals[k], lty = 2, col = "darkgrey")
    abline(h = params$py, lty = 2, col = "darkgrey")
    
    if(k == 1) legend("bottomleft", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                                               paste0("Mid Unc (", uncertainty_nu[2], ")"),
                                               paste0("Low Unc (", uncertainty_nu[3], ")")), 
                      col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n")
  }
  mtext(paste0("Avg. Subsidies Paths (X0 = ", X0_titles[v], ")"), outer = TRUE, cex = 1.2, font = 2)
  dev.off()
  
  ### Plot for Avg. Quantities (q) ====
  png(paste0("avg_quantities_paths_", X0_tags[v], ".png"), width = 800, height = 600)
  par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
  for (k in 1:4) {
    plot(time_seq_Nt, mc_low[[k]]$q_avg, type = "l", col = "darkgreen", lwd = 2,
         main = paste("theta_bar =", th_vals[k]),
         xlab = "Time (Years)", ylab = "Quantities",
         ylim = range(c(mc_low[[k]]$q_avg, mc_mid[[k]]$q_avg, mc_high[[k]]$q_avg, params$py, th_vals[k])))
    
    lines(time_seq_Nt, mc_mid[[k]]$q_avg, col = "black", lwd = 2)
    lines(time_seq_Nt, mc_high[[k]]$q_avg, col = "orange", lwd = 2)
    
    abline(h = params$q_max, lty = 2, col = "darkgrey")
    abline(h = 1.96, lty = 2, col = "darkgrey")
    
    if(k == 1) legend("topleft", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                                            paste0("Mid Unc (", uncertainty_nu[2], ")"),
                                            paste0("Low Unc (", uncertainty_nu[3], ")")), 
                      col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n")
  }
  mtext(paste0("Avg. Quantities Paths (X0 = ", X0_titles[v], ")"), outer = TRUE, cex = 1.2, font = 2)
  dev.off()
  
  ### Plot for Avg. Investment Rates (gamma) ====
  png(paste0("avg_investment_rates_paths_", X0_tags[v], ".png"), width = 800, height = 600)
  par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
  for (k in 1:4) {
    plot(time_seq_Nt, mc_low[[k]]$gamma_avg, type = "l", col = "darkgreen", lwd = 2,
         main = paste("theta_bar =", th_vals[k]),
         xlab = "Time (Years)", ylab = "Investment Rates",
         ylim = range(c(mc_low[[k]]$gamma_avg, mc_mid[[k]]$gamma_avg, 
                        mc_high[[k]]$gamma_avg, params$py, th_vals[k])))
    
    lines(time_seq_Nt, mc_mid[[k]]$gamma_avg, col = "black", lwd = 2)
    lines(time_seq_Nt, mc_high[[k]]$gamma_avg, col = "orange", lwd = 2)
    
    abline(h = 0, lty = 2, col = "darkgrey")
    
    if(k == 1) legend("bottomleft", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                                               paste0("Mid Unc (", uncertainty_nu[2], ")"),
                                               paste0("Low Unc (", uncertainty_nu[3], ")")), 
                      col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n")
  }
  mtext(paste0("Avg. Investment Rates Paths (X0 = ", X0_titles[v], ")"), outer = TRUE, cex = 1.2, font = 2)
  dev.off()
}


## ==============================================================================
# PHASE 4.2.1: Process plots combined (Rows = X0, Columns = theta_bar) ==========
## ==============================================================================
time_seq <- (0:Nt) * dt
time_seq_Nt <- (1:Nt) * dt

# Rows are X0, Columns are theta_bar
n_rows <- length(X0_vals)
n_cols <- length(th_vals)

# Define A4 dimensions in inches (Landscape orientation)
a4_width  <- 11.69
a4_height <- 8.27
resolution <- 300 

### Panel Plot for Avg. Investment Paths ====
png("avg_investment_paths_combined.png", width = a4_width, height = a4_height, units = "in", res = resolution)

par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))

col_ylim <- list()
for (k in 1:n_cols) {
  row_vals <- c()
  for(v in 1:n_rows) {
    row_vals <- c(row_vals, 
                  mc_results[[v]]$low[[k]]$S_avg, 
                  mc_results[[v]]$mid[[k]]$S_avg, 
                  mc_results[[v]]$high[[k]]$S_avg)
  }
  col_ylim[[k]] <- range(row_vals, na.rm = TRUE) # Use double brackets
}

for (v in 1:n_rows) {
  for (k in 1:n_cols) {
    mc_low  <- mc_results[[v]]$low
    mc_mid  <- mc_results[[v]]$mid
    mc_high <- mc_results[[v]]$high
    
    # Initialize empty plot for background
    plot(time_seq, mc_low[[k]]$S_avg, type = "n", 
         xlab = "Time (Years)", ylab = "Investment",
         ylim = col_ylim[[k]])
    
    # Add grid
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    
    # Draw data lines on top of grid
    lines(time_seq, mc_low[[k]]$S_avg, col = "darkgreen", lwd = 2)
    lines(time_seq, mc_mid[[k]]$S_avg, col = "black", lwd = 2)
    lines(time_seq, mc_high[[k]]$S_avg, col = "orange", lwd = 2)
    
    # Headers
    if (v == 1) { # Top row gets the column headers (theta_bar)
      mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    if (k == 1) { # First column gets the row headers (X0)
      mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}
mtext("Avg. Investment (X) over Time", outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                            paste0("Mid Unc (", uncertainty_nu[2], ")"),
                            paste0("Low Unc (", uncertainty_nu[3], ")")), 
       col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
dev.off()


### Panel Plot for Avg. Subsidies (theta) ====
png("avg_subsidies_paths_combined.png", width = a4_width, height = a4_height, units = "in", res = resolution)
par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))

col_ylim <- list()
for (k in 1:n_cols) {
  row_vals <- c()
  for(v in 1:n_rows) {
    row_vals <- c(row_vals, 
                  mc_results[[v]]$low[[k]]$theta_avg, 
                  mc_results[[v]]$mid[[k]]$theta_avg, 
                  mc_results[[v]]$high[[k]]$theta_avg)
  }
  col_ylim[[k]] <- range(row_vals, na.rm = TRUE) # Use double brackets
}

for (v in 1:n_rows) {
  for (k in 1:n_cols) {
    mc_low  <- mc_results[[v]]$low
    mc_mid  <- mc_results[[v]]$mid
    mc_high <- mc_results[[v]]$high
    
    # Initialize empty plot
    plot(time_seq_Nt, mc_low[[k]]$theta_avg, type = "n", 
         xlab = "Time (Years)", ylab = "Subsidies",
         ylim = col_ylim[[k]])
    
    # Add grid and background reference lines
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    abline(h = th_vals[k], lty = 2, col = "darkgrey")
    abline(h = params$py, lty = 2, col = "darkgrey")
    
    # Draw data lines on top
    lines(time_seq_Nt, mc_low[[k]]$theta_avg, col = "darkgreen", lwd = 2)
    lines(time_seq_Nt, mc_mid[[k]]$theta_avg, col = "black", lwd = 2)
    lines(time_seq_Nt, mc_high[[k]]$theta_avg, col = "orange", lwd = 2)
    
    if (v == 1) {
      mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    if (k == 1) {
      mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}
mtext("Avg. Subsidies Paths", outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)

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
                  mc_results[[v]]$low[[k]]$q_avg, 
                  mc_results[[v]]$mid[[k]]$q_avg, 
                  mc_results[[v]]$high[[k]]$q_avg)
  }
  col_ylim[[k]] <- range(row_vals, na.rm = TRUE) # Use double brackets
}

for (v in 1:n_rows) {
  for (k in 1:n_cols) {
    mc_low  <- mc_results[[v]]$low
    mc_mid  <- mc_results[[v]]$mid
    mc_high <- mc_results[[v]]$high
    
    # Initialize empty plot
    plot(time_seq_Nt, mc_low[[k]]$q_avg, type = "n", 
         xlab = "Time (Years)", ylab = "Quantities",
         ylim = col_ylim[[k]])
    
    # Add grid and background reference lines
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    abline(h = params$q_max, lty = 2, col = "darkgrey")
    abline(h = 1.96, lty = 2, col = "darkgrey")
    
    # Draw data lines on top
    lines(time_seq_Nt, mc_low[[k]]$q_avg, col = "darkgreen", lwd = 2)
    lines(time_seq_Nt, mc_mid[[k]]$q_avg, col = "black", lwd = 2)
    lines(time_seq_Nt, mc_high[[k]]$q_avg, col = "orange", lwd = 2)
    
    if (v == 1) {
      mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    if (k == 1) {
      mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}
mtext("Avg. Quantities Paths", outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                            paste0("Mid Unc (", uncertainty_nu[2], ")"),
                            paste0("Low Unc (", uncertainty_nu[3], ")")), 
       col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
dev.off()


### Panel Plot for Avg. Investment Rates (gamma) ====
png("avg_investment_rates_paths_combined.png", width = a4_width, height = a4_height, units = "in", res = resolution)
par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))

col_ylim <- list()
for (k in 1:n_cols) {
  row_vals <- c()
  for(v in 1:n_rows) {
    row_vals <- c(row_vals, 
                  mc_results[[v]]$low[[k]]$gamma_avg, 
                  mc_results[[v]]$mid[[k]]$gamma_avg, 
                  mc_results[[v]]$high[[k]]$gamma_avg)
  }
  col_ylim[[k]] <- range(row_vals, na.rm = TRUE) # Use double brackets
}

for (v in 1:n_rows) {
  for (k in 1:n_cols) {
    mc_low  <- mc_results[[v]]$low
    mc_mid  <- mc_results[[v]]$mid
    mc_high <- mc_results[[v]]$high
    
    # Initialize empty plot
    plot(time_seq_Nt, mc_low[[k]]$gamma_avg, type = "n", 
         xlab = "Time (Years)", ylab = "Investment Rates",
         ylim = col_ylim[[k]])
    
    # Add grid and background reference lines
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    abline(h = 0, lty = 2, col = "darkgrey")
    
    # Draw data lines on top
    lines(time_seq_Nt, mc_low[[k]]$gamma_avg, col = "darkgreen", lwd = 2)
    lines(time_seq_Nt, mc_mid[[k]]$gamma_avg, col = "black", lwd = 2)
    lines(time_seq_Nt, mc_high[[k]]$gamma_avg, col = "orange", lwd = 2)
    
    if (v == 1) {
      mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    if (k == 1) {
      mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}
mtext("Avg. Investment Rates Paths", outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)

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

# Pre-calculate global ylim
col_ylim <- list()
for (k in 1:n_cols) {
  col_vals <- c(mc_results[[X0_state]]$low[[k]]$q_avg, 
                mc_results[[X0_state]]$mid[[k]]$q_avg, 
                mc_results[[X0_state]]$high[[k]]$q_avg)
  
  lims <- range(col_vals, na.rm = TRUE)
  col_ylim[[k]] <- c(0, max(lims[2], params$q_max)) 
}

for (r in 1:n_rows_nu) {
  for (k in 1:n_cols) {
    
    if (r == 1) curr_data <- mc_results[[X0_state]]$high[[k]]
    if (r == 2) curr_data <- mc_results[[X0_state]]$mid[[k]]
    if (r == 3) curr_data <- mc_results[[X0_state]]$low[[k]]
    
    q_val <- curr_data$q_avg
    
    # CALCULATE Pg_x
    x_val <- curr_data$S_avg[-1] 
    Pg_x <- params$p_g * pmax(x_val - params$x_bar, 0)
    
    # Prevent overlap bugs
    q_val <- pmax(q_val, Pg_x)
    
    # 2. Draw Plot
    plot(time_seq_Nt, q_val, type = "n", 
         xlab = "Time (Years)", ylab = "Quantities",
         ylim = col_ylim[[k]])
    
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    abline(h = params$q_max, lty = 2, col = "darkgrey")
    abline(h = 1.96, lty = 2, col = "darkgrey")
    
    # 3. Draw Green Area
    polygon(x = c(time_seq_Nt, rev(time_seq_Nt)), 
            y = c(rep(0, length(time_seq_Nt)), rev(Pg_x)), 
            col = "lightgreen", border = NA)
    
    # 4. Draw Brown Area
    polygon(x = c(time_seq_Nt, rev(time_seq_Nt)), 
            y = c(Pg_x, rev(q_val)), 
            col = "tan", border = NA)
    
    # Add Borders
    lines(time_seq_Nt, Pg_x, col = "darkgreen", lwd = 2)
    lines(time_seq_Nt, q_val, col = "saddlebrown", lwd = 2)
    
    # 5. Margins and Titles
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
mtext(paste0("Avg. Quantities Paths (Initial State X0 = ", X0_titles[X0_state], ")"), 
      outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)

# Outer Legend
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", 
       legend = c(expression("P"["g"]*"(x) (Green Area)"), expression("q - P"["g"]*"(x) (Brown Area)")), 
       fill = c("lightgreen", "tan"), 
       border = c("darkgreen", "saddlebrown"), 
       bty = "n", horiz = TRUE, cex = 1.2)
dev.off()


## ==============================================================================
# PHASE 5: Value plots ==========================================================
## ==============================================================================

color_vec <- c("darkblue", "blue", "lightblue", "lightgrey", "grey")
step_values <- c(0.25, 0.5, 0.75, 1) * Nt

plot_value_grid <- function(res_list, nu_val) {
  png(paste0("game_value_nu_", nu_val, ".png"), width = 800, height = 600)
  par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
  for(k in 1:4) {
    U_val <- matrix(unlist(res_list[[k]]$value), nrow = Nt + 1, ncol = Nx, byrow = TRUE)
    
    plot(X, U_val[1,], type = "l", col = color_vec[1], lwd = 2,
         xlab = "Investment X", ylab = "Game Value",
         main = paste("theta_bar =", th_vals[k]),
         ylim = range(U_val))
    
    for (i in 1:4){
      lines(X, U_val[step_values[i], ], col = color_vec[i+1], lwd = 2)
    }
    
    if (k == 1) {
      legend("bottomright", legend = c("t=0", paste("t=", step_values)), 
             col = color_vec, lwd = 2, bty = "n", cex=0.8)
    }
  }
  mtext(paste0("Game value over X (nu = ", nu_val, ")"), outer = TRUE, cex = 1.2, font = 2)
  dev.off()
}

plot_value_grid(res_high_unc, uncertainty_nu[1])
plot_value_grid(res_mid_unc, uncertainty_nu[2])
plot_value_grid(res_low_unc, uncertainty_nu[3])


## ==============================================================================
# PHASE 5.1: Panel Value plots ==================================================
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

plot_derivative_grid <- function(res_list, nu_val, line_color) {
  png(paste0("marginal_value_benchmark_nu_", nu_val, ".png"), width = 800, height = 600)
  par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
  for(k in 1:4) {
    Ux <- diff(res_list[[k]]$value[[1]]) / dx
    
    plot(X[-Nx], Ux, type = "l", col = line_color, lwd = 2,
         xlab = "Investment Space (X)", ylab = "U_x",
         main = paste("theta_bar =", th_vals[k]))
    abline(v = params$x_bar, lty = 3, col = "red")
  }
  mtext(paste0("First Derivative (Marginal Value) Benchmark, nu = ", nu_val), outer = TRUE, cex = 1.2, font = 2)
  dev.off()
}

plot_derivative_grid(res_high_unc, uncertainty_nu[1], "purple")
plot_derivative_grid(res_mid_unc, uncertainty_nu[2], "black")
plot_derivative_grid(res_low_unc, uncertainty_nu[3], "darkblue")

## ==============================================================================
# PHASE 7: SENSITIVITY ANALYSIS (DELTA & SIGMA) (Sequential) ====================
## ==============================================================================

base_nu <- uncertainty_nu[2]
delta_vals <- c(0.01, params$delta, 0.05) 
sigma_vals <- c(0.01, params$sigma, 1)   

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
    sigma_high = lapply(1:4, function(k) { p <- params; p$sigma <- sigma_vals[3]; simulate_paths(res_sigma_high[[k]], p, Npaths, Nt, dx, dt, W_sigma_high, X0_current) })
  )
}

### PLOTS: DELTA & SIGMA SENSITIVITY ====
for (v in 1:length(X0_vals)) {
  mc_delta_low  <- mc_sens_results[[v]]$delta_low
  mc_delta_high <- mc_sens_results[[v]]$delta_high
  mc_sigma_low  <- mc_sens_results[[v]]$sigma_low
  mc_sigma_high <- mc_sens_results[[v]]$sigma_high
  
  mc_mid <- mc_results[[v]]$mid
  
  # DELTA SENSITIVITY PLOT
  png(paste0("sens_delta_avg_investment_", X0_tags[v], ".png"), width = 800, height = 600)
  par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
  for (k in 1:4) {
    plot(time_seq, mc_delta_high[[k]]$S_avg, type = "l", col = "red", lwd = 2,
         main = paste("theta_bar =", th_vals[k]),
         xlab = "Time (Years)", ylab = "Investment",
         ylim = range(c(mc_delta_low[[k]]$S_avg, mc_mid[[k]]$S_avg, mc_delta_high[[k]]$S_avg)))
    
    lines(time_seq, mc_mid[[k]]$S_avg, col = "black", lwd = 2) 
    lines(time_seq, mc_delta_low[[k]]$S_avg, col = "blue", lwd = 2)
    
    if(k == 1) legend("bottomright", legend = c(paste0("High Delta (", delta_vals[3], ")"), 
                                                paste0("Base Delta (", delta_vals[2], ")"),
                                                paste0("Low Delta (", delta_vals[1], ")")), 
                      col = c("red", "black", "blue"), lwd = 2, bty = "n")
  }
  mtext(paste0("Sensitivity to Depreciation (Delta): Avg. Investment (X0 = ", X0_titles[v], ")"), outer = TRUE, cex = 1.2, font = 2)
  dev.off()
  
  # SIGMA SENSITIVITY PLOT
  png(paste0("sens_sigma_avg_investment_", X0_tags[v], ".png"), width = 800, height = 600)
  par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
  for (k in 1:4) {
    plot(time_seq, mc_sigma_high[[k]]$S_avg, type = "l", col = "purple", lwd = 2,
         main = paste("theta_bar =", th_vals[k]),
         xlab = "Time (Years)", ylab = "Investment",
         ylim = range(c(mc_sigma_low[[k]]$S_avg, mc_mid[[k]]$S_avg, mc_sigma_high[[k]]$S_avg)))
    
    lines(time_seq, mc_mid[[k]]$S_avg, col = "black", lwd = 2) 
    lines(time_seq, mc_sigma_low[[k]]$S_avg, col = "cyan", lwd = 2)
    
    if(k == 1) legend("bottomright", legend = c(paste0("High Sigma (", sigma_vals[3], ")"), 
                                                paste0("Base Sigma (", sigma_vals[2], ")"),
                                                paste0("Low Sigma (", sigma_vals[1], ")")), 
                      col = c("purple", "black", "cyan"), lwd = 2, bty = "n")
  }
  mtext(paste0("Sensitivity to Volatility (Sigma): Avg. Investment (X0 = ", X0_titles[v], ")"), outer = TRUE, cex = 1.2, font = 2)
  dev.off()
}

par(mfrow = c(1, 1))

## ==============================================================================
# PHASE 7.1: SENSITIVITY ANALYSIS (DELTA & SIGMA) Combined Plot ================
## ==============================================================================

### COMBINED PLOTS: DELTA & SIGMA SENSITIVITY Combined Panel Plots ====
time_seq <- (0:Nt) * dt

n_rows <- length(X0_vals)
n_cols <- length(th_vals)

a4_width  <- 11.69
a4_height <- 8.27
resolution <- 300 

#### DELTA SENSITIVITY PLOT ====
png("sens_delta_avg_investment_combined_nu_1.png", width = a4_width, height = a4_height, units = "in", res = resolution)
par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))

# Calculate limits by row (for each X0 across all theta_bar columns)
row_ylim <- list()
for (v in 1:n_rows) {
  vals <- c()
  for (k in 1:n_cols) {
    vals <- c(vals, 
              mc_sens_results[[v]]$delta_low[[k]]$S_avg, 
              mc_results[[v]]$mid[[k]]$S_avg, 
              mc_sens_results[[v]]$delta_high[[k]]$S_avg)
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
         xlab = "Time (Years)", ylab = "Investment",
         ylim = row_ylim[[v]])
    
    # Add grid
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    
    # Draw lines
    lines(time_seq, mc_delta_low[[k]]$S_avg, col = "blue", lwd = 2)
    lines(time_seq, mc_mid[[k]]$S_avg, col = "black", lwd = 2) 
    lines(time_seq, mc_delta_high[[k]]$S_avg, col = "red", lwd = 2)
    
    # Headers
    if (v == 1) {
      mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    if (k == 1) {
      mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}
mtext("Sensitivity to Depreciation (Delta): Avg. Investment", outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)

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

# Calculate limits by row (for each X0 across all theta_bar columns)
row_ylim <- list()
for (v in 1:n_rows) {
  vals <- c()
  for (k in 1:n_cols) {
    vals <- c(vals, 
              mc_sens_results[[v]]$sigma_low[[k]]$S_avg, 
              mc_results[[v]]$mid[[k]]$S_avg, 
              mc_sens_results[[v]]$sigma_high[[k]]$S_avg)
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
         xlab = "Time (Years)", ylab = "Investment",
         ylim = row_ylim[[v]])
    
    # Add grid
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    
    # Draw lines
    lines(time_seq, mc_sigma_low[[k]]$S_avg, col = "cyan", lwd = 2)
    lines(time_seq, mc_mid[[k]]$S_avg, col = "black", lwd = 2) 
    lines(time_seq, mc_sigma_high[[k]]$S_avg, col = "purple", lwd = 2)
    
    # Headers
    if (v == 1) {
      mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    if (k == 1) {
      mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}
mtext("Sensitivity to Volatility (Sigma): Avg. Investment", outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", legend = c(paste0("High Sigma (", sigma_vals[3], ")"), 
                            paste0("Base Sigma (", sigma_vals[2], ")"),
                            paste0("Low Sigma (", sigma_vals[1], ")")), 
       col = c("purple", "black", "cyan"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
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


#res_higher_tdisc <- lapply(th_vals, function(th) simulate_fixed_p(th, uncertainty_nu[1], params))
#params$r <- 0.01
#res_lower_r  <- lapply(th_vals, function(th) simulate_fixed_p(th, uncertainty_nu[2], params))
#params$r <- prev_r

### Saddle point processes ====

# 1. Panel Plot for Optimal Subsidy (theta)
png("optimal_subsidy_theta_new_terminal_discount.png", width = 800, height = 600)
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
       xlab = "Investment (x)", ylab = "Subsidy (theta)",
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
mtext("Optimal Subsidy (theta) across scenarios", outer = TRUE, cex = 1.2, font = 2)
dev.off()

# 2. Panel Plot for Optimal Production (q)
png("optimal_production_q_new_terminal_discount.png", width = 800, height = 600)
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
       xlab = "Investment (x)", ylab = "Production (q)",
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
mtext("Optimal Production (q) across scenarios", outer = TRUE, cex = 1.2, font = 2)
dev.off()


### Phase Plots (theta vs q and q vs theta) ====

# 1. Optimal Subsidy (theta) given Production (q)
png("optimal_subsidy_given_q_new_terminal_discount.png", width = 800, height = 600)
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
       xlab = "Optimal Production (q)", ylab = "Optimal Subsidy (theta)",
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
mtext("Optimal Subsidy (theta) given Production (q)", outer = TRUE, cex = 1.2, font = 2)
dev.off()

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

# Define A4 dimensions in inches (Landscape orientation)
a4_width  <- 11.69
a4_height <- 8.27
resolution <- 300 

### Panel Plot for Avg. Investment Paths ====
png("avg_investment_paths_combined_new_terminal_discount.png", width = a4_width, height = a4_height, units = "in", res = resolution)

par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))

col_ylim <- list()
for (k in 1:n_cols) {
  row_vals <- c()
  for(v in 1:n_rows) {
    row_vals <- c(row_vals, 
                  mc_results[[v]]$low[[k]]$S_avg, 
                  mc_results[[v]]$mid[[k]]$S_avg, 
                  mc_results[[v]]$high[[k]]$S_avg)
  }
  col_ylim[[k]] <- range(row_vals, na.rm = TRUE) # Use double brackets
}

for (v in 1:n_rows) {
  for (k in 1:n_cols) {
    mc_low  <- mc_results[[v]]$low
    mc_mid  <- mc_results[[v]]$mid
    mc_high <- mc_results[[v]]$high
    
    # Initialize empty plot for background
    plot(time_seq, mc_low[[k]]$S_avg, type = "n", 
         xlab = "Time (Years)", ylab = "Investment",
         ylim = col_ylim[[k]])
    
    # Add grid
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    
    # Draw data lines on top of grid
    lines(time_seq, mc_low[[k]]$S_avg, col = "darkgreen", lwd = 2)
    lines(time_seq, mc_mid[[k]]$S_avg, col = "black", lwd = 2)
    lines(time_seq, mc_high[[k]]$S_avg, col = "orange", lwd = 2)
    
    # Headers
    if (v == 1) { # Top row gets the column headers (theta_bar)
      mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    if (k == 1) { # First column gets the row headers (X0)
      mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}
mtext(paste0("Avg. Investment (X) over Time at ", 100*(1- params$terminal_discount), "% Terminal Depreciation"),
      outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                            paste0("Mid Unc (", uncertainty_nu[2], ")"),
                            paste0("Low Unc (", uncertainty_nu[3], ")")), 
       col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
dev.off()


### Panel Plot for Avg. Subsidies (theta) ====
png("avg_subsidies_paths_combined_new_terminal_discount.png", width = a4_width, height = a4_height, units = "in", res = resolution)
par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))

col_ylim <- list()
for (k in 1:n_cols) {
  row_vals <- c()
  for(v in 1:n_rows) {
    row_vals <- c(row_vals, 
                  mc_results[[v]]$low[[k]]$theta_avg, 
                  mc_results[[v]]$mid[[k]]$theta_avg, 
                  mc_results[[v]]$high[[k]]$theta_avg)
  }
  col_ylim[[k]] <- range(row_vals, na.rm = TRUE) # Use double brackets
}

for (v in 1:n_rows) {
  for (k in 1:n_cols) {
    mc_low  <- mc_results[[v]]$low
    mc_mid  <- mc_results[[v]]$mid
    mc_high <- mc_results[[v]]$high
    
    # Initialize empty plot
    plot(time_seq_Nt, mc_low[[k]]$theta_avg, type = "n", 
         xlab = "Time (Years)", ylab = "Subsidies",
         ylim = col_ylim[[k]])
    
    # Add grid and background reference lines
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    abline(h = th_vals[k], lty = 2, col = "darkgrey")
    abline(h = params$py, lty = 2, col = "darkgrey")
    
    # Draw data lines on top
    lines(time_seq_Nt, mc_low[[k]]$theta_avg, col = "darkgreen", lwd = 2)
    lines(time_seq_Nt, mc_mid[[k]]$theta_avg, col = "black", lwd = 2)
    lines(time_seq_Nt, mc_high[[k]]$theta_avg, col = "orange", lwd = 2)
    
    if (v == 1) {
      mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    if (k == 1) {
      mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}
mtext(paste0("Avg. Subsidies Paths at ", 100*(1- params$terminal_discount), "% Terminal Depreciation"), 
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
                  mc_results[[v]]$low[[k]]$q_avg, 
                  mc_results[[v]]$mid[[k]]$q_avg, 
                  mc_results[[v]]$high[[k]]$q_avg)
  }
  col_ylim[[k]] <- range(row_vals, na.rm = TRUE) # Use double brackets
}

for (v in 1:n_rows) {
  for (k in 1:n_cols) {
    mc_low  <- mc_results[[v]]$low
    mc_mid  <- mc_results[[v]]$mid
    mc_high <- mc_results[[v]]$high
    
    # Initialize empty plot
    plot(time_seq_Nt, mc_low[[k]]$q_avg, type = "n", 
         xlab = "Time (Years)", ylab = "Quantities",
         ylim = col_ylim[[k]])
    
    # Add grid and background reference lines
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    abline(h = params$q_max, lty = 2, col = "darkgrey")
    abline(h = 1.96, lty = 2, col = "darkgrey")
    
    # Draw data lines on top
    lines(time_seq_Nt, mc_low[[k]]$q_avg, col = "darkgreen", lwd = 2)
    lines(time_seq_Nt, mc_mid[[k]]$q_avg, col = "black", lwd = 2)
    lines(time_seq_Nt, mc_high[[k]]$q_avg, col = "orange", lwd = 2)
    
    if (v == 1) {
      mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    if (k == 1) {
      mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}
mtext(paste0("Avg. Quantities Paths at ", 100*(1- params$terminal_discount), "% Terminal Depreciation"),
      outer = TRUE, side = 3, cex = 1.5, font = 2, line = 1)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", legend = c(paste0("High Unc (", uncertainty_nu[1], ")"), 
                            paste0("Mid Unc (", uncertainty_nu[2], ")"),
                            paste0("Low Unc (", uncertainty_nu[3], ")")), 
       col = c("orange", "black", "darkgreen"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.2)
dev.off()


### Panel Plot for Avg. Investment Rates (gamma) ====
png("avg_investment_rates_paths_combined_new_terminal_discount.png", width = a4_width, height = a4_height, units = "in", res = resolution)
par(mfrow = c(n_rows, n_cols), oma = c(4, 2, 4, 1), mar = c(4, 6, 3, 1))

col_ylim <- list()
for (k in 1:n_cols) {
  row_vals <- c()
  for(v in 1:n_rows) {
    row_vals <- c(row_vals, 
                  mc_results[[v]]$low[[k]]$gamma_avg, 
                  mc_results[[v]]$mid[[k]]$gamma_avg, 
                  mc_results[[v]]$high[[k]]$gamma_avg)
  }
  col_ylim[[k]] <- range(row_vals, na.rm = TRUE) # Use double brackets
}

for (v in 1:n_rows) {
  for (k in 1:n_cols) {
    mc_low  <- mc_results[[v]]$low
    mc_mid  <- mc_results[[v]]$mid
    mc_high <- mc_results[[v]]$high
    
    # Initialize empty plot
    plot(time_seq_Nt, mc_low[[k]]$gamma_avg, type = "n", 
         xlab = "Time (Years)", ylab = "Investment Rates",
         ylim = col_ylim[[k]])
    
    # Add grid and background reference lines
    grid(col = "lightgray", lty = "dotted", lwd = 1)
    abline(h = 0, lty = 2, col = "darkgrey")
    
    # Draw data lines on top
    lines(time_seq_Nt, mc_low[[k]]$gamma_avg, col = "darkgreen", lwd = 2)
    lines(time_seq_Nt, mc_mid[[k]]$gamma_avg, col = "black", lwd = 2)
    lines(time_seq_Nt, mc_high[[k]]$gamma_avg, col = "orange", lwd = 2)
    
    if (v == 1) {
      mtext(paste0("theta_bar = ", th_vals[k]), side = 3, line = 1, cex = 1.1, font = 2)
    }
    if (k == 1) {
      mtext(paste0("X0 = ", X0_titles[v]), side = 2, line = 4.5, cex = 1.1, font = 2, las = 0)
    }
  }
}
mtext(paste0("Avg. Investment Rates Paths at ", 100*(1- params$terminal_discount), "% Terminal Depreciation"),
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
