# ==============================================================================
# Zebrafish PFAS mixture PBTK model
# Deterministic PBTK simulation and Shapley/Owen process attribution
# Yi Lu, Huaijun Xie, Xin He, Huimin Zhao
# ==============================================================================

library(deSolve)

# ------------------------------------------------------------------------------
# 1. General settings
# ------------------------------------------------------------------------------

PFAS <- c("PFOS", "PFHxS", "PFOA", "PFBS")
N_PFAS <- length(PFAS)
TIMES <- seq(0, 21, by = 0.1)

solver_method <- "lsodes"
solver_rtol <- 1e-8
solver_atol <- 1e-10

KT <- function(T, TR, TA) exp(TA / TR - TA / T)

trapz_auc <- function(time, conc) {
  ok <- is.finite(time) & is.finite(conc)
  time <- time[ok]; conc <- conc[ok]
  o <- order(time); time <- time[o]; conc <- conc[o]
  sum(diff(time) * (head(conc, -1) + tail(conc, -1)) / 2)
}

select_state <- function(single, mixture, use_mixture) {
  if (isTRUE(use_mixture)) mixture else single
}

subsets <- function(x) {
  out <- list(character(0))
  if (length(x)) for (k in seq_along(x)) out <- c(out, combn(x, k, simplify = FALSE))
  out
}

# ------------------------------------------------------------------------------
# 2. Experimentally constrained toxicokinetic parameters
# ------------------------------------------------------------------------------

PVF <- 0.685

fu_single <- c(PFOS = 0.0395, PFHxS = 0.0825, PFOA = 0.0895, PFBS = 0.1090)
fu_mix    <- c(PFOS = 0.0460, PFHxS = 0.1010, PFOA = 0.1090, PFBS = 0.1430)

free_single <- fu_single / PVF
free_mix    <- fu_mix / PVF

# Exact unrounded apparent non-fecal clearance values used in the model.
CLNF_single <- c(
  PFOS  = 0.010189099327616679,
  PFHxS = 0.022188040621311186,
  PFOA  = 0.025355212606299953,
  PFBS  = 0.043039623976775486
)

CLNF_mix <- c(
  PFOS  = 0.013225126377632823,
  PFHxS = 0.030702554769637560,
  PFOA  = 0.033963217304149420,
  PFBS  = 0.064352592551082580
)

CLNFu_single <- CLNF_single / free_single
CLNFu_mix    <- CLNF_mix / free_mix

CLLG_single <- c(
  PFOS  = 0.167000046588436,
  PFHxS = 0.133000420520232,
  PFOA  = 0.133000174217354,
  PFBS  = 0.0849978281561792
)

CLLG_mix <- c(
  PFOS  = 0.152000006320296,
  PFHxS = 0.111000199847804,
  PFOA  = 0.109999456168975,
  PFBS  = 0.0650002504244171
)

# Measured mixture-exposure water concentrations (µg/mL).
water_mix <- c(PFOS = 11.4e-3, PFHxS = 10.6e-3, PFOA = 9.8e-3, PFBS = 10.3e-3)

stopifnot(
  max(abs(CLNFu_single * free_single - CLNF_single)) < 1e-12,
  max(abs(CLNFu_mix * free_mix - CLNF_mix)) < 1e-12
)

# ------------------------------------------------------------------------------
# 3. Physiological and chemical parameters
# ------------------------------------------------------------------------------

phys <- c(
  BW_card_ref = 0.5,
  BW_VO2_ref = 0.4,
  TA = 3000,
  TR_Fcard = 299.5,
  TR_VO2 = 300.15,
  V_water = 1e12,
  F_card_ref = 41.328,
  VO2_ref = 9.84,
  OEE = 0.71,
  Sat = 0.90,
  Fbile = 0.0245,
  sc_blood = 0.0222,
  sc_liv = 0.019,
  sc_GB = 0.003,
  sc_gon = 0.04,
  sc_git = 0.129,
  sc_muscle = 0.2176,
  liv_frac = 0.01942,
  gon_frac = 0.01252,
  git_frac = 0.174,
  muscle_frac = 0.53689,
  delta_Kow = 3.1,
  TC = 26,
  BW = 0.8,
  exposure_end = 14
)

phys <- c(
  phys,
  sc_rest = 1 - sum(phys[c("sc_blood", "sc_liv", "sc_GB", "sc_gon", "sc_git", "sc_muscle")]),
  rest_frac = 1 - sum(phys[c("liv_frac", "gon_frac", "git_frac", "muscle_frac")])
)

chem <- list(
  ion_type = c(1, 1, 1, 1),
  pKa      = c(-3.27, -3.34, 0.50, -3.57),
  logKow   = c(5.77, 4.70, 4.65, 3.90),
  Plivb    = c(0.80, 0.75, 0.74, 0.81),
  Pgonb    = c(0.25, 0.62, 0.80, 0.90),
  Pgitb    = c(0.75, 0.70, 0.70, 0.45),
  Pmuscleb = c(0.05, 0.06, 0.08, 0.09),
  Prestb   = c(0.4865, 0.4773, 0.5101, 0.3852),
  Rloss    = rep(0.03, 4)
)

# ------------------------------------------------------------------------------
# 4. State vector
# ------------------------------------------------------------------------------

blocks <- c(
  "A_absorb_gill", "A_NF", "A_feces", "A_GB", "A_blood",
  "A_liv", "A_gon", "A_git", "A_muscle", "A_rest", "A_water"
)

block_index <- function(block) {
  k <- match(block, blocks)
  ((k - 1L) * N_PFAS + 1L):(k * N_PFAS)
}

state_names <- unlist(lapply(blocks, function(x) paste0(x, "_", PFAS)))

initial_state <- function(Cw) {
  y <- setNames(rep(0, length(blocks) * N_PFAS), state_names)
  y[block_index("A_water")] <- unname(Cw) * phys["V_water"]
  y
}

Y0 <- initial_state(water_mix)

end_exposure <- function(t, y, parms, ...) {
  y[block_index("A_water")] <- 0
  y
}

events <- list(func = end_exposure, time = unname(phys["exposure_end"]))

# ------------------------------------------------------------------------------
# 5. PBTK model
# ------------------------------------------------------------------------------

pbtk_model <- function(t, y, parms, tk) {
  
  # State variables
  A_GB     <- unname(y[block_index("A_GB")])
  A_blood  <- unname(y[block_index("A_blood")])
  A_liv    <- unname(y[block_index("A_liv")])
  A_gon    <- unname(y[block_index("A_gon")])
  A_git    <- unname(y[block_index("A_git")])
  A_muscle <- unname(y[block_index("A_muscle")])
  A_rest   <- unname(y[block_index("A_rest")])
  A_water  <- unname(y[block_index("A_water")])
  
  # Volumes
  Vblood  <- phys["sc_blood"]  * phys["BW"]
  Vliv    <- phys["sc_liv"]    * phys["BW"]
  VGB     <- phys["sc_GB"]     * phys["BW"]
  Vgon    <- phys["sc_gon"]    * phys["BW"]
  Vgit    <- phys["sc_git"]    * phys["BW"]
  Vmuscle <- phys["sc_muscle"] * phys["BW"]
  Vrest   <- phys["sc_rest"]   * phys["BW"]
  
  # Concentrations
  Cb      <- A_blood  / Vblood
  Cliv    <- A_liv    / Vliv
  CGB     <- A_GB     / VGB
  Cgon    <- A_gon    / Vgon
  Cgit    <- A_git    / Vgit
  Cmuscle <- A_muscle / Vmuscle
  Crest   <- A_rest   / Vrest
  
  # Cardiac output and tissue blood flows
  Tk <- phys["TC"] + 273.15
  Fcard <- phys["F_card_ref"] * KT(Tk, phys["TR_Fcard"], phys["TA"]) *
    (phys["BW"] / phys["BW_card_ref"])^(-0.1) * phys["BW"]
  
  Fliv    <- phys["liv_frac"]    * Fcard
  Fgon    <- phys["gon_frac"]    * Fcard
  Fgit    <- phys["git_frac"]    * Fcard
  Fmuscle <- phys["muscle_frac"] * Fcard
  Frest   <- phys["rest_frac"]   * Fcard
  
  # Gill uptake
  VO2 <- phys["VO2_ref"] * KT(Tk, phys["TR_VO2"], phys["TA"]) *
    (phys["BW"] / phys["BW_VO2_ref"])^(-0.1) * phys["BW"]
  
  Co2w <- ((-0.24 * phys["TC"] + 14.04) * phys["Sat"]) / 1e3
  Ywater <- VO2 / (phys["OEE"] * Co2w) / (1000^0.25 * phys["BW"]^0.75)
  
  logKow_ion <- chem$logKow - phys["delta_Kow"]
  fn <- 1 / (1 + 10^(chem$ion_type * (7.4 - chem$pKa)))
  Dow <- fn * 10^chem$logKow + (1 - fn) * 10^logKow_ion
  
  PBW <- 0.008 * 0.3 * Dow + 0.007 * 2 * Dow^0.94 +
    0.134 * 2.9 * Dow^0.63 + 0.851
  
  Yblood <- Fcard * PBW / (1000^0.25 * phys["BW"]^0.75)
  
  kx <- ((phys["BW"] / 1000)^0.75) /
    (2.8e-3 + 68 / Dow + 1 / Ywater + 1 / Yblood) * 1000
  
  F_bile <- phys["Fbile"] * phys["BW"]^0.75
  
  # Uptake and elimination fluxes
  J_abs <- kx * A_water / phys["V_water"]
  J_NF  <- tk$CLNFu * tk$free_NF * Cb
  J_fec <- chem$Rloss * F_bile * CGB
  
  # Tissue exchange
  Ceq_muscle <- Cmuscle / chem$Pmuscleb
  Ceq_rest   <- Crest / chem$Prestb
  
  dA_muscle <- tk$free_dist * Fmuscle * (Cb - Ceq_muscle)
  dA_rest   <- tk$free_dist * Frest   * (Cb - Ceq_rest)
  dA_gon    <- tk$free_dist * Fgon    * (Cb - Cgon / chem$Pgonb)
  
  dA_liv <- tk$free_dist * (
    Fliv * Cb +
      Fgit * Cgit / chem$Pgitb -
      (Fliv + Fgit) * Cliv / chem$Plivb
  ) - (tk$CLLG + F_bile) * Cliv
  
  dA_GB <- (tk$CLLG + F_bile) * Cliv - F_bile * CGB
  
  dA_git <- tk$free_dist * Fgit * (Cb - Cgit / chem$Pgitb) +
    (1 - chem$Rloss) * F_bile * CGB
  
  dA_blood <- J_abs - J_NF -
    tk$free_dist * Cb * (Fliv + Fgon + Fgit + Fmuscle + Frest) +
    tk$free_dist * (
      Fmuscle * Ceq_muscle +
        Frest * Ceq_rest +
        Fgon * Cgon / chem$Pgonb +
        (Fliv + Fgit) * Cliv / chem$Plivb
    )
  
  dA_water <- J_NF + J_fec - J_abs
  
  # Cumulative bookkeeping states retained from the original model.
  deriv <- c(
    J_abs, J_NF, J_fec,
    dA_GB, dA_blood, dA_liv, dA_gon, dA_git, dA_muscle, dA_rest, dA_water
  )
  
  list(deriv, C_blood = Cb)
}

# ------------------------------------------------------------------------------
# 6. Counterfactual simulations
# ------------------------------------------------------------------------------

run_scenario <- function(F1, F2, H, N) {
  
  tk <- list(
    free_dist = select_state(free_single, free_mix, F1),
    free_NF   = select_state(free_single, free_mix, F2),
    CLLG      = select_state(CLLG_single, CLLG_mix, H),
    CLNFu     = select_state(CLNFu_single, CLNFu_mix, N)
  )
  
  as.data.frame(ode(
    y = Y0,
    times = TIMES,
    func = pbtk_model,
    parms = NULL,
    tk = tk,
    method = solver_method,
    rtol = solver_rtol,
    atol = solver_atol,
    events = events
  ))
}

states <- expand.grid(
  F1 = c(FALSE, TRUE),
  F2 = c(FALSE, TRUE),
  H  = c(FALSE, TRUE),
  N  = c(FALSE, TRUE),
  KEEP.OUT.ATTRS = FALSE
)

AUC <- matrix(
  NA_real_,
  nrow = nrow(states),
  ncol = N_PFAS,
  dimnames = list(NULL, paste0("AUC_", PFAS))
)

blood_cols <- paste0("C_blood", seq_len(N_PFAS))

for (i in seq_len(nrow(states))) {
  s <- states[i, ]
  sim <- run_scenario(s$F1, s$F2, s$H, s$N)
  
  for (j in seq_len(N_PFAS)) {
    AUC[i, j] <- trapz_auc(sim$time, sim[[blood_cols[j]]])
  }
}

states <- cbind(states, AUC)

# ------------------------------------------------------------------------------
# 7. Counterfactual value retrieval
# ------------------------------------------------------------------------------

get_value <- function(tab, state, metric) {
  keep <- rep(TRUE, nrow(tab))
  for (nm in names(state)) keep <- keep & tab[[nm]] == state[[nm]]
  stopifnot(sum(keep) == 1L)
  tab[[metric]][keep]
}

# ------------------------------------------------------------------------------
# 8. First-level F/H/N Shapley attribution
# ------------------------------------------------------------------------------

shapley_FHN <- function(tab, metric) {
  
  players <- c("F", "H", "N")
  m <- length(players)
  
  baseline <- get_value(tab, setNames(rep(FALSE, m), players), metric)
  
  v <- function(S) {
    baseline - get_value(tab, setNames(players %in% S, players), metric)
  }
  
  phi <- setNames(numeric(m), players)
  
  for (player in players) {
    others <- setdiff(players, player)
    
    for (S in subsets(others)) {
      k <- length(S)
      w <- factorial(k) * factorial(m - k - 1L) / factorial(m)
      phi[player] <- phi[player] + w * (v(c(S, player)) - v(S))
    }
  }
  
  list(
    baseline = baseline,
    full = get_value(tab, setNames(rep(TRUE, m), players), metric),
    total = v(players),
    phi = phi
  )
}

# ------------------------------------------------------------------------------
# 9. Nested Owen decomposition of F into F1 and F2
# ------------------------------------------------------------------------------

owen_F <- function(tab, metric) {
  
  children <- c("F1", "F2")
  other_groups <- c("H", "N")
  
  baseline <- get_value(
    tab,
    c(F1 = FALSE, F2 = FALSE, H = FALSE, N = FALSE),
    metric
  )
  
  v <- function(state) baseline - get_value(tab, state, metric)
  omega <- setNames(numeric(2), children)
  
  for (child in children) {
    other_child <- setdiff(children, child)
    
    for (G in subsets(other_groups)) {
      g <- length(G)
      wg <- factorial(g) * factorial(3L - g - 1L) / factorial(3L)
      
      for (C in subsets(other_child)) {
        c_n <- length(C)
        wc <- factorial(c_n) * factorial(2L - c_n - 1L) / factorial(2L)
        
        before <- c(
          F1 = "F1" %in% C,
          F2 = "F2" %in% C,
          H  = "H" %in% G,
          N  = "N" %in% G
        )
        
        after <- before
        after[child] <- TRUE
        
        omega[child] <- omega[child] + wg * wc * (v(after) - v(before))
      }
    }
  }
  
  omega
}

# ------------------------------------------------------------------------------
# 10. Attribution results
# ------------------------------------------------------------------------------

# First-level F means F1 and F2 are switched jointly.
primary <- states[states$F1 == states$F2, ]
primary$F <- primary$F1
primary <- primary[, c("F", "H", "N", paste0("AUC_", PFAS))]

results <- vector("list", N_PFAS)

for (j in seq_along(PFAS)) {
  
  pfas <- PFAS[j]
  metric <- paste0("AUC_", pfas)
  
  shp <- shapley_FHN(primary, metric)
  own <- owen_F(states, metric)
  
  scale <- 100 / shp$baseline
  
  total <- shp$total * scale
  F  <- unname(shp$phi["F"] * scale)
  H  <- unname(shp$phi["H"] * scale)
  N  <- unname(shp$phi["N"] * scale)
  F1 <- unname(own["F1"] * scale)
  F2 <- unname(own["F2"] * scale)
  
  stopifnot(
    abs(F + H + N - total) < 1e-6,
    abs(F1 + F2 - F) < 1e-6
  )
  
  results[[j]] <- data.frame(
    PFAS = pfas,
    AUC_reduction_pct = total,
    F_pp = F,
    H_pp = H,
    N_pp = N,
    F1_pp = F1,
    F2_pp = F2
  )
}

results <- do.call(rbind, results)
row.names(results) <- NULL

print(results, digits = 4)
