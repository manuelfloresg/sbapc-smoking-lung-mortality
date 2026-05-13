BAPC_VERBOSE <- FALSE                         # verbose engine logging
# =============================================================
# 0A) PERILLAS GENERALES
# =============================================================

# ---------- Rutas de insumos ----------
BAPC_PROJECT_ROOT <- if (exists("BAPC_PATHS")) BAPC_PATHS$project_root else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
BAPC_DATA_ROOT    <- normalizePath(file.path(BAPC_PROJECT_ROOT, ".."), winslash = "/", mustWork = FALSE)

.first_existing_path <- function(candidates) {
  candidates <- unique(Filter(function(x) !is.null(x) && length(x) && nzchar(x), as.list(candidates)))
  hits <- Filter(file.exists, candidates)
  if (length(hits)) normalizePath(hits[[1]], winslash = "/", mustWork = FALSE) else NA_character_
}

.resolve_input_path <- function(option_name, env_name, candidates, required = TRUE) {
  opt <- tryCatch(getOption(option_name), error = function(e) NULL)
  env <- Sys.getenv(env_name, unset = "")
  cand <- c(opt, if (nzchar(env)) env else NULL, candidates)
  path <- .first_existing_path(cand)
  if (isTRUE(required) && (is.na(path) || !nzchar(path))) {
    stop(sprintf(
      "No pude resolver %s. Definilo con options(%s = '...') o con la variable de entorno %s.",
      option_name, option_name, env_name
    ))
  }
  path
}

PATH_MORT_CSV <- .resolve_input_path(
  option_name = "BAPC_PATH_MORT_CSV",
  env_name = "BAPC_PATH_MORT_CSV",
  candidates = c(
    file.path(BAPC_DATA_ROOT, "Mortalidad", "muertes_suavizadas_cancer.csv")
  )
)                                # csv de muertes suavizadas por causa
PATH_POP_DTA  <- .resolve_input_path(
  option_name = "BAPC_PATH_POP_DTA",
  env_name = "BAPC_PATH_POP_DTA",
  candidates = c(
    file.path(BAPC_DATA_ROOT, "Base de datos", "Proyecciones población", "poblacion_1950_2070_empalmada.dta"),
    file.path(BAPC_DATA_ROOT, "Base de datos", "Proyecciones poblacion", "poblacion_1950_2070_empalmada.dta")
  )
)  # proyecciones de población
PATH_PREV_DTA <- .resolve_input_path(
  option_name = "BAPC_PATH_PREV_DTA",
  env_name = "BAPC_PATH_PREV_DTA",
  candidates = c(
    file.path(BAPC_DATA_ROOT, "Base de datos", "base_completa.dta")
  )
)                                         # base de prevalencia
PATH_INC_CSV  <- .resolve_input_path(
  option_name = "BAPC_PATH_INC_CSV",
  env_name = "BAPC_PATH_INC_CSV",
  candidates = c(
    file.path(BAPC_DATA_ROOT, "Resultados", "incidencia_suavizada_1998_2022.csv")
  )
)                           # incidencia suavizada histórica
# ---------- Rutas de salida ----------
# Nota: los runners modernos crean subcarpetas específicas por corrida y escenario
# debajo de results/. Estos defaults quedan como fallback/legado para helpers sueltos.
BASE_RESULTS_DIR <- if (exists("BAPC_PATHS")) BAPC_PATHS$results else "results"
MASTER_RESULTS_DIR      <- file.path(BASE_RESULTS_DIR, "master")
AGGREGATE_RESULTS_DIR   <- file.path(BASE_RESULTS_DIR, "aggregate")
DIAGNOSTICS_RESULTS_DIR <- file.path(BASE_RESULTS_DIR, "diagnostics")
PLOTS_TOTAL_DIR         <- file.path(BASE_RESULTS_DIR, "plots_total")

# ---------- Ventanas APC (entrenamiento) ----------
PERIOD_M_MIN <- 1998; PERIOD_M_MAX <- 2022  # años usados para el ajuste de mortalidad
AGE_P_MIN    <- 20;   AGE_P_MAX    <- 65    # edades de prevalencia
AGE_I_MIN    <- 35;   AGE_I_MAX    <- 89    # edades de incidencia
AGE_M_MIN    <- 35;   AGE_M_MAX    <- 89    # edades usadas en el ajuste de mortalidad

# Maximum memory for former smoker risk reversal (years)
QUIT_A_I_MAX <- 50

# ---------- Horizonte ----------
PROJ_TO <- 2070                               # internal technical ceiling; published outputs are truncated by the endogenous horizon policy


# Endogenous projection-horizon policy (Lexis/border reliability)
HORIZON_SUPPORT_FLOOR_CREDIBLE <- 0.50
HORIZON_SUPPORT_FLOOR_CAUTION  <- 0.33
HORIZON_SUPPORT_FLOOR_MAX      <- 0.20
HORIZON_EDGE_SHARE_CREDIBLE    <- 0.30
HORIZON_EDGE_SHARE_CAUTION     <- 0.50
HORIZON_EDGE_SHARE_MAX         <- 0.80

# ---------- Regla post-fit para beta_P (canal PREV -> Incidencia) ----------
BETA_P_POSTFIT_RULE <- "floor0"

# ---------- Diagnóstico estructural en simulación ----------
SIM_DIAG_STRUCTURAL <- FALSE
SIM_DIAG_STRUCTURAL_PRINT <- FALSE
SIM_EVAL_PROFILE_DEFAULT <- "paper"
EMIT_PREV_DIAG_CONSOLE <- TRUE

# ---- Helpers técnicos ----
.safe_num <- function(x) {
  if (is.null(x)) return(numeric(0))
  suppressWarnings(as.numeric(x))
}

# ---- Helper: normalizar texto (minúsculas + sin tildes) ----
norm_txt <- function(x) {
  x |>
    stringr::str_trim() |>
    stringr::str_to_lower() |>
    stringi::stri_trans_general("Latin-ASCII")
}

# === Escenarios: niveles, etiquetas y colores (fuente única de verdad) ===
scenario_levels <- c("freeze","up1pc","down1pc","down3pc","quit")
scenario_labels <- c(
  quit    = "Full cessation",
  up1pc   = "Up 1% per year",
  down1pc = "Down 1% per year",
  freeze  = "Frozen at 2022",
  down3pc = "Down 3% per year"
)
scenario_colors <- c(
  "Frozen at 2022" = "#d62728",  # red
  "Up 1% per year" = "#9467bd",   # purple
  "Down 1% per year" = "#ff7f0e",   # orange
  "Down 3% per year" = "#2ca02c",   # green
  "Full cessation" = "#1f77b4"   # blue
)
# Orden deseado de la leyenda (nombres mostrados)
ESC_LEVELS <- unname(scenario_labels[scenario_levels])

# === Scenarios (English labels for paper-ready plots) ===
scenario_labels_en <- scenario_labels

# Keep the SAME hex colors, but keyed by English display labels
scenario_hex_by_code <- c(
  freeze  = "#d62728",  # red
  up1pc   = "#9467bd",  # purple
  down1pc = "#ff7f0e",  # orange
  down3pc = "#2ca02c",  # green
  quit    = "#1f77b4"   # blue
)

ESC_LEVELS_EN    <- unname(scenario_labels_en[scenario_levels])
scenario_colors_en <- setNames(unname(scenario_hex_by_code[names(scenario_labels_en)]), unname(scenario_labels_en[names(scenario_labels_en)]))

SCENARIOS_REAL_9SITES <- c("freeze","down1pc","down3pc","quit")
SCENARIOS_METHOD_LUNG <- c("freeze","up1pc","down1pc","quit")

normalize_prev_scenario_name <- function(x) {
  x <- as.character(x)[1]
  allowed <- c("freeze","up1pc","down1pc","down3pc","quit")
  if (!x %in% allowed) stop("Invalid prevalence scenario: ", x)
  x
}

make_prev_config <- function(
  scenario = NULL,
  axis = NULL,
  annual_rate = PREV_ANNUAL_RATE,
  annual_rate_down3 = PREV_ANNUAL_RATE_DOWN3,
  base_year = PREV_BASE_YEAR,
  backbone = PREV_BACKBONE,
  quit_mode = QUIT_MODE,
  quit_floor_sd = QUIT_FLOOR_SD,
  quit_floor_sd_M = NA_real_,
  quit_floor_sd_F = NA_real_,
  quit_half_life = QUIT_HALF_LIFE,
  quit_ramp_years = QUIT_RAMP_YEARS,
  prev_base_M = PREV_BASE_M,
  prev_base_F = PREV_BASE_F,
  prev_base_default = PREV_BASE_DEFAULT
) {
  if (is.null(scenario)) scenario <- PREV_SCENARIO
  if (is.null(axis)) axis <- PREV_SCENARIO_AXIS
  list(
    scenario = normalize_prev_scenario_name(as.character(scenario)[1]),
    axis = match.arg(as.character(axis)[1], c("period", "cohort_expo")),
    annual_rate = annual_rate,
    annual_rate_down3 = annual_rate_down3,
    base_year = base_year,
    backbone = match.arg(backbone, c("freeze", "forecast")),
    quit_mode = match.arg(quit_mode, c("decay", "floor", "ramp", "none")),
    quit_floor_sd = quit_floor_sd,
    quit_floor_sd_M = quit_floor_sd_M,
    quit_floor_sd_F = quit_floor_sd_F,
    quit_half_life = quit_half_life,
    quit_ramp_years = quit_ramp_years,
    prev_base_M = prev_base_M,
    prev_base_F = prev_base_F,
    prev_base_default = prev_base_default
  )
}

get_prev_config <- function(scenario = NULL, axis = NULL, ...) {
  make_prev_config(scenario = scenario, axis = axis, ...)
}

sex_labels_en <- c(M = "Males", F = "Females")

cause_labels_en <- c(
  oralphar  = "Oral cavity and pharynx",
  cervix    = "Cervix",
  stomach   = "Stomach",
  esophagus = "Esophagus",
  larynx    = "Larynx",
  pancreas  = "Pancreas",
  lung      = "Lung",
  kidney    = "Kidney",
  bladder   = "Bladder"
)

get_cause_label_en <- function(cause_id, fallback = NULL) {
  if (!is.null(cause_labels_en[[cause_id]])) return(unname(cause_labels_en[[cause_id]]))
  if (!is.null(fallback)) return(fallback)
  return(cause_id)
}



# ---- Registro de las 8 causas ----
causes <- tibble::tribble(
  ~cause_id, ~label,                    ~mort_regex,                ~inc_regex,                 ~AGE_P_MIN, ~AGE_P_MAX, ~AGE_I_MIN, ~AGE_I_MAX, ~AGE_M_MIN, ~AGE_M_MAX, ~L_I_MAX_YEARS, ~MORT_SHOCK_YEARS, ~DOWNWEIGHT_F,
  "oralphar", "Oral cavity and pharynx", "cavidad oral|faringe",     "cavidad oral|faringe",            20,        65,          35,         89,          35,         89,           3L,      integer(0),       integer(0),
  "cervix",   "Cervix",                  "cuello.*utero",            "cuello.*utero",                   20,        65,          25,         89,          25,         89,           3L,      integer(0),       integer(0),
  "stomach",  "Stomach",                 "estomago",                 "estomago",                        20,        65,          35,         89,          35,         89,           3L,      integer(0),       integer(0),
  "esophagus",  "Esophagus",               "esofago",                 "esofago",                         20,        65,          35,         89,          35,         89,           3L,      integer(0),       integer(0),
  "larynx",   "Larynx",                  "laringe",                  "laringe",                         20,        65,          35,         89,          35,         89,           3L,      integer(0),       integer(0),
  "pancreas", "Pancreas",                "pancreas",                 "pancreas",                        20,        65,          35,         89,          35,         89,           3L,      integer(0),       integer(0),
  "lung",     "Lung",                    "pulmon",                   "pulmon",                          20,        65,          35,         89,          35,         89,           3L,      integer(0),       integer(0),
  "kidney",   "Kidney",                  "rinon",                    "rinon",                           20,        65,          35,         89,          35,         89,           3L,      integer(0),       integer(0),
  "bladder",  "Bladder",                 "vejiga",                   "vejiga",                          20,        65,          35,         89,          35,         89,           3L,      integer(0),       integer(0)
)

# ---- FAP actual por cáncer y sexo (Sandoya & Bianco, 2011)----
FAP_BY_CAUSE_SEX <- tibble::tribble(
  ~cause_id,  ~sex, ~FAP_now,
  "oralphar", "M",  0.75,
  "oralphar", "F",  0.50,
  "cervix",   "M",  NA_real_,  # no aplica; se usará FAP_DEFAULT si queda NA
  "cervix",   "F",  0.13,
  "stomach",  "M",  0.28,
  "stomach",  "F",  0.13,
  "esophagus", "M",  0.72,
  "esophagus", "F",  0.61,
  "larynx",   "M",  0.84,
  "larynx",   "F",  0.75,
  "pancreas", "M",  0.25,
  "pancreas", "F",  0.27,
  "lung",     "M",  0.89,
  "lung",     "F",  0.74,
  "kidney",   "M",  0.39,
  "kidney",   "F",  0.06,
  "bladder",  "M",  0.48,
  "bladder",  "F",  0.31
)

# Valor por defecto si falta una FAP o es NA explícito (podés cambiarlo)
FAP_DEFAULT <- 0.30

# ---- Probabilidades post-diagnóstico para el canal INC -> MORT ----
# Base homogénea proxy construida desde REDECAN (España, 2013--2017),
# derivada de supervivencia neta a 1, 3 y 5 años.
# Las probabilidades son incondicionales sobre el total de diagnosticados.
MORT_POSTDX_DEATH_TABLE <- tibble::tribble(
  ~cause_id,   ~sex, ~p_0_1, ~p_1_3, ~p_3_5, ~p_le_5,
  "oralphar",  "M",  0.298,  0.194,  0.067,  0.559,
  "oralphar",  "F",  0.202,  0.138,  0.066,  0.406,
  "esophagus", "M",  0.565,  0.242,  0.049,  0.856,
  "esophagus", "F",  0.566,  0.220,  0.039,  0.825,
  "stomach",   "M",  0.464,  0.208,  0.060,  0.732,
  "stomach",   "F",  0.426,  0.172,  0.036,  0.634,
  "pancreas",  "M",  0.672,  0.183,  0.042,  0.897,
  "pancreas",  "F",  0.633,  0.205,  0.041,  0.879,
  "larynx",    "M",  0.143,  0.160,  0.071,  0.374,
  "larynx",    "F",  0.151,  0.130,  0.071,  0.352,
  "lung",      "M",  0.588,  0.196,  0.053,  0.837,
  "lung",      "F",  0.500,  0.202,  0.063,  0.765,
  "cervix",    "F",  0.141,  0.135,  0.052,  0.328,
  "kidney",    "M",  0.174,  0.091,  0.056,  0.321,
  "kidney",    "F",  0.196,  0.077,  0.041,  0.314,
  "bladder",   "M",  0.113,  0.085,  0.049,  0.247,
  "bladder",   "F",  0.137,  0.075,  0.013,  0.225,
  "sim",       "M",  0.588,  0.196,  0.053,  0.837,
  "sim",       "F",  0.500,  0.202,  0.063,  0.765,
  "spec_linear", "M",  0.588,  0.196,  0.053,  0.837,
  "spec_linear", "F",  0.500,  0.202,  0.063,  0.765,
  "misspec_tanh", "M",  0.588,  0.196,  0.053,  0.837,
  "misspec_tanh", "F",  0.500,  0.202,  0.063,  0.765
)

MORT_POSTDX_KERNEL_MODE   <- "midyear_uniform" # diagnóstico medio a mitad de año + densidad uniforme dentro de cada tramo
MORT_POSTDX_USE_AGE_SHIFT <- TRUE              # una muerte a edad a en t mira incidencia a-k en t-k
MORT_POSTDX_CLIP_MIN_AGE  <- TRUE              # si a-k cae fuera del rango, usar la mínima edad disponible
MORT_POSTDX_CLIP_MIN_YEAR <- TRUE              # si t-k cae antes del inicio, usar el primer año disponible

# =============================================================
# 0B) PERILLAS APC COMUNES (modelos y priors)
# =============================================================

# -------- Modelos APC por componente (rw1 / rw2) --------
PREV_AGE_MODEL <- "rw2"; PREV_PER_MODEL <- "rw1"; PREV_COH_MODEL <- "rw1"  # "rw1" más diente
INC_AGE_MODEL  <- "rw2"; INC_PER_MODEL  <- "rw1"; INC_COH_MODEL  <- "rw1"  # "rw1" más diente
MORT_AGE_MODEL <- "rw2"; MORT_PER_MODEL <- "rw1"; MORT_COH_MODEL <- "rw1"  # "rw1" más diente

# ---------- Hiperparámetros PC-prior (P(sigma > u) = alpha) ----------
# Prevalencia
PREV_AGE_PC_U <- 0.35; PREV_AGE_PC_A <- 0.01   # edad en P (↓u/α = más suavizado; atenúa “mariposa” de bordes)
PREV_PER_PC_U <- 1; PREV_PER_PC_A <- 0.01    # período en P (↓ = curva más lisa)
PREV_COH_PC_U <- 1; PREV_COH_PC_A <- 0.01    # cohorte en P (↓ = más liso; sube señal si ↑)
# Incidencia
INC_AGE_PC_U  <- 0.35; INC_AGE_PC_A  <- 0.01    # edad en I
INC_PER_PC_U  <- 1;  INC_PER_PC_A  <- 0.01    # período en I (↑u/α = más diente; ↓ = más liso)
INC_COH_PC_U  <- 1; INC_COH_PC_A  <- 0.01    # cohorte en I
# Mortalidad
MORT_AGE_PC_U <- 0.35; MORT_AGE_PC_A <- 0.01    # edad en M
MORT_PER_PC_U <- 1; MORT_PER_PC_A <- 0.01    # período en M (menor u/α = más suavizado de tendencia)
MORT_COH_PC_U <- 1; MORT_COH_PC_A <- 0.01    # cohorte en M
# Shock opcional de período en M
MORT_SHOCK_PC_U <- 1.0; MORT_SHOCK_PC_A <- 0.01 # prior para iid de “shock” (↑u/α = más variabilidad permitida)

# -------- Cohorte ponderada (anti-esquinas) --------
USE_WEIGHTED_COHORT <- FALSE                    # cohort weighting apagado por consistencia entre niveles

# -------- Edge-weighting histórico (edad + período) --------
EDGE_WEIGHTING_ON <- FALSE
EDGE_WEIGHT_GEOMETRY <- "additive_mean"
EDGE_WEIGHT_K_AGE <- 5L
EDGE_WEIGHT_K_PERIOD <- 5L
EDGE_WEIGHT_STRENGTH <- 0.5
EDGE_WEIGHT_MIN <- 0.25

# =============================================================
# 0C) PREVALENCIA (nivel P)
# =============================================================

# --- Proyección de γ^P (para construir el índice P→I) ---
COHORT_FC_METHOD <- "damped_trend"             # "freeze","arima","trend","damped_trend"
COHORT_FC_WINDOW <- 5L                         # ventana reciente para pendiente local
COHORT_FC_DAMPING <- 0.8                       # amortiguación geométrica de la pendiente (0<δ<=1)

GAMMAP_METHOD <- COHORT_FC_METHOD              # método para extender γ^P a futuro
TREND_TYPE    <- "level"                       # si método=="trend": "level" o "trend"
PREV_BACKBONE <- "freeze"                      # "freeze" o "forecast"
PREV_TREND_DEGREE <- 1                         # 0: sin tendencia explícita en P; 1: lineal; 2: cuadrática
PREV_TREND_PRIOR_SD <- 1                       # sd del prior para coef(s) de tendencia en P

# --- Escenarios de PREVALENCIA a futuro (aplican desde 2023) ---
PREV_SCENARIO      <- "freeze"                  # opciones: "down1pc","freeze","down3pc","quit"
PREV_ANNUAL_RATE   <- 0.01                      # 1% anual 
PREV_ANNUAL_RATE_DOWN3 <- 0.03                   # 3% anual — escenario "down3pc"
PREV_BASE_YEAR     <- 2022                      # a partir de periodo > 2022 aplica el escenario
PREV_SCENARIO_AXIS <- "period"                  # opciones: "period" o "cohort_expo"

# Prevalencia base (año PREV_BASE_YEAR), por sexo
PREV_BASE_M <- 0.235   # <-- reemplazar por Uruguay 2022
PREV_BASE_F <- 0.135   # <-- reemplazar por Uruguay 2022
PREV_BASE_DEFAULT <- 0.20

# Escenario "quit": modo de transición hacia el piso
QUIT_MODE       <- "decay"   # "floor" (actual), "decay" (exponencial), "ramp" (lineal)
QUIT_FLOOR_SD   <- -2.5        # piso en z_prev (en sigmas), sigue valiendo para todos los modos
QUIT_FLOOR_SD_M <- NA_real_    # opcional: piso específico para M
QUIT_FLOOR_SD_F <- NA_real_    # opcional: piso específico para F
QUIT_HALF_LIFE  <- 2         # solo si QUIT_MODE="decay": semivida (años) hacia el piso
QUIT_RAMP_YEARS <- 3         # solo si QUIT_MODE="ramp": años hasta llegar al piso

# --- Enganche P → I (índice a incidencia) ---
W_I  <- 0.7                                      # peso del índice en I (↑ = canal P→I más fuerte)
PREV_W_I    <- W_I                             # alias para consistencia con código previo

# ---- Nuevo canal PREV -> INC (v1; infraestructura, sin activar aún el rewire del motor)
PREV_INC_CHANNEL_MODE <- "stock_former"          # ruta principal nueva PREV -> INC
PREV_INC_MAX_QUIT_YEARS <- NA_integer_           # si NA, usar máximo disponible de la risk-reversion schedule
PREV_BACKCAST_MODE <- "freeze_first_period"     # backcast del efecto de período en PREV
PREV_BACKCAST_COHORT_MODE <- "freeze_oldest"    # cohortes PREV no observadas hacia atrás
PREV_POST65_MODE <- "carry_states"              # edades > AGE_P_MAX: transportar estados por cohorte
BETA_MODE <- "fixed_rr_offset"                 # "estimate","prior_ols","offset","fixed_rr_offset"

# --- Riesgos relativos de incidencia atribuibles a tabaquismo ---
# Fuente preferida para la rama fixed_rr_offset: Pichon-Riviere et al. (2013),
# con valores diferenciados por sexo para los 9 sitios actualmente modelados.
INC_RR_DEFAULT <- 4.0
INC_RR_BY_CAUSE <- c(
  oralphar = 4.5,
  cervix   = 1.59,
  stomach  = 1.6,
  esophagus= 5.0,
  larynx   = 10.0,
  pancreas = 2.0,
  lung     = 20.0,
  kidney   = 2.0,
  bladder  = 3.0,
  sim      = 6.0,
  sim_spec_linear = 6.0,
  spec_linear = 20.0,
  misspec_tanh = 20.0
)
INC_RR_BY_CAUSE_SEX <- list(
  oralphar = c(M = 10.89, F = 5.08),
  esophagus = c(M = 6.76, F = 7.75),
  stomach = c(M = 1.96, F = 1.36),
  pancreas = c(M = 2.31, F = 2.25),
  larynx = c(M = 14.60, F = 13.02),
  lung = c(M = 23.26, F = 12.69),
  cervix = c(M = NA_real_, F = 1.59),
  kidney = c(M = 2.72, F = 1.29),
  bladder = c(M = 3.27, F = 2.22),
  sim = c(M = 6.0, F = 6.0),
  sim_spec_linear = c(M = 6.0, F = 6.0),
  spec_linear = c(M = 23.26, F = 12.69),
  misspec_tanh = c(M = 23.26, F = 12.69)
)
SIM_RR_I_DEFAULT <- 4.0

# --- Link Uncertainty (Jitter) ---
# Set SD > 0 to sample Relative Risks from a Log-Normal distribution before each run.
# This propagates literature uncertainty into the modeling pipeline.
BAPC_JITTER_RR_SD <- 0.05   # Default 5% uncertainty on the log-scale

.normalize_rr_lookup_sex <- function(sex, default = NA_character_) {
  sx <- suppressWarnings(as.character(sex)[1])
  if (!length(sx) || is.na(sx) || !nzchar(sx)) return(default)
  sx <- toupper(substr(trimws(sx), 1L, 1L))
  if (!sx %in% c("M", "F")) return(default)
  sx
}

get_inc_rr_by_cause <- function(cause_id, default = INC_RR_DEFAULT) {
  cause_id <- suppressWarnings(as.character(cause_id)[1])
  if (!length(cause_id) || is.na(cause_id) || !nzchar(cause_id)) cause_id <- NA_character_
  rr <- suppressWarnings(as.numeric(INC_RR_BY_CAUSE[cause_id]))[1]
  if ((!is.finite(rr) || rr <= 1) && !is.na(cause_id) && !is.null(INC_RR_BY_CAUSE_SEX[[cause_id]])) {
    rr <- suppressWarnings(max(as.numeric(INC_RR_BY_CAUSE_SEX[[cause_id]]), na.rm = TRUE))
  }
  if (!is.finite(rr) || rr <= 1) rr <- suppressWarnings(as.numeric(default))[1]
  if (!is.finite(rr) || rr <= 1) rr <- 2.0
  rr
}

get_inc_rr_by_cause_sex <- function(cause_id, sex, default = INC_RR_DEFAULT) {
  cause_id <- suppressWarnings(as.character(cause_id)[1])
  if (!length(cause_id) || is.na(cause_id) || !nzchar(cause_id)) cause_id <- NA_character_
  sx <- .normalize_rr_lookup_sex(sex)
  rr <- NA_real_
  rr_map <- tryCatch(INC_RR_BY_CAUSE_SEX[[cause_id]], error = function(e) NULL)
  if (!is.null(rr_map) && !is.na(sx)) {
    rr <- suppressWarnings(as.numeric(rr_map[sx]))[1]
  }
  if (!is.finite(rr) || rr <= 1) {
    rr <- get_inc_rr_by_cause(cause_id = cause_id, default = default)
  }
  if (!is.finite(rr) || rr <= 1) rr <- suppressWarnings(as.numeric(default))[1]
  if (!is.finite(rr) || rr <= 1) rr <- 2.0
  rr
}

# --- Efecto de prevalencia en I (I|P) ---
SD_THETA_IP <- 0.4                               # sd del prior para el coeficiente de z_prev (↓ = empuja hacia 0; ↑ = más libre)

# =============================================================
# 0D) INCIDENCIA (nivel I)
# =============================================================

# --- Rugosidad extra en PERÍODO de I ---
INC_PER_EXTRA_IID <- FALSE                     # TRUE agrega iid de período (más diente local)
INC_PER_IID_PC_U  <- 1.2                       # prior u del iid (↑u/α = más variación permitida)
INC_PER_IID_PC_A  <- 0.1                       # prior α del iid

# --- Tendencias en I (para escenarios; el ajuste histórico puede ir sin tendencia) ---
INC_TREND_DEGREE   <- 0                          # 0: sin tendencia; 1: lineal
INC_TREND_ON       <- "none"                     # dónde actúa la tendencia: "period", "cohort" o "none"
INC_TREND_SCENARIO <- "freeze"                   # "freeze" (se congela tras 2022), "continue" (sigue pendiente), "delta" (ajuste exógeno)
DELTA_INC          <- 0                          # si "delta": cambio anual adicional (p.ej. -0.02 = -2%/año)
SD_BETA_I          <- 1                          # sd del prior de coef(s) de tendencia en I (↓ = más ancla; ↑ = más libertad)

# --- Pronóstico de componente APC  ---
INC_COEF_FC_TARGET <- "cohort"                   # "none","cohort","period"
INC_COEF_FC_METHOD <- COHORT_FC_METHOD           # método para el pronóstico APC explícito
INC_COEF_FC_LOCK_MODE <- "none"                  # legado: bloquea el offset completo coef_fc_offset_I
INC_COEF_FC_LOCK_VALUE <- 0                      # valor fijo usado cuando INC_COEF_FC_LOCK_MODE == "fixed"
INC_COEF_FC_RECENTER_LOCK_MODE <- "none"         # "none","zero","fixed": bloquea sólo el componente de recenterización
INC_COEF_FC_RECENTER_LOCK_VALUE <- 0             # valor fijo usado cuando INC_COEF_FC_RECENTER_LOCK_MODE == "fixed"
INC_COEF_FC_POSTHOC_LOCK_MODE <- "none"          # "none","zero","fixed": ajusta post hoc el offset futuro manteniendo fija la señal APC ya estimada
INC_COEF_FC_POSTHOC_LOCK_VALUE <- 0              # valor fijo usado cuando INC_COEF_FC_POSTHOC_LOCK_MODE == "fixed"

# =============================================================
# 0E) MORTALIDAD (nivel M)
# =============================================================

# --- Proyección futura de cohorte en M ---
MORT_COHORT_FC_ON <- TRUE                      # TRUE: reemplaza el clamp/freeze futuro por ajuste post hoc de cohorte

# --- Tendencia de “progreso técnico” en M (aplica a BAPC y anclado) ---
MORT_TREND_DEGREE   <- 1                       # 0: sin tendencia; 1: lineal
MORT_TREND_SCENARIO <- "freeze"                # "freeze","continue","delta" para las ramas informadas (residuo APC)
MORT_BAPC_TREND_SCENARIO <- "continue"         # convención para la rama autónoma de mortalidad (línea punteada)
DELTA_TECH          <- -0.02                   # si "delta": mejora anual adicional (negativo = reduce mortalidad)
SD_BETA_M           <- 2                       # sd del prior para coef(s) de tendencia en M (↓ = más ancla)
MORT_TREND_PRIOR_MEAN_M <- -0.01               # media prior para H
MORT_TREND_PRIOR_SD_M   <-  0.08               # sd prior (más chico = más ancla)
MORT_TREND_PRIOR_MEAN_F <- -0.01               # media prior para M
MORT_TREND_PRIOR_SD_F   <-  0.005              # sd prior (ligeramente más fuerte en F)

# --- Intervenciones en años ---
MORT_PERIOD_SHOCK_YEARS  <- integer(0)         # años con shock de período (p.ej. c(2020) para COVID/registro) / integer(0) no incluye shocks
MORT_DOWNWEIGHT_YEARS_F  <- integer(0)         # años a penalizar p.e. c(2014:2018)
MORT_DOWNWEIGHT_WEIGHT_F <- 0.3                # peso relativo (1 = sin cambio)

# --- Enlace Incidencia → Mortalidad ---
MORT_I_LINK_MODE   <- "external_kernel"        # enlace externo distribuido por años post-diagnóstico
SD_BETA_I_MORT     <- 0.1                      # legado: ya no se usa en el camino principal

# --- Suavización del enganche historia/proyección ---
ANCHOR_PSEUDO_W    <- 0.01                     # peso de las pseudo-observaciones de continuidad (0–1). ↑ valor = “enganche” más rígido entre 2022 y 2023. 0 apaga el truco de continuidad.
BRIDGE_INC_YEARS   <- 0                        # legado: el nuevo offset usa directamente la incidencia proyectada, sin bridge adicional
MORT_ANNUAL_BRIDGE  <- FALSE                    # si FALSE, no reescala las series anuales de mortalidad para empalmar 2022; usar versiones raw como canónicas

# --- Legado del rezago I -> M (ya no calibra el camino principal) ---
L_I_DEFAULT       <- 20L
L_I_MAX_YEARS     <- 3L
L_I_EXCLUDE_YEARS <- c(2020, 2021)
DA_I               <- 0

# --- Cohorte residual en el anclado (suavidad del “enganche”) ---
SD_COHORT_RESID <- 0.02                        # sd del residuo de cohorte en el anclado (↓ = anclaje más fuerte/continuo)

if (!exists("SD_BETA", inherits = TRUE))   SD_BETA   <- 1   # probá 0.7–1.0 para afinar
if (!exists("SD_COHFIX", inherits = TRUE)) SD_COHFIX <- 0.01  # probá 0.005–0.02

# =========================
# Backtesting / LOO deltas
# =========================
BT_ENABLE        <- TRUE
BT_HOLDOUT_YEARS <- 8L   # ej.: deja fuera los últimos 10 años observados


# =============================================================
# 0F) PERILLAS GLOBALES / LEGADO
# =============================================================

SD_BETA_FIXED <- .5                            # sd prior genérica para coeficientes fijos cuando se estiman (↓ = más ancla)
# Legacy aliases (kept for backward compatibility during transition)
if (exists("INC_INCLUDE_TREND"))  INC_TREND_DEGREE  <- if (isTRUE(INC_INCLUDE_TREND)) 1 else 0
if (exists("MORT_INCLUDE_TREND")) MORT_TREND_DEGREE <- if (isTRUE(MORT_INCLUDE_TREND)) 1 else 0
if (exists("MORT_TREND_SCENARIO_BAPC")) MORT_BAPC_TREND_SCENARIO <- MORT_TREND_SCENARIO_BAPC
EXPOSE_SELECTED_LAGS_TO_ENV <- TRUE
