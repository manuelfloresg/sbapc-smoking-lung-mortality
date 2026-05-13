# 3) Carga de datos y preprocesamiento
# =============================================================

# ---------------------------
# 3A) Mortalidad y población
# ---------------------------
mort_hist <- readr::read_csv(PATH_MORT_CSV, show_col_types = FALSE) %>%
  mutate(
    causa2_norm = norm_txt(causa2),       # <= normalizado para regex robusto
    period = as.integer(año),
    age    = as.integer(edad),
    sex    = factor(ifelse(mujer == 1, "F", "M"), levels = c("M","F")),
    deaths = as.numeric(muertes_suavizadas),
    cohort = period - age,
    cause  = causa2_norm                  # <= *** NUEVO: garantizamos 'cause' ***
  ) %>%
  select(period, age, sex, cause, deaths, cohort) %>%   # <= usamos 'cause'
  arrange(sex, period, age)

# === 3A) Mortalidad y población — CARGA DE POBLACIÓN (REEMPLAZO COMPLETO) ===
pop_raw <- haven::read_dta(PATH_POP_DTA)

nm <- names(pop_raw)
year_col <- intersect(c("period","year","anio","ano","año"), nm)[1]
age_col  <- intersect(c("edad","age"), nm)[1]
sex_col  <- intersect(c("sexo","sex","mujer"), nm)[1]
pop_col  <- intersect(c("poblacion","population","exposure","pob","N","n"), nm)[1]

if (is.na(year_col)) stop("No encuentro columna de año en la base de población.")
if (is.na(age_col))  stop("No encuentro columna de edad en la base de población.")
if (is.na(sex_col))  stop("No encuentro columna de sexo en la base de población.")
if (is.na(pop_col))  stop("No encuentro columna de población/exposición en la base de población.")

# Año máximo disponible en el archivo original (aún con su nombre original)
MAX_POP_YEAR <- max(pop_raw[[year_col]], na.rm = TRUE)

# Si PROJ_TO todavía no existe, usamos MAX_POP_YEAR como tope
proj_to_lim <- if (exists("PROJ_TO", inherits = FALSE)) min(get("PROJ_TO"), MAX_POP_YEAR) else MAX_POP_YEAR

pop_all <- pop_raw %>%
  dplyr::transmute(
    period   = as.integer(.data[[year_col]]),
    age      = as.integer(.data[[age_col]]),
    sex      = {
      if (identical(sex_col, "mujer")) {
        factor(ifelse(as.numeric(.data[[sex_col]]) == 1, "F", "M"), levels = c("M","F"))
      } else {
        factor(ifelse(as.numeric(.data[[sex_col]]) == 2, "F", "M"), levels = c("M","F"))
      }
    },
    exposure = as.numeric(.data[[pop_col]])
  ) %>%
  dplyr::filter(dplyr::between(period, 1952L, proj_to_lim)) %>%
  dplyr::arrange(sex, period, age)

if (!"exposure.pop" %in% names(pop_all)) {
  pop_all <- dplyr::mutate(pop_all, exposure.pop = exposure)
}

pop_future <- pop_all %>%
  dplyr::filter(period >= PERIOD_M_MAX + 1L) %>%
  dplyr::select(age, period, sex, exposure)


# Loader general de mortalidad por causa (usa regex normalizado)
load_mortality_by_cause <- function(path_csv, cause_regex, pop_all_tbl, cause_id = NA_character_) {
  readr::read_csv(path_csv, show_col_types = FALSE) %>%
    dplyr::mutate(causa2_norm = norm_txt(causa2)) %>%
    dplyr::filter(stringr::str_detect(causa2_norm, cause_regex)) %>%
    dplyr::transmute(
      period = as.integer(año),
      age    = as.integer(edad),
      sex    = factor(ifelse(mujer == 1, "F", "M"), levels = c("M","F")),
      cause  = ifelse(is.na(cause_id), cause_regex, cause_id),
      deaths = as.numeric(muertes_suavizadas)
    ) %>%
    dplyr::mutate(cohort = period - age) %>%
    dplyr::left_join(pop_all_tbl, by = c("period","age","sex")) %>%
    { if (any(is.na(.$exposure))) stop("Exposures faltantes en mort_hist para ", cause_regex); . } %>%
    dplyr::arrange(sex, period, age)
}

# ---------------------------
# 3B) Prevalencia (constructor desde micro)
# ---------------------------
build_prev_from_micro_df <- function(micro_df,
                                     sex_sel = c("M","F"),
                                     period_min = 1998, period_max = 2022,
                                     age_min = 15, age_max = 65,
                                     min_neff = 5) {
  sex_sel <- match.arg(sex_sel)
  req <- c("period","age","cohort","sex","fuma","w","d_act","d_12m","d_30d")
  miss <- setdiff(req, names(micro_df))
  if (length(miss)) stop("Faltan columnas en micro_df: ", paste(miss, collapse = ", "))

  micro <- micro_df %>%
    transmute(
      period = as.integer(period),
      age    = as.integer(age),
      cohort = as.integer(cohort),
      sex    = factor(as.character(sex), levels = c("M","F")),
      fuma   = as.numeric(fuma),
      w      = as.numeric(w),
      d_act  = as.integer(d_act),
      d_12m  = as.integer(d_12m),
      d_30d  = as.integer(d_30d)
    ) %>%
    filter(period >= period_min, period <= period_max,
           age >= age_min, age <= age_max,
           sex == sex_sel)

  inst_count_ok <- with(micro, d_act + d_12m + d_30d)
  if (!all(inst_count_ok %in% c(0, 1))) stop("Más de una dummy de instrumento = 1")

  micro <- micro %>%
    filter(inst_count_ok == 1) %>%
    mutate(inst = case_when(
      d_act == 1 ~ "act",
      d_12m == 1 ~ "12m",
      d_30d == 1 ~ "30d",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(inst)) %>%
    mutate(inst = factor(inst, levels = c("act", "12m", "30d")))

  micro %>%
    group_by(age, period, cohort, sex, inst) %>%
    summarise(
      w_sum   = sum(w, na.rm = TRUE),
      w2_sum  = sum(w^2, na.rm = TRUE),
      smoke_w = sum(w * fuma, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      p_hat  = ifelse(w_sum > 0, smoke_w / w_sum, NA_real_),
      neff_k = ifelse(w_sum > 0, (w_sum^2) / pmax(w2_sum, 1e-12), 0),
      neff   = pmax(0, round(neff_k)),
      y_eff  = pmin(neff, pmax(0, round(p_hat * neff)))
    ) %>%
    filter(neff >= min_neff) %>%
    select(age, period, cohort, sex, inst, y_eff, neff) %>%
    arrange(period, age, inst)
}

build_prev_from_micro <- function(path_dta,
                                  sex_sel = c("M","F"),
                                  period_min = 1998, period_max = 2022,
                                  age_min = 15, age_max = 65,
                                  min_neff = 5) {
  sex_sel <- match.arg(sex_sel)
  micro <- read_dta(path_dta) %>%
    transmute(
      period = as.integer(año),
      age    = as.integer(edad),
      cohort = as.integer(coh),
      sex    = factor(ifelse(mujer == 1, "F", "M"), levels = c("M", "F")),
      fuma   = as.numeric(fuma),
      w      = as.numeric(expansor),
      d_act  = as.integer(d_act),
      d_12m  = as.integer(d_12m),
      d_30d  = as.integer(d_30d)
    )

  build_prev_from_micro_df(
    micro_df   = micro,
    sex_sel    = sex_sel,
    period_min = period_min,
    period_max = period_max,
    age_min    = age_min,
    age_max    = age_max,
    min_neff   = min_neff
  )
}
message("build_prev_from_micro() and build_prev_from_micro_df() loaded.")

diag_prev <- function(df) {
  message("Filas: ", nrow(df)); stopifnot(nrow(df) > 0)
  df <- df %>% mutate(neff = as.integer(neff), y_eff = as.integer(y_eff))
  if (any(is.na(df$neff) | is.na(df$y_eff))) stop("NA en neff o y_eff")
  if (any(df$neff < 1)) stop("neff < 1")
  if (any(df$y_eff < 0 | df$y_eff > df$neff)) stop("y_eff fuera de [0, neff]")
  n_age <- n_distinct(df$age); n_per <- n_distinct(df$period); n_coh <- n_distinct(df$cohort)
  message("unique ages=", n_age, " periods=", n_per, " cohorts=", n_coh)
  if (n_age < 3) stop("Muy pocas edades para RW2(age)")
  if (n_per < 3) stop("Muy pocos períodos para RW2(period)")
  if (n_coh < 2) stop("Muy pocas cohortes para RW1(cohort)")
  invisible(df)
}

# ---------------------------
# 3C) Incidencia (cáncer pulmón) desde CSV
# ---------------------------
load_incidence_lung <- function(path_csv, pop_all_tbl) {
  inc_raw <- readr::read_csv(path_csv, show_col_types = FALSE)
  
  cause_col <- intersect(c("causa","causa2","diagnostico","dx"), names(inc_raw))[1]
  year_col  <- intersect(c("año","anio","ano","year","periodo"), names(inc_raw))[1]
  age_col   <- intersect(c("edad","age"), names(inc_raw))[1]
  sex_col   <- intersect(c("mujer","sexo","sex"), names(inc_raw))[1]
  case_col  <- intersect(c("casos_pred","casos_suavizados","casos","y","n","incid_suavizada","incid_suavizadas","conteo"),
                         names(inc_raw))[1]
  if (any(is.na(c(cause_col, year_col, age_col, sex_col, case_col)))) {
    stop(sprintf("Incidencia: faltan columnas clave. Disponibles: %s", paste(names(inc_raw), collapse=", ")))
  }
  
  inc_raw %>%
    mutate(causa_norm = str_to_lower(str_trim(.data[[cause_col]]))) %>%
    filter(str_detect(causa_norm, "pulm")) %>%
    transmute(
      period = as.integer(.data[[year_col]]),
      age    = as.integer(.data[[age_col]]),
      sex    = case_when(
        sex_col == "mujer" ~ ifelse(as.numeric(.data[[sex_col]]) == 1, "F", "M"),
        TRUE ~ ifelse(str_to_lower(as.character(.data[[sex_col]])) %in% c("f","female","mujer","1"), "F", "M")
      ),
      sex = factor(sex, levels = c("M","F")),
      cases = as.numeric(.data[[case_col]])
    ) %>%
    mutate(cohort = period - age) %>%
    left_join(pop_all_tbl, by = c("period","age","sex")) %>%
    filter(is.finite(exposure), exposure > 0)
}

# Loader general de incidencia por causa (usa regex normalizado)
load_incidence_by_cause <- function(path_csv, cause_regex, pop_all_tbl) {
  inc_raw <- readr::read_csv(path_csv, show_col_types = FALSE)
  
  cause_col <- intersect(c("causa2","causa","diagnostico","dx"), names(inc_raw))[1]
  year_col  <- intersect(c("año","anio","ano","year","periodo"), names(inc_raw))[1]
  age_col   <- intersect(c("edad","age"), names(inc_raw))[1]
  sex_col   <- intersect(c("mujer","sexo","sex"), names(inc_raw))[1]
  case_col  <- intersect(c("casos_pred","casos_suavizados","casos","y","n",
                           "incid_suavizada","incid_suavizadas","conteo"), names(inc_raw))[1]
  
  inc_raw %>%
    dplyr::mutate(causa_norm = norm_txt(.data[[cause_col]])) %>%
    dplyr::filter(stringr::str_detect(causa_norm, cause_regex)) %>%
    dplyr::transmute(
      period = as.integer(.data[[year_col]]),
      age    = as.integer(.data[[age_col]]),
      sex    = {
        if (sex_col == "mujer")
          ifelse(as.numeric(.data[[sex_col]]) == 1, "F", "M")
        else
          ifelse(stringr::str_to_lower(as.character(.data[[sex_col]])) %in% c("f","female","mujer","1"), "F", "M")
      },
      sex = factor(sex, levels = c("M","F")),
      cases = as.numeric(.data[[case_col]])
    ) %>%
    dplyr::mutate(cohort = period - age) %>%
    dplyr::left_join(pop_all_tbl, by = c("period","age","sex")) %>%
    dplyr::filter(is.finite(exposure), exposure > 0) %>%
    dplyr::arrange(sex, period, age)
}


# =============================================================


# ---------------------------
# 3D) Contrato de inputs estandarizado
# ---------------------------
make_bapc_inputs <- function(mort_hist_tbl,
                             pop_all_tbl,
                             inc_hist_tbl = NULL,
                             prev_path = PATH_PREV_DTA,
                             prev_data = NULL,
                             metadata = list()) {
  out <- list(
    mort_hist_tbl = mort_hist_tbl,
    pop_all_tbl   = pop_all_tbl,
    inc_hist_tbl  = inc_hist_tbl,
    prev_path     = prev_path,
    prev_data     = prev_data,
    metadata      = metadata
  )
  class(out) <- c("bapc_inputs", class(out))
  validate_bapc_inputs(out)
  out
}

validate_bapc_inputs <- function(inputs) {
  stopifnot(is.list(inputs))
  req <- c("mort_hist_tbl", "pop_all_tbl")
  miss <- req[!req %in% names(inputs)]
  if (length(miss)) stop("Faltan componentes en inputs: ", paste(miss, collapse = ", "))
  if (!is.data.frame(inputs$mort_hist_tbl) || nrow(inputs$mort_hist_tbl) == 0) stop("inputs$mort_hist_tbl debe ser un data.frame no vacío.")
  if (!is.data.frame(inputs$pop_all_tbl) || nrow(inputs$pop_all_tbl) == 0) stop("inputs$pop_all_tbl debe ser un data.frame no vacío.")
  invisible(inputs)
}

build_inputs_real_cause <- function(cfg_row,
                                    path_mort_csv = PATH_MORT_CSV,
                                    path_inc_csv  = PATH_INC_CSV,
                                    pop_all_tbl   = pop_all,
                                    prev_path     = PATH_PREV_DTA) {
  cfg_row <- tibble::as_tibble(cfg_row)
  stopifnot(nrow(cfg_row) == 1)
  mort_c <- load_mortality_by_cause(path_mort_csv, cfg_row$mort_regex[[1]], pop_all_tbl, cause_id = cfg_row$cause_id[[1]])
  inc_c  <- load_incidence_by_cause(path_inc_csv,  cfg_row$inc_regex[[1]],  pop_all_tbl)
  make_bapc_inputs(
    mort_hist_tbl = mort_c,
    pop_all_tbl   = pop_all_tbl,
    inc_hist_tbl  = inc_c,
    prev_path     = prev_path,
    metadata = list(
      cause_id = cfg_row$cause_id[[1]],
      label    = cfg_row$label[[1]],
      source   = "real"
    )
  )
}
