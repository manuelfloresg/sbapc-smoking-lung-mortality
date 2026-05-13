
q_q <- readRDS('results/20260507_0153_STABLE_V14/raw_data/res_spec_linear_s4_quit.rds')
nms <- names(q_q$resM$inc_fit$rates_all_full)
cat('z_prev present? ', 'z_prev' %in% nms, '\n')
cat('offset_prev_rr present? ', 'offset_prev_rr' %in% nms, '\n')
cat('z_prev.x present? ', 'z_prev.x' %in% nms, '\n')
cat('z_prev.y present? ', 'z_prev.y' %in% nms, '\n')
