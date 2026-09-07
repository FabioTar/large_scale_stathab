
dev.off()
rm(list = setdiff(ls(), c("riverATLAS_eu", "wGAP_driver", "wGAP", "results_df",
                          "JDB_test")))

library(sf)
library(fst)
library(dplyr)
library(terra)
library(data.table)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)



wGAP_driver = st_read("C:/Users/ftarena/Desktop/DAM DATABASES/WaterGAP/dryver_net_eu.shp") # waterGAP dryver network shapefile
wGAP = read_fst("C:/Users/ftarena/Desktop/DAM DATABASES/WaterGAP/watergap_dis_net_eu.fst") # waterGAP discharges table

JDB = read.csv("C:/Users/ftarena/Desktop/DAM DATABASES/R_analysis/Databases/JDB_damtype12.csv", header = TRUE) # joined database
results_df = read.csv("C:/Users/ftarena/Desktop/DAM DATABASES/R_analysis/Other documents/AlgoHH_depth_new.csv", header = TRUE) # result of the bypass algorithm run (table associating each dam to the bypassed segments downstream of it)










# LIMIT 100 km      (Impose maximum bypassed length as 100km)

# Rename original columns
results_df$segments_or = results_df$segments
results_df$n_segments_or = results_df$n_segments
results_df$outlet_or = results_df$outlet

# Lookup table for segment lengths (in km)
seg_length = setNames(wGAP_driver$LENGTH_GEO / 1000,
                      wGAP_driver$DRYVER_RIV)

# Prepare new columns
results_df$`100km` = "no"

max_len = 100

for(i in seq_len(nrow(results_df))){
  
  cat("\r", i, "/", nrow(results_df),
      "- dam_id:", results_df$dam_id[i])
  flush.console()
  
  segs = unlist(strsplit(results_df$segments_or[i], ";"))
  segs_num = as.numeric(segs)
  
  # Get lengths in km
  lens = seg_length[as.character(segs_num)]
  
  # Remove NA lengths if any
  valid = !is.na(lens)
  
  segs_num = segs_num[valid]
  lens = lens[valid]
  
  total_len = sum(lens)
  
  # If total length exceeds threshold
  if(total_len > max_len){
    
    cum_len = cumsum(lens)
    
    keep_idx = which(cum_len <= max_len)
    
    # Add first segment exceeding threshold
    keep_idx = which(cum_len <= max_len)
    
    segs_new = segs_num[keep_idx]
    
    results_df$segments[i] = paste(segs_new, collapse = ";")
    results_df$n_segments[i] = length(segs_new)
    results_df$outlet[i] = tail(segs_new, 1)
    results_df$`100km`[i] = "yes"
  }
}









rm(list = setdiff(ls(), c("riverATLAS_eu", "wGAP_driver", "wGAP", "results_df",
                          "JDB")))











# Script with STOCK


# Functions 



# STO regulation function (compute release from storage dams with Hanasaki algorithm)

compute_Qreg = function(flow_ts, C, Imean) { # needs: flow time series, storage capacity, mean flow
  
  
  if (length(flow_ts) < 5) {
    Qreg = rep(mean(flow_ts, na.rm = TRUE), length(flow_ts))
    return(Qreg)
  }
  
  flow_ts = ifelse(is.na(flow_ts), 0, flow_ts)
  Imean = ifelse(is.na(Imean), mean(flow_ts, na.rm = TRUE), Imean)
  
  DOR = C / (Imean * 365 * 24 * 3600)
  
  r0 = Imean
  Qeco = 0.1 * Imean
  
  is_recharge = ifelse(is.na(flow_ts), FALSE, flow_ts > Imean)
  n = length(flow_ts)
  
  for (i in 3:(n-2)) {
    
    w = is_recharge[(i-2):(i+2)]
    
    if (all(!is.na(w[c(1,2,4,5)])) &&
        all(w[c(1,2,4,5)] == TRUE) &&
        !is.na(w[3]) && w[3] == FALSE) {
      
      is_recharge[i] = TRUE
    }
    
    if (all(!is.na(w[c(1,2,4,5)])) &&
        all(w[c(1,2,4,5)] == FALSE) &&
        !is.na(w[3]) && w[3] == TRUE) {
      
      is_recharge[i] = FALSE
    }
  }
  
  change = c(TRUE, diff(is_recharge) != 0)
  groups = cumsum(change)
  
  group_type = tapply(is_recharge, groups, function(x) x[1])
  group_lengths = tapply(flow_ts, groups, length)
  
  release_groups = which(group_type == FALSE)
  
  preceding_recharge = sapply(release_groups, function(g) {
    if (g > 1 && group_type[g-1]) group_lengths[g-1] else 0
  })
  
  if (length(release_groups) == 0 || all(preceding_recharge == 0)) {
    start_month = 1
  } else {
    best_group = release_groups[which.max(preceding_recharge)]
    idx = which(groups == best_group)
    start_month = if (length(idx) > 0) months[idx[1]] else 1
  }
  
  S = numeric(length(flow_ts) + 1)
  S[1] = 0.5 * C
  
  Qreg = numeric(length(flow_ts))
  krls = NA
  
  for (m in 1:length(flow_ts)) {
    
    if (m == 1 || months[m] == start_month) {
      krls = S[m] / (0.85 * C)
    }
    
    I_m = flow_ts[m]
    
    if (DOR >= 0.5) {
      r = krls * r0
    } else {
      alpha = (DOR / 0.5)^2
      r = alpha * krls * r0 + (1 - alpha) * I_m
    }
    
    Qreg[m] = max(r, Qeco)
    
    Vin = I_m * sec_month[m]
    Vout = Qreg[m] * sec_month[m]
    
    S[m+1] = S[m] + Vin - Vout
    
    if (S[m+1] < 0) S[m+1] = 0
    if (S[m+1] > C) S[m+1] = C
  }
  
  Qreg
}


# Get all downstream segments (function to get all downstream segments following the network until end)

downstream_lookup = split(wGAP_driver$DRYVER_RIV,
                          wGAP_driver$from_node)

seg_to_node = setNames(
  as.character(wGAP_driver$to_node),
  wGAP_driver$DRYVER_RIV
)


get_downstream_all = function(seg) {
  
  visited = character()
  stack = seg
  
  while (length(stack) > 0) {
    
    current_seg = stack[1]
    stack = stack[-1]
    
    if (!current_seg %in% visited) {
      
      visited = c(visited, current_seg)
      
      to_node = seg_to_node[[current_seg]]
      
      nxt = downstream_lookup[[to_node]]
      
      if (!is.null(nxt)) {
        stack = c(stack, nxt)
      }
    }
  }
  
  visited[-1]
}










# time setup
months = 1:12
days_in_month = c(31,28,31,30,31,30,31,31,30,31,30,31)
sec_month = days_in_month * 24 * 3600



# natural flow input --> computations for random year (chosen 2018 for now)
wGAP_2018 = wGAP %>% select(DRYVER_RIV, contains("_2018"))
wGAP_2018 = wGAP_2018 %>% mutate(DRYVER_RIV = as.character(DRYVER_RIV))







# Modification to results_df structure for next computation
seg_map = results_df %>%
  mutate(segments = str_split(segments, ";")) %>%
  unnest(segments) %>%
  group_by(dam_id) %>%
  mutate(seg_order = row_number(),
         nseg = n()) %>%
  ungroup() %>%
  mutate(segments = as.character(segments))

seg_map = seg_map %>% left_join(JDB, by = c("dam_id" = "ID"))






JDB = JDB %>% mutate(DRYVER_RIV = as.character(DRYVER_RIV)) # make dryver id as a character

JDB = JDB %>% left_join(
  wGAP_2018 %>%
    mutate(Imean = rowMeans(select(., -DRYVER_RIV), na.rm = TRUE)) %>% # Add Imean column (mean flow) to JDB (for compute_Qreg function)
    select(DRYVER_RIV, Imean),
  by = "DRYVER_RIV"
)



# initial state
Q_state = split(wGAP_2018[, -1], wGAP_2018$DRYVER_RIV) # copy of wGAP dataframe to be updated in the loop
Qstock_state = lapply(Q_state, function(x) rep(0, length(x))) # keep track of upstream stock



dam_list = JDB %>% arrange(desc(dam_order)) # dam ordering on barrage order high to low (for computation)






dam_storage = list() # container for plotting



# main loop
for (i in seq_len(nrow(dam_list))) { # following dam_list (so in barrage order)
  
  print(i)
  
  
  # Define things for current dam
  dam = dam_list[i, ]
  dam_id2 = dam$ID # "dam_id2" because "dam_id" already exists (results_df and seg_map)
  dam_type = dam$dam_type4
  
  segs = seg_map %>%
    filter(dam_id == dam_id2) %>% # keep only impacted segments for current dam
    arrange(seg_order)
  
  if (nrow(segs) == 0) {
    next
  }
  
  seg_list = segs$segments
  nseg = length(seg_list)
  
  seg1 = seg_list[1] # first segment DRYVER_RIV (id)
  
  Qin = as.numeric(Q_state[[seg1]]) # natural Q of the dam reach modified by the previous loops 
                                    # (i.e. by the already calculated upstream dams)
                                    # -Qstock (everywhere ds) and -Qturb (if bypassed)
  
  
  
  
  # Define Qbase = what's in the dam segment = Qin (ROR), Qreg (STO)
  if (dam_type == "STO") {
    
    Qreg = compute_Qreg(Qin, dam$res_volume3 * 1e6, dam$Imean) # Hanasaki --> Qreg = Q regulation = reservoir release
    Qbase = Qreg # STO dam segment --> Q = reservoir release = Qin - Qstock
    Qstock = Qin - Qreg
    
  } else {
    
    Qbase = Qin # ROR dam segment --> Q = natural flow 
    Qstock = rep(0, length(Qin)) # if ROR --> Qstock = 0s
    
  }
  

  Qnat = as.numeric(wGAP_2018[wGAP_2018$DRYVER_RIV == seg1, paste0("Q_", 1:12, "_2018")]) # original wGAP (for eflow)
  
  Qeflow = 0.1 * mean(Qnat, na.rm = TRUE) # defined on unmodified natural flow
  Qcut = 0 # cut in flow (debit armement)
  
  
  
  
  # Design Flow Qd
  if (dam_type == "ROR") {
    Qd = sort.int(Qnat, decreasing = TRUE)[2]
  } else {
    Qd = 1.5 * mean(Qnat, na.rm = TRUE)
  }
  
  
  
  # Definition Qturb (Turbined Flow) on Qbase
  if (nseg == 1) { # case dams with one segment (Qturb = 0, no bypass)
    
    Qturb = rep(0, length(Qbase))
    
  } else { # cases impacted segments > 2
    
    Qturb = rep(0, length(Qbase)) # vector initialisation
    
    for (m in 1:length(Qbase)) {
      
      if (Qbase[m] <= (Qeflow + Qcut)) {
        
        Qturb[m] = 0
        
      } else if (!is.na(Qbase[m])) {
        
        if (Qbase[m] >= (Qd + Qeflow)) {
          Qturb[m] = Qd
        } else {
          Qturb[m] = Qbase[m] - Qeflow
        }
        
      }
      
    }
    
  }
  
  
  
  
 
  
  
  
  
  # BYPASSED REACHES (-Qturb)
  
  for (j in seq_along(seg_list)) { # in order of bypassed segment for current dam (from us to ds)
    
    seg = seg_list[j] # current segment ID
    Qseg = as.numeric(Q_state[[seg]])
    
    Q_new = Qseg - Qstock # take out stocked Q
    Qstock_state[[seg]] = Qstock_state[[seg]] + Qstock # keep track of Qstock accumulation
    
    if (j < nseg) {           
      Q_new = Q_new - Qturb   # j = nseg (last segment = outlet) --> excluded because Qturb is given back 
    }
    
    Q_state[[seg]] = Q_new 
  }
  
  
  
  
  
  
  
  # All downstream propagation of current Qstock 
  # (from outlet to end of network)
  
    outlet_seg = seg_list[nseg] # outlet id
    downstream_segs = get_downstream_all(outlet_seg) # all segments ids downstream of outlet

    for (down_seg in downstream_segs) {

      Q_state[[down_seg]] = Q_state[[down_seg]] - Qstock # take out current Qstock from previous situation (previous loops)
      Qstock_state[[down_seg]] = Qstock_state[[down_seg]] + Qstock # keep track of Qstock accumulation

    }
    
    
    
    
    # For following plots (save dam by dam info)
    dam_storage[[as.character(dam_id2)]] = list(
      seg_list2 = seg_list,
      Qturb2 = Qturb,
      Qreg2 = if (dam_type == "STO") Qreg else rep(NA, length(Qbase)),
      Qbase2 = Qbase
    )
    
    
    
    
  
}





# New Modifed WATERGAP
Q_final = bind_cols(
  tibble(DRYVER_RIV = names(Q_state)),
  as.data.frame(do.call(rbind, Q_state))
)


















# Plot impacted reaches per dam (single dam)

dam_id2 = "5010"   # change

dam_info = dam_storage[[dam_id2]]
dam_type = dam_list$dam_type4[dam_list$ID == dam_id2]

seg_list = dam_info$seg_list2
Qturb = dam_info$Qturb2
Qreg = dam_info$Qreg2


dev.off()

df_plot_all = data.frame()

# loop over segments
for (seg in seg_list) {
  
  Qorig = wGAP_2018 %>%
    filter(DRYVER_RIV == seg) %>%
    select(starts_with("Q_")) %>%
    unlist() %>%
    as.numeric()
  
  Qfinal = as.numeric(Q_state[[seg]])
  Qbase = dam_info$Qbase2
  
  df_plot = data.frame(
    seg = seg,
    month = 1:12,
    Qorig = Qorig,
    Qfinal = Qfinal,
    Qturb = Qturb,
    Qstock = as.numeric(Qstock_state[[seg]]),
    Qbase = Qbase
  )
  
  # df_plot$Qin = df_plot$Qfinal + df_plot$Qstock
    
  df_plot_all = rbind(df_plot_all, df_plot)
  
  
  
  
  
  p = ggplot(df_plot, aes(x = month)) +
    geom_line(aes(y = Qorig, color = "WaterGAP"), linewidth = 1.1) +
    geom_line(aes(y = Qbase, color = "WaterGAP - Stock"), linewidth = 1.1) +
    geom_line(aes(y = Qfinal, color = "Final"), linewidth = 1.1) +
    geom_line(aes(y = Qturb, color = "Turbined"), linewidth = 1, linetype = "dashed") +
    # geom_line(aes(y = Qstock, color = "Ups Storage"), linewidth = 1,linetype = "dashed") +
    
    {
      if (dam_type == "STO") {
        geom_line(
          aes(y = Qreg, color = "Hanasaki Release"),
          linewidth = 1.1,
          linetype = "dotted"
        )
      } else {
        NULL
      }
    } +
    
    scale_color_manual(values = c(
      "WaterGAP" = "black",
      Final = "red",
      Turbined = "orange",
      # "Ups Storage" = "purple",
      "Hanasaki Release" = "green",
      "WaterGAP - Stock" = "blue"
    )) +
    
    scale_y_log10() +
    
    labs(
      title = paste0("Dam ", dam_id2, " [", dam_type, "] - Segment ", seg),
      x = "Time [month]",
      y = "Discharge [m3/s]",
      color = ""
    ) +
    
    theme_minimal()
  
  print(p)
}


# Filter on segment id
test = df_plot_all %>%
  filter(seg == 4908052)













# Plot first reach of each dam

dev.off()

df_plot_all = data.frame()

for (dam_id2 in names(dam_storage)[5000:5010]) {
  
  dam_info = dam_storage[[dam_id2]]
  dam_type = dam_list$dam_type4[dam_list$ID == dam_id2]
  dam_order = dam_list$dam_order[dam_list$ID == dam_id2]
  
  seg = dam_info$seg_list2[1]   # first reach only
  
  Qturb = dam_info$Qturb2
  Qreg = dam_info$Qreg2
  Qbase = dam_info$Qbase2
  
  Qorig = wGAP_2018 %>%
    filter(DRYVER_RIV == seg) %>%
    select(starts_with("Q_")) %>%
    unlist() %>%
    as.numeric()
  
  Qfinal = as.numeric(Q_state[[seg]])
  
  df_plot = data.frame(
    dam_id = dam_id2,
    seg = seg,
    month = 1:12,
    Qorig = Qorig,
    Qfinal = Qfinal,
    Qturb = Qturb,
    Qbase = Qbase,
    Qstock = as.numeric(Qstock_state[[seg]])
  )
  
  
  # add to global dataset
  df_plot_all = rbind(df_plot_all, df_plot)
  
  # plot immediately
  p = ggplot(df_plot, aes(x = month)) +
    geom_line(aes(y = Qorig, color = "WaterGAP"), linewidth = 1.1) +
    geom_line(aes(y = Qbase, color = "WaterGAP - Qstock"), linewidth = 1.1) +
    geom_line(aes(y = Qfinal, color = "Final"), linewidth = 1.1) +
    geom_line(aes(y = Qturb, color = "Turbined"),
              linewidth = 1,
              linetype = "dashed") +
    {
      if (dam_type == "STO") {
        geom_line(aes(y = Qreg, color = "Hanasaki Release"),
                  linewidth = 1.1,
                  linetype = "dotted")
      } else {
        NULL
      }
    } +
    scale_color_manual(values = c(
      "WaterGAP" = "black",
      Final = "red",
      Turbined = "orange",
      "Hanasaki Release" = "green",
      "WaterGAP - Qstock" = "blue"
    )) +
    scale_y_log10() +
    labs(
      title = paste0(
        "Dam ", dam_id2,
        " (order ", dam_order, ")",
        " [", dam_type, "] - First reach ", seg
      ),
      x = "Time [month]",
      y = "Discharge [m3/s]",
      color = ""
    ) +
    theme_minimal()
  
  print(p)
}




















