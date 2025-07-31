
# Filter data
filter_df <- function(cycle_data, rename_list, df_cols, sep) {
  purrr::map(cycle_data, ~ .x %>%
               dplyr::rename(!!!rename_list) %>%
               separate(col = V1, into = df_cols, sep = sep) %>%
               mutate(transcriptId = sub("\\.\\d+$", "", transcriptId)) %>%
               filter(!grepl("_", T8)))
}

filter_raw_df <- function(cycle_data, rename_list, df_cols, sep) {
  purrr::map(cycle_data, ~ .x %>%
               dplyr::rename(!!!rename_list) %>%
               separate(col = V1, into = df_cols, sep = sep) %>%
               mutate(transcriptId = sub("\\.\\d+$", "", transcriptId)) %>%
               dplyr::select(UMI, U8)
  )
}

load_data <- function(path, pattern=".tsv"){
  files_list <- list.files(path = path, recursive = T, pattern = pattern, full.names = TRUE)
  
  data.pcr <- list()
  for (path in files_list){
    file_name <- strsplit(path, split ="/")[[1]]
    tool <- file_name[[length(file_name)-1]]
    file_name <- as.character(file_name[length(file_name)])
    
    cycle <- file_name %>%
      strsplit(split="_") %>% unlist() %>%
      magrittr::extract(length(.)) %>%
      strsplit(split=".tsv") %>% unlist()
    
    if (!is.element(cycle, names(data.pcr))){
      data.pcr[[cycle]] <- NULL
    }
    df <- read.csv(path, sep = "\t", header = F)
    data.pcr[[cycle]][[tool]] <- df
  }
  return(data.pcr)
}

# Function to calculate accuracy 
accuracy_umi <- function(tool_data, col) {
  tool_data %>%
    dplyr::mutate(match = .data[[col[1]]] == .data[[col[2]]]) %>%
    dplyr::summarise(accuracy = sum(match, na.rm = TRUE) / n(), .groups = 'drop') %>%
    pull(accuracy)
}

# Process each tool with all matching conditions and return a data frame
process_tool <- function(tool_name, tool_data) {
  purrr::map2_dfr(names(mcs.li), mcs.li, ~data.frame(
    tool = tool_name,
    accuracy = .x,
    val = accuracy_umi(tool_data, .y) 
  ))
}

# Process each cycle, applying process_tool to each tool within the PCR cycle
process_cycle <- function(cycle_name, cycle_data) {
  tools_df <- bind_rows(
    purrr::map2_dfr(names(cycle_data), cycle_data, process_tool)
  )
  tools_df$cycle <- cycle_name
  tools_df
}

# Function to calculate Precision, recall for each class
prec_recall_fn <- function(levs, true_labels, predicted_labels){
  registerDoParallel(cores = max(4,detectCores() - 3))
  foreach(label = levels(levs), .combine = rbind) %dopar% { 
    tp <- sum(true_labels == label & predicted_labels == label, na.rm=TRUE)
    fp <- sum(true_labels != label & predicted_labels == label, na.rm=TRUE)
    fn <- sum(true_labels == label & predicted_labels != label, na.rm=TRUE)
    
    precision <- tp / (tp + fp)
    recall <- tp / (tp + fn)
    
    return(c(label, precision, recall))
  }
}

# Function to calculate precision and recall for a given tool's data
calculate_prec_recall <- function(true_labels, predicted_labels) {
  
  levs = c(true_labels, predicted_labels)
  
  results_df <- prec_recall_fn(levs, true_labels, predicted_labels)
  results_df <- as.data.frame(results_df)
  colnames(results_df) <- c("Label", "Precision", "Recall")
  
  macro_precision <-  mean(as.numeric(results_df$Precision), na.rm = TRUE)
  macro_recall <-  mean(as.numeric(results_df$Recall), na.rm = TRUE)
  
  return(data.frame(Macro_Precision = macro_precision, Macro_Recall = macro_recall))
}

# Function to process each tool within a cycle
prec_recc_by_tool <- function(tool_name, tool_data) {
  true_labels <- factor(tool_data$transcriptId)
  predicted_labels <- factor(tool_data$T8)
  levs = c(true_labels, predicted_labels)
  
  prec_recall_df <- calculate_prec_recall(levs, true_labels, predicted_labels)
  prec_recall_df$tool <- tool_name  
  return(prec_recall_df)
}

list.to.df <- function(prec.li){
  df_list <- list()
  for (cycle_name in names(prec.li)) {
    cycle_data <- prec.li[[cycle_name]]
    temp_df <- do.call(rbind, lapply(cycle_data, function(x) {
      data.frame(precision = x$Macro_Precision, recall = x$Macro_Recall)
    }))
    temp_df$tool <- rownames(temp_df)
    temp_df$cycle <- cycle_name
    df_list[[cycle_name]] <- temp_df}
  return(df_list)
}
