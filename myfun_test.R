#LOADING RAW DATA

loading_data <- function(file_path){
  raw_data <- data.table::fread(file_path) %>% #Carica file
    tibble::as_tibble(.name_repair=janitor::make_clean_names) #Modifica file: il nome delle colonne senza spazi e con minuscolo
  return(raw_data)
}

#MAKE EXPERIMENTAL DESIGN
exp_design <- function(data, pattern_interest){
  
  exp_des <- data %>% 
    dplyr::select(gene_names, dplyr::starts_with(pattern_interest)) %>%
    tidyr::pivot_longer(!gene_names, names_to="key", values_to= "intensity") %>% #Raggruppa tutte le lfq in unica colonna
    dplyr::distinct(key) %>%
    dplyr::mutate(label = stringr::str_remove(key, "lfq_intensity_bc_")) %>%
    dplyr::mutate(condition = stringr::str_remove(label, "_[^_]*$")) %>%
    dplyr::mutate(replicate = stringr::str_remove(label, ".*_"))
  
  return(exp_des)
}

#DATA PREPROCESSING
pre_process <- function(data, pattern_interest){
  
  data_pre <- data %>% 
    dplyr::mutate(dplyr::across(dplyr::starts_with(pattern_interest), ~ log2(.))) %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with(pattern_interest), ~ dplyr::na_if(.,-Inf)))
  
  return(data_pre)
}

#DATA WRANGLING
data_wrangling <- function(data, 
                           pep_filter, 
                           pep_thr, 
                           rev=TRUE,
                           cont=TRUE,
                           oibs=TRUE){
  
  data_wrang <- data %>% #DATA STANDARDIZED
    dplyr::select(protein_i_ds, gene_names, id) %>%
    dplyr::mutate(gene_names = stringr::str_extract(gene_names, "[^;]*")) %>%
    ## every protein groups now have only 1 gene name associated to it
    dplyr::rename(unique_gene_names = gene_names) %>%
    janitor::get_dupes(unique_gene_names) %>%
    dplyr::mutate(unique_gene_names = dplyr::case_when(
      unique_gene_names != "" ~ paste0(
        unique_gene_names, "__",
        stringr::str_extract(protein_i_ds, "[^;]*")),
      TRUE ~ stringr::str_extract(protein_i_ds, "[^;]*"))) %>%
    dplyr::select(unique_gene_names, id) %>%
    dplyr::right_join(data, by = "id") %>%
    dplyr::mutate(gene_names = dplyr::case_when(unique_gene_names != "" ~ unique_gene_names,
                                                TRUE ~ gene_names)) %>%
    dplyr::select(-unique_gene_names) %>%
    dplyr::mutate(gene_names = dplyr::if_else(gene_names == "",
                                              stringr::str_extract(protein_i_ds, "[^;]*"),
                                              gene_names)) %>%
    dplyr::mutate(gene_names = stringr::str_extract(gene_names, "[^;]*")) %>% 
    dplyr::select(gene_names,
                  dplyr::all_of(exp_des$key),
                  peptides,
                  razor_unique_peptides,
                  unique_peptides,
                  reverse,
                  potential_contaminant,
                  only_identified_by_site) %>% 
    tidyr::pivot_longer(!c(gene_names,
                           peptides,
                           razor_unique_peptides,
                           unique_peptides,
                           reverse,
                           potential_contaminant,
                           only_identified_by_site),
                        names_to = "key",
                        values_to = "raw_intensity") %>% 
    dplyr::inner_join(., exp_des, by = "key") %>%   #aggiunge righe e colonne che matchano tra expdesign e data
    dplyr::mutate(bin_intensity = dplyr::if_else(is.na(raw_intensity), 0, 1)) %>%  #Nuova colonna. 1 se valore esite, 0 se NA
    dplyr::select(-key) %>% 
    {if(rev)dplyr::filter(., !reverse == "+") else .} %>% #DATA WRANGLING
    {if(cont)dplyr::filter(., !potential_contaminant == "+") else .} %>%
    {if(oibs)dplyr::filter(., !only_identified_by_site == "+") else .} %>% 
    ## filter on peptides:
    {if(pep_filter == "peptides"){dplyr::filter(., peptides >= pep_thr)}
      else if (pep_filter == "unique") {dplyr::filter(., unique_peptides >= pep_thr)}
      else {dplyr::filter(., razor_unique_peptides >= pep_thr)}}
  
  return(data_wrang)
  
}

#DATA FILTERING
data_filtered <- function(data,
                          valid_val_filter,
                          valid_val_thr) {
  data_filt<- data %>% 
    {if(valid_val_filter == "total")dplyr::group_by(., gene_names)
      else dplyr::group_by(., gene_names, condition)} %>% 
    dplyr::mutate(miss_val = dplyr::n() - sum(bin_intensity)) %>% 
    dplyr::mutate(n_size = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(gene_names) %>%
    ## rage compreso tra 0 e 100% espresso in valori tra 0 e 1
    {if(valid_val_filter == "alog") dplyr::filter(., any(miss_val <= round(n_size * (1 - valid_val_thr), 0)))
      else dplyr::filter(., all(miss_val <= round(n_size * (1 - valid_val_thr), 0)))} %>%
    dplyr::ungroup() %>%
    dplyr::select(gene_names, label, condition, replicate, bin_intensity, raw_intensity) %>% 
    dplyr::rename(intensity = raw_intensity)
  return(data_filt)
}

#DATA IMPUTATION
data_imputed <- function(data,
                         shift, 
                         scale, 
                         unique_visual = FALSE){
  
  imputed_data <- data %>%
    dplyr::group_by(gene_names, condition) %>%
    dplyr::mutate(for_mean_imp = dplyr::if_else((sum(bin_intensity) / dplyr::n()) >= 0.75, TRUE, FALSE)) %>%
    dplyr::mutate(mean_grp = mean(intensity, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(imp_intensity = dplyr::case_when(
      bin_intensity == 0 & for_mean_imp ~ mean_grp,
      TRUE ~ as.numeric(intensity))) %>%
    dplyr::mutate(intensity = imp_intensity) %>% 
    dplyr::select(-c(for_mean_imp, mean_grp, imp_intensity))%>%
    dplyr::group_by(label) %>%
    # Define statistic to generate the random distribution relative to sample
    dplyr::mutate(mean = mean(intensity, na.rm = TRUE),
                  sd = sd(intensity, na.rm = TRUE),
                  n = sum(!is.na(intensity)),
                  total = nrow(data) - n) %>%
    dplyr::ungroup() %>%
    # Impute missing values by random draws from a distribution
    # which is left-shifted by parameter 'shift' * sd and scaled by parameter 'scale' * sd.
    dplyr::mutate(imp_intensity = dplyr::case_when(is.na(intensity) ~ rnorm(total,
                                                                            mean = (mean - shift * sd), 
                                                                            sd = sd * scale),
                                                   TRUE ~ intensity)) %>%
    dplyr::mutate(intensity = imp_intensity) %>%
    dplyr::select(-c(mean, sd, n, total, imp_intensity)) %>%
    dplyr::group_by(condition)
  
  return(imputed_data) 
}

#COLOR PALETTE (8 options, from A to h)
define_colors = function(){
  n_of_color <- max(exp_des %>% dplyr::count(replicate) %>% dplyr::pull(n))
  color_palette <- viridis::viridis(n = n_of_color , direction = -1, end = 0.70, begin = 0.30)
}

#BARPLOT COUNT
barplot_count <- function (data,
                           color="black",
                           width_barplot=0.5){
  bar <- data %>% 
    dplyr::group_by(label) %>%
    dplyr::summarise(counts = sum(bin_intensity)) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(., exp_des, by = "label") %>%
    dplyr::mutate(replicate = as.factor(replicate)) %>%
    dplyr::group_by(condition) %>% 
    
    ggplot2::ggplot(aes(x=label, y=counts, fill=condition))+
    geom_bar(stat="identity", width = width_barplot, color=color)+
    scale_fill_manual(values=color_palette)+
    theme_cuc()+
    theme(axis.text.x = element_text(angle = 30))+
    labs(title= "Protein counts", x="Label", y="Counts")+
    geom_text(aes(label=counts), vjust=-0.2, size=4)
  
  return(bar)
}

#BARPLOT COVERAGE
barplot_cover <- function (data,
                           fillcolor="#21908CFF",
                           color="black",
                           width_barplot=0.5){
  bar <- data %>% 
    dplyr::group_by(gene_names) %>%
    dplyr::summarise(counts = sum(bin_intensity)) %>%
    dplyr::ungroup() %>%
    dplyr::count(counts) %>% 
    dplyr::rename(occurrence = n) %>% 
    
    ggplot2::ggplot(aes(x=counts, y=occurrence))+
    geom_bar(stat="identity", width = width_barplot, color=color, fill=fillcolor)+
    theme_cuc()+
    scale_fill_manual(values= c(fillcolor))+
    labs(title= "Protein coverage", x="Counts", y= "Occurence")+
    geom_text(aes(label=counts), vjust=-0.2, size=4)+
    scale_x_continuous(breaks = 4:15 )
  
  return(bar)
}

#BOXPLOT
boxplot <- function (data,
                     width_boxplot=0.5,
                     color="black"){
  box <- ggplot2::ggplot(data, aes(x=label, y=intensity, fill=condition))+
    geom_boxplot(width=width_boxplot, color=color)+
    scale_fill_manual(values=color_palette)+
    theme_cuc()+
    theme(axis.text.x = element_text(angle = 30))+
    labs(title="Normalized data distribution", x="Label", y="Intensity")
  
  return(box)
}

#DENSITYPLOT
densityplot <- function (data,
                         color2="#481567FF",
                         color1="#21908CFF"){
  den <- data %>% 
    dplyr::group_by(gene_names) %>%
    dplyr::summarise(mean = mean(intensity, na.rm = TRUE),
                     missval = any(is.na(intensity))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(missing_value = dplyr::if_else(missval, "Missing", "Valid")) %>%
    dplyr::mutate(missing_value = factor(missing_value, levels = c("Valid", "Missing"))) %>%
    dplyr::group_by(missing_value) %>%
    
    ggplot2::ggplot(aes(x=mean, color=missing_value))+
    geom_density(linewidth=0.8, linetype='solid')+
    scale_color_manual(values = c(color1, color2))+
    theme_cuc()+
    labs(title="Density plot", x="Log2Intensity", y="Density")
  
  return(den)
}

#BARPLOT MISSING VALUES
barplot_missval <- function (data,
                             color="black",
                             color1="#481567FF",
                             color2="#21908CFF",
                             width_barplot=0.8){
  bar <- data %>% 
    dplyr::group_by(label) %>% 
    dplyr::mutate(bin_intensity = dplyr::if_else(bin_intensity == 1, "Valid", "Missing")) %>%
    dplyr::count(bin_intensity) %>% 
    dplyr::mutate(bin_intensity = as.factor(bin_intensity)) %>% 
    dplyr::rename(Data= bin_intensity) %>% 
    
    
    ggplot2::ggplot(aes(x=label, y=n, fill=Data))+
    geom_bar(stat="identity", width = width_barplot, color=color)+
    scale_fill_manual(values = c(color1, color2 ))+
    theme_cuc()+
    theme(axis.text.x = element_text(angle = 30))+
    labs(title= "Plot valid and missing data", x="Samples", y= "Counts")+
    geom_text(aes(label=n), position= position_stack(vjust = 0.5), size=3, show.legend = FALSE, color="white")
  
  return(bar)
}


#DENSITYPLOT IMPUTATION
den_imput <- function (data,
                       color1="#43BF71FF",
                       color2="#21908CFF",
                       color3="#35608DFF"){
  den <- data %>%
    dplyr::group_by(gene_names, condition) %>%
    dplyr::mutate(for_mean_imp = dplyr::if_else((sum(bin_intensity) / dplyr::n()) >= 0.75, TRUE, FALSE)) %>%
    dplyr::mutate(mean_grp = mean(intensity, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(imp_intensity = dplyr::case_when(
      bin_intensity == 0 & for_mean_imp ~ mean_grp,
      TRUE ~ as.numeric(intensity))) %>%
    dplyr::mutate(intensity = imp_intensity) %>% 
    dplyr::select(-c(for_mean_imp, mean_grp, imp_intensity))%>%
    dplyr::group_by(label) %>%
    # Define statistic to generate the random distribution relative to sample
    dplyr::mutate(mean = mean(intensity, na.rm = TRUE),
                  sd = sd(intensity, na.rm = TRUE),
                  n = sum(!is.na(intensity)),
                  total = nrow(data) - n) %>%
    dplyr::ungroup() %>%
    # Impute missing values by random draws from a distribution
    # which is left-shifted by parameter 'shift' * sd and scaled by parameter 'scale' * sd.
    dplyr::mutate(imp_intensity = dplyr::case_when(is.na(intensity) ~ rnorm(total,
                                                                            mean = (mean - 1.8 * sd), 
                                                                            sd = sd * 0.3),
                                                   TRUE ~ intensity)) %>%
    dplyr::mutate(intensity = imp_intensity) %>%
    dplyr::select(-c(mean, sd, n, total, imp_intensity)) %>%
    dplyr::group_by(condition) %>% 
    
    ggplot2::ggplot(aes(x=intensity, color=condition))+
    geom_density(linewidth=0.8, linetype='solid')+
    scale_color_manual(values = c(color1, color2, color3 ))+
    theme_cuc()+
    labs(title="Imputation plot", x="Log2Intensity", y="Density")
  
  return(den)
}

#HEATMAP
htmap <- function(data){
  map <- data %>%
    dplyr::select(gene_names, label, intensity) %>%
    tidyr::pivot_wider(names_from = label, values_from = intensity) %>%
    dplyr::filter(dplyr::if_all(.cols = dplyr::everything(), .fns = ~ !is.na(.x))) %>%
    tibble::column_to_rownames("gene_names") %>%
    cor() %>%
    round(digits = 2) %>% 
    Heatmap(name="Correlation",
            col = color_palette,
            column_title = "Correlation plot",
            column_title_gp=gpar(fontsize=18),
            cluster_rows = FALSE,
            cluster_columns = FALSE)
  
  return(map)
}

#PCA
#Make matrix
mat <- function(data) {
  mat <- data %>% 
    dplyr::select(gene_names, label, intensity) %>%
    tidyr::pivot_wider(id_cols = "gene_names",
                       names_from = "label",
                       values_from = "intensity") %>%
    tibble::column_to_rownames("gene_names") %>%
    as.matrix()
  return(mat)
}

#PCAPLOT
pca_plot <- function(data,
                     color1="#43BF71FF",
                     color2="#21908CFF",
                     color3="#35608DFF"){
  
  pca <- ggplot2::ggplot(data, aes(x=x, y=y),group=condition)+
  geom_point(size=3, shape=19, aes(color=condition))+
  theme_cuc()+
  geom_hline(yintercept = 0, linetype="longdash")+
  geom_vline(xintercept = 0, linetype="longdash")+
  scale_color_manual(values = c(color1, color2, color3))+
  labs(title="PCA", subtitle = "Principal component analysis", x="PC1", y="PC2")+
  geom_text(aes(label=replicate), size=3, position = "dodge", hjust=1.5)
  
  return(pca)
}


#DATASET_SIGNIFICANT
significant <- function(data, test){
  data <- data %>% 
    dplyr::select(gene_names, starts_with(test)) %>% 
    rename_at(vars(matches(test)), ~ str_remove(., paste0(test, "_"))) %>% 
    dplyr::filter(., !significant==FALSE) %>% 
    dplyr::mutate(regulation = dplyr::if_else(fold_change > 0, "Up", "Down")) %>% 
    dplyr::select(gene_names, fold_change, p_val, p_adj, regulation) %>%
    dplyr::mutate(gene_names = stringr::str_extract(gene_names, "[^;_]*")) %>% 
    dplyr::distinct(gene_names, .keep_all = TRUE) 
  return(data)
}


#FILTER_UP_DOWN_REGULATED
sig_up_down <- function(data, remove){
  data <- data %>% 
    dplyr::filter(., !regulation == remove)
  return(data)
}


#NODES
nodes <- function(data, nodesize=5){
  data <- data %>% 
    dplyr::pull(gene_names) %>% 
    rba_string_map_ids(species = 9606, echo_query = TRUE) %>%
    dplyr::select(stringId, preferredName, annotation ) %>% 
    dplyr::rename(name=preferredName) %>% 
    dplyr::left_join(sig, by= c("name"="gene_names")) %>% 
    dplyr::select(-c(stringId, annotation, fold_change)) %>% 
    dplyr::rename(value =p_adj, size = p_val, grp = regulation) %>% 
    dplyr::mutate(value = -log10(value)) %>% 
    dplyr::mutate(size = -log10(size)*nodesize) %>% 
    dplyr::mutate(grp = as.factor(grp))
  
  
  return(data)
}

#EDGES
edges <- function(data, score=0.4, edgesize=5){
  data <- data %>% 
    dplyr::pull(name) %>%
    rba_string_interactions_network(species = 9606) %>% 
    dplyr::mutate(Names= paste0(preferredName_A, "_", preferredName_B)) %>% 
    dplyr::distinct(Names, .keep_all = TRUE) %>% 
    dplyr::rename(FROM = preferredName_A, TO=preferredName_B) %>%
    dplyr::mutate(totscore = escore+dscore) %>% 
    dplyr::filter(!totscore== "0") %>%
    dplyr::select(FROM,TO, escore, dscore) %>%
    dplyr::mutate(Score1 = (escore-0.041)*(1-0.041)) %>% 
    dplyr::mutate(Score2 = (dscore-0.041)*(1-0.041)) %>% 
    dplyr::mutate(Score_combin = 1-(1-Score1)*(1-Score2)) %>% 
    dplyr::mutate(tot_score= Score_combin+0.041*(1-Score_combin)) %>% 
    dplyr::filter(tot_score>score) %>% #Filtro sullo score
    dplyr::mutate(tot_score = tot_score*edgesize) %>% 
    dplyr::select(FROM,TO, tot_score) %>% 
    dplyr::rename(Source=FROM, Target=TO, Value=tot_score)
  
  return(data)
}

#FILTER_NODES_DEGREE_ZERO
filter_nodes <- function(datanodes, data_edges){
  data <- datanodes %>% 
    left_join(data_edges, by=c("name"="Source")) %>% 
    left_join(data_edges, by=c("name"="Target")) %>% 
    dplyr::mutate(Source = dplyr::if_else(is.na(Source), 0, 1)) %>% 
    dplyr::mutate(Target = dplyr::if_else(is.na(Target), 0, 1)) %>% 
    dplyr::mutate(Exists = Source+Target) %>% 
    dplyr::filter(!Exists== "0") %>% 
    dplyr::select(name, size, value, grp) %>% 
    dplyr::distinct(name, .keep_all = TRUE)
  return(data)
}

#PLOT_NETWORK
plot_net <- function(datanode, dataedge, animation=FALSE, layout="force"){
  p <- echarts4r::e_charts(animation=FALSE) %>% 
    echarts4r::e_graph(roam=TRUE, 
                       force= list( initLayout = layout, 
                                    repulsion=100,
                                    edgeLength=30,
                                    layoutAnimation=animation),
                       itemStyle=list(opacity=0.65),
                       lineStyle=list(curveness=0.2), 
                       emphasis=list(focus="adjacency", 
                                     lineStyle=list(width=10))) %>% 
    echarts4r::e_graph_nodes(nodes=datanode, names = name,  value=value, size=size, category=grp) %>% 
    echarts4r::e_graph_edges(edges=dataedge, source=Source, target=Target, value=Value, size = Value) %>%
    echarts4r::e_color(c("blue", "red")) %>% 
    echarts4r::e_labels(datanode$name, font_size=4)%>%
    echarts4r:: e_title("Network", "Up & Down Regulated") %>%
    echarts4r::e_tooltip()
  
  return(p)
}

#PLOT_ALL
plot_all <- function(data, test, remove, score=0.4, animation=FALSE, layout="force", no_edge=TRUE){
  sig <- significant(data, test)
  up_down <- sig_up_down(data=sig, remove)
  nodes <- nodes(data=up_down)
  edges <- edges(data= nodes, score)
  if (no_edge) {nodes<- filter_nodes (datanodes =nodes, data_edges = edges)}
  p <- plot_net(datanode=nodes, dataedge=edges, animation, layout)
  
  return(p)
}


### DATABASE SEARCH

#SOURCE_DB
source_db <- function (data){
  data <- data %>% 
    dplyr::select(Source, Target, Score) %>% 
    dplyr::pull(Source) %>%
    bitr(fromType="UNIPROT", toType="SYMBOL", OrgDb="org.Hs.eg.db") %>% 
    dplyr::rename("Source"="SYMBOL") %>% 
    dplyr::left_join(unique_df, by= c("UNIPROT"="Source")) %>% 
    dplyr::select(Source, Target, Score)
  
  return(data) 
}