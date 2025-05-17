# 代码需要保证统一
# 只能是 sce.all.int

# ------------------------------ 1. 仅保留神经相关标记 -------------------------------
neural_markers_list <- list(
  Oligodendrocytes = c('OLIG1','OLIG2','SOX10'),
  Astrocytes = c('GFAP','AQP4','SLC1A3'),
  Microglia = c('P2RY12','TMEM119','CX3CR1'),
  Vascular = c('CLDN5','FLT1','PECAM1'),
  Progenitors = c('ASCL1','DCX','SOX2'),
  Excitatory_Neurons = c('SLC17A7','SLC17A6'),
  Inhibitory_Neurons = c('GAD1','GAD2')
)

# 添加到 markers_list 变量中
markers_list <- list(neural_markers_list = neural_markers_list)

# ------------------------------ 2. 基因名称适配函数 -------------------------------
# 自动处理人/小鼠基因名大小写
format_gene_names <- function(genes, species) {
  if (species == 'human') {
    return(toupper(genes))
  } else if (species == 'mouse') {
    return(stringr::str_to_title(genes))
  } else {
    stop("Unsupported species. Use 'human' or 'mouse'.")
  }
}

# ------------------------------ 3. 核心可视化部分 -------------------------------
p_umap <- DimPlot(sce.all.int, reduction = "umap", raster = F, label = T, repel = T) 

if (exists('sp')) {
  # 检查物种是否定义，若未定义则默认人类
  if (!exists('sp') || !(sp %in% c('human', 'mouse'))) {
    warning("Species (sp) not defined. Defaulting to 'human'.")
    sp <- 'human'
  }
  
  # 动态过滤不存在基因
  valid_genes <- lapply(neural_markers_list, function(genes) {
    formatted_genes <- format_gene_names(genes, sp)
    existing_genes <- formatted_genes[formatted_genes %in% rownames(sce.all.int)]
    if (length(existing_genes) == 0) {
      warning(paste("No valid genes found for group:", names(genes)))
      return(NULL)
    }
    return(existing_genes)
  })
  
  # 移除空组
  valid_genes <- valid_genes[!sapply(valid_genes, is.null)]
  
  if (length(valid_genes) > 0) {
    # 生成DotPlot
    p_dot <- DotPlot(sce.all.int, features = valid_genes, 
                    assay = 'RNA') + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # 保存结果
    w <- length(unique(unlist(valid_genes))) / 5 + 6
    ggsave('neural_markers_dotplot.pdf', plot = p_dot, width = w)
  } else {
    warning("No valid genes to plot.")
  }
  
  # 合并UMAP和DotPlot
  if (exists('p_dot')) {
    combined_plot <- p_umap + p_dot + plot_layout(ncol = 2)
    ggsave('combined_umap_dotplot.pdf', combined_plot, width = 16, height = 8)
  }
}

# ------------------------------ 4. 保留基础QC图 -------------------------------
if ("percent_mito" %in% colnames(sce.all.int@meta.data)) {
  qc_plots <- VlnPlot(sce.all.int, 
                      features = c("nFeature_RNA", "nCount_RNA", "percent_mito"),
                      pt.size = 0.1, 
                      ncol = 3)
  ggsave("qc_metrics.pdf", qc_plots, width = 12, height = 5)
}