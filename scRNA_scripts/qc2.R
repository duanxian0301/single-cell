basic_qc <- function(input_sce) {
  # 计算线粒体基因比例
  mito_genes <- rownames(input_sce)[grep("^MT-", rownames(input_sce), ignore.case = TRUE)]
  print(mito_genes)
  input_sce <- PercentageFeatureSet(input_sce, features = mito_genes, col.name = "percent_mito")
  
  # 计算核糖体基因比例
  ribo_genes <- rownames(input_sce)[grep("^Rp[sl]", rownames(input_sce), ignore.case = TRUE)]
  print(ribo_genes)
  input_sce <- PercentageFeatureSet(input_sce, features = ribo_genes, col.name = "percent_ribo")
  
  # 计算红血细胞基因比例
  Hb_genes <- rownames(input_sce)[grep("^Hb[^(p)]", rownames(input_sce), ignore.case = TRUE)]
  print(Hb_genes)
  input_sce <- PercentageFeatureSet(input_sce, features = Hb_genes, col.name = "percent_hb")
  
  # 可视化 QC 指标
  feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
  VlnPlot(input_sce, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3)
  
  # 合并过滤条件
  filtered_cells <- WhichCells(
    input_sce,
    expression = nFeature_RNA > 200 &         # 基因数大于 200
                 nFeature_RNA < 5000 &        # 基因数小于 5000
                 percent_mito < 30 &          # 线粒体基因比例小于 30%
                 percent_hb < 5               # 红血细胞基因比例小于 5%
  )
  
  # 应用过滤
  input_sce.filt <- subset(input_sce, cells = filtered_cells)
  
  # 检查过滤后的样本分布
  table(input_sce.filt$orig.ident)
  
  return(input_sce.filt)
}