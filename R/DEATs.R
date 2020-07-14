# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

## Function to map distribution of differentially expressed genes
DE_distr <- function(res) {
  significant <- sum(res$padj < 0.01,na.rm = TRUE)
  hist_col <- c("red", rep("grey", 99))
  hist_labs <- c(significant, rep("", 99))
  hist(res$padj,breaks = 100, col = hist_col, ylim = c(0,2000), labels = hist_labs,
       main = "Adjusted P-value Distribution", xlab="Adjusted P-value")
}

##Function to filter expression matrix for only significant genes
sig_df_filter <- function(rld, rld_df, res){
  rownames(rld_df) <- unlist(lapply(rownames(rld_df), endstrip))
  idx <- match(rownames(rld_df), genemap$ensembl_gene_id)
  rld_df$entrez <- genemap$entrezgene_id[idx]
  rld_df$mgi_symbol <- genemap$mgi_symbol[idx]
  print("did this connect to github")

  res$entrez <- genemap$entrezgene_id[idx]
  res$mgi_symbol <- genemap$mgi_symbol[idx]

  sig_genes <- which(res$padj < 0.01) ##p <0.05 not consistent between HFD individuals

  rownames(rld) <- rld_df$mgi_symbol
  sig_df <- assay(rld)[sig_genes,]
  sig_df <- if (length(grep("Gm", rownames(sig_df))) != 0){sig_df[-grep("Gm", rownames(sig_df)),]}else{sig_df}
  return(sig_df)
}

## function to trim chip labels to create proper ensembl values
endstrip <- function(x) {strsplit(x,split = ".", 1)[[1]][1]}
