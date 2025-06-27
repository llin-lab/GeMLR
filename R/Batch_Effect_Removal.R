library(tidyverse)
repo <- "~/Desktop/SENNET/"
setwd(repo)

# Import Olink cleaned data
Olink_df <- read_csv("TNK_Olink_data.csv")
metadata <- read_csv("TNK_metadata.csv")
Olink_df$is_overlapping <- FALSE

## Import Olink control data and clean, align the column in need
Olink_control1 <- readxl::read_excel("./UM1-HV_NPX 052825.xlsx")
Olink_control1 %>% dplyr::filter(Assay == "Sample Control" | Assay == "TCS -1" | Assay == "TCS-2") -> Olink_control1
Olink_control1$subject_id <- paste("sample control", 1:nrow(Olink_control1))
colnames(Olink_control1) <- tolower(gsub(" ", "_", gsub("-", "_", colnames(Olink_control1))))
colnames(Olink_control1)[colnames(Olink_control1) == "mic_a/b"] <- "mic_a_b"
Olink_control1$Organ <- "Lung"
Olink_control1$is_overlapping <- TRUE

Olink_control2 <- readxl::read_excel("./U54 Murdock_NPX 052825.xlsx")
Olink_control2 %>% dplyr::filter(Assay == "Sample Control" | Assay == "Sample Control2") -> Olink_control2
Olink_control2$subject_id <- paste("sample control", 1:nrow(Olink_control2))
colnames(Olink_control2) <- tolower(gsub(" ", "_", gsub("-", "_", colnames(Olink_control2))))
colnames(Olink_control2)[colnames(Olink_control2) == "mic_a/b"] <- "mic_a_b"
Olink_control2$Organ <- "Colon"
Olink_control2$is_overlapping <- TRUE


shared_columns <- intersect(colnames(Olink_control1), colnames(Olink_df))
data_mat <- rbind(Olink_df[, shared_columns], Olink_control1[, shared_columns], Olink_control2[, shared_columns])
data_mat[,1:92] <- lapply(data_mat[,1:92], function(x) {
  if (is.character(x) || is.factor(x)) suppressWarnings(as.numeric(as.character(x))) else x
})

# Separate metadata and expression data
meta_cols <- c("subject_id", "Organ", "is_overlapping")
meta_info <- data_mat[, meta_cols]
expr_mat <- as.matrix(data_mat[, !colnames(data_mat) %in% meta_cols])

# Ensure batch_group is a factor
meta_info$Organ <- factor(meta_info$Organ)

library(limma)

# Subset overlapping samples
overlap_idx <- which(meta_info$is_overlapping == TRUE)
expr_overlap <- t(expr_mat[overlap_idx, ])  # Limma expects features x samples
meta_overlap <- meta_info[overlap_idx, ]

# Design matrix with batch group
design <- model.matrix(~ 0 + Organ, data = meta_overlap)
colnames(design) <- levels(meta_overlap$Organ)

# Blocking by subject_id
block <- factor(meta_overlap$subject_id)

# Estimate correlation across replicates
corfit <- duplicateCorrelation(expr_overlap, design, block = block)

# Fit model
fit <- lmFit(expr_overlap, design, block = block, correlation = corfit$consensus)

# Estimate contrast: Lung - Colon
contrast <- makeContrasts(Lung - Colon, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

# Estimated batch effects (per marker)
batch_effect <- fit2$coefficients[, 1]

# Apply correction to full data
expr_mat_corrected <- t(expr_mat)  # Now features x samples
for (i in 1:nrow(meta_info)) {
  if (meta_info$Organ[i] == "Lung") {
    expr_mat_corrected[, i] <- expr_mat_corrected[, i] - batch_effect
  }
}
expr_mat_corrected <- t(expr_mat_corrected)  # Back to samples x markers
expr_mat_corrected <- as.data.frame(expr_mat_corrected)

expr_mat_corrected$subject_id <- data_mat$subject_id
expr_mat_corrected$Organ <- data_mat$Organ
expr_mat_corrected$is_overlapping <- data_mat$is_overlapping

set.seed(259)
result_limma <- prcomp(expr_mat_corrected[,1:92])

Olink_pca_limma <- as.data.frame(result_limma$x)
Olink_pca_limma$Organ <- as.factor(expr_mat_corrected$Organ)

Olink_pca_limma %>% ggplot(aes(x = PC1, y = PC2, color = Organ)) +
  geom_point(size = 4) +
  theme_classic() +
  guides(fill = "none", color = guide_legend(title = "", override.aes = list(size = 3))) +
  theme(axis.title.x = element_text(size = 26, face = "bold"), axis.title.y = element_text(size = 26, face = "bold"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        legend.position = c(0.9, 0.9),           # Position inside plot (x, y in [0,1])
        legend.background = element_rect(fill = "white", color = "white"),
        legend.text = element_text(size = 22),   # Adjust font size of legend labels
        legend.title = element_text(size = 22)) +
  ggtitle("PCA of Olink data wrt batch")
ggsave("PCA of Olink full wrt batch after limma remove batch effect by batch with overlapped samples.png", width = 10, height = 8, device = "png", dpi = 1000)

write_csv(expr_mat_corrected, "expr_mat_corrected_by_limma_with_overlapped_samples.csv")


kmeans_result <- kmeans(result_limma$x, centers = 2)
Olink_pca_limma$cluster <- as.factor(kmeans_result$cluster)
Olink_pca_limma %>% ggplot(aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 4) +
  theme_classic() +
  guides(fill = "none", color = guide_legend(title = "", override.aes = list(size = 3))) +
  theme(axis.title.x = element_text(size = 26, face = "bold"), axis.title.y = element_text(size = 26, face = "bold"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        legend.position = c(0.9, 0.9),           # Position inside plot (x, y in [0,1])
        legend.background = element_rect(fill = "white", color = "white"),
        legend.text = element_text(size = 22),   # Adjust font size of legend labels
        legend.title = element_text(size = 22)) +
  ggtitle("PCA of Olink data wrt batch")
ggsave("PCA of Olink full wrt batch after limma remove batch effect by cluster with overlapped samples.png", width = 10, height = 8, device = "png", dpi = 1000)









