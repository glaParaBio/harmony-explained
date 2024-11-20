library(data.table)
library(tidyverse)
library(ggthemes)
library(ggrepel)
library(harmony)
library(patchwork)
library(tidyr)

## Useful util functions

cosine_normalize <- function(X, margin) {
    if (margin == 1) {
        res <- sweep(as.matrix(X), 1, sqrt(rowSums(X ^ 2)), '/')
        row.names(res) <- row.names(X)
        colnames(res) <- colnames(X)        
    } else {
        res <- sweep(as.matrix(X), 2, sqrt(colSums(X ^ 2)), '/')
        row.names(res) <- row.names(X)
        colnames(res) <- colnames(X)
    }
    return(res)
}

onehot <- function(vals) {
    t(model.matrix(~0 + as.factor(vals)))
}


colors_use <- c(`jurkat` = rgb(129, 15, 124, maxColorValue=255),
                `t293` = rgb(208, 158, 45, maxColorValue=255),
                `half` = rgb(0, 109, 44, maxColorValue=255))


do_scatter <- function(umap_use, meta_data, label_name, no_guides = TRUE, do_labels = TRUE, nice_names, 
                       palette_use = colors_use,
                       pt_size = 4, point_size = .5, base_size = 10, do_points = TRUE, do_density = FALSE, h = 4, w = 8) {
    umap_use <- umap_use[, 1:2]
    colnames(umap_use) <- c('X1', 'X2')
    plt_df <- umap_use %>% data.frame() %>% 
        cbind(meta_data) %>% 
        dplyr::sample_frac(1L) 
    plt_df$given_name <- plt_df[[label_name]]
    
    if (!missing(nice_names)) {
        plt_df %<>%
            dplyr::inner_join(nice_names, by = "given_name") %>% 
            subset(nice_name != "" & !is.na(nice_name))

        plt_df[[label_name]] <- plt_df$nice_name        
    }
        
    plt <- plt_df %>% 
        ggplot(aes(X1, X2, colour = .data[[label_name]], fill = .data[[label_name]])) + 
        theme_tufte(base_size = base_size) + 
        theme(panel.background = element_rect(fill = NA, color = "black")) + 
        guides(color = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16, size = 4)), alpha = FALSE) +
        scale_color_manual(values = palette_use) + 
        scale_fill_manual(values = palette_use) +    
        theme(plot.title = element_text(hjust = .5)) + 
        labs(x = "UMAP 1", y = "UMAP 2") 
    
    if (do_points) 
        plt <- plt + geom_point(size = 0.2)
    if (do_density) 
        plt <- plt + geom_density_2d()    
        

    if (no_guides)
        plt <- plt + guides("none")
    
    if (do_labels) 
        plt <- plt + geom_label_repel(data = data.table(plt_df)[, .(X1 = mean(X1), X2 = mean(X2)), by = label_name], label.size = NA,
                                      aes(label = .data[[label_name]]), color = "white", size = pt_size, alpha = 1, segment.size = 0) + 
        guides(col = FALSE, fill = FALSE)
    return(plt)
}
