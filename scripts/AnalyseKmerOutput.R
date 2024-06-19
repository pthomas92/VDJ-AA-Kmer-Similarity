library(tidyverse); library(dendextend); library(RColorBrewer)

mode = 'Sanger'
similarity = c('hamming', 
               'kidera')
plot_heats = F

p = paste('../mAb feedback/datasets/VH-VK/Correct reclustering/20230519_CDR3_IMGT_3mer_complete/',
          mode, sep = '')

heatmaps <- list()

for(sim in similarity){
  
  mats = lapply(paste(p, list.files(p, pattern = paste('[^mAb]_', sim, sep = '')), sep = '/'),
                read.csv,
                sep = '\t')
  
  names(mats) <- sapply(strsplit(list.files(p, pattern = paste('[^mAb]_', sim, sep = '')), '_'), function(x) return(x[1]))
  
  hclust_ = lapply(mats, function(x){
    
    return(hclust(as.dist(x[,2:ncol(x)]), method = 'complete'))
    
  })
  
  cophens_ = lapply(hclust_, cophenetic)
  
  names(cophens_) = names(mats)
  cophens_ = do.call('rbind', lapply(cophens_, function(x){
    
    m = as.matrix(x)
    
    xy = t(combn(colnames(m), 2))
    
    data.frame(xy, COPHENETIC_DISTANCE= m[xy])
    
  })) 
  
  if(mode == 'Sanger'){
    
    cophens_ = cophens_ %>% 
      select(COPHENETIC_DISTANCE) %>% rownames_to_column('ID') %>%
      mutate(SAMPLE_ID = gsub('\\..*?$', '', ID))
    
    agg_ <- list()
    
    for(i in 1:100){
      
      agg_[[i]] = aggregate(COPHENETIC_DISTANCE ~  SAMPLE_ID,
                            cophens_ %>% 
                              slice_sample(n = 20),
                            mean)
      
    }
    
    agg_ = do.call('rbind', agg_)
    
    print(sim)
    
    test_ = t.test(agg_ %>% filter(SAMPLE_ID == 'gp120') %>% pull(COPHENETIC_DISTANCE),
                   agg_ %>% filter(SAMPLE_ID == 'gp120 + Rapamycin') %>% pull(COPHENETIC_DISTANCE))
    
    if(test_$p.value < 0.05){
      
      plot_bar = T
      
    } else {
      
      plot_bar = F
      
    }
    
    plot_ = 
      ggplot(data = agg_,
             aes(x = SAMPLE_ID, y = COPHENETIC_DISTANCE, colour = SAMPLE_ID))+
      geom_boxplot()+
      scale_colour_manual(values = cbbPalette[-2])+
      ylim(NA, max(agg_$COPHENETIC_DISTANCE) * 1.1)+
      theme_bw()+
      labs(x = 'Immunisation',
           y = 'Cophenetic distance',
           colour = 'Immunisation')+
      theme(text = element_text(size = 15),
            legend.position = 'bottom')
    
    if(plot_bar == T){
      
      print(
        plot_ + 
          annotate('segment',
                   x = 1, xend = 2,
                   y = max(agg_$COPHENETIC_DISTANCE) * 1.05,
                   yend = max(agg_$COPHENETIC_DISTANCE) * 1.05)+
          annotate('segment',
                   x = c(1, 2), xend = c(1, 2),
                   y = max(agg_$COPHENETIC_DISTANCE) * 1.05,
                   yend = max(agg_$COPHENETIC_DISTANCE) * 1.035)
      )
      
    }
    
    print(test_)
     
  } else if(mode == 'Bulk'){
    
    cophens_ = cophens_ %>% 
      select(COPHENETIC_DISTANCE) %>% rownames_to_column('ID') %>%
      mutate(MOUSE_ID = gsub('\\..*?$', '', ID)) %>% 
      mutate(SAMPLE_ID = ifelse(grepl('59[0|1]', MOUSE_ID), 'gp120', 
                                ifelse(grepl('59[2|3]',
                                             MOUSE_ID), 'gp120 + mAb', 'gp120 + Rapamycin')))
    
    agg_ = aggregate(COPHENETIC_DISTANCE ~ MOUSE_ID + SAMPLE_ID, cophens_, mean)
    
    test_ = t.test(agg_ %>% filter(SAMPLE_ID == 'gp120') %>% pull(COPHENETIC_DISTANCE),
                   agg_ %>% filter(SAMPLE_ID == 'gp120 + Rapamycin') %>% pull(COPHENETIC_DISTANCE))
    
    print(test_)
    
    if(test_$p.value < 0.05){
      
      plot_bar = T
      
    } else {
      
      plot_bar = F
      
    }
    
    plot_ = ggplot(data = agg_, aes(x = SAMPLE_ID, y = COPHENETIC_DISTANCE, colour = SAMPLE_ID))+
      geom_boxplot()+
      scale_colour_manual(values = cbbPalette[-2])+
      theme_bw()+
      ylim(NA, max(agg_$COPHENETIC_DISTANCE) * 1.1)+
      labs(x = 'Immunisation',
           y = 'Cophenetic distance',
           colour = 'Immunisation')+
      theme(text = element_text(size = 15))
    
    if(plot_bar == T){
      plot_ = 
        plot_ + 
          annotate('segment',
                   x = 1, xend = 2,
                   y = max(agg_$COPHENETIC_DISTANCE) * 1.05,
                   yend = max(agg_$COPHENETIC_DISTANCE) * 1.05)+
          annotate('segment',
                   x = c(1, 2), xend = c(1, 2),
                   y = max(agg_$COPHENETIC_DISTANCE) * 1.05,
                   yend = max(agg_$COPHENETIC_DISTANCE) * 1.035)
    
    }
    
      print(plot_)
  
  }
  
  if(sim == 'kidera'){
    
    max_dist = as.numeric(unlist(sapply(mats, function(x) x[,2:ncol(x)])))
    max_dist = ceiling(max_dist[which.max(max_dist)])
    
    c_ = colorRampPalette(brewer.pal(n = 8, name = "Blues"))(100)
    b_ = seq(0, max_dist, by = max_dist/100)
    
  } else {

    n = 100
    c_ = colorRampPalette(brewer.pal(n = 8, name = "Blues"))(n)
    b_ = seq(0, 1, 1 / n)
    
  }
  
  if(plot_heats == T){
  
    if(mode == 'Bulk'){
    
      heatmaps[[sim]] <- 
        lapply(mats[c(6, 15)], function(x){
          hm = pheatmap::pheatmap(x[,2:ncol(x)],
                                    show_rownames = F,
                                    show_colnames = F,
                                    border_color = NA,
                                    color = c_,
                                    breaks = b_)
          
          return(hm)
          
        })
      
    } else {
      
      heatmaps[[sim]] <- 
        lapply(mats, function(x){
          hm = pheatmap::pheatmap(x[,2:ncol(x)],
                                  show_rownames = F,
                                  show_colnames = F,
                                  border_color = NA,
                                  color = c_,
                                  breaks = b_)
          
          return(hm)
          
        })
      
    }
    
  }

}


