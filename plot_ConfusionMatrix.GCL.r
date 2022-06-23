plot_ConfusionMatrix.GCL <- function(confusion_matrix, plot = c("group_group", "pop_group", "pop_pop")[1], group_names = NULL, pop_names = NULL, high_color = "#132B43", low_color = "#56B1F7", text_color = "yellow"){  
  
  ####################################################################################################################################################################################################################################################################
  #
  # This function takes the output object from ConfusionMatrices.GCL() and plots a confusion matrix in a heatmap using ggplot2::geom_tile().
  # 
  # Input parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  # confusion_matrix - a list of confusion matrix tibble(s) produced by ConfusionMatrices.GCL()
  #
  # plot - a chactater vector of lengh == 1 indicating the type of matrix you want to plot.
  #
  # group_names - an optional vector of group names used to keep groups in order; if NULL, the groups will be ordered alphabetically
  #
  # pop_names - an optional vector of population names used to keep pops in order; if NULL, the pops will be ordered alphabetically
  #
  # high_color - the color used for high "mean scale likelihood values by ggplot2:: scale_fill_gradient(); either supply hexadecimal or color name from colors()
  #
  # low_color - the color used for low values by ggplot2:: scale_fill_gradient(); either supply hexadecimal or color name from colors()
  # 
  # text_color - the color used by ggplot2::geom_text() to color the Mean Scaled Likelihood values on the plot.
  #
  # Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  # A heatmap of mean scaled likelihood
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  # baseline <- read_csv("V:/Analysis/5_Coastwide/Chum/NPen2WA_Chum_baseline/rubias/baseline/NPen2Wa_Chum_227pops_91loci_base.csv") 
  # 
  # ConfusionMatrices_out <- ConfusionMatrices.GCL(reference = baseline , gen_start_col = 5, output = c("group_group", "pop_group", "pop_pop"))
  # 
  # groups <- baseline$repunit %>% unique()
  # pops <- baseline$collection %>% unique()
  #
  # plot_ConfusionMatrix.GCL(confusion_matrix = ConfusionMatrices_out, group_names = groups, pop_names = NULL, plot = c("group_group", "pop_group", "pop_pop")[1], high_color = "darkblue", low_color = "lightblue", text_color = "yellow")
  #
  ####################################################################################################################################################################################################################################################################

  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse)  # Install packages, if not in library and then load them

  if(!plot %in% names(confusion_matrix)){
    
    stop(paste0("The supplied confusion_matrix list object does not contain a tibble named '", plot, "'"))
    
  }
  
  if(length(plot) > 1){stop("Only one type of matrix can be plotted at a time with this function, make sure length(plot)==1")}
  
  cm <- confusion_matrix[[plot]] # subset confusion_matrix
  
  # Group by group
  
  if(plot == "group_group"){
    
    if(is.null(group_names)){group_names <- cm$group %>% unique()}
    
    grcheck <- dplyr::setdiff(group_names, cm$group %>% unique()) 
    
    if(length(grcheck)>=1){
      
      stop(paste0("The following reporting groups in group_names do not exist in confusion_matrix: ", paste0(grcheck, collapse = ", ")))
      
    }
    
    hm <- cm %>% 
      dplyr::mutate(group = factor(x = group, levels = group_names), inferred_group = factor(inferred_group, levels = group_names)) %>% 
      ggplot2::ggplot(aes(x = group, y = inferred_group, z = mean_group_group_scaled_like, fill = mean_group_group_scaled_like)) +
      ggplot2::geom_tile()+
      ggplot2::scale_fill_gradient(high = high_color, low = low_color, name = "Mean Scaled Likelihood") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::geom_text(aes(x = group, y = inferred_group, label = round(mean_group_group_scaled_like, digits = 3)), color = text_color, show.legend = FALSE, size = 2.5)+
      ggplot2::ggtitle(label = "Group by Group confusion matrix", subtitle = paste0(length(group_names)," reporting groups"))+
      ggplot2::xlab("Reporting Group") +
      ggplot2::ylab("Inferred Reporting Group")
   
  }
  
  # Pop by group
  
  if(plot == "pop_group"){
    
    if(is.null(group_names)){group_names <- cm$inferred_group %>% unique()}
    
    grcheck <- dplyr::setdiff(group_names, cm$inferred_group %>% unique()) 
    
    if(length(grcheck)>=1){
      
      stop(paste0("The following reporting groups in group_names do not exist in confusion_matrix: ", paste0(grcheck, collapse = ", ")))
      
    }
    
    if(is.null(pop_names)){pop_names <- cm$pop %>% unique()}
    
    popcheck <- dplyr::setdiff(pop_names, cm$pop %>% unique()) 
    
    if(length(popcheck)>=1){
      
      stop(paste0("The following populations in pop_names do not exist in confusion_matrix: ", paste0(popcheck, collapse = ", ")))
      
    }
    
    hm <- cm %>% 
      dplyr::mutate(pop = factor(x = pop, levels = pop_names), inferred_group = factor(inferred_group, levels = group_names)) %>% 
      ggplot2::ggplot(aes(x = pop, y = inferred_group, z = mean_pop_group_scaled_like, fill = mean_pop_group_scaled_like)) +
      ggplot2::geom_tile()+
      ggplot2::scale_fill_gradient(high = high_color, low = low_color, name = "Mean Scaled Likelihood") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::geom_text(aes(x = pop, y = inferred_group, label = round(mean_pop_group_scaled_like, digits = 3)), color = text_color, show.legend = FALSE, size = 2.5)+
      ggplot2::ggtitle(label = "Pop by Group confusion matrix", subtitle = paste0(length(group_names)," reporting groups and ", length(pop_names), "populations"))+
      ggplot2::xlab("Population") +
      ggplot2::ylab("Inferred Reporting Group")
    
  }
  
  # Pop by pop
  
  if(plot == "pop_pop"){
  
    if(is.null(pop_names)){pop_names <- cm$pop %>% unique()}
    
    popcheck <- dplyr::setdiff(pop_names, cm$pop %>% unique()) 
    
    if(length(popcheck)>=1){
      
      stop(paste0("The following populations in pop_names do not exist in confusion_matrix: ", paste0(grcheck, collapse = ", ")))
      
    }
    
    hm <- cm %>% 
      dplyr::mutate(pop = factor(pop, levels = pop_names), inferred_pop = factor(inferred_pop, levels = pop_names)) %>% 
      ggplot2::ggplot(aes(x = pop, y = inferred_pop, z = mean_pop_pop_scaled_like, fill = mean_pop_pop_scaled_like)) +
      ggplot2::geom_tile()+
      ggplot2::scale_fill_gradient(high = high_color, low = low_color, name = "Mean Scaled Likelihood") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::geom_text(aes(x = pop, y = inferred_pop, label = round(mean_pop_pop_scaled_like, digits = 3)), color = text_color, show.legend = FALSE, size = 2.5)+
      ggplot2::ggtitle(label = "Pop by Pop confusion matrix", subtitle = paste0(length(pop_names), "populations"))+
      ggplot2::xlab("Population") +
      ggplot2::ylab("Inferred Population")
    
  }
  
  print(hm) # Print the plot
  
}
  
  
  