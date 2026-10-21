#' build functions to generate plots for this project
#' functional enrichment plots
#' SNP finemapping plots
#' gene finemapping plots
#' 
#' 

# theme setup for publication
theme_pub <- function(base_size = 8, base_family = "Helvetica") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_line(color = "grey90", size = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      axis.line = element_line(color = "black", size = 0.4),
      axis.ticks = element_line(color = "black", size = 0.4),
      axis.text = element_text(size = base_size, color = "black"),
      axis.title = element_text(size = base_size + 1, color = "black"),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 2),
      legend.key = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0),
      plot.margin = margin(1, 1, 1, 1)
    )
}


enrichment_bubble_plot <- function(data, use_tauP = NULL, 
                                   xTextAngle = FALSE,
                                   fdr_method = "bonferroni"){
  if(!(is.null(use_tauP))){
    data <- data %>% mutate(p_value= as.numeric(tau_star_p))
  }else{
    data <- data %>% mutate(p_value = Enrichment_p)
  }
  ggplot(
    data %>% mutate(FDR = p.adjust(p_value, method = fdr_method),
                    Enrichment = ifelse(Enrichment > 0, Enrichment, NA)),
    aes(
      x = trait,
      y = Category,
      colour = -log10(FDR),
      size = Enrichment
      #label = ifelse(p.adjust(p_value, method = fdr_method) <= 0.05, "*", "")
    )
  ) +
    scale_colour_viridis_b(alpha = 0.9) + 
    # scale_color_gradientn(
    #   colours = c("#998ec3","#f7f7f7","#f1a340"),
    #   values = rescale(c(0, 1, 3))) +
    scale_size_continuous() + 
    geom_point() + xlab("") + ylab("") + 
    theme(axis.text.x = element_text(size = 16, angle = 45))
    #geom_text(size = 7, nudge_y = 0.1, color = "brown") + 
    # theme_ArchR(baseSize = 16, legendTextSize = 16, 
    #             xText90 = xTextAngle,
    #             legendPosition = "right")
}

snp_enrichment_plot<-function(output.tbl, log.based=FALSE, y.order = NULL, intercept=FALSE,
                              y.label = NULL,split.data="Category", tolabel="pval", label.size=7,
                              bar.width=0.1, text.size=12, dot.size=4, nudgeX = 1, textVjust = 0.8, 
                              trait=""){
  ## strip off column names
  stopifnot(dim(output.tbl)[2] > 3)

  if(is.null(y.label)){
    y.label=waiver()
  }
  if(intercept==TRUE){
    output.tbl<-output.tbl[-1,]
  }
  if(is.null(y.order)){
    y.order=waiver()
  }else{
    ## order the terms of the enrichment or change the names shown in y-axis
    output.tbl$Category<-factor(output.tbl$Category, levels=y.order)
  }
  
  if(log.based){
    ##Torus results
    if(tolabel=="h2g"){
      p<-ggplot(output.tbl, aes(y=Category, x=log2(exp(estimate)), color=Category, label=h2g))
      #p<-ggplot(output.tbl, aes(y=Category, x=exp(estimate), color=Category, label=h2g))
    
      p<-p + geom_text(vjust=-0.4, hjust=0.1, nudge_x=nudgeX, size=3)
      
    }else{
      p<-ggplot(output.tbl, aes(y=Category, x=log2(exp(estimate)), colour=Category))#get0(split.data)))
    }
    # show estimate as a point and the condifence interval as an error bar
    p <- p + geom_point(size=2) + geom_errorbar(aes(xmax=log2(exp(high)), 
                                                    xmin=log2(exp(low))),
                                                    width=bar.width) +
             theme(legend.position = "none") +
             xlab(sprintf("%s variant enrichment(log2)", trait))  
  }else{
    ## S-LDSC results
    # if tolabel=TRUE, label the plot with percent of heritability explained.
    output.tbl$Category<-factor(output.tbl$Category)
    labels<- lapply(levels(output.tbl$Category), function(i) {
      prop <- round(unique(as.numeric(output.tbl[output.tbl$Category==i, "Prop._SNPs"])) * 100, 1)
      sprintf("%s (%s%%)", i, as.character(prop))
    })
    output.tbl$Category<-factor(output.tbl$Category, labels = labels)
    
    if(tolabel=="h2g"){

      #output.tbl["h2g"]<-paste0(round(output.tbl$Prop._h2*100,1), "%")
      output.tbl <-
        output.tbl %>%
        mutate(h2g = sprintf(
          "%s (%s)",
          paste0(round(output.tbl$Prop._h2 * 100), "%"),
          format(Enrichment_p, scientific = TRUE, digits = 1)
        ))
      p <-
        ggplot(output.tbl,
               aes(
                 x = Enrichment,
                 y = Category,
                 color = get0(split.data),
                 label = h2g
               ))
      p <-
        p + geom_text(
          vjust = textVjust,
          hjust = 0.5,
          nudge_x = nudgeX,
          size = label.size
        ) +
        theme(legend.position = "none")
      
      }else if(tolabel=="pval"){
        p<-ggplot(output.tbl, aes(x=Enrichment, y=Category, 
                                  color = get0(split.data),
                                  label = format(Enrichment_p, scientific = TRUE, digits = 2)))
        p<-p + geom_text(vjust=textVjust, hjust=0.5, size=label.size) + 
          theme(legend.position = "none")
        
      }else if(tolabel=="tau_p"){
        p<-ggplot(output.tbl, aes(x=Enrichment, y=Category, 
                                  color = get0(split.data),
                                  # p-value of tau indicates cell-type specific contributions
                                  label = sprintf("%s(%s)", tau_star_p, round(tau_star,1))))
                                    # format(2*pnorm(-abs(Coefficient_z.score)), 
                                    #              scientific = TRUE, digits = 1)))
        p<-p + geom_text(vjust=textVjust, hjust=0.5, size=label.size) + 
          theme(legend.position = "none")
        
      }else if(tolabel=="pval.sig"){
        p<-ggplot(output.tbl, aes(x=Enrichment, y=Category, color=ifsignif)) + 
                 scale_colour_discrete(name="", labels=c('p<0.05','o.w.'))
      }else{
        p<-ggplot(output.tbl, aes(x=Enrichment, y=Category, color=Category)) +
          theme(legend.position = "none")
        }
    
      p <- p + geom_point(size = dot.size) + 
    # p <- p + geom_col(aes(x = Enrichment,
    #                       y = Category,
    #                       fill = get0(split.data)),
    #                       width=0.5) +
      geom_errorbar(aes(xmin=Enrichment-Enrichment_std_error,
                        xmax=Enrichment+Enrichment_std_error,
                        width=bar.width,
                        )) + 
      geom_vline(xintercept = 1, linetype="dashed") + 
      xlab(sprintf("%s variant enrichment", trait))
      #xlab("")
    
    
    #p<- ggplot(output.tbl, aes(y=term, x=estimate, colour=get0(split.data))) +
      # show estimate as a point and the condifence interval as an error bar
      #geom_point(size=2) + geom_linerange(aes(xmax=high, xmin=low)) +
      #xlab(sprintf("%s variant enrichment", trait))
  }
  
  
   ## customize axis
  q<- p + #xlab(sprintf("%s variant enrichment(log2)", trait)) +
          scale_y_discrete(labels=y.label) + ylab("") +
          #change legend label
          labs(color=split.data) + 
          # add vertical line for no enrichment
          # geom_vline(xintercept=0, linetype="dashed", color="black") +  # LDSC
          # geom_vline(xintercept=0, linetype="dashed", color="black") +  # Torus
          # customize theme
          theme(axis.text=element_text(size=text.size), 
                axis.title=element_text(size=text.size),
                strip.text.x=element_text(size = 20, face="bold"),
                panel.background = element_rect(fill="white"), 
                panel.spacing.x = unit(2, "lines"),
                panel.spacing.y = unit(1, "lines"),
                panel.grid.major.y=element_blank(),
                panel.grid.major.x=element_line(colour="grey", linetype = "dotted"),
                panel.border = element_rect(fill = NA, 
                                            colour = "black"))
  return(q)
}

plot_torus_enrichment<-function(df, text_angle=45, vjust=+0.7){
    ggplot(df, aes(x=term, y=log2(exp(est)), fill= highlight)) +
            geom_bar(stat="identity") +
            geom_errorbar(aes(ymin=log2(exp(low)),
                             ymax=log2(exp(high)))) + 
            ylab("log2 Enrichment") + xlab("") + 
            theme(axis.text.x=element_text(angle=text_angle, vjust = +0.7)) + 
            scale_fill_manual(values=c("FALSE"="grey", "TRUE"="orange"), guide=F) 
  }


plot_ldsc_enrichment_old<-function(df, tolabel=TRUE, tohighlight=NULL){
  if(is.null(tohighlight)){
    tohighlight=waiver()
    df["fill"]="TRUE"
  }else{
    #df$Category[(dim(df)[1]-length(tohighlight)+1):dim(df)[1]]<-tohighlight
    df["fill"]<-df$Category %in% tohighlight
  }
  # if tolabel=TRUE, label the plot with percent of heritability explained.
  if(tolabel){
    df["h2g"]<-paste0(round(df$Prop._h2*100,2), "%")
    p1<-ggplot(df, aes(x=Category, y=Enrichment, color=Category, label=h2g))
    p1<-p1 + geom_text(hjust=0, nudge_x=0.05) 
  }else{
    p1<-ggplot(df, aes(x=Category, y=Enrichment, color=Category))
  }
  
  p1<- p1 + geom_point() + 
    geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error,
                      ymax=Enrichment+Enrichment_std_error,
                      width=.1)) + 
    theme(axis.text.x=element_text(angle=90), legend.position = "none") +
    geom_hline(yintercept = 1, linetype="dashed")
  #scale_fill_manual(values=c("FALSE"="grey", "TRUE"="orange"), guide=F) 
  return(p1)
}

plot_ldsc_enrichment<-function(df, tolabel=TRUE, y.order=NULL){
  
  if(is.null(y.order)){
    y.order=waiver()
  }else{
    ## order the terms of the enrichment or change the names shown in y-axis
    df$Category<-factor(output.tbl$Category, levels=y.order)
  }
  
  # if tolabel=TRUE, label the plot with enrichment p-values
  if(tolabel){
    df["h2g"]<-paste0(round(df$Prop._h2*100,1), "%")
    p1<-ggplot(df, aes(x=Category, y=Enrichment, color=Category, label=h2g))
    p1<-p1 + geom_text(hjust=0, nudge_x=0.05) 
  }else{
    p1<-ggplot(df, aes(x=Category, y=Enrichment, color=Category))
  }
  
  p1<- p1 + geom_point() + 
    geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error,
                      ymax=Enrichment+Enrichment_std_error,
                      width=.1)) + 
    theme(axis.text.x=element_text(angle=45,vjust=+0.7), legend.position = "none") +
    geom_hline(yintercept = 1, linetype="dashed")
  #scale_fill_manual(values=c("FALSE"="grey", "TRUE"="orange"), guide=F) 
  return(p1)
}

plot_cs_size<-function(tbl){
  #plot the distribution of credible sizes for the fine-mapped GWAS loci
  # Numbers on the stacked barplot denote the number of independent loci with
  # the corresponding credible size 
  sizes <- tbl$cs_size
  pips <- tbl$susie_pip
  sizes_breaks <- cut(sizes, breaks = c(-1, 0, 1, 5, 10, 20, 50, 200))
  sizes_breaks <-
    factor(
      sizes_breaks,
      levels = c(
        "(0,1]",
        "(1,5]",
        "(5,10]",
        "(10,20]",
        "(20,50]",
        "(50,200]",
        "(-1,0]"
      )
    )
  tbl["sizes_breaks"]<-sizes_breaks
  
  loci.by.breaks <- lapply(names(table(sizes_breaks)), function(i) {
    return(unique(tbl[which(sizes_breaks == i), "locus"]))
  })
  num.loci<-unlist(lapply(loci.by.breaks, length))
  # For the the break of (-1,0], we need to subtract all the loci 
  # covered in other breaks from the total number of loci.
  num.unresolved<-length(unique(tbl$locus))-length(unique(unlist(loci.by.breaks[1:(length(loci.by.breaks)-1)])))
  num.loci[length(num.loci)]<-num.unresolved
  gginput <- data.frame(table(sizes_breaks), num.loci)
  
  p<-ggplot(data = gginput, aes(x = factor(1), y = num.loci, fill = sizes_breaks)) +
    geom_bar(position = position_stack(reverse = T),
             stat = "identity",
             width = 0.2) + coord_flip() +
    scale_fill_brewer(
      palette = 2,
      direction = -1,
      name = "",
      labels = c("1", "2-5", "6-10", "11-20", "21-50", "50+", "not resolved")
    )  +
    geom_text(aes(label = num.loci), position = position_stack(vjust = 0.5, reverse =T)) +
    xlab("") + ylab("") + ggtitle("Number of variants in each signal") +
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_text(size = 10)
    ) 
  return(p)
}
