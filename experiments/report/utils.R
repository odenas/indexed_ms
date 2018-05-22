

read_many_csv <- function(flist, add_bpath=TRUE){
  message(sprintf(flist[1]))
  do.call(rbind, lapply(flist, function(f) {
    if(add_bpath)
      read_csv(f) %>% mutate(b_path=basename(f))
    else
      read_csv(f)
  }))
}

performance_boxplot <- function(p, tit){
  p + geom_boxplot(outlier.shape = 2, outlier.size = 0.5, outlier.colour = 'black') +
    geom_jitter(alpha=0.2, width=0.1) +
    scale_colour_brewer(palette = "Dark2") +
    theme(#axis.text.x = element_text(angle = 0, hjust = 1),
      legend.position="bottom") +
    labs(title = tit) + ylab("Time(ms)") + xlab("Input type")
}

prange_plot <- function(dt){
  qup <- function(x) quantile(x, 3/4)
  qdown <- function(x) quantile(x, 1/4)

  ggplot(dt, aes(inp_type, speedup)) +
    geom_pointrange(stat = "summary", fatten = 1, fun.y = median, fun.ymin = qdown, fun.ymax = qup) +
    geom_jitter(alpha=0.1, width = 0.1) +
    geom_hline(yintercept = 1.0, color='blue') +
    facet_grid(~alp) +
    #ylim(-0.1, 5) +
    ylab("Speedup") + xlab("Input type")
}


