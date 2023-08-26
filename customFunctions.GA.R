#' Thin summary statistics for faster plotting
#' 
#' This function makes qq-plot generation much faster. Careful: It must be applied only for visualization, not for quantile distribution generation and $\lambda$ calculation
#'
#' @param res.df \code{data.frame/tibble} containing GWAS summary statistics (essential is the p.value column) 
#' @param p.col \code{character} with the variable name in the set
#' @param log.p \code{logical} indicating whether p-values are as crude numbers or transformed as $-\log_{10}(P-value)$
#' @param p.thresh Values of p-value up to which you want to thin the points. E.g. Default is 2, so all $0 \leqslant -\log_{10}(P-value) \leqslant 2$ will be thinned. If `log.p = TRUE`, then you should provide the crude p-value, which will then be thinned in the range $p_thresh \leqslant P-value \leqslant 1$
#'
#' @return A \code{data.frame} object with fewer dots to plot

#' @export
#'
#' @examples
#' 
#' summary_stats <- data.table::fread("/path/to/regenie/summstats.regenie.gz", nThreads = 8)
#' dim(summary_stats)
#' 
#' summars_stats_thinned <- thin_gwas_summStats(summary_stats, p.col = "LOG10P", log.p = T, p-thresh = 2)

thin_gwas_summStats <- function(res.df, p.col, log.p = T, p.thresh = 2){
  if(!log.p){
    res.df <- mutate(res.df, 
                     !!sym(p.col) := -log(!!sym(p.col), 10)
                     )
    p.thresh <- -log(p.thresh, 10)
  }
  
  data0 <- filter(
    res.df,
    !!sym(p.col) <= p.thresh) %>% 
    mutate(p.cat = cut(!!sym(p.col), seq(0,p.thresh, 0.1)))
  
  max_p.n <-  data0 %>% 
    group_by(p.cat) %>%
    tally() %>% 
    .[[nrow(.),2]]
  
  res.df_under <- data0 %>% 
    group_by(p.cat) %>% 
    sample_n(max_p.n, replace = FALSE) %>% 
    ungroup() %>% 
    select(-p.cat)
    
  
  res_thinned <- rbind.data.frame(res.df_under, res.df %>% filter(!!sym(p.col) > p.thresh))
  
  return(res_thinned)
  
  
}


#' QQ plot faster
#'
#'  QQ plot generation with a faster plotting thanks to the thinning strategy 
#'  (used only in the plotting step). The rest of the function was designed by 
#'  Michele Filosi and it is part of the gitlab repository regenie_pipeline
#'
#' @param df \code{data.frame} of the summary statistics with Minor Allele Frequency (MAF) and P-value columns
#' @param mafbreaks \code{numeric} vector with the breaks to `cut` the MAF column 
#' @param mafcol \code{character}, the column containing MAF data
#' @param pvalcol \code{character}, the column containing P-value data
#' @param log \code{logical} telling whether the P-value column is in Log10 base
#' @param thin.p \code{logical} telling whether the data frame should be thinned or not. Default is `TRUE`. If set to FALSE, it is simply Michele Filosi's function
#' @param title \code{character}, The title of the plot. Default is `QQPlot`
#'
#' @return
#' @export
#'
#' @examples
myqqplot_mafbreak.GA <- function(
    df,
    mafbreaks=c(0, 0.0005, 0.005, 0.05, 1),
    mafcol = "AF_Allele2",
    pvalcol = "p.value",
    log=FALSE,
    thin.p = TRUE,
    title="QQPLOT"
){
  
  df$MAFBRK <- cut(df %>% pull({{ mafcol }}), breaks=mafbreaks)
  
  
  df <- df %>%
    group_by(MAFBRK)  %>%
    mutate(
      across(
        any_of(pvalcol),
        ~create_nulldist(.x, log=log),
        .names="EXP")
    )
  overall.lmb <- compute.inflation(df %>% pull({{ pvalcol }}), log=log)
  
  # mylambdas <- df %>% group_by(MAFBRK) %>%
  #   summarise(lambda = compute.inflation(-log10(p.value))$lambda,
  #             EXP=0, p.value=-log10(min(p.value)))
  
  if (!log){
    df <- df %>% mutate_at(vars({{pvalcol}}), ~log10)
  }
  
  if(thin.p){
    qqp <- ggplot(df %>% thin_gwas_summStats(p.col = pvalcol, log.p = TRUE, p.thresh = 2.5), aes_string("EXP", pvalcol, color="MAFBRK")) + geom_abline()  
  }else{
    qqp <- ggplot(df, aes_string("EXP", pvalcol, color="MAFBRK")) + geom_abline()
  }
  
  qqp <- qqp + geom_point() + theme_bw()
  qqp <- qqp + labs(x = "Expected -log(P.Value)", y = "Observed -log(P.Value)", color="MAF intervals") +
    theme(panel.border = element_blank(), legend.position = "right")
  qqp <- qqp + ggtitle(title)
  
  xannot <- max(df$EXP) - (max(df$EXP)*0.1)
  qqp <- qqp + annotate(
    "text",
    x = xannot,
    y = 2,
    label =
      paste0(
        "lambda: ",
        round(overall.lmb$lambda, 2)
      )
  )
  return(qqp)
}


##---- Myqqplot ----
myqqplot <- function(
    x,
    gcb = FALSE,
    log=FALSE
)
{
  n <- length(x)
  if (log){
    observed <- sort(x = x, decreasing= TRUE)
    ix <- order(x, decreasing = TRUE)
  } else {
    observed <- -log10(sort(x, decreasing = FALSE))
    ix <- order(x, decreasing = FALSE)
  }
  expected <- -log10(1:n / n)
  
  lmbd <- compute.inflation(x, log=log)$lambda
  
  xmin <- 0
  # if (sum(varres$FDR[ix]<0.05, na.rm=TRUE)) xmin <-
  # min(expected[varres$FDR[ix] < 0.05])
  
  ## QQplot
  qqp <- ggplot(data.frame(EXP = expected, OBS = observed), aes(EXP, OBS)) +
    geom_point() +
    geom_abline() +
    theme_bw()
  qqp <- qqp + labs(x = "Expected -log(P.Value)", y = "Observed -log(P.Value)") +
    theme(panel.border = element_blank(), legend.position = "right") # + guides(color=gcb)
  # qqp <- qqp + scale_color_manual(values=c('black', 'red'), labels =
  # c('FDR>0.05', 'FDR<=0.05'), name='')
  if (xmin > 0) {
    qqp <- qqp + geom_vline(xintercept = xmin, col = "grey60", linetype = "dashed")
  }
  
  qqp <- qqp + annotate("text", x = 0.5, y = max(observed), label = paste0(
    "lambda: ",
    round(lmbd, 2)
  ))
  
  return(qqp)
}