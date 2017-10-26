# adapted from https://github.com/joey711/phyloseq/issues/143

library('phyloseq')
data('GlobalPatterns')



set.seed(42)

calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

rarefaction_curve_data <- calculate_rarefaction_curves(ps, c('Observed', 'Shannon'), rep(c(1, 10, 100, 1000, 1:100 * 10000), each = 10))
summary(rarefaction_curve_data)

rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, 
                                        c('Depth', 'Sample', 'Measure'), summarise, 
                                        Alpha_diversity_mean = mean(Alpha_diversity), 
                                        Alpha_diversity_sd = sd(Alpha_diversity))

rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(psdata)), by.x = 'Sample', by.y = 'row.names')

library('ggplot2')

ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    colour = SampleType,
    group = Sample
  )
) + geom_line(
) + geom_pointrange(
) + facet_wrap(
  facets = ~ Measure,
  scales = 'free_y'
)