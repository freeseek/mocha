#!/usr/bin/env Rscript
###
#  The MIT License
#
#  Copyright (C) 2017-2019 Giulio Genovese
#
#  Author: Giulio Genovese <giulio.genovese@gmail.com>
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
###

for (x in c('argparse', 'ggplot2')) {
  if ( x %in% .packages(all.available = TRUE) ) {
    suppressPackageStartupMessages(library(x, character.only = TRUE))
  } else {
    cat(paste0('package ', x, ' is not available. To install ', x, ' run:\n'))
    cat('Rscript -e \'dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)\' \\\n')
    cat(paste0('        -e \'install.packages(pkgs = "', x, '", lib = Sys.getenv("R_LIBS_USER"), repos = "https://cran.rstudio.com")\'\n'))
    q()
  }
}

parser <- ArgumentParser(description='Plot MoChA summary statistics (version 2019-08-28)')
parser$add_argument('--stats', metavar = '<file.tsv>', type = 'character', required = TRUE, help = 'input MoChA stats file')
parser$add_argument('--calls', metavar = '<file.tsv>', type = 'character', required = TRUE, help = 'input MoChA calls file')
parser$add_argument('--pdf', metavar = '<file.pdf>', type = 'character', help = 'output PDF file')
parser$add_argument('--fontsize', metavar = '<integer>', type = 'integer', default = 16, help = 'font size [16]')
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

df_stats <- read.table(args$stats, sep = '\t', header = TRUE)
df_calls <- read.table(args$calls, sep = '\t', header = TRUE)
df_calls$CHROM <- as.factor(gsub('^chr', '', gsub('^chrM', 'MT', df_calls$CHROM)))
ord <- order(as.numeric(gsub('MT', '26', gsub('Y', '24', gsub('X', '23', levels(df_calls$CHROM))))))
df_calls$CHROM <- factor(df_calls$CHROM, levels(df_calls$CHROM)[ord])

if ('LRR_MEDIAN' %in% colnames(df_stats)) {
  model <- 'array'
  df_stats$X_NONPAR_ADJ_LRR_MEDIAN <- df_stats$X_NONPAR_LRR_MEDIAN - df_stats$LRR_MEDIAN
  df_stats$Y_NONPAR_ADJ_LRR_MEDIAN <- df_stats$Y_NONPAR_LRR_MEDIAN - df_stats$LRR_MEDIAN
  df_stats$MT_ADJ_LRR_MEDIAN <- df_stats$MT_LRR_MEDIAN - df_stats$LRR_MEDIAN
} else {
  model <- 'wgs'
  df_stats$X_NONPAR_ADJ_LRR_MEDIAN <- log(df_stats$X_NONPAR_COV_MEDIAN) - log(df_stats$COV_MEDIAN)
  df_stats$Y_NONPAR_ADJ_LRR_MEDIAN <- log(df_stats$Y_NONPAR_COV_MEDIAN) - log(df_stats$COV_MEDIAN)
  df_stats$MT_ADJ_LRR_MEDIAN <- log(df_stats$MT_COV_MEDIAN) - log(df_stats$COV_MEDIAN)
}

df_calls$SV <- factor(df_calls$CHROM, levels = c('0-2 Mbp', '2-10 Mbp', '10-50 Mbp', '50-250 Mbp'))
df_calls$SV[df_calls$LENGTH <= 250e6] <- '50-250 Mbp'
df_calls$SV[df_calls$LENGTH <= 50e6] <- '10-50 Mbp'
df_calls$SV[df_calls$LENGTH <= 10e6] <- '2-10 Mbp'
df_calls$SV[df_calls$LENGTH <= 2e6] <- '0-2 Mbp'
df_calls$BDEV[is.na(df_calls$BDEV)] <- -.05

if (!is.null(args$pdf)) {
  pdf(args$pdf)
} else {
  png(args$png, width = 7, height = 7, units = 'in', res = 150)
}

idx <- !( df_calls$CHROM %in% c('X', 'Y', 'MT') )
p <- ggplot(df_calls[idx,], aes(x=BDEV, y=REL_COV, color=TYPE, shape=TYPE)) +
  geom_hline(yintercept = c(1.0, 2.0, 3.0), color = 'gray', size = .5, linetype = 'dashed') +
  geom_segment(aes(x = 0.0, y = 2.0, xend = 1.0/6.0, yend = 3.0), color = 'gray', size = .5, linetype = 'dashed') +
  geom_segment(aes(x = 0.0, y = 2.0, xend = 1.0/6.0, yend = 1.5), color = 'gray', size = .5, linetype = 'dashed') +
  geom_point(size = 1, alpha = 1/2) +
  scale_color_manual('', values = c('CN-LOH' = 'orange', 'CNP Deletion' = 'lightblue', 'CNP Duplication' = 'violetred', 'Deletion' = 'blue', 'Duplication' = 'red', 'Undetermined' = 'gray50')) +
  scale_shape_manual('', values = 0:5) +
  theme_bw(base_size = args$fontsize) +
  theme(strip.background = element_rect(color=NA, fill=NA), legend.position = 'bottom', legend.box = 'horizontal') +
  facet_wrap(~SV)
print(p + scale_x_continuous('BAF deviation (Bdev)', breaks = c(-0.05, 0.0, 0.05, 0.1, 0.15, 0.2)) + scale_y_continuous(expression(paste('Relative coverage (rescaled ', 2^LRR, ')')), breaks = c(0, 1, 2, 3, 4)) + coord_cartesian(xlim = c(-0.05, 0.2), ylim = c(0, 4)))
print(p + scale_x_continuous('BAF deviation (Bdev)', breaks = c(0.0, 0.01, 0.02, 0.03, 0.04, 0.05)) + scale_y_continuous(expression(paste('Relative coverage (rescaled ', 2^LRR, ')')), breaks = c(1.8, 1.9, 2.0, 2.1, 2.2)) + coord_cartesian(xlim = c(0.0, 0.05), ylim = c(1.8, 2.2)))

idx <- !( df_calls$CHROM %in% c('X', 'Y', 'MT') ) & !( df_calls$TYPE %in% c('CNP Deletion', 'CNP Duplication') ) & (df_calls$LENGTH > 5e5 | df_calls$BDEV<.12 | df_calls$REL_COV<2.5 & df_calls$REL_COV>1.5) # non-germline events
p <- ggplot(df_calls[idx,], aes(x=CHROM, fill=TYPE)) +
  geom_bar(stat = 'count', color = 'black') +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual('', values = c('CN-LOH' = 'orange', 'Deletion' = 'blue', 'Duplication' = 'red', 'Undetermined' = 'gray50')) +
  theme_bw(base_size = args$fontsize) +
  theme(legend.position = 'bottom', legend.box = 'horizontal')
print(p)

p <- ggplot(df_stats, aes(x=X_NONPAR_ADJ_LRR_MEDIAN, y=X_NONPAR_NHETS, color=SEX)) +
  geom_point(size = .5, alpha = 1/2) +
  scale_x_continuous('X nonPAR median LRR (autosome corrected)') +
  scale_y_continuous('X nonPAR number of heterozygous sites') +
  scale_color_manual('', values = c('M' = 'blue', 'F' = 'orchid', 'U' = 'gray'), labels = c('M' = 'Male', 'F' = 'Female', 'U' = 'Undetermined')) +
  theme_bw(base_size = args$fontsize) +
  theme(legend.position = 'bottom', legend.box = 'horizontal')
print(p)

if ('LRR_AUTO' %in% colnames(df_stats)) { col_x <- 'LRR_AUTO'; lbl_x <- 'GC-adjusted LRR auto-correlation'
} else if ('COV_AUTO' %in% colnames(df_stats)) { col_x <- 'COV_AUTO'; lbl_x <- 'GC-adjusted coverage auto-correlation' }
p <- ggplot(df_stats, aes_string(x = col_x, y = 'BAF_AUTO', color = 'SEX')) +
  geom_point(size = .5, alpha = 1/2) +
  scale_x_continuous(lbl_x) +
  scale_y_continuous('BAF auto-correlation') +
  scale_color_manual('', values = c('M' = 'blue', 'F' = 'orchid', 'U' = 'gray'), labels = c('M' = 'Male', 'F' = 'Female', 'U' = 'Undetermined')) +
  theme_bw(base_size = args$fontsize) +
  theme(legend.position = 'bottom', legend.box = 'horizontal')
print(p)

if ('BAF_SD' %in% colnames(df_stats)) { col_x <- 'BAF_SD'; lbl_x <- 'Standard deviation BAF'
} else if ('BAF_CORR' %in% colnames(df_stats)) { col_x <- 'BAF_CORR'; lbl_x <- 'Beta-binomial correlation BAF' }
if ('LRR_SD' %in% colnames(df_stats)) { col_y <- 'LRR_SD'; lbl_y <- 'Standard deviation LRR'
} else if ('COV_SD' %in% colnames(df_stats)) { col_y <- 'COV_SD'; lbl_y <- 'Standard deviation coverage' }
p <- ggplot(df_stats, aes_string(x = col_x, y = col_y, color = 'SEX')) +
  geom_point(size = .5, alpha = 1/2) +
  scale_x_continuous(lbl_x) +
  scale_y_continuous(lbl_y) +
  scale_color_manual('', values = c('M' = 'blue', 'F' = 'orchid', 'U' = 'gray'), labels = c('M' = 'Male', 'F' = 'Female', 'U' = 'Undetermined')) +
  theme_bw(base_size = args$fontsize) +
  theme(legend.position = 'bottom', legend.box = 'horizontal')
print(p)

invisible(dev.off())
