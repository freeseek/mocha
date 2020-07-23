#!/usr/bin/env Rscript
###
#  The MIT License
#
#  Copyright (C) 2017-2020 Giulio Genovese
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

library(optparse)
library(ggplot2)
options(bitmapType = 'cairo')

parser <- OptionParser('usage: summary_plot.R [options] --stats <file.tsv> --calls <file.tsv> --pdf <file.pdf>')
parser <- add_option(parser, c('--stats'), type = 'character', help = 'input MoChA stats file', metavar = '<file.tsv>')
parser <- add_option(parser, c('--calls'), type = 'character', help = 'input MoChA calls file', metavar = '<file.tsv>')
parser <- add_option(parser, c('--pdf'), type = 'character', help = 'output PDF file', metavar = '<file.pdf>')
parser <- add_option(parser, c('--width'), type = 'integer', default = 7, help = 'inches width of the output file [7]', metavar = '<integer>')
parser <- add_option(parser, c('--height'), type = 'integer', default = 7, help = 'inches height of the output file [7]', metavar = '<integer>')
parser <- add_option(parser, c('--fontsize'), type = 'integer', default = 16, help = 'font size [16]', metavar = '<integer>')
args <- parse_args(parser, commandArgs(trailingOnly = TRUE))

if (is.null(args$stats)) {print_help(parser); stop('option --stats is required')}
if (is.null(args$calls)) {print_help(parser); stop('option --calls is required')}
if (is.null(args$pdf)) {print_help(parser); stop('option --pdf is required')}

df_stats <- read.table(args$stats, sep = '\t', header = TRUE)
df_calls <- read.table(args$calls, sep = '\t', header = TRUE)
df_calls$chrom <- as.factor(gsub('^chr', '', gsub('^chrM', 'MT', df_calls$chrom)))
ord <- order(as.numeric(gsub('MT', '26', gsub('Y', '24', gsub('X', '23', levels(df_calls$chrom))))))
df_calls$chrom <- factor(df_calls$chrom, levels(df_calls$chrom)[ord])

if ('lrr_median' %in% colnames(df_stats)) {
  model <- 'array'
  df_stats$x_nonpar_adj_lrr_median <- df_stats$x_nonpar_lrr_median - df_stats$lrr_median
  df_stats$y_nonpar_adj_lrr_median <- df_stats$y_nonpar_lrr_median - df_stats$lrr_median
  df_stats$mt_adj_lrr_median <- df_stats$mt_lrr_median - df_stats$lrr_median
} else {
  model <- 'wgs'
  df_stats$x_nonpar_adj_lrr_median <- log(df_stats$x_nonpar_cov_median) - log(df_stats$cov_median)
  df_stats$y_nonpar_adj_lrr_median <- log(df_stats$y_nonpar_cov_median) - log(df_stats$cov_median)
  df_stats$mt_adj_lrr_median <- log(df_stats$mt_cov_median) - log(df_stats$cov_median)
}

df_calls$sv <- factor(df_calls$chrom, levels = c('0-2 Mbp', '2-10 Mbp', '10-50 Mbp', '50-250 Mbp'))
df_calls$sv[df_calls$length <= 250e6] <- '50-250 Mbp'
df_calls$sv[df_calls$length <= 50e6] <- '10-50 Mbp'
df_calls$sv[df_calls$length <= 10e6] <- '2-10 Mbp'
df_calls$sv[df_calls$length <= 2e6] <- '0-2 Mbp'
df_calls$bdev[is.na(df_calls$bdev)] <- -.05

pdf(args$pdf, width = args$width, height = args$height)

idx <- !( df_calls$chrom %in% c('X', 'Y', 'MT') )
p <- ggplot(df_calls[idx,], aes(x=bdev, y=rel_cov, color=type, shape=type)) +
  geom_hline(yintercept = c(1.0, 2.0, 3.0), color = 'gray', size = .5, linetype = 'dashed') +
  geom_segment(aes(x = 0.0, y = 2.0, xend = 1.0/6.0, yend = 3.0), color = 'gray', size = .5, linetype = 'dashed') +
  geom_segment(aes(x = 0.0, y = 2.0, xend = 1.0/6.0, yend = 1.5), color = 'gray', size = .5, linetype = 'dashed') +
  geom_point(size = 1, alpha = 1/2) +
  scale_color_manual('', values = c('CN-LOH' = 'orange', 'CNP Deletion' = 'lightblue', 'CNP Duplication' = 'violetred', 'Deletion' = 'blue', 'Duplication' = 'red', 'Undetermined' = 'gray50')) +
  scale_shape_manual('', values = 0:5) +
  theme_bw(base_size = args$fontsize) +
  theme(strip.background = element_rect(color=NA, fill=NA), legend.position = 'bottom', legend.box = 'horizontal') +
  facet_wrap(~sv)
print(p + scale_x_continuous('BAF deviation (Bdev)', breaks = c(-0.05, 0.0, 0.05, 0.1, 0.15, 0.2)) + scale_y_continuous(expression(paste('Relative coverage (rescaled ', 2^LRR, ')')), breaks = c(0, 1, 2, 3, 4)) + coord_cartesian(xlim = c(-0.05, 0.2), ylim = c(0, 4)))
print(p + scale_x_continuous('BAF deviation (Bdev)', breaks = c(0.0, 0.01, 0.02, 0.03, 0.04, 0.05)) + scale_y_continuous(expression(paste('Relative coverage (rescaled ', 2^LRR, ')')), breaks = c(1.8, 1.9, 2.0, 2.1, 2.2)) + coord_cartesian(xlim = c(0.0, 0.05), ylim = c(1.8, 2.2)))

idx <- !( df_calls$chrom %in% c('X', 'Y', 'MT') ) & !( df_calls$type %in% c('CNP Deletion', 'CNP Duplication') ) & (df_calls$length > 5e5 | df_calls$bdev<.12 | df_calls$rel_cov<2.5 & df_calls$rel_cov>1.5) # non-germline events
p <- ggplot(df_calls[idx,], aes(x=chrom, fill=type)) +
  geom_bar(stat = 'count', color = 'black') +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual('', values = c('CN-LOH' = 'orange', 'Deletion' = 'blue', 'Duplication' = 'red', 'Undetermined' = 'gray50')) +
  theme_bw(base_size = args$fontsize) +
  theme(legend.position = 'bottom', legend.box = 'horizontal')
print(p)

p <- ggplot(df_stats, aes(x=x_nonpar_adj_lrr_median, y=x_nonpar_nhets, color=mocha_gender)) +
  geom_point(size = .5, alpha = 1/2) +
  scale_x_continuous('X nonPAR median LRR (autosome corrected)') +
  scale_y_continuous('X nonPAR number of heterozygous sites') +
  scale_color_manual('', values = c('M' = 'blue', 'F' = 'orchid', 'U' = 'gray'), labels = c('M' = 'Male', 'F' = 'Female', 'U' = 'Undetermined')) +
  theme_bw(base_size = args$fontsize) +
  theme(legend.position = 'bottom', legend.box = 'horizontal')
print(p)

if ('lrr_auto' %in% colnames(df_stats)) { col_x <- 'lrr_auto'; lbl_x <- 'GC-adjusted LRR auto-correlation'
} else if ('cov_auto' %in% colnames(df_stats)) { col_x <- 'cov_auto'; lbl_x <- 'GC-adjusted coverage auto-correlation' }
p <- ggplot(df_stats, aes_string(x = col_x, y = 'baf_auto', color = 'mocha_gender')) +
  geom_point(size = .5, alpha = 1/2) +
  scale_x_continuous(lbl_x) +
  scale_y_continuous('BAF auto-correlation') +
  scale_color_manual('', values = c('M' = 'blue', 'F' = 'orchid', 'U' = 'gray'), labels = c('M' = 'Male', 'F' = 'Female', 'U' = 'Undetermined')) +
  theme_bw(base_size = args$fontsize) +
  theme(legend.position = 'bottom', legend.box = 'horizontal')
print(p)

if ('baf_sd' %in% colnames(df_stats)) { col_x <- 'baf_sd'; lbl_x <- 'Standard deviation BAF'
} else if ('baf_corr' %in% colnames(df_stats)) { col_x <- 'baf_corr'; lbl_x <- 'Beta-binomial correlation BAF' }
if ('lrr_sd' %in% colnames(df_stats)) { col_y <- 'lrr_sd'; lbl_y <- 'Standard deviation LRR'
} else if ('cov_sd' %in% colnames(df_stats)) { col_y <- 'cov_sd'; lbl_y <- 'Standard deviation coverage' }
p <- ggplot(df_stats, aes_string(x = col_x, y = col_y, color = 'mocha_gender')) +
  geom_point(size = .5, alpha = 1/2) +
  scale_x_continuous(lbl_x) +
  scale_y_continuous(lbl_y) +
  scale_color_manual('', values = c('M' = 'blue', 'F' = 'orchid', 'U' = 'gray'), labels = c('M' = 'Male', 'F' = 'Female', 'U' = 'Undetermined')) +
  theme_bw(base_size = args$fontsize) +
  theme(legend.position = 'bottom', legend.box = 'horizontal')
print(p)

invisible(dev.off())
