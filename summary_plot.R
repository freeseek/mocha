#!/usr/bin/env Rscript
###
#  The MIT License
#
#  Copyright (C) 2017-2022 Giulio Genovese
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

summary_plot_version <- '2022-12-21'

library(optparse)
library(ggplot2)
options(bitmapType = 'cairo')

parser <- OptionParser('usage: summary_plot.R [options] --stats <file.tsv> --calls <file.tsv> --pdf <file.pdf>')
parser <- add_option(parser, c('--n-chrs'), type = 'integer', default = 23, help = 'number of chromosomes including X [23]', metavar = '<integer>')
parser <- add_option(parser, c('--stats'), type = 'character', help = 'input MoChA stats file', metavar = '<file.tsv>')
parser <- add_option(parser, c('--calls'), type = 'character', help = 'input MoChA calls file', metavar = '<file.tsv>')
parser <- add_option(parser, c('--call-rate-thr'), type = 'double', default = 0.97, help = 'minimum call rate threshold [0.97]', metavar = '<float>')
parser <- add_option(parser, c('--baf-auto-thr'), type = 'double', default = 0.03, help = 'maximum BAF autocorrelation threshold [0.03]', metavar = '<float>')
parser <- add_option(parser, c('--pdf'), type = 'character', help = 'output PDF file', metavar = '<file.pdf>')
parser <- add_option(parser, c('--width'), type = 'double', default = 7.0, help = 'inches width of the output file [7.0]', metavar = '<float>')
parser <- add_option(parser, c('--height'), type = 'double', default = 7.0, help = 'inches height of the output file [7.0]', metavar = '<float>')
parser <- add_option(parser, c('--fontsize'), type = 'integer', default = 16, help = 'font size [16]', metavar = '<integer>')
args <- parse_args(parser, commandArgs(trailingOnly = TRUE), convert_hyphens_to_underscores = TRUE)

write(paste('summary_plot.R', summary_plot_version, 'https://github.com/freeseek/mocha'), stderr())

if (is.null(args$stats)) {print_help(parser); stop('option --stats is required')}
if (is.null(args$calls)) {print_help(parser); stop('option --calls is required')}
if (is.null(args$pdf)) {print_help(parser); stop('option --pdf is required')}

df_stats <- read.table(args$stats, sep = '\t', header = TRUE)
df_calls <- read.table(args$calls, sep = '\t', header = TRUE)
df_calls$chrom <- as.factor(gsub('^chr', '', gsub('^chrM', 'MT', df_calls$chrom)))
ord <- order(as.numeric(gsub('MT', args$n_chrs + 3, gsub('Y', args$n_chrs + 1, gsub('X', args$n_chrs, levels(df_calls$chrom))))))
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

df_calls$sv <- factor(df_calls$chrom, levels = c('<0.5 Mbp', '0.5-5 Mbp', '5-50 Mbp', '>50 Mbp'))
df_calls$sv[df_calls$length <= 5e8] <- '>50 Mbp'
df_calls$sv[df_calls$length <= 5e7] <- '5-50 Mbp'
df_calls$sv[df_calls$length <= 5e6] <- '0.5-5 Mbp'
df_calls$sv[df_calls$length <= 5e5] <- '<0.5 Mbp'

# avoid plotting BDEV for calls with not enough heterozygous sites
df_calls$bdev[df_calls$n_hets < 5] <- NaN
df_calls$bdev[is.na(df_calls$bdev)] <- 0

pdf(args$pdf, width = args$width, height = args$height)

idx <- !( df_calls$sample_id %in% df_stats$sample_id[df_stats$call_rate < args$call_rate_thr | df_stats$baf_auto > args$baf_auto_thr] |
          df_calls$chrom %in% c('X', 'Y', 'MT') | grepl('^CNP', df_calls$type) )

if (sum(idx) > 0) {
  p <- ggplot(df_calls[idx,], aes(x = bdev, y = rel_cov, color = type)) +
    geom_hline(yintercept = c(1.0, 2.0, 3.0), color = 'gray', size = .5, linetype = 'dashed') +
    geom_segment(aes(x = 0.0, y = 2.0, xend = 1.0/6.0, yend = 3.0), color = 'gray', size = .5, linetype = 'dashed') +
    geom_segment(aes(x = 0.0, y = 2.0, xend = 1.0/6.0, yend = 1.5), color = 'gray', size = .5, linetype = 'dashed') +
    geom_point(shape = 20, size = .5, alpha = 1/2) +
    scale_color_manual('', values = c('CN-LOH' = 'orange', 'Loss' = 'blue', 'Gain' = 'red', 'Undetermined' = 'gray50')) +
    theme_bw(base_size = args$fontsize) +
    theme(strip.background = element_rect(color = NA, fill = NA), legend.position = 'bottom', legend.box = 'horizontal') +
    facet_wrap(~sv)
  print(p + scale_x_continuous('BAF deviation (Bdev)', breaks = c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30)) + scale_y_continuous(expression(paste('Relative coverage (rescaled ', 2^LRR, ')')), breaks = c(0, 1, 2, 3, 4)) + coord_cartesian(xlim = c(0.00, 0.30), ylim = c(0, 4)))
  print(p + scale_x_continuous('BAF deviation (Bdev)', breaks = c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06)) + scale_y_continuous(expression(paste('Relative coverage (rescaled ', 2^LRR, ')')), breaks = c(1.8, 1.9, 2.0, 2.1, 2.2)) + coord_cartesian(xlim = c(0.00, 0.06), ylim = c(1.8, 2.2)))
}

if (sum(idx & df_calls$type == 'CN-LOH') > 0) {
  p <- ggplot(df_calls[idx & df_calls$type == 'CN-LOH',], aes(x = bdev, y = n50_hets, color = type)) +
    geom_point(shape = 20, size = .5, alpha = 1/2) +
    scale_x_continuous('BAF deviation (Bdev)') +
    scale_y_log10('N50 heterozygous sites distance') +
    scale_color_manual('', values = c('CN-LOH' = 'orange')) +
    theme_bw(base_size = args$fontsize) +
    theme(strip.background = element_rect(color = NA, fill = NA), legend.position = 'bottom', legend.box = 'horizontal') +
    facet_wrap(~sv)
  print(p)
}

p <- ggplot(df_stats, aes(x = 1 - call_rate, y = baf_auto, color = computed_gender)) +
  geom_vline(xintercept = 1 - args$call_rate_thr, color = 'black', size = .5, alpha = 1/2) +
  geom_hline(yintercept = args$baf_auto_thr, color = 'black', size = .5, alpha = 1/2) +
  geom_point(data = df_stats[df_stats$call_rate < args$call_rate_thr | df_stats$baf_auto > args$baf_auto_thr,], color = 'black', shape = 1, size = .5, alpha = 1/2) +
  geom_point(shape = 20, size = .5, alpha = 1/2) +
  scale_x_log10('1 - call rate') +
  scale_y_continuous('BAF auto-correlation') +
  scale_color_manual('', values = c('M' = 'blue', 'F' = 'orchid', 'K' = 'orange', 'U' = 'gray'), labels = c('M' = 'Male', 'F' = 'Female', 'K' = 'Klinefelter', 'U' = 'Undetermined')) +
  theme_bw(base_size = args$fontsize) +
  theme(legend.position = 'bottom', legend.box = 'horizontal')
print(p)

p <- ggplot(df_stats, aes(x = x_nonpar_adj_lrr_median, y = y_nonpar_adj_lrr_median, color = computed_gender)) +
  geom_point(data = df_stats[df_stats$call_rate < args$call_rate_thr | df_stats$baf_auto > args$baf_auto_thr,], color = 'black', shape = 1, size = .5, alpha = 1/2) +
  geom_point(shape = 20, size = .5, alpha = 1/2) +
  scale_x_continuous('X nonPAR median LRR (autosome corrected)') +
  scale_y_continuous('Y nonPAR median LRR (autosome corrected)') +
  scale_color_manual('', values = c('M' = 'blue', 'F' = 'orchid', 'K' = 'orange', 'U' = 'gray'), labels = c('M' = 'Male', 'F' = 'Female', 'K' = 'Klinefelter', 'U' = 'Undetermined')) +
  theme_bw(base_size = args$fontsize) +
  theme(legend.position = 'bottom', legend.box = 'horizontal')
print(p)

if ('baf_sd' %in% colnames(df_stats)) { col_x <- 'baf_sd'; lbl_x <- 'Standard deviation BAF'
} else if ('baf_corr' %in% colnames(df_stats)) { col_x <- 'baf_corr'; lbl_x <- 'Beta-binomial correlation BAF' }
if ('lrr_sd' %in% colnames(df_stats)) { col_y <- 'lrr_sd'; lbl_y <- 'Standard deviation LRR'
} else if ('cov_sd' %in% colnames(df_stats)) { col_y <- 'cov_sd'; lbl_y <- 'Standard deviation coverage' }
p <- ggplot(df_stats, aes_string(x = col_x, y = col_y, color = 'computed_gender')) +
  geom_point(data = df_stats[df_stats$call_rate < args$call_rate_thr | df_stats$baf_auto > args$baf_auto_thr,], color = 'black', shape = 1, size = .5, alpha = 1/2) +
  geom_point(shape = 20, size = .5, alpha = 1/2) +
  scale_x_continuous(lbl_x) +
  scale_y_continuous(lbl_y) +
  scale_color_manual('', values = c('M' = 'blue', 'F' = 'orchid', 'K' = 'orange', 'U' = 'gray'), labels = c('M' = 'Male', 'F' = 'Female', 'K' = 'Klinefelter', 'U' = 'Undetermined')) +
  theme_bw(base_size = args$fontsize) +
  theme(legend.position = 'bottom', legend.box = 'horizontal')
print(p)

idx_mlox <- df_calls$computed_gender == 'F' & df_calls$chrom == 'X' & df_calls$length > 1e8
if (sum(idx_mlox) > 0) {
  df_merge <- merge(df_stats[df_stats$computed_gender == 'F' & df_stats$call_rate >= args$call_rate_thr & df_stats$baf_auto <= args$baf_auto_thr, c('sample_id', 'x_nonpar_adj_lrr_median')],
                    df_calls[idx_mlox, c('sample_id', 'bdev')], all.x = TRUE)
  df_merge$mlox <- !is.na(df_merge$bdev)
  df_merge$bdev[is.na(df_merge$bdev)] <- 0
  p <- ggplot(df_merge, aes(x = bdev, y = x_nonpar_adj_lrr_median, color = mlox)) +
    geom_point(shape = 20, size = .5, alpha = 1/2) +
    scale_x_continuous('X BAF deviation') +
    scale_y_continuous('X nonPAR median LRR (autosome corrected)') +
    scale_color_manual(paste(sum(df_merge$mlox), 'mLOX'), values = c('FALSE' = 'gray', 'TRUE' = 'red'), labels = c('FALSE' = 'no', 'TRUE' = 'yes')) +
    theme_bw(base_size = args$fontsize) +
    theme(legend.position = 'bottom', legend.box = 'horizontal')
  print(p)
}

idx_mloy <- df_calls$computed_gender == 'M' & df_calls$chrom == 'X' & df_calls$length > 2e6
if (sum(idx_mloy) > 0) {
  df_merge <- merge(df_stats[df_stats$computed_gender == 'M' & df_stats$call_rate >= args$call_rate_thr & df_stats$baf_auto <= args$baf_auto_thr, c('sample_id', 'y_nonpar_adj_lrr_median')],
                    df_calls[idx_mloy, c('sample_id', 'bdev')], all.x = TRUE)
  df_merge$mloy <- !is.na(df_merge$bdev)
  df_merge$bdev[is.na(df_merge$bdev)] <- 0
  p <- ggplot(df_merge, aes(x = bdev, y = y_nonpar_adj_lrr_median, color = mloy)) +
    geom_point(shape = 20, size = .5, alpha = 1/2) +
    scale_x_continuous('PAR1 BAF deviation') +
    scale_y_continuous('Y nonPAR median LRR (autosome corrected)') +
    scale_color_manual(paste(sum(df_merge$mloy), 'mLOY'), values = c('FALSE' = 'gray', 'TRUE' = 'red'), labels = c('FALSE' = 'no', 'TRUE' = 'yes')) +
    theme_bw(base_size = args$fontsize) +
    theme(legend.position = 'bottom', legend.box = 'horizontal')
  print(p)
}

if (sum(idx_mlox | idx_mloy) > 0) {
  p <- ggplot(df_calls[idx_mlox | idx_mloy,], aes(x = bdev)) +
    geom_histogram(binwidth = .001, color = 'black', fill = 'transparent') +
    scale_x_continuous('mLOX/mLOY BAF deviation', limits = c(0.0, 0.05), expand = c(0, 0)) +
    theme_bw(base_size = args$fontsize) +
    facet_grid(computed_gender ~ ., scales = 'free_y')
  print(p)
}

invisible(dev.off())
