#!/usr/bin/env Rscript
###
#  The MIT License
#
#  Copyright (C) 2017-2018 Giulio Genovese
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

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
out_pdf <- args[1]
stats <- args[2]
mocha <- args[3]

df_stats <- read.table(stats, sep = '\t', header = TRUE)
df <- read.table(mocha, sep = '\t', header = TRUE)

if ('LRR_MEDIAN' %in% colnames(df_stats)) {
  model <- 'array'
  df_stats$X_NONPAR_ADJ_LRR_MEDIAN <- df_stats$X_NONPAR_LRR_MEDIAN - df_stats$LRR_MEDIAN
  df_stats$Y_NONPAR_ADJ_LRR_MEDIAN <- df_stats$Y_NONPAR_LRR_MEDIAN - df_stats$LRR_MEDIAN
  df_stats$MT_ADJ_LRR_MEDIAN <- df_stats$MT_LRR_MEDIAN - df_stats$LRR_MEDIAN
} else {
  model <- 'wgs'
  # df$LDEV <- log(df$REL_COV/2)
  df_stats$X_NONPAR_ADJ_LRR_MEDIAN <- log(df_stats$X_NONPAR_COV_MEDIAN) - log(df_stats$COV_MEDIAN)
  df_stats$Y_NONPAR_ADJ_LRR_MEDIAN <- log(df_stats$Y_NONPAR_COV_MEDIAN) - log(df_stats$COV_MEDIAN)
  df_stats$MT_ADJ_LRR_MEDIAN <- log(df_stats$MT_COV_MEDIAN) - log(df_stats$COV_MEDIAN)
}

idx <- !( df$CHROM %in% c('X', 'Y', 'MT') ) & (df$REL_COV > 0.25)
df$SV <- factor(df$CHROM, levels = c('0-2 Mbp', '2-10 Mbp', '10-50 Mbp', '50-250 Mbp'))
df$SV[df$LENGTH < 250e6] <- '50-250 Mbp'
df$SV[df$LENGTH < 50e6] <- '10-50 Mbp'
df$SV[df$LENGTH < 10e6] <- '2-10 Mbp'
df$SV[df$LENGTH <= 2e6] <- '0-2 Mbp'

bdev <- .3;
df$BDEV[is.na(df$BDEV)] <- bdev
cnf_dup <- 1 / (0.5 - bdev);
cnf_del <- 1 / (0.5 + bdev);
fs <- 16

p1 <- ggplot(df[idx,], aes(x=BDEV, y=REL_COV, color=TYPE, shape=TYPE)) +
  geom_hline(yintercept = c(1.0, 2.0, 3.0), color = 'gray', size = .5, linetype = 'dashed') +
  geom_segment(aes(x = 0, y = 2.0, xend = bdev, yend = cnf_dup), color = 'gray', size = .5, linetype = 'dashed') +
  geom_segment(aes(x = 0, y = 2.0, xend = bdev, yend = cnf_del), color = 'gray', size = .5, linetype = 'dashed') +
  geom_point(size=1, alpha=1/2) +
  scale_color_manual('', breaks = c("CNN-LOH", "CNP Deletion", "CNP Duplication", "Deletion", "Duplication", "Undetermined"), values=c("green", "orange", "cyan", "red", "blue", "gray")) +
  scale_shape_manual('', values=0:5) +
  scale_x_continuous('Bdev', breaks = c(0, .1, .2, .3)) +
  scale_y_log10('Rel_Cov', breaks = c(1, 2, 3, 4)) +
  theme_bw(base_size = fs) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~SV)

p2 <- ggplot(df[idx,], aes(x=BDEV, y=REL_COV, color=CHROM, shape=CHROM)) +
  geom_hline(yintercept = c(1.0, 2.0, 3.0), color = 'gray', size = .5, linetype = 'dashed') +
  geom_segment(aes(x = 0, y = 2.0, xend = bdev, yend = cnf_dup), color = 'gray', size = .5, linetype = 'dashed') +
  geom_segment(aes(x = 0, y = 2.0, xend = bdev, yend = cnf_del), color = 'gray', size = .5, linetype = 'dashed') +
  geom_point(size=1, alpha=1/2) +
  scale_color_discrete('Chromosome') +
  scale_shape_manual('Chromosome', values=0:23) +
  scale_x_continuous('Bdev', breaks = c(0, .1, .2, .3)) +
  scale_y_log10('Rel_Cov', breaks = c(1, 2, 3, 4)) +
  theme_bw(base_size = fs) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~SV)

p3 <- ggplot(df[idx,], aes(x=BDEV, y=REL_COV, color=(NFLIPS!=-1), shape=(NFLIPS!=-1))) +
  geom_hline(yintercept = c(1.0, 2.0, 3.0), color = 'gray', size = .5, linetype = 'dashed') +
  geom_segment(aes(x = 0, y = 2.0, xend = bdev, yend = cnf_dup), color = 'gray', size = .5, linetype = 'dashed') +
  geom_segment(aes(x = 0, y = 2.0, xend = bdev, yend = cnf_del), color = 'gray', size = .5, linetype = 'dashed') +
  geom_point(size=1, alpha=1/2) +
  scale_color_discrete('Model', labels = c('LRR+BAF', 'BAF+phase')) +
  scale_shape_manual('Model', values=0:1, labels = c('LRR+BAF', 'BAF+phase')) +
  scale_x_continuous('Bdev', breaks = c(0,.1,.2,.3)) +
  scale_y_log10('Rel_Cov', breaks = c(1,2,3,4)) +
  theme_bw(base_size = fs) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~SV)

idx <- df$LENGTH > 2e6 | df$BDEV<.12 | df$REL_COV<2.5 # non-germline events
df <- as.data.frame.matrix(table(df[idx, c('CHROM', 'TYPE')]))
df$CHROM <- rownames(df)
idx <- !( df$CHROM %in% c('X', 'Y', 'MT') )
p4 <- ggplot(df[idx,], aes(x=Deletion, y=Duplication, color=CHROM, shape=CHROM, label=CHROM)) +
  geom_point() +
  geom_text(color = 'black', hjust=-.1, vjust=-.1) +
  scale_color_discrete(guide = FALSE) +
  scale_shape_manual(values=0:23, guide=FALSE) +
  theme_bw(base_size = fs)

p5 <- ggplot(df_stats, aes(x=X_NONPAR_ADJ_LRR_MEDIAN, y=X_NONPAR_NHETS, color=SEX)) +
  geom_point(alpha=1/3) +
  scale_x_continuous('X nonPAR median LRR (autosome corrected)') +
  scale_y_continuous('X nonPAR number of heterozygous sites') +
  scale_color_discrete(guide = FALSE) +
  theme_bw(base_size = fs)
p6 <- ggplot(df_stats, aes(x=X_NONPAR_ADJ_LRR_MEDIAN, y=Y_NONPAR_ADJ_LRR_MEDIAN, color=SEX)) +
  geom_point(alpha=1/3) +
  scale_x_continuous('X nonPAR median LRR (autosome corrected)') +
  scale_y_continuous('Y nonPAR median LRR (autosome corrected)') +
  scale_color_discrete(guide = FALSE) +
  theme_bw(base_size = fs)
p7 <- ggplot(df_stats, aes(x=X_NONPAR_ADJ_LRR_MEDIAN, y=MT_ADJ_LRR_MEDIAN, color=SEX)) +
  geom_point(alpha=1/3) +
  scale_x_continuous('X nonPAR median LRR (autosome corrected)') +
  scale_y_continuous('MT median LRR (autosome corrected)') +
  scale_color_discrete(guide = FALSE) +
  theme_bw(base_size = fs)
p8 <- ggplot(df_stats, aes(x=Y_NONPAR_ADJ_LRR_MEDIAN, y=MT_ADJ_LRR_MEDIAN, color=SEX)) +
  geom_point(alpha=1/3) +
  scale_x_continuous('Y nonPAR median LRR (autosome corrected)') +
  scale_y_continuous('MT median LRR (autosome corrected)') +
  scale_color_discrete(guide = FALSE) +
  theme_bw(base_size = fs)
if (model == 'wgs') {
  p9 <- ggplot(df_stats, aes(x=BAF_CORR, y=BAF_CONC)) +
    geom_point(alpha=1/3) +
    scale_x_continuous('Beta-binomial correlation BAF') +
    scale_y_continuous('BAF concordance') +
    theme_bw(base_size = fs)
  p10 <- ggplot(df_stats, aes(x=REL_ESS, y=COV_AUTO)) +
    geom_point(alpha=1/3) +
    scale_x_continuous('Relative coverage variance explained by GC') +
    scale_y_continuous('GC-adjusted coverage auto-correlation') +
    theme_bw(base_size = fs)
  p11 <- ggplot(df_stats, aes(x=COV_SD, y=COV_AUTO)) +
    geom_point(alpha=1/3) +
    scale_x_continuous('Standard deviation coverage') +
    scale_y_continuous('GC-adjusted coverage auto-correlation') +
    theme_bw(base_size = fs)
} else {
  p9 <- ggplot(df_stats, aes(x=BAF_SD, y=BAF_CONC)) +
    geom_point(alpha=1/3) +
    scale_x_continuous('Standard deviation BAF') +
    scale_y_continuous('BAF concordance') +
    theme_bw(base_size = fs)
  p10 <- ggplot(df_stats, aes(x=REL_ESS, y=LRR_AUTO)) +
    geom_point(alpha=1/3) +
    scale_x_continuous('Relative LRR variance explained by GC') +
    scale_y_continuous('GC-adjusted LRR auto-correlation') +
    theme_bw(base_size = fs)
  p11 <- ggplot(df_stats, aes(x=LRR_SD, y=LRR_AUTO)) +
    geom_point(alpha=1/3) +
    scale_x_continuous('Standard deviation LRR') +
    scale_y_continuous('GC-adjusted LRR auto-correlation') +
    theme_bw(base_size = fs)
}

pdf(out_pdf)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
print(p8)
print(p9)
print(p10)
print(p11)
invisible(dev.off())
