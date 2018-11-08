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

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))

args <- commandArgs(trailingOnly = TRUE)
model <- args[1]
rules <- args[2]
out_pdf <- args[3]
vcf <- args[4]
if (length(args)>4) {
  maxlen <- as.numeric(args[5])
} else {
  maxlen <- 0
}

title <- sub('.mocha$', '', sub('.vcf$', '', sub('.gz$', '', sub('.bcf$', '', basename(vcf)))))

if ( model == 'array' ) {
  fmt <- '"%CHROM\t%POS[\t%GT\t%BAF\t%LRR\t%Ldev\t%Bdev\t%Bdev_Phase]\\n"'
  names <- c('CHROM', 'POS', 'GT', 'BAF', 'LRR', 'LDEV', 'BDEV', 'BDEV_Phase')
} else if ( model == 'wgs' ) {
  fmt <- '"%CHROM\t%POS[\t%GT\t%AD{0}\t%AD{1}\t%Ldev\t%Bdev\t%Bdev_Phase]\\n"'
  names <- c('CHROM', 'POS', 'GT', 'AD0', 'AD1', 'LDEV', 'BDEV', 'BDEV_Phase')
} else {
  stop("Model parameter needs to be \"array\" or \"wgs\"")
}

chrs <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y')
if ( rules == 'GRCh37' ) {
  chrlen <- c(249251621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63026520, 48129895, 51305566, 155270560, 59373566)
} else if ( rules == 'GRCh38' ) {
  chrlen <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
} else {
  stop("Rules parameter needs to be \"GRCh37\" or \"GRCh38\"")
}
names(chrlen) <- chrs

# load main table from VCF file
cmd <- paste('bcftools query --format', fmt, vcf)
write(paste('Command:', cmd), stderr())
df <- setNames(fread(cmd, sep = '\t', header = FALSE, na.strings = '.', data.table = FALSE), names)

# fills in variables of interest
df$CHROM <- as.factor(gsub('^chr', '', gsub('^chrM', 'MT', df$CHROM)))
ord <- order(as.numeric(gsub('MT', '26', gsub('Y', '24', gsub('X', '23', levels(df$CHROM))))))
df$CHROM <- factor(df$CHROM, levels(df$CHROM)[ord])
df$PHASE = (df$GT == '0|1') - (df$GT == '1|0')
if ( model == 'wgs' ) {
  df$BAF <- ( df$AD1 + 0.5 ) / ( df$AD0 + df$AD1 + 1.0 )
  df$LRR <- log2( df$AD0 + df$AD1 ) / log2( median( df$AD0 + df$AD1 ) )
}
df$BAF[unname(sapply(df$GT, (function(x) substr(x,1,1)==substr(x,3,3))))] <- NaN
df$GT <- NULL
df$pBAF <- (df$BAF - 0.5) * df$PHASE + 0.5
df$pBAF[df$PHASE==0] <- NaN

pdf(out_pdf)

fs <- 8
ggplot(df[!is.na(df$BAF),], aes(x = POS/1e6, y = BAF, color = as.factor(BDEV))) +
  geom_bin2d(binwidth = c(1, .05), alpha = 1/4, show.legend = FALSE) +
  scale_x_continuous('Mbp position', limits = c(0, chrlen['1']/1e6), expand = c(0, 0)) +
  scale_y_continuous('BAF', limits = c(0, 1), breaks = c()) +
  scale_color_discrete(guide = FALSE) +
  ggtitle(title) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_grid(CHROM ~ .) +
  theme_bw(base_size = fs) +
  theme(plot.title = element_text(hjust = 0.5))

if (any(df$LDEV!=0, na.rm=TRUE)) {
  boundaries <- which(diff(c(0,df$LDEV,0))!=0)
  mocha_start <- boundaries[which(df$LDEV[boundaries]!=0)]
  mocha_end <- boundaries[which(df$LDEV[boundaries]!=0)+1]-1

  for (i in 1:length(mocha_start)) {
    chrom <- df$CHROM[mocha_start[i]]
    start_pos <- df$POS[mocha_start[i]]
    end_pos <- df$POS[mocha_end[i]]
    bdev <- df$BDEV[mocha_start[i]]
    lrr_median <- median(df$LRR[mocha_start[i]:mocha_end[i]])
    if ( ( end_pos - start_pos < 1e5 ) || ( end_pos - start_pos < maxlen ) || ( chrom == 'X') ) next # skip small events and chromosome X events
    write(paste0('Plotting call: ', chrom, ':', start_pos, '-', end_pos, ' ', bdev, ' ', lrr_median), stderr())

    p <- list()

    title <- paste0('chr', chrom, ':', format(start_pos, scientific = FALSE, big.mark=","), '-', format(end_pos, scientific = FALSE, big.mark=","))
    use <- !is.na(df$BAF) & df$CHROM==chrom & df$POS>=start_pos & df$POS<=end_pos
    p[['hist']] <- ggplot(df[use, ], aes(x = BAF)) +
      geom_histogram(binwidth = .01, color = 'black', fill = 'white') +
      geom_vline(xintercept = c(1/3, 1/2, 2/3), alpha=1/2, linetype = 'dashed', color = 'gray') +
      scale_x_continuous('', limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous('count',  expand = c(0, 0)) +
      ggtitle(title) +
      theme_bw(base_size = fs) +
      theme(plot.title = element_text(hjust = 0.5))

    left <- max(4 * start_pos - 3 * end_pos, 0)
    right <- min(4 * end_pos - 3 * start_pos, chrlen[chrom])
    idx <- df$CHROM == chrom & df$POS >= left & df$POS <= right
    for (s in c(1e5, 1e6)) {
      df_scaled <- data.frame(POS = unname(tapply(df$POS[idx], round(df$POS[idx]/s), median, na.rm=TRUE)),
                              pBAF = unname(tapply(df$pBAF[idx], round(df$POS[idx]/s), median, na.rm=TRUE)),
                              LRR = unname(tapply(df$LRR[idx], round(df$POS[idx]/s), median, na.rm=TRUE)))
      df_scaled$CNV <- df_scaled$POS >= start_pos & df_scaled$POS <= end_pos

      title <- paste0('Scaling of ', format(s, scientific = FALSE, big.mark=","), 'bp')
      use <- !is.na(df_scaled$LRR)
      p[[paste0('lrr', as.character(s))]] <- ggplot(df_scaled[use,], aes(x = POS/1e6, y = LRR, color = CNV)) +
        geom_vline(xintercept = c(start_pos/1e6, end_pos/1e6), alpha=1/2, linetype = 'dashed', color = 'gray') +
        geom_point(size = 2, shape = 'x') +
        scale_x_continuous('', limits = c(left, right)/1e6, expand = c(0, 0)) +
        scale_y_continuous('LRR') +
        ggtitle(title) +
        scale_color_manual(guide = FALSE, values = c('TRUE' = 'red', 'FALSE' = 'gray')) +
        theme_bw(base_size = fs) +
        theme(plot.title = element_text(hjust = 0.5))

      use <- !is.na(df_scaled$pBAF)
      p[[paste0('baf', as.character(s))]] <- ggplot(df_scaled[use,], aes(x = POS/1e6, y = pBAF, color = CNV)) +
        geom_hline(yintercept = c(1/3, 1/2, 2/3), alpha=1/2, linetype = 'dashed', color = 'gray') +
        geom_vline(xintercept = c(start_pos/1e6, end_pos/1e6), alpha=1/2, linetype = 'dashed', color = 'gray') +
        geom_point(size = 2, shape = 'x') +
        scale_x_continuous('', limits = c(left, right)/1e6, expand = c(0, 0)) +
        scale_y_continuous('pBAF') +
        scale_color_manual(guide = FALSE, values = c('TRUE' = 'red', 'FALSE' = 'black')) +
        theme_bw(base_size = fs)

      use <- !is.na(df_scaled$pBAF) & df_scaled$CNV
      p[[paste0('hist', as.character(s))]] <- ggplot(df_scaled[use, ], aes(x = pBAF)) +
        geom_histogram(binwidth = .01, color = 'black', fill = 'white') +
        geom_vline(xintercept = c(1/3, 1/2, 2/3), alpha=1/2, linetype = 'dashed', color = 'black') +
        scale_x_continuous('', limits = c(0, 1), expand = c(0, 0)) +
        scale_y_continuous('count',  expand = c(0, 0)) +
        theme_bw(base_size = fs)
    }
    do.call(grid.arrange, c(p, nrow = length(p), ncol = 1))
  }
}

invisible(dev.off())
