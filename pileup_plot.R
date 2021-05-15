#!/usr/bin/env Rscript
###
#  The MIT License
#
#  Copyright (C) 2017-2021 Giulio Genovese
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

pileup_plot_version <- '2021-05-14'

library(optparse)
library(ggplot2)
options(bitmapType = 'cairo')

parser <- OptionParser('usage: pileup_plot.R [options] --cytoband <cytoband.txt.gz> --stats <file.tsv> --calls <file.tsv> --pdf <file.pdf>')
parser <- add_option(parser, c('--cytoband'), type = 'character', help = 'cytoband file', metavar = '<cytoband.txt.gz>')
parser <- add_option(parser, c('--stats'), type = 'character', help = 'input MoChA stats file', metavar = '<file.tsv>')
parser <- add_option(parser, c('--calls'), type = 'character', help = 'input MoChA calls file', metavar = '<file.tsv>')
parser <- add_option(parser, c('--loss-cumulative'), action = 'store_true', default = FALSE, help = 'whether to plot cumulative losses')
parser <- add_option(parser, c('--cn-loh-cumulative'), action = 'store_true', default = FALSE, help = 'whether to plot cumulative CN-LOHs')
parser <- add_option(parser, c('--gain-cumulative'), action = 'store_true', default = FALSE, help = 'whether to plot cumulative gains')
parser <- add_option(parser, c('--call-rate-thr'), type = 'double', default = 0.97, help = 'minimum call rate threshold [0.97]', metavar = '<float>')
parser <- add_option(parser, c('--baf-auto-thr'), type = 'double', default = 0.03, help = 'maximum BAF autocorrelation threshold [0.03]', metavar = '<float>')
parser <- add_option(parser, c('--pdf'), type = 'character', help = 'output PDF file', metavar = '<file.pdf>')
parser <- add_option(parser, c('--width'), type = 'integer', default = 7, help = 'inches width of the output file [7]', metavar = '<integer>')
parser <- add_option(parser, c('--height'), type = 'integer', default = 7, help = 'inches height of the output file [7]', metavar = '<integer>')
parser <- add_option(parser, c('--fontsize'), type = 'integer', default = 16, help = 'font size [16]', metavar = '<integer>')
args <- parse_args(parser, commandArgs(trailingOnly = TRUE), convert_hyphens_to_underscores = TRUE)

write(paste('pileup_plot.R', pileup_plot_version, 'https://github.com/freeseek/mocha'), stderr())

if (is.null(args$cytoband)) {print_help(parser); stop('option --cytoband is required')}
if (is.null(args$stats)) {print_help(parser); stop('option --stats is required')}
if (is.null(args$calls)) {print_help(parser); stop('option --calls is required')}
if (is.null(args$pdf)) {print_help(parser); stop('option --pdf is required')}

df_stats <- read.table(args$stats, sep = '\t', header = TRUE)
xcl_smpls <- df_stats$sample_id[df_stats$call_rate < args$call_rate_thr & df_stats$baf_auto > args$baf_auto_thr]
df_calls <- read.table(args$calls, sep = '\t', header = TRUE)
df_calls$chrom <- as.factor(gsub('^chr', '', gsub('^chrM', 'MT', df_calls$chrom)))
ord <- order(as.numeric(gsub('MT', '26', gsub('Y', '24', gsub('X', '23', levels(df_calls$chrom))))))
df_calls$chrom <- factor(df_calls$chrom, levels(df_calls$chrom)[ord])
beg_pos <- names(df_calls)[grepl('^beg_', names(df_calls))]
end_pos <- names(df_calls)[grepl('^end_', names(df_calls))]

df_cyto <- setNames(read.table(args$cytoband, sep = '\t', header = FALSE), c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain'))
df_cyto$chrom <- gsub('chr', '', df_cyto$chrom)
chrlen <- tapply(df_cyto$chromEnd, df_cyto$chrom, max)
idx <- df_cyto$gieStain %in% c('acen', 'gvar', 'stalk')
cen_beg <- tapply(df_cyto$chromEnd[idx], df_cyto$chrom[idx], min)
cen_end <- tapply(df_cyto$chromEnd[idx], df_cyto$chrom[idx], max)
chrs <- unique(df_cyto$chrom)
modified_chrs <- gsub('MT', '26', gsub('Y', '24', gsub('X', '23', chrs)))
ord <- order(suppressWarnings(as.numeric(modified_chrs)))
chrs <- chrs[ord]

df_cen <- rbind(cbind(setNames(df_cyto[df_cyto$gieStain == 'acen' & substr(df_cyto$name, 1, 3) == 'p11', c('chrom', 'name', 'chromStart')], c('chrom', 'name', 'x')), y = -1),
                cbind(setNames(df_cyto[df_cyto$gieStain == 'acen' & substr(df_cyto$name, 1, 3) == 'p11', c('chrom', 'name', 'chromEnd')], c('chrom', 'name', 'x')), y = -1/2),
                cbind(setNames(df_cyto[df_cyto$gieStain == 'acen' & substr(df_cyto$name, 1, 3) == 'p11', c('chrom', 'name', 'chromStart')], c('chrom', 'name', 'x')), y = 0),
                cbind(setNames(df_cyto[df_cyto$gieStain == 'acen' & substr(df_cyto$name, 1, 3) == 'q11', c('chrom', 'name', 'chromEnd')], c('chrom', 'name', 'x')), y = -1),
                cbind(setNames(df_cyto[df_cyto$gieStain == 'acen' & substr(df_cyto$name, 1, 3) == 'q11', c('chrom', 'name', 'chromStart')], c('chrom', 'name', 'x')), y = -1/2),
                cbind(setNames(df_cyto[df_cyto$gieStain == 'acen' & substr(df_cyto$name, 1, 3) == 'q11', c('chrom', 'name', 'chromEnd')], c('chrom', 'name', 'x')), y = 0))

pdf(args$pdf, width = args$width, height = args$height)

for (chr in c(1:22, 'X')) {
  if (chr == 'X') {
    idx <- !(df_calls$sample_id %in% xcl_smpls) & df_calls$chrom == chr & df_calls$computed_gender == 'F' & df_calls$type %in% c('Undetermined', 'CN-LOH', 'Loss', 'Gain')
  } else {
    idx <- !(df_calls$sample_id %in% xcl_smpls) & df_calls$chrom == chr & df_calls$type %in% c('Undetermined', 'CN-LOH', 'Loss', 'Gain')
  }
  if ( sum(idx) == 0 ) next
  df_calls$index[idx] <- rank(df_calls[idx, beg_pos] - df_calls[idx, 'length']/chrlen[chr] + 1e9 * (3 * (df_calls$type[idx]=='Loss') + 2 * (df_calls$type[idx]=='CN-LOH') + (df_calls$type[idx]=='Gain')) , ties.method = 'first') - .5

  p <- ggplot(data = df_calls[idx,], aes_string(x = paste0(beg_pos, '/1e6'), y = 'index', color = 'type')) +
    geom_segment(aes_string(x = paste0(beg_pos, '/1e6'), xend = paste0(end_pos, '/1e6'), y = 'index', yend = 'index')) +
    theme_bw() +
    scale_x_continuous(paste('Chromosome', chr, '(Mbp position)')) +
    scale_y_continuous(NULL, breaks = NULL) +
    scale_color_manual(NULL, values = c('Undetermined' = 'gray', 'CN-LOH' = 'orange', 'Loss' = 'blue', 'Gain' = 'red'),
                             labels = c('Undetermined' = paste0('Undetermined (n=', sum(df_calls$type[idx] == 'Undetermined'), ')'),
                                        'CN-LOH' = paste0('CN-LOH (n=', sum(df_calls$type[idx] == 'CN-LOH'), ')'),
                                        'Loss' = paste0('Loss (n=', sum(df_calls$type[idx] == 'Loss'), ')'),
                                        'Gain' = paste0('Gain (n=', sum(df_calls$type[idx] == 'Gain'), ')'))) +
    theme(legend.position = 'bottom', legend.box = 'horizontal')

  loss_cumulative <- args$loss_cumulative && sum(idx & df_calls$type == 'Loss') > 0
  cn_loh_cumulative <- args$cn_loh_cumulative && sum(idx & df_calls$type == 'CN-LOH') > 0
  gain_cumulative <- args$gain_cumulative && sum(idx & df_calls$type == 'Gain') > 0
  if (loss_cumulative || cn_loh_cumulative || gain_cumulative) {
    if (loss_cumulative) {
      begs <- as.data.frame(table(df_calls[idx & df_calls$type == 'Loss', beg_pos]))
      begs$Var1 <- as.numeric(as.character(begs$Var1))
      ends <- as.data.frame(table(df_calls[idx & df_calls$type == 'Loss', end_pos]))
      ends$Var1 <- as.numeric(as.character(ends$Var1))
      df <- merge(begs, ends, by = 'Var1', all = TRUE)
      df[is.na(df)] <- 0
      df_loss <- data.frame(x = rep(df$Var1, each = 2), y = c(0, head(rep(cumsum(df$Freq.x-df$Freq.y), each = 2), -1)))
      p <- p + geom_line(data = df_loss, aes(x = x/1e6, y = sum(idx) / 20 * (21 + 4 * y/max(y))), color = 'blue')
    }

    if (cn_loh_cumulative) {
      begs <- as.data.frame(table(df_calls[idx & df_calls$type == 'CN-LOH', beg_pos]))
      begs$Var1 <- as.numeric(as.character(begs$Var1))
      ends <- as.data.frame(table(df_calls[idx & df_calls$type == 'CN-LOH', end_pos]))
      ends$Var1 <- as.numeric(as.character(ends$Var1))
      df <- merge(begs, ends, by = 'Var1', all = TRUE)
      df[is.na(df)] <- 0
      df_cnloh <- data.frame(x = rep(df$Var1, each = 2), y = c(0, head(rep(cumsum(df$Freq.x-df$Freq.y), each = 2), -1)))
      p <- p + geom_line(data = df_cnloh, aes(x = x/1e6, y = sum(idx) / 20 * (21 + 4 * y/max(y))), color = 'orange')
    }

    if (gain_cumulative) {
      begs <- as.data.frame(table(df_calls[idx & df_calls$type == 'Gain', beg_pos]))
      begs$Var1 <- as.numeric(as.character(begs$Var1))
      ends <- as.data.frame(table(df_calls[idx & df_calls$type == 'Gain', end_pos]))
      ends$Var1 <- as.numeric(as.character(ends$Var1))
      df <- merge(begs, ends, by = 'Var1', all = TRUE)
      df[is.na(df)] <- 0
      df_cnloh <- data.frame(x = rep(df$Var1, each = 2), y = c(0, head(rep(cumsum(df$Freq.x-df$Freq.y), each = 2), -1)))
      p <- p + geom_line(data = df_cnloh, aes(x = x/1e6, y = sum(idx) / 20 * (21 + 4 * y/max(y))), color = 'red')
    }

    p <- p  +
      geom_rect(data = df_cyto[df_cyto$chrom == chr & df_cyto$gieStain != 'acen',], aes(x = NULL, y = NULL, xmin = chromStart/1e6, xmax = chromEnd/1e6, fill = gieStain, shape = NULL), ymin = sum(idx), ymax = sum(idx) * 21 / 20, color = 'black', size = 1/4, show.legend = FALSE) +
      geom_polygon(data = df_cen[df_cen$chrom == chr,], aes(x = x/1e6, y = sum(idx) - sum(idx)/20 * y, shape = NULL, group = name), color = 'black', fill = 'red', size = 1/8) +
      scale_fill_manual(values = c('gneg' = 'white', 'gpos25' = 'lightgray', 'gpos50' = 'gray50', 'gpos75' = 'darkgray', 'gpos100' = 'black', 'gvar' = 'lightblue', 'stalk' = 'slategrey')) +
      coord_cartesian(xlim = c(0, chrlen[chr] / 1e6), ylim = c(0, sum(idx) * 25 / 20), expand = FALSE)
  } else {
    p <- p  +
      geom_rect(data = df_cyto[df_cyto$chrom == chr & df_cyto$gieStain != 'acen',], aes(x = NULL, y = NULL, xmin = chromStart/1e6, xmax = chromEnd/1e6, fill = gieStain, shape = NULL), ymin = -sum(idx)/20, ymax = 0, color = 'black', size = 1/4, show.legend = FALSE) +
      geom_polygon(data = df_cen[df_cen$chrom == chr,], aes(x = x/1e6, y = sum(idx)/20 * y, shape = NULL, group = name), color = 'black', fill = 'red', size = 1/8) +
      scale_fill_manual(values = c('gneg' = 'white', 'gpos25' = 'lightgray', 'gpos50' = 'gray50', 'gpos75' = 'darkgray', 'gpos100' = 'black', 'gvar' = 'lightblue', 'stalk' = 'slategrey')) +
      coord_cartesian(xlim = c(0, chrlen[chr] / 1e6), ylim = c(-sum(idx)/20, sum(idx)), expand = FALSE)
  }

  print(p)
}

invisible(dev.off())
