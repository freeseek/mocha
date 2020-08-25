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

pileup_plot_version <- '2020-08-25'

library(optparse)
library(ggplot2)
options(bitmapType = 'cairo')

parser <- OptionParser('usage: pileup_plot.R [options] --cytoband <cytoband.txt.gz> --stats <file.tsv> --calls <file.tsv> --pdf <file.pdf>')
parser <- add_option(parser, c('--cytoband'), type = 'character', help = 'cytoband file', metavar = '<cytoband.txt.gz>')
parser <- add_option(parser, c('--stats'), type = 'character', help = 'input MoChA stats file', metavar = '<file.tsv>')
parser <- add_option(parser, c('--calls'), type = 'character', help = 'input MoChA calls file', metavar = '<file.tsv>')
parser <- add_option(parser, c('--pdf'), type = 'character', help = 'output PDF file', metavar = '<file.pdf>')
parser <- add_option(parser, c('--width'), type = 'integer', default = 7, help = 'inches width of the output file [7]', metavar = '<integer>')
parser <- add_option(parser, c('--height'), type = 'integer', default = 7, help = 'inches height of the output file [7]', metavar = '<integer>')
parser <- add_option(parser, c('--fontsize'), type = 'integer', default = 16, help = 'font size [16]', metavar = '<integer>')
args <- parse_args(parser, commandArgs(trailingOnly = TRUE))

write(paste('pileup_plot.R', pileup_plot_version, 'https://github.com/freeseek/mocha'), stderr())

if (is.null(args$cytoband)) {print_help(parser); stop('option --cytoband is required')}
if (is.null(args$stats)) {print_help(parser); stop('option --stats is required')}
if (is.null(args$calls)) {print_help(parser); stop('option --calls is required')}
if (is.null(args$pdf)) {print_help(parser); stop('option --pdf is required')}

df_stats <- read.table(args$stats, sep = '\t', header = TRUE)
xcl_smpls <- df_stats$sample_id[df_stats$call_rate < .97 & df_stats$baf_conc > .51]
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
  df_calls$index[idx] <- rank(df_calls[idx, beg_pos] + 1e9 * (3 * (df_calls$type[idx]=='Loss') + 2 * (df_calls$type[idx]=='CN-LOH') + (df_calls$type[idx]=='Gain')) , ties.method = 'first') - .5
  p <- ggplot(data=df_calls[idx,], aes_string(x=paste0(beg_pos, '/1e6'), y='index', color='type')) +
    geom_segment(aes_string(x=paste0(beg_pos, '/1e6'), xend=paste0(end_pos, '/1e6'), y='index', yend='index')) +
    theme_bw() +
    scale_x_continuous(paste('Chromosome', chr, '(Mbp position)')) +
    scale_y_continuous(NULL, breaks = NULL) +
    scale_color_manual(NULL, values = c('Undetermined' = 'gray', 'CN-LOH' = 'orange', 'Loss' = 'blue', 'Gain' = 'red'),
                             labels = c('Undetermined' = paste0('Undetermined (n=', sum(df_calls$type[idx] == 'Undetermined'), ')'),
                                        'CN-LOH' = paste0('CN-LOH (n=', sum(df_calls$type[idx] == 'CN-LOH'), ')'),
                                        'Loss' = paste0('Loss (n=', sum(df_calls$type[idx] == 'Loss'), ')'),
                                        'Gain' = paste0('Gain (n=', sum(df_calls$type[idx] == 'Gain'), ')'))) +
    theme(legend.position='bottom', legend.box = 'horizontal')
  p <- p  +
    geom_rect(data = df_cyto[df_cyto$chrom == chr & df_cyto$gieStain != 'acen',], aes(x = NULL, y = NULL, xmin = chromStart/1e6, xmax = chromEnd/1e6, fill = gieStain, shape = NULL), ymin = sum(idx), ymax = sum(idx) * 21 / 20, color = 'black', size = 1/4, show.legend = FALSE) +
    geom_polygon(data = df_cen[df_cen$chrom == chr,], aes(x = x/1e6, y = sum(idx) - sum(idx)/20 * y, shape = NULL, group = name), color = 'black', fill = 'red', size = 1/8) +
    scale_fill_manual(values = c('gneg' = 'white', 'gpos25' = 'lightgray', 'gpos50' = 'gray50', 'gpos75' = 'darkgray', 'gpos100' = 'black', 'gvar' = 'lightblue', 'stalk' = 'slategrey')) +
    coord_cartesian(xlim = c(0, chrlen[chr] / 1e6), ylim = c(0, sum(idx) * 21 / 20), expand = FALSE)
  print(p)
}

invisible(dev.off())
