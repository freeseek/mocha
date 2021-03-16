#!/usr/bin/env Rscript
###
#  The MIT License
#
#  Copyright (C) 2021 Giulio Genovese
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

shift_plot_version <- '2021-03-15'

library(optparse)
library(data.table)
library(ggplot2)
options(bitmapType = 'cairo')

parser <- OptionParser('usage: shift_plot.R [options] --vcf <file.vcf> --region <region>')
parser <- add_option(parser, c('--cytoband'), type = 'character', help = 'cytoband file', metavar = '<cytoband.txt.gz>')
parser <- add_option(parser, c('--vcf'), type = 'character', help = 'input VCF file', metavar = '<file.vcf>')
parser <- add_option(parser, c('--region'), type = 'character', help = 'region to plot', metavar = '<region>')
parser <- add_option(parser, c('--min-phred'), type = 'integer', default = 10, help = 'minimum phred score [10]', metavar = '<integer>')
parser <- add_option(parser, c('--cyto-ratio'), type = 'integer', default = 25, help = 'plot to cytoband ratio [25]', metavar = '<integer>')
parser <- add_option(parser, c('--pdf'), type = 'character', help = 'output PDF file', metavar = '<file.pdf>')
parser <- add_option(parser, c('--png'), type = 'character', help = 'output PNG file', metavar = '<file.png>')
parser <- add_option(parser, c('--width'), type = 'integer', default = 7, help = 'inches width of the output file [7]', metavar = '<integer>')
parser <- add_option(parser, c('--height'), type = 'integer', default = 7, help = 'inches height of the output file [7]', metavar = '<integer>')
parser <- add_option(parser, c('--fontsize'), type = 'integer', default = 12, help = 'font size [12]', metavar = '<integer>')

args <- parse_args(parser, commandArgs(trailingOnly = TRUE), convert_hyphens_to_underscores = TRUE)

write(paste('shift_plot.R', shift_plot_version, 'https://github.com/freeseek/mocha'), stderr())

if (is.null(args$vcf)) {print_help(parser); stop('option --vcf is required')}
if (is.null(args$region)) {print_help(parser); stop('option --region is required')}
if (is.null(args$pdf) && is.null(args$png)) {print_help(parser); stop('either --pdf or --png is required')}
if (!is.null(args$pdf) && !is.null(args$png)) {print_help(parser); stop('cannot use --pdf and --png at the same time')}
if (!is.null(args$png) && !capabilities('png')) {print_help(parser); stop('unable to start device PNG: no png support in this version of R\nyou need to reinstall R with support for PNG to use the --png option')}

if (grepl(':', args$region) && grepl('-', args$region)) {
  left <- as.numeric(gsub('^.*:', '', gsub('-.*$', '', args$region)))
  right <- as.numeric(gsub('^.*-', '', args$region))
}

if (!is.null(args$cytoband)) {
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
  df_chrs <- data.frame(chrlen = chrlen[chrs], cen_beg = cen_beg[chrs], cen_end = cen_end[chrs], CHROM = chrs)
}

fmt <- '"%CHROM\\t%POS\\t%AS{0}\\t%AS{1}\\n"'
names <- c('chrom', 'pos', 'as0', 'as1')
cmd <- paste('bcftools query --format', fmt, '--regions', args$region, args$vcf)

write(paste('Command:', cmd), stderr())
if (packageVersion("data.table") < '1.11.6') {
  df <- setNames(fread(cmd, sep = '\t', header = FALSE, na.strings = '.', colClasses = list(character = c(1)), data.table = FALSE), names)
} else {
  df <- setNames(fread(cmd = cmd, sep = '\t', header = FALSE, na.strings = '.', colClasses = list(character = c(1)), data.table = FALSE), names)
}

chrom <- gsub('chr', '', unique(df$chrom))
if (length(chrom) > 1) stop('option --region must provide a region with one chromosome only')

df$pbinom_as = -10 * (pbinom(pmin(df$as0, df$as1), df$as0 + df$as1, .5, log.p = TRUE) + log(2)) / log(10)

idx <- df$pbinom_as > args$min_phred
p <- ggplot(df[idx, ], aes(x=pos/1e6, y=pbinom_as/10)) +
  geom_point(alpha = 1/2, size = 1/2) +
  scale_x_continuous('Mbp position', expand = c(.01,.01)) +
  scale_y_continuous('-log10(p)', expand = c(.01,.01)) +
  theme_bw(base_size = args$fontsize)
if (grepl(':', args$region) && grepl('-', args$region)) p <- p + coord_cartesian(xlim = c(left, right)/1e6)
if (!is.null(args$cytoband)) {
  cyto_height <- (max(df$pbinom_as) - args$min_phred) / args$cyto_ratio
  p <- p  +
    geom_rect(data = df_cyto[df_cyto$chrom == chrom & df_cyto$gieStain != 'acen',], aes(x = NULL, y = NULL, xmin = chromStart/1e6, xmax = chromEnd/1e6, fill = gieStain, shape = NULL), ymin = args$min_phred/10 - cyto_height/10, ymax = args$min_phred/10, color = 'black', size = 1/4, show.legend = FALSE) +
    geom_polygon(data = df_cen[df_cen$chrom == chrom,], aes(x = x/1e6, y = args$min_phred/10 + cyto_height/10 * y, shape = NULL, group = name), color = 'black', fill = 'red', size = 1/8) +
    scale_fill_manual(values = c('gneg' = 'white', 'gpos25' = 'lightgray', 'gpos50' = 'gray50', 'gpos75' = 'darkgray', 'gpos100' = 'black', 'gvar' = 'lightblue', 'stalk' = 'slategrey'))
}

if (!is.null(args$pdf)) {
  pdf(args$pdf, width = args$width, height = args$height)
} else {
  png(args$png, width = args$width, height = args$height, units = 'in', res = 150)
}
print(p)
invisible(dev.off())
