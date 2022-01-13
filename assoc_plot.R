#!/usr/bin/env Rscript
###
#  The MIT License
#
#  Copyright (C) 2021-2022 Giulio Genovese
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

assoc_plot_version <- '2022-01-12'

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
if (capabilities()[['cairo']]) options(bitmapType = 'cairo')

parser <- OptionParser('usage: assoc_plot.R [options] --genome <GRCh37|GRCh38>|--cytoband <cytoband.txt.gz> --vcf|--tbx <file>')
parser <- add_option(parser, c('--genome'), type = 'character', help = 'genome assembly (e.g. GRCh38)', metavar = '<assembly>')
parser <- add_option(parser, c('--cytoband'), type = 'character', help = 'cytoband file', metavar = '<cytoband.txt.gz>')
parser <- add_option(parser, c('--vcf'), type = 'character', help = 'input VCF file', metavar = '<file.vcf>')
parser <- add_option(parser, c('--csq'), action = 'store_true', default = FALSE, help = 'whether coding variant should be flagged as red')
parser <- add_option(parser, c('--tbx'), type = 'character', help = 'input REGENIE/PLINK summary statistics', metavar = '<file.gz>')
parser <- add_option(parser, c('--region'), type = 'character', help = 'region to plot', metavar = '<region>')
parser <- add_option(parser, c('--min-a1freq'), type = 'double', help = 'minimum minor allele frequency', metavar = '<floatr>')
parser <- add_option(parser, c('--min-phred'), type = 'integer', help = 'minimum phred score [20]', metavar = '<integer>')
parser <- add_option(parser, c('--cyto-ratio'), type = 'integer', default = 25, help = 'plot to cytoband ratio [25]', metavar = '<integer>')
parser <- add_option(parser, c('--max-height'), type = 'integer', help = 'phred score ceiling', metavar = '<integer>')
parser <- add_option(parser, c('--spacing'), type = 'integer', default = 20, help = 'spacing between chromosomes [10]', metavar = '<integer>')
parser <- add_option(parser, c('--pdf'), type = 'character', help = 'output PDF file', metavar = '<file.pdf>')
parser <- add_option(parser, c('--png'), type = 'character', help = 'output PNG file', metavar = '<file.png>')
parser <- add_option(parser, c('--width'), type = 'double', default = 7.0, help = 'inches width of the output file [7.0]', metavar = '<float>')
parser <- add_option(parser, c('--height'), type = 'double', help = 'inches height of the output file [3.5/7.0]', metavar = '<float>')
parser <- add_option(parser, c('--fontsize'), type = 'integer', default = 12, help = 'font size [12]', metavar = '<integer>')

args <- parse_args(parser, commandArgs(trailingOnly = TRUE), convert_hyphens_to_underscores = TRUE)

write(paste('assoc_plot.R', assoc_plot_version, 'https://github.com/freeseek/mocha'), stderr())

if (is.null(args$genome) && is.null(args$cytoband)) {print_help(parser); stop('either --genome or --cytoband is required')}
if (!is.null(args$genome) && !is.null(args$cytoband)) {print_help(parser); stop('cannot use --genome and --cytoband at the same time')}
if (!is.null(args$genome) && args$genome != 'GRCh37' && args$genome != 'GRCh38') {print_help(parser); stop('--genome accepts only GRCh37 or GRCh38')}
if (is.null(args$vcf) && is.null(args$tbx)) {print_help(parser); stop('either --vcf or --tbx is required')}
if (!is.null(args$vcf) && !is.null(args$tbx)) {print_help(parser); stop('either --vcf or --tbx is required')}
if (is.null(args$vcf) && args$csq)  {print_help(parser); stop('--csq requires --vcf')}
if (is.null(args$pdf) && is.null(args$png)) {print_help(parser); stop('either --pdf or --png is required')}
if (!is.null(args$pdf) && !is.null(args$png)) {print_help(parser); stop('cannot use --pdf and --png at the same time')}
if (!is.null(args$png) && !capabilities('png')) {print_help(parser); stop('unable to start device PNG: no png support in this version of R\nyou need to reinstall R with support for PNG to use the --png option')}

if (is.null(args$min_phred)) {
  if (is.null(args$region)) {
    args$min_phred <- 20
  } else {
    args$min_phred <- 0
  }
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

  chrlen <- chrlen[c(1:22,'X')]
} else if ( args$genome == 'GRCh37' ) {
  chrlen <- setNames(c(249251621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63026520, 48129895, 51305566, 155270560), c(1:22,'X'))
} else if ( args$genome == 'GRCh38' ) {
  chrlen <- setNames(c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895), c(1:22,'X'))
}

if ( !is.null(args$vcf) ) {
  if (args$csq) {
    fmt <- '"%CHROM\\t%POS\\t%AS{0}\\t%AS{1}\\t%Consequence\\n"'
    names <- c('chrom', 'pos', 'as0', 'as1', 'consequence')
    if (is.null(args$region)) {
      cmd <- paste('bcftools +split-vep --output-type u --columns Consequence', args$vcf, '| bcftools query --format', fmt)
    } else {
      cmd <- paste('bcftools +split-vep --output-type u --columns Consequence --regions', args$region, args$vcf, '| bcftools query --format', fmt)
    }
  } else {
    fmt <- '"%CHROM\\t%POS\\t%AS{0}\\t%AS{1}\\n"'
    names <- c('chrom', 'pos', 'as0', 'as1')
    if (is.null(args$region)) {
      cmd <- paste('bcftools query --format', fmt, args$vcf)
    } else {
      cmd <- paste('bcftools query --format', fmt, '--regions', args$region, args$vcf)
    }
  }
} else {
  if ( !is.null(args$region) ) {
    cmd <- paste('tabix --print-header', args$tbx, strsplit(args$region,','))
  } else {
    cmd <- paste('zcat', args$tbx)
  }
  if (!is.null(args$min_a1freq)) {
    filter <- paste0('$a1freq>', args$min_a1freq, ' && $a1freq<', 1-args$min_a1freq, ' && $log10p!="NA" && $log10p>', args$min_phred / 10.0)
  } else {
    filter <- paste0('$log10p!="NA" && $log10p>', args$min_phred / 10.0)
  }
  cmd <- paste(cmd, '| awk \'NR==1 {for (i=1; i<=NF; i++) f[$i] = i; if ("CHROM" in f) chrom=f["CHROM"]; else chrom=f["#CHROM"]; if ("GENPOS" in f) pos=f["GENPOS"]; else pos=f["POS"]; if ("A1FREQ" in f) a1freq=f["A1FREQ"]; else a1freq=f["A1_FREQ"]; if ("LOG10P" in f) log10p=f["LOG10P"]; else log10p=f["LOG10_P"]} NR==1 || NR>1 &&', filter , '{print $chrom"\\t"$pos"\\t"$log10p}\'')
  names <- c('chrom', 'pos', 'log10p')
}

write(paste('Command:', cmd), stderr())
if (packageVersion('data.table') < '1.11.6') {
  df <- setNames(fread(cmd, sep = '\t', header = TRUE, na.strings = '.', colClasses = list(character = c(1)), data.table = FALSE), names)
} else {
  df <- setNames(fread(cmd = cmd, sep = '\t', header = TRUE, na.strings = '.', colClasses = list(character = c(1)), data.table = FALSE), names)
}

if (!is.null(args$vcf)) {
  df$log10p = -(log(2) + pbinom(pmin(df$as0, df$as1), df$as0 + df$as1, .5, log.p = TRUE)) / log(10)
  df <- df[df$log10p > args$min_phred / 10.0,]
}

df$chrom <- as.factor(gsub('^chr', '', gsub('^chrM', 'MT', df$chrom)))
ord <- order(as.numeric(gsub('MT', '26', gsub('Y', '24', gsub('X', '23', levels(df$chrom))))))
df$chrom <- factor(df$chrom, levels(df$chrom)[ord])

if (is.null(args$height)) {
  if (length(unique(df$chrom)) == 1) {
    args$height <- 7.0
  } else {
    args$height <- 3.5
  }
}

if (!is.null(args$max_height)) {
  max_height <- max(args$max_height, -10.0 * log10(5e-8))
} else {
  max_height <- max(10.0 * max(df$log10p), -10.0 * log10(5e-8))
}

if (length(unique(df$chrom)) == 1) {
  p <- ggplot(df, aes(x = pos/1e6, y = log10p)) +
    geom_hline(yintercept = -log10(5e-8), color = 'gray', lty = 'longdash') +
    geom_point(size = 1/2) +
    scale_x_continuous(paste('Chromosome', unique(df$chrom), '(Mbp position)'), expand = c(.01,.01)) +
    scale_y_continuous('-log10(p-value)', expand = c(.01,.01)) +
    theme_bw(base_size = args$fontsize)
} else {
  df$chrompos <- cumsum(c(0,args$spacing * 1e6 + chrlen))[as.numeric(df$chrom)] + df$pos
  p <- ggplot(df, aes(x = chrompos/1e6, y = log10p, color = as.numeric(chrom)%%2 == 1)) +
    geom_hline(yintercept = -log10(5e-8), color = 'gray', lty = 'longdash') +
    geom_point(size = 1/2) +
    scale_x_continuous(NULL, breaks = (cumsum(args$spacing * 1e6 + chrlen) - chrlen/2 - args$spacing * 1e6) / 1e6, labels = names(chrlen), expand = c(.01,.01)) +
    scale_y_continuous('-log10(p-value)', expand = c(.01,.01)) +
    scale_color_manual(guide = FALSE, values = c('FALSE' = 'dodgerblue', 'TRUE' = 'gray')) +
    theme_bw(base_size = args$fontsize)
}

if (!is.null(args$cytoband) && length(unique(df$chrom)) == 1) {
  cyto_height <- (max_height - args$min_phred) / args$cyto_ratio
  p <- p  +
    geom_rect(data = df_cyto[df_cyto$chrom == unique(df$chrom) & df_cyto$gieStain != 'acen',], aes(x = NULL, y = NULL, xmin = chromStart/1e6, xmax = chromEnd/1e6, fill = gieStain, shape = NULL), ymin = args$min_phred/10 - cyto_height/10, ymax = args$min_phred/10, color = 'black', size = 1/4, show.legend = FALSE) +
    geom_polygon(data = df_cen[df_cen$chrom == unique(df$chrom),], aes(x = x/1e6, y = args$min_phred/10 + cyto_height/10 * y, shape = NULL, group = name), color = 'black', fill = 'red', size = 1/8) +
    scale_fill_manual(values = c('gneg' = 'white', 'gpos25' = 'lightgray', 'gpos50' = 'gray50', 'gpos75' = 'darkgray', 'gpos100' = 'black', 'gvar' = 'lightblue', 'stalk' = 'slategrey'))
} else {
  cyto_height <- 0
}

xlim <- FALSE
if ( !is.null(args$region) ) {
  if (grepl(':', args$region) && grepl('-', args$region)) {
    left <- as.numeric(gsub('^.*:', '', gsub('-.*$', '', args$region)))
    right <- as.numeric(gsub('^.*-', '', args$region))
    xlim <- TRUE
  }
}

if (!is.null(args$max_height) && xlim) {
  p <- p + coord_cartesian(xlim = c(left, right)/1e6, ylim = c(args$min_phred - cyto_height, max_height)/10)
} else if (xlim) {
  p <- p + coord_cartesian(xlim = c(left, right)/1e6)
} else if (!is.null(args$max_height)) {
  p <- p + coord_cartesian(ylim = c(args$min_phred - cyto_height, max_height)/10)
}

if (!is.null(args$pdf)) {
  pdf(args$pdf, width = args$width, height = args$height)
} else {
  png(args$png, width = args$width, height = args$height, units = 'in', res = 150)
}
print(p)
invisible(dev.off())
