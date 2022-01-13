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
#  FITNESS FOR A PARTICULAR PURposE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
###

mocha_plot_version <- '2022-01-12'

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
if (capabilities()[['cairo']]) options(bitmapType = 'cairo')

parser <- OptionParser('usage: mocha_plot.R [options] --genome <GRCh37|GRCh38>|--cytoband <cytoband.txt.gz> --vcf <file.vcf> --samples <list>')
parser <- add_option(parser, c('--genome'), type = 'character', help = 'genome assembly (e.g. GRCh38)', metavar = '<assembly>')
parser <- add_option(parser, c('--cytoband'), type = 'character', help = 'cytoband file', metavar = '<cytoband.txt.gz>')
parser <- add_option(parser, c('--wgs'), action = 'store_true', default = FALSE, help = 'whether the input VCF file contains WGS data')
parser <- add_option(parser, c('--mocha'), action = 'store_true', default = FALSE, help = 'whether the input VCF file contains Ldev/Bdev data')
parser <- add_option(parser, c('--no-adjust'), action = 'store_true', default = FALSE, help = 'for array data whether BAF and LRR should not be adjusted')
parser <- add_option(parser, c('--stats'), type = 'character', help = 'input MoChA stats file', metavar = '<file.tsv>')
parser <- add_option(parser, c('--vcf'), type = 'character', help = 'input VCF file', metavar = '<file.vcf>')
parser <- add_option(parser, c('--exclude'), type = 'character', help = 'regions to exclude listed in a file', metavar = '<file.bed>')
parser <- add_option(parser, c('--pdf'), type = 'character', help = 'output PDF file', metavar = '<file.pdf>')
parser <- add_option(parser, c('--png'), type = 'character', help = 'output PNG file', metavar = '<file.png>')
parser <- add_option(parser, c('--width'), type = 'double', default = 7.0, help = 'inches width of the output file [7.0]', metavar = '<float>')
parser <- add_option(parser, c('--height'), type = 'double', default = 7.0, help = 'inches height of the output file [7.0]', metavar = '<float>')
parser <- add_option(parser, c('--samples'), type = 'character', help = 'comma-separated list of samples to plot', metavar = '<list>')
parser <- add_option(parser, c('--regions'), type = 'character', help = 'comma-separated list of regions to plot [all]', metavar = '<list>')
parser <- add_option(parser, c('--fontsize'), type = 'integer', default = 12, help = 'font size [12]', metavar = '<integer>')
parser <- add_option(parser, c('--roll'), type = 'integer', default = 20, help = 'width of the rolling window [20]', metavar = '<integer>')
parser <- add_option(parser, c('--min-depth'), type = 'integer', default = 10, help = 'minimum depth coverage (WGS data only) [10]', metavar = '<integer>')
parser <- add_option(parser, c('--clump'), type = 'integer', default = 10, help = 'width of the clumping window (WGS data only) [10]', metavar = '<integer>')
args <- parse_args(parser, commandArgs(trailingOnly = TRUE), convert_hyphens_to_underscores = TRUE)

write(paste('mocha_plot.R', mocha_plot_version, 'https://github.com/freeseek/mocha'), stderr())

if (is.null(args$vcf)) {print_help(parser); stop('option --vcf is required')}
if (is.null(args$samples)) {print_help(parser); stop('option --samples is required')}
if (is.null(args$genome) && is.null(args$cytoband)) {print_help(parser); stop('either --genome or --cytoband is required')}
if (!is.null(args$genome) && !is.null(args$cytoband)) {print_help(parser); stop('cannot use --genome and --cytoband at the same time')}
if (!is.null(args$genome) && args$genome != 'GRCh37' && args$genome != 'GRCh38') {print_help(parser); stop('--genome accepts only GRCh37 or GRCh38')}
if (is.null(args$pdf) && is.null(args$png)) {print_help(parser); stop('either --pdf or --png is required')}
if (!is.null(args$pdf) && !is.null(args$png)) {print_help(parser); stop('cannot use --pdf and --png at the same time')}
regions <- unlist(strsplit(args$regions, ','))
if (!is.null(args$png) && length(regions) > 1) {print_help(parser); stop('cannot print a png file with more than one image')}
if (!is.null(args$png) && !capabilities('png')) {print_help(parser); stop('unable to start device PNG: no png support in this version of R\nyou need to reinstall R with support for PNG to use the --png option')}

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
} else if (!is.null(args$genome)) {
  if ( args$genome == 'GRCh37' ) {
    chrlen <- c(249251621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63026520, 48129895, 51305566, 155270560, 59373566)
    cen_beg <- c(121535434, 92326171, 90504854, 49660117, 46405641, 58830166, 58054331, 43838887, 47367679, 39254935, 51644205, 34856694, 0, 0, 0, 35335801, 22263006, 15460898, 24681782, 26369569, 0, 0, 58632012, 10104553)
    cen_end <- c(142535434, 95326171, 93504854, 52660117, 49405641, 61830166, 61054331, 46838887, 65367679, 42254935, 54644205, 37856694, 19000000, 19000000, 20000000, 46335801, 25263006, 18460898, 27681782, 29369569, 14288129, 16000000, 61632012, 13104553)
  } else if ( args$genome == 'GRCh38' ) {
    chrlen <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
    cen_beg <- c(122026459, 92188145, 90772458, 49712061, 46485900, 58553888, 58169653, 44033744, 43389635, 39686682, 51078348, 34769407, 0, 0, 0, 36311158, 22813679, 15460899, 24498980, 26436232, 0, 0, 58605579, 10316944)
    cen_end <- c(143184587, 94090557, 93655574, 51743951, 50059807, 59829934, 61528020, 45877265, 60518558, 41593521, 54425074, 37185252, 18051248, 18173523, 19725254, 46280682, 26616164, 20861206, 27190874, 30038348, 12915808, 15054318, 62412542, 10544039)
  }
  chrs <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y')
  names(chrlen) <- chrs
  names(cen_beg) <- chrs
  names(cen_end) <- chrs
  df_chrs <- data.frame(chrlen = chrlen[chrs], cen_beg = cen_beg[chrs], cen_end = cen_end[chrs], CHROM = chrs)
}

# load main table from VCF file
fmt <- '"[%CHROM\\t%POS\\t%REF\\t%ALT\\t%SAMPLE\\t%GT'
names <- c('chrom', 'pos', 'ref', 'alt', 'sample_id', 'gt')
if (!args$wgs) {
  if (args$mocha && !args$no_adjust)
  {
    fmt <- paste0(fmt, paste0('\\t%INFO/ADJ_COEFF{', 0:8, '}', collapse = ''))
    names <- c(names, c(outer(c('BAF0', 'BAF1', 'LRR0'), c('AA_', 'AB_', 'BB_'), FUN=function(x, y) paste0(y, x))))
  }
  fmt <- paste0(fmt, '\\t%INFO/GC\\t%INFO/ALLELE_A\\t%INFO/ALLELE_B\\t%BAF\\t%LRR')
  names <- c(names, c('gc', 'allele_a', 'allele_b', 'BAF', 'LRR'))
} else {
  fmt <- paste0(fmt, '\\t%AD{0}\\t%AD{1}')
  names <- c(names, c('ad0', 'ad1'))
}
if (args$mocha) {
  fmt <- paste0(fmt, '\\t%Ldev\\t%Bdev')
  names <- c(names, 'ldev', 'bdev')
}
fmt <- paste0(fmt, '\\n]"')
cmd <- paste('bcftools query --format', fmt, args$vcf, '--samples', args$samples)

contigs <- unlist(lapply(regions[regions!='all'], function(x) unlist(strsplit(x, ':'))[1]))
chroms <- gsub('^chr', '', gsub('^chrM', 'MT', contigs))
begs <- as.numeric(unlist(lapply(regions[regions != 'all'], function(x) unlist(strsplit(unlist(strsplit(x, ':'))[2], '-'))[1])))
ends <- as.numeric(unlist(lapply(regions[regions != 'all'], function(x) unlist(strsplit(unlist(strsplit(x, ':'))[2], '-'))[2])))
lefts <- round(pmax(1.5 * begs - .5 * ends, 0, na.rm = TRUE))
rights <- round(pmin(1.5 * ends - .5 * begs, chrlen[chroms], na.rm = TRUE))
if (!('all' %in% regions)) cmd <- paste(cmd, '--regions', paste0(contigs, ':', lefts, '-', rights, collapse = ','))

if (!is.null(args$exclude)) cmd <- paste0(cmd, ' --targets-file ^', args$exclude)

write(paste('Command:', cmd), stderr())
if (packageVersion('data.table') < '1.11.6') {
  df <- setNames(fread(cmd, sep = '\t', header = FALSE, na.strings = '.', colClasses = list(character = c(1,3:6)), data.table = FALSE), names)
} else {
  df <- setNames(fread(cmd = cmd, sep = '\t', header = FALSE, na.strings = '.', colClasses = list(character = c(1,3:6)), data.table = FALSE), names)
}

if (args$wgs) {
  nonref <- grepl(',', df$alt)
  if (any(nonref)) df <- df[!nonref,]
  df$allele_a <- 0
  df$allele_b <- 1
}

allele_0 <- lapply(df$gt, function(x) suppressWarnings(as.numeric(substr(x,1,1))))
phased <- lapply(df$gt, function(x) substr(x,2,2)) == '|'
allele_1 <- lapply(df$gt, function(x) suppressWarnings(as.numeric(substr(x,3,3))))
df$gts <- 'NC'
df$gts[allele_0 == df$allele_a & allele_1 == df$allele_a] <- 'AA'
df$gts[allele_0 == df$allele_a & allele_1 == df$allele_b |
       allele_0 == df$allele_b & allele_1 == df$allele_a] <- 'AB'
df$gts[allele_0 == df$allele_b & allele_1 == df$allele_b] <- 'BB'
df$phase <- 0
df$phase[phased & allele_0 == df$allele_a & allele_1 == df$allele_b] <- 1
df$phase[phased & allele_0 == df$allele_b & allele_1 == df$allele_a] <- -1
if (!args$wgs && args$mocha && !args$no_adjust) {
  for (gt in c('AA', 'AB', 'BB')) {
    idx <- df$gts == gt
    df$BAF[idx] <- df$BAF[idx] - df[idx, paste0(gt, '_BAF1')] * df$LRR[idx] - df[idx, paste0(gt, '_BAF0')]
    df$LRR[idx] <- df$LRR[idx] - df[idx, paste0(gt, '_LRR0')]
  }
}
if (!args$wgs && !is.null(args$stats)) {
  df_stats <- read.table(args$stats, sep = '\t', header = TRUE)
  lrr_gc_order <- sum(grepl('^lrr_gc_[0-9]', names(df_stats))) - 1
  df <- merge(df, df_stats[, c('sample_id', paste0('lrr_gc_', 0:lrr_gc_order))])
  for (i in 0:lrr_gc_order) {
    df$LRR <- df$LRR - as.numeric(df$gc)^i * df[, paste0('lrr_gc_', i)]
  }
}
df$sample_id <- factor(df$sample_id, levels = strsplit(args$samples, ',')[[1]])

# fills in variables of interest
df$chrom <- as.factor(gsub('^chr', '', gsub('^chrM', 'MT', df$chrom)))
ord <- order(as.numeric(gsub('MT', '26', gsub('Y', '24', gsub('X', '23', levels(df$chrom))))))
df$chrom <- factor(df$chrom, levels(df$chrom)[ord])
if (!args$wgs) {
  df$eLRR <- exp(df$LRR)
  df$pBAF <- NaN
  idx <- df$phase == 1
  df$pBAF[idx] <- df$BAF[idx]
  idx <- df$phase == -1
  df$pBAF[idx] <- 1 - df$BAF[idx]
  cov_var <- 'eLRR'
  plot_vars <- c('eLRR', 'BAF', 'pBAF')
  df_horiz <- data.frame(variable = factor(c('eLRR', 'pBAF'), levels = plot_vars), value = c(1.0, 0.5))
} else {
  df$logDP <- log(df$ad0 + df$ad1)
  # subset to variants that are heterozygous
  df <- df[df$gts == 'AB' & df$logDP > log(args$min_depth),]
  df$BAF <- NaN
  df$BAF <- df$ad1 / (df$ad0 + df$ad1)
  df$hd0 <- NaN
  df$hd1 <- NaN
  df$hd0[df$phase == 1] <- df$ad0[df$phase == 1]
  df$hd1[df$phase == 1] <- df$ad1[df$phase == 1]
  df$hd0[df$phase == -1] <- df$ad1[df$phase == -1]
  df$hd1[df$phase == -1] <- df$ad0[df$phase == -1]
  if (args$clump == 1)
  {
    df$pBAF <- NaN
    df$pBAF[df$phase == 1] <- df$BAF[df$phase == 1]
    df$pBAF[df$phase == -1] <- 1 - df$BAF[df$phase == -1]
  }
  cov_var <- 'logDP'
  plot_vars <- c('logDP', 'pBAF')
  df_horiz <- data.frame(variable = factor('pBAF', levels = plot_vars), value = 0.5)
}
if (!is.null(args$cytoband)) {
  df_cyto$variable <- factor('pBAF', levels = plot_vars)
  df_cen$variable <- factor('pBAF', levels = plot_vars)
}
if (args$mocha) {
  df$color <- as.factor(1e3 * df$ldev + df$bdev)
  df$color[df$ldev == 0 & df$bdev == 0 | is.na(df$color)] <- NA
} else {
  df$color <- NA
}

if (!is.null(args$pdf)) {
  pdf(args$pdf, width = args$width, height = args$height)
} else {
  png(args$png, width = args$width, height = args$height, units = 'in', res = 150)
}

if ('all' %in% regions) {
  write(paste0('Plotting region: all (', sum(!is.na(df$BAF) & df$chrom %in% chrs), ' heterozygous sites)'), stderr())
  p <- ggplot(df[!is.na(df$BAF) & df$chrom %in% chrs,], aes(x = pos/1e6, y = BAF, color = color))
  if (!is.null(args$cytoband) | !is.null(args$genome)) {
    p <- p + geom_vline(data = df_chrs, aes(xintercept = chrlen/1e6), color = 'black', size = 1, alpha = 1/2)
  }
  p <- p + geom_rect(data = df_chrs, mapping = aes(x = NULL, y = NULL, xmin = cen_beg/1e6, xmax = cen_end/1e6), color = 'transparent', fill = 'gray', ymin = 0, ymax = 1, alpha = 1/2) +
    scale_x_continuous('Mbp position') +
    scale_y_continuous('B Allele Frequency (BAF)', breaks = NULL) +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    scale_color_discrete(guide = FALSE) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    facet_grid(chrom ~ sample_id) +
    theme_bw(base_size = args$fontsize) +
    theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5))
  if (!args$mocha || all(is.na(df$color))) {
    if (!is.null(args$png)) {
      p <- p + geom_point(alpha = 1/4, size = 1/16, show.legend = FALSE, color = 'gray50')
    } else {
      p <- p + geom_bin2d(binwidth = c(1, .05), alpha = 1/4, show.legend = FALSE, color = 'gray50')
    }
  } else {
    if (!is.null(args$png)) {
      p <- p + geom_point(alpha = 1/4, size = 1/16, show.legend = FALSE)
    } else {
      p <- p + geom_bin2d(binwidth = c(1, .05), alpha = 1/4, show.legend = FALSE)
    }
  }
  print(p)
  regions <- regions[regions!='all']
}

if (length(regions)>0) {
  for (i in 1:length(regions)) {
    write(paste('Subsetting region:', regions[i]), stderr())
    idx <- df$chrom == chroms[i] & df$pos >= lefts[i] & df$pos <= rights[i]
    if (!args$wgs) {
      df_melt <- rbind(cbind(setNames(df[idx, c('pos', 'color', 'gts', 'sample_id', 'eLRR')], c('pos', 'color', 'gts', 'sample_id', 'value')), variable = 'eLRR'),
                       cbind(setNames(df[idx, c('pos', 'color', 'gts', 'sample_id', 'BAF')], c('pos', 'color', 'gts', 'sample_id', 'value')), variable = 'BAF'),
                       cbind(setNames(df[idx, c('pos', 'color', 'gts', 'sample_id', 'pBAF')], c('pos', 'color', 'gts', 'sample_id', 'value')), variable = 'pBAF'))
      df_melt$variable <- factor(df_melt$variable, c('eLRR', 'BAF', 'pBAF'))
      df_melt$value[df_melt$variable == 'eLRR' & df_melt$value > 2 |
                    df_melt$variable == 'BAF' & abs(df_melt$value - 0.5) > 0.6 |
                    df_melt$variable == 'pBAF' & abs(df_melt$value - 0.5) > 0.35] <- NA
    } else {
      if (args$clump == 1) {
        df_melt <- rbind(cbind(setNames(df[idx, c('pos', 'color', 'gts', 'sample_id', 'logDP')], c('pos', 'color', 'gts', 'sample_id', 'value')), variable = 'logDP'),
                         cbind(setNames(df[idx, c('pos', 'color', 'gts', 'sample_id', 'BAF')], c('pos', 'color', 'gts', 'sample_id', 'value')), variable = 'BAF'),
                         cbind(setNames(df[idx, c('pos', 'color', 'gts', 'sample_id', 'pBAF')], c('pos', 'color', 'gts', 'sample_id', 'value')), variable = 'pBAF'))
        df_melt$variable <- factor(df_melt$variable, c('logDP', 'BAF', 'pBAF'))
        median_logdp <- median(df_melt$value[df_melt$variable == 'logDP'], na.rm = TRUE)
        df_melt$value[df_melt$variable == 'logDP' & df_melt$value > 2 * median_logdp] <- NA
      } else {
        l <- list()
        for (sm in unique(df$sample_id)) {
          sm_idx <- idx & df$sample_id == sm
          index <- round(.5 + (1:sum(sm_idx)) / sum(sm_idx) * ceiling(sum(sm_idx) / args$clump))
          l[[sm]] <- data.frame(pos = unname(tapply(df$pos[sm_idx], index, median, na.rm = TRUE)),
                                color = unname(tapply(as.numeric(df$color[sm_idx]), index, max, -1, na.rm = TRUE)),
                                sample_id = sm,
                                logDP = unname(tapply(df$logDP[sm_idx], index, mean, na.rm = TRUE)),
                                hd0 = unname(tapply(df$hd0[sm_idx], index, sum, na.rm = TRUE)),
                                hd1 = unname(tapply(df$hd1[sm_idx], index, sum, na.rm = TRUE)))
        }
        tmp <- rbindlist(l)
        tmp$color <- as.factor(tmp$color)
        tmp$color[tmp$color == -1] <- NaN
        tmp$pBAF = tmp$hd1 / (tmp$hd0 + tmp$hd1)
        df_melt <- reshape2::melt(tmp[, c('pos', 'color', 'sample_id', 'logDP', 'pBAF')], id.vars = c('pos', 'color', 'sample_id'))
      }
    }
    df_melt$sample_id <- factor(df_melt$sample_id, levels = strsplit(args$samples, ',')[[1]])

    idx <- df_melt$variable == cov_var & !is.na(df_melt$value)
    if (sum(idx) >= args$roll) df_melt$smooth[idx] <- filter(df_melt$value[idx], rep(1 / args$roll, args$roll))
    idx <- df_melt$variable == 'pBAF' & !is.na(df_melt$value)
    if (sum(idx) >= args$roll) df_melt$smooth[idx] <- filter(df_melt$value[idx], rep(1 / args$roll, args$roll))

    write(paste('Plotting region:', regions[i]), stderr())
    if ('gts' %in% names(df_melt)) {
      p <- ggplot(df_melt[!is.na(df_melt$value),], aes(x = pos/1e6, y = value, color = color, shape = gts)) +
        scale_shape_manual(guide = FALSE, values = c('AA' = 3, 'AB' = 8, 'BB' = 4, 'NC' = 1))
    } else {
      p <- ggplot(df_melt[!is.na(df_melt$value),], aes(x = pos/1e6, y = value, color = color))
    }
    if (!is.na(begs[i]) || !is.na(ends[i])) {
      p <- p + geom_vline(xintercept = c(begs[i]/1e6, ends[i]/1e6), alpha = 1/2, linetype = 'dashed', color = 'gray')
    }
    p <- p + scale_x_continuous(paste('Chromosome', chroms[i], '(Mbp position)')) +
      scale_y_continuous(NULL) +
      scale_color_discrete(guide = FALSE) +
      theme_bw(base_size = args$fontsize) +
      theme(legend.position = 'bottom', legend.box = 'horizontal', strip.background = element_blank()) +
      facet_grid(variable ~ sample_id, scales = 'free_y') +
      coord_cartesian(xlim = c(lefts[i], rights[i])/1e6, expand = FALSE)
    if (all(is.na(df_melt$color))) {
      p <- p + geom_point(size = 1/4, color = 'gray50')
    } else {
      p <- p + geom_point(size = 1/4)
    }
    p <- p + geom_hline(data = df_horiz, aes(yintercept = value), linetype = 'dashed', color = 'black', size = 1/4)
    if ('smooth' %in% names(df_melt)) {
      p <- p + geom_line(data = df_melt[!is.na(df_melt$smooth),], aes(y = smooth, shape = NULL), color = 'blue', size = 1/4)
    }
    if (!is.null(args$cytoband)) {
      bottom <- floor(min(0.35, df_melt$value[df_melt$variable == 'pBAF'], na.rm = TRUE) * 20) / 20
      p <- p  +
        geom_rect(data = df_cyto[df_cyto$chrom == chroms[i] & df_cyto$gieStain != 'acen',], aes(x = NULL, y = NULL, xmin = chromStart/1e6, xmax = chromEnd/1e6, fill = gieStain, shape = NULL), ymin = bottom -0.05, ymax = bottom, color = 'black', size = 1/4, show.legend = FALSE) +
        geom_polygon(data = df_cen[df_cen$chrom == chroms[i],], aes(x = x/1e6, y = bottom + 0.05 * y, shape = NULL, group = name), color = 'black', fill = 'red', size = 1/8) +
        scale_fill_manual(values = c('gneg' = 'white', 'gpos25' = 'lightgray', 'gpos50' = 'gray50', 'gpos75' = 'darkgray', 'gpos100' = 'black', 'gvar' = 'lightblue', 'stalk' = 'slategrey'))
    }
    print(p)
  }
}

invisible(dev.off())
