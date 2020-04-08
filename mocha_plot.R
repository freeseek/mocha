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
library(data.table)
library(ggplot2)

parser <- OptionParser('usage: mocha_plot.R [options] --rules <GRCh37|GRCh38>|--cytoband <cytoband.txt.gz> --vcf <file.vcf> --samples <list>')
parser <- add_option(parser, c('--rules'), type = 'character', help = 'genome assembly (e.g. GRCh38)', metavar = '<assembly>')
parser <- add_option(parser, c('--cytoband'), type = 'character', help = 'cytoband file', metavar = '<cytoband.txt.gz>')
parser <- add_option(parser, c('--wgs'), action = 'store_true', default = FALSE, help = 'whether the input VCF file contains WGS data')
parser <- add_option(parser, c('--mocha'), action = 'store_true', default = FALSE, help = 'whether the input VCF file contains Ldev/Bdev data')
parser <- add_option(parser, c('--no-adjust'), action = 'store_true', default = FALSE, help = 'for array data whether BAF and LRR should not be adjusted')
parser <- add_option(parser, c('--vcf'), type = 'character', help = 'input VCF file', metavar = '<file.vcf>')
parser <- add_option(parser, c('--exclude'), type = 'character', help = 'regions to exclude listed in a file', metavar = '<file.bed>')
parser <- add_option(parser, c('--pdf'), type = 'character', help = 'output PDF file', metavar = '<file.pdf>')
parser <- add_option(parser, c('--png'), type = 'character', help = 'output PNG file', metavar = '<file.png>')
parser <- add_option(parser, c('--width'), type = 'integer', default = 7, help = 'inches width of the output file [7]', metavar = '<integer>')
parser <- add_option(parser, c('--height'), type = 'integer', default = 7, help = 'inches height of the output file [7]', metavar = '<integer>')
parser <- add_option(parser, c('--samples'), type = 'character', help = 'comma-separated list of samples to plot', metavar = '<list>')
parser <- add_option(parser, c('--regions'), type = 'character', help = 'comma-separated list of regions to plot [all]', metavar = '<list>')
parser <- add_option(parser, c('--fontsize'), type = 'integer', default = 12, help = 'font size [12]', metavar = '<integer>')
parser <- add_option(parser, c('--roll'), type = 'integer', default = 20, help = 'width of the rolling window [20]', metavar = '<integer>')
parser <- add_option(parser, c('--min-depth'), type = 'integer', default = 10, help = 'minimum depth coverage (WGS data only) [10]', metavar = '<integer>')
parser <- add_option(parser, c('--clump'), type = 'integer', default = 10, help = 'width of the clumping window (WGS data only) [10]', metavar = '<integer>')
args <- parse_args(parser, commandArgs(trailingOnly = TRUE), convert_hyphens_to_underscores = TRUE)

if (is.null(args$vcf)) {print_help(parser); stop('option --vcf is required')}
if (is.null(args$samples)) {print_help(parser); stop('option --samples is required')}
if (!is.null(args$rules) && !is.null(args$cytoband)) {print_help(parser); stop('either --rules or --cytoband is required')}
if (!is.null(args$rules) && !is.null(args$cytoband)) {print_help(parser); stop('cannot use --rules and --cytoband at the same time')}
if (!is.null(args$rules) && args$rules != 'GRCh37' && args$rules != 'GRCh38') {print_help(parser); stop('--rules accepts only GRCh37 or GRCh38')}
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

  df_cen <- rbind(melt(setNames(df_cyto[df_cyto$gieStain == 'acen' & substr(df_cyto$name, 1, 3) == 'p11', c('chrom', 'name', 'chromStart', 'chromEnd', 'chromStart')], c('chrom', 'name', -1, -.5, 0)), id = c('chrom', 'name')),
                  melt(setNames(df_cyto[df_cyto$gieStain == 'acen' & substr(df_cyto$name, 1, 3) == 'q11', c('chrom', 'name', 'chromEnd', 'chromStart', 'chromEnd')], c('chrom', 'name', -1, -.5, 0)), id = c('chrom', 'name')))
  df_cen$y <- as.numeric(levels(df_cen$variable)[df_cen$variable])
  df_chrs <- data.frame(chrlen = chrlen[chrs], cen_beg = cen_beg[chrs], cen_end = cen_end[chrs], CHROM = chrs)
} else if (!is.null(args$rules)) {
  if ( args$rules == 'GRCh37' ) {
    chrlen <- c(249251621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63026520, 48129895, 51305566, 155270560, 59373566)
    cen_beg <- c(121535434, 92326171, 90504854, 49660117, 46405641, 58830166, 58054331, 43838887, 47367679, 39254935, 51644205, 34856694, 0, 0, 0, 35335801, 22263006, 15460898, 24681782, 26369569, 0, 0, 58632012, 10104553)
    cen_end <- c(142535434, 95326171, 93504854, 52660117, 49405641, 61830166, 61054331, 46838887, 65367679, 42254935, 54644205, 37856694, 19000000, 19000000, 20000000, 46335801, 25263006, 18460898, 27681782, 29369569, 14288129, 16000000, 61632012, 13104553)
  } else if ( args$rules == 'GRCh38' ) {
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
names <- c('CHROM', 'POS', 'REF', 'ALT', 'SAMPLE', 'GT')
if (!args$wgs) {
  if (!args$no_adjust)
  {
    fmt <- paste0(fmt, '\\t%INFO/ADJUST_BAF_LRR{0}\\t%INFO/ADJUST_BAF_LRR{1}\\t%INFO/ADJUST_BAF_LRR{2}\\t%INFO/ADJUST_BAF_LRR{3}')
    names <- c(names, c('BAF_AB', 'LRR_AA', 'LRR_AB', 'LRR_BB'))
  }
  fmt <- paste0(fmt, '\\t%INFO/ALLELE_A\\t%INFO/ALLELE_B\\t%BAF\\t%LRR')
  names <- c(names, c('ALLELE_A', 'ALLELE_B', 'BAF', 'LRR'))
} else {
  fmt <- paste0(fmt, '\\t%AD{0}\\t%AD{1}')
  names <- c(names, c('AD0', 'AD1'))
}
if (args$mocha) {
  fmt <- paste0(fmt, '\\t%Ldev\\t%Bdev')
  names <- c(names, 'LDEV', 'BDEV')
}
fmt <- paste0(fmt, '\\n]"')
cmd <- paste('bcftools query --format', fmt, args$vcf, '--samples', args$samples)

contigs <- unlist(lapply(regions[regions!='all'], function(x) unlist(strsplit(x, ':'))[1]))
chroms <- gsub('^chr', '',gsub('^chrM', 'MT', contigs))
begs <- as.numeric(unlist(lapply(regions[regions != 'all'], function(x) unlist(strsplit(unlist(strsplit(x, ':'))[2], '-'))[1])))
ends <- as.numeric(unlist(lapply(regions[regions != 'all'], function(x) unlist(strsplit(unlist(strsplit(x, ':'))[2], '-'))[2])))
lefts <- round(pmax(1.5 * begs - .5 * ends, 0))
if (!is.null(args$cytoband) | !is.null(args$rules)) {
  rights <- round(pmin(1.5 * ends - .5 * begs, chrlen[chroms]))
} else {
  rights <- round(1.5 * ends - .5 * begs)
}
if (!('all' %in% regions)) cmd <- paste(cmd, '--regions', paste0(contigs, ':', lefts, '-', rights, collapse = ','))

if (!is.null(args$exclude)) cmd <- paste0(cmd, ' --targets-file ^', args$exclude)

write(paste('Command:', cmd), stderr())
if (packageVersion("data.table") < '1.11.6') {
  df <- setNames(fread(cmd, sep = '\t', header = FALSE, na.strings = '.', colClasses = list(character = c(1,3:6)), data.table = FALSE), names)
} else {
  df <- setNames(fread(cmd = cmd, sep = '\t', header = FALSE, na.strings = '.', colClasses = list(character = c(1,3:6)), data.table = FALSE), names)
}

allele_0 <- lapply(df$GT, function(x) suppressWarnings(as.numeric(substr(x,1,1))))
phased <- lapply(df$GT, function(x) substr(x,2,2)) == '|'
allele_1 <- lapply(df$GT, function(x) suppressWarnings(as.numeric(substr(x,3,3))))
if (args$wgs) {
  nonref <- grepl(',', df$ALT)
  if (any(nonref)) df <- df[!nonref,]
  df$ALLELE_A <- 0
  df$ALLELE_B <- 1
}
df$GTS <- 'NC'
df$GTS[allele_0 == df$ALLELE_A & allele_1 == df$ALLELE_A] <- 'AA'
df$GTS[allele_0 == df$ALLELE_A & allele_1 == df$ALLELE_B | 
       allele_0 == df$ALLELE_B & allele_1 == df$ALLELE_A] <- 'AB'
df$GTS[allele_0 == df$ALLELE_B & allele_1 == df$ALLELE_B] <- 'BB'
df$PHASE <- 0
df$PHASE[phased & allele_0 == df$ALLELE_A & allele_1 == df$ALLELE_B] <- 1
df$PHASE[phased & allele_0 == df$ALLELE_B & allele_1 == df$ALLELE_A] <- -1
if ( !args$wgs && !args$no_adjust )
{
  idx <- !is.na(df$BAF) & (df$GTS == 'NC' & (df$BAF > df$BAF_AB) & (df$BAF - df$BAF_AB < 1) | df$GTS == 'AB')
  df$BAF[idx] <- df$BAF[idx] - df$BAF_AB[idx]

  idx <- !is.na(df$LRR) & df$GTS == 'AB'
  df$LRR[idx] <- df$LRR[idx] - df$LRR_AB[idx]
  
  idx <- !is.na(df$LRR) & df$GTS == 'AA'
  df$LRR[idx] <- df$LRR[idx] - df$LRR_AA[idx]

  idx <- !is.na(df$LRR) & df$GTS == 'BB'
  df$LRR[idx] <- df$LRR[idx] - df$LRR_BB[idx]
}

# fills in variables of interest
df$CHROM <- as.factor(gsub('^chr', '', gsub('^chrM', 'MT', df$CHROM)))
ord <- order(as.numeric(gsub('MT', '26', gsub('Y', '24', gsub('X', '23', levels(df$CHROM))))))
df$CHROM <- factor(df$CHROM, levels(df$CHROM)[ord])
if (!args$wgs) {
  df$eLRR <- exp(df$LRR)
  df$pBAF <- NaN
  idx <- df$PHASE == 1
  df$pBAF[idx] <- df$BAF[idx]
  idx <- df$PHASE == -1
  df$pBAF[idx] <- 1 - df$BAF[idx]
  cov_var <- 'eLRR'
  plot_vars <- c('eLRR', 'BAF', 'pBAF')
  df_horiz <- data.frame(variable = factor(c('eLRR', 'pBAF'), levels = plot_vars), value = c(1.0, 0.5))
} else {
  df$DP <- df$AD0 + df$AD1
  # subset to variants that are heterozygous
  df <- df[df$GTS == 'AB' & df$DP > args$min_depth,]
  df$BAF <- NaN
  df$BAF <- df$AD1/ df$DP
  df$HD0 <- NaN
  df$HD1 <- NaN
  df$HD0[df$PHASE == 1] <- df$AD0[df$PHASE == 1]
  df$HD1[df$PHASE == 1] <- df$AD1[df$PHASE == 1]
  df$HD0[df$PHASE == -1] <- df$AD1[df$PHASE == -1]
  df$HD1[df$PHASE == -1] <- df$AD0[df$PHASE == -1]
  if (args$clump == 1)
  {
    df$pBAF <- NaN
    df$pBAF[df$PHASE == 1] <- df$BAF[df$PHASE == 1]
    df$pBAF[df$PHASE == -1] <- 1 - df$BAF[df$PHASE == -1]
  }
  cov_var <- 'DP'
  plot_vars <- c('DP', 'pBAF')
  df_horiz <- data.frame(variable = factor('pBAF', levels = plot_vars), value = 0.5)
}
if (!is.null(args$cytoband)) {
  df_cyto$variable <- factor('pBAF', levels = plot_vars)
  df_cen$variable <- factor('pBAF', levels = plot_vars)
}
if (args$mocha) {
  df$COLOR <- as.factor(1e3 * df$LDEV + df$BDEV)
  df$COLOR[df$LDEV == 0 & df$BDEV == 0 | is.na(df$COLOR)] <- NA
} else {
  df$COLOR <- NA
}

if (!is.null(args$pdf)) {
  pdf(args$pdf, width = args$width, height = args$height)
} else {
  png(args$png, width = args$width, height = args$height, units = 'in', res = 150)
}

if ('all' %in% regions) {
  write(paste0('Plotting region: all (', sum(!is.na(df$BAF) & df$CHROM %in% chrs), ' heterozygous sites)'), stderr())
  p <- ggplot(df[!is.na(df$BAF) & df$CHROM %in% chrs,], aes(x = POS/1e6, y = BAF, color = COLOR))
  if (!is.null(args$cytoband) | !is.null(args$rules)) {
    p <- p + geom_vline(data = df_chrs, aes(xintercept = chrlen/1e6), color = 'black', size = 1, alpha = 1/2)
  }
  p <- p + geom_rect(data = df_chrs, mapping = aes(x = NULL, y = NULL, xmin = cen_beg/1e6, xmax = cen_end/1e6), color = 'transparent', fill = 'gray', ymin = 0, ymax = 1, alpha = 1/2) +
    scale_x_continuous('Mbp position') +
    scale_y_continuous('B Allele Frequency (BAF)', breaks = NULL) +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    scale_color_discrete(guide = FALSE) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    facet_grid(CHROM ~ SAMPLE) +
    theme_bw(base_size = args$fontsize) +
    theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5))
  if (!args$mocha || all(is.na(df$COLOR))) {
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
    idx <- df$CHROM == chroms[i] & df$POS >= lefts[i] & df$POS <= rights[i]
    if (!args$wgs) {
      df_melt <- melt(df[idx, c('POS', 'COLOR', 'GTS', 'SAMPLE', 'eLRR', 'BAF', 'pBAF')], id.vars = c('POS', 'COLOR', 'GTS', 'SAMPLE'))
      df_melt$value[df_melt$variable == 'eLRR' & df_melt$value > 2 | df_melt$variable == 'pBAF' & abs(df_melt$value - 0.5) > 0.35] <- NA
    } else {
      if (args$clump == 1) {
        df_melt <- melt(df[idx, c('POS', 'COLOR', 'GTS', 'SAMPLE', 'DP', 'BAF', 'pBAF')], id.vars = c('POS', 'COLOR', 'GTS', 'SAMPLE'))
        medianDP <- median(df_melt$value[df_melt$variable == 'DP'], na.rm = TRUE)
        df_melt$value[df_melt$variable == 'DP' & df_melt$value > 2 * medianDP] <- NA
      } else {
        l <- list()
        for (sm in unique(df$SAMPLE)) {
          sm_idx <- idx & df$SAMPLE == sm
          index <- round(.5 + (1:sum(sm_idx)) / sum(sm_idx) * ceiling(sum(sm_idx) / args$clump))
          l[[sm]] <- data.frame(POS = unname(tapply(df$POS[sm_idx], index, median, na.rm = TRUE)),
                                COLOR = unname(tapply(as.numeric(df$COLOR[sm_idx]), index, max, -1, na.rm = TRUE)),
                                SAMPLE = sm,
                                DP = unname(tapply(df$DP[sm_idx], index, mean, na.rm = TRUE)),
                                HD0 = unname(tapply(df$HD0[sm_idx], index, sum, na.rm = TRUE)),
                                HD1 = unname(tapply(df$HD1[sm_idx], index, sum, na.rm = TRUE)))
        }
        tmp <- rbindlist(l)
        tmp$COLOR <- as.factor(tmp$COLOR)
        tmp$COLOR[tmp$COLOR == -1] <- NaN
        tmp$pBAF = tmp$HD1 / (tmp$HD0 + tmp$HD1)
        df_melt <- melt(tmp[, c('POS', 'COLOR', 'SAMPLE', 'DP', 'pBAF')], id.vars = c('POS', 'COLOR', 'SAMPLE'))
      }
    }

    idx <- df_melt$variable == cov_var & !is.na(df_melt$value)
    if (sum(idx) >= args$roll) df_melt$smooth[idx] <- filter(df_melt$value[idx], rep(1 / args$roll, args$roll))
    idx <- df_melt$variable == 'pBAF' & !is.na(df_melt$value)
    if (sum(idx) >= args$roll) df_melt$smooth[idx] <- filter(df_melt$value[idx], rep(1 / args$roll, args$roll))

    write(paste('Plotting region:', regions[i]), stderr())
    if ('GTS' %in% names(df_melt)) {
      p <- ggplot(df_melt[!is.na(df_melt$value),], aes(x = POS/1e6, y = value, color = COLOR, shape = GTS)) +
        scale_shape_manual(guide = FALSE, values = c('AA' = 3, 'AB' = 8, 'BB' = 4, 'NC' = 1))
    } else {
      p <- ggplot(df_melt[!is.na(df_melt$value),], aes(x = POS/1e6, y = value, color = COLOR))
    }
    p <- p + geom_vline(xintercept = c(begs[i]/1e6, ends[i]/1e6), alpha = 1/2, linetype = 'dashed', color = 'gray') +
      scale_x_continuous(paste('Chromosome', chroms[i], '(Mbp position)')) +
      scale_y_continuous(NULL) +
      scale_color_discrete(guide = FALSE) +
      theme_bw(base_size = args$fontsize) +
      theme(legend.position = 'bottom', legend.box = 'horizontal', strip.background = element_blank()) +
      facet_grid(variable ~ SAMPLE, scales = 'free_y') +
      coord_cartesian(xlim = c(lefts[i], rights[i])/1e6, expand = FALSE)
    if (all(is.na(df_melt$COLOR))) {
      p <- p + geom_point(size = 1/4, color = 'gray50')
    } else {
      p <- p + geom_point(size = 1/4)
    }
    p <- p + geom_hline(data = df_horiz, aes(yintercept = value), linetype = 'dashed', color = 'black', size = 1/4) +
      geom_line(data = df_melt[!is.na(df_melt$smooth),], aes(y = smooth, shape = NULL), color = 'blue', size = 1/4)
    if (!is.null(args$cytoband)) {
      bottom <- floor(min(0.35, df_melt$value[df_melt$variable == 'pBAF'], na.rm = TRUE) * 20) / 20
      p <- p  +
        geom_rect(data = df_cyto[df_cyto$chrom == chroms[i] & df_cyto$gieStain != 'acen',], aes(x = NULL, y = NULL, xmin = chromStart/1e6, xmax = chromEnd/1e6, fill = gieStain, shape = NULL), ymin = bottom -0.05, ymax = bottom, color = 'black', size = 1/4, show.legend = FALSE) +
        geom_polygon(data = df_cen[df_cen$chrom == chroms[i],], aes(x = value/1e6, y = bottom + 0.05 * y, shape = NULL, group = name), color = 'black', fill = 'red', size = 1/8) +
        scale_fill_manual(values = c('gneg' = 'white', 'gpos25' = 'lightgray', 'gpos50' = 'gray50', 'gpos75' = 'darkgray', 'gpos100' = 'black', 'gvar' = 'lightblue', 'stalk' = 'slategrey'))
    }
    print(p)
  }
}

invisible(dev.off())
