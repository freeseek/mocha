#!/usr/bin/env Rscript
###
#  The MIT License
#
#  Copyright (C) 2017-2019 Giulio Genovese
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

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))

parser <- ArgumentParser(description='Plot MoChA calls from VCF file (version 2019-02-15)')
parser$add_argument('--rules', metavar = '<assembly>', type = 'character', required = TRUE, help = 'genome assembly (e.g. GRCh38)')
parser$add_argument('--vcf', metavar = '<file.vcf>', type = 'character', required = TRUE, help = 'input VCF file')
parser$add_argument('--pdf', metavar = '<file.pdf>', type = 'character', help = 'output PDF file')
parser$add_argument('--png', metavar = '<file.png>', type = 'character', help = 'output PNG file')
parser$add_argument('--samples', metavar = '<id>', type = 'character', required = TRUE, help = 'samples to extract')
parser$add_argument('--regions', metavar = '<file.pdf>', type = 'character', default = 'all', help = 'comma-separated list of regions to plot [all]')
parser$add_argument('--fontsize', metavar = '<integer>', type = 'integer', default = 12, help = 'font size [12]')
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
# args <- parser$parse_args(c('--rules', 'GRCh37', '--vcf', '/tmp/csmd1_dlg2/cirm.core.gsa.mocha.bcf', '--pdf', '/tmp/cazzo.pdf', '--samples', 'CW70255C,RP-1468_CW70255', '--regions', '8:2792875-4852328'))
if (is.null(args$pdf) && is.null(args$png)) stop('either --pdf or --png is required')
if (!is.null(args$pdf) && !is.null(args$png)) stop('cannot use --pdf and --png at the same time')
regions <- unlist(strsplit(args$regions, ','))

chrs <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y')
if ( args$rules == 'GRCh37' ) {
  chrlen <- c(249251621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63026520, 48129895, 51305566, 155270560, 59373566)
  names(chrlen) <- chrs
} else if ( args$rules == 'GRCh38' ) {
  chrlen <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
  names(chrlen) <- paste0('chr', chrs)
} else {
  stop("Rules parameter needs to be \"GRCh37\" or \"GRCh38\"")
}

# load main table from VCF file
fmt <- '"[%CHROM\\t%POS\\t%REF\\t%ALT\\t%SAMPLE\\t%GT\\t%BAF\\t%LRR\\t%Ldev\\t%Bdev\\t%Bdev_Phase\\n]"'
names <- c('CHROM', 'POS', 'REF', 'ALT', 'SAMPLE', 'GT', 'BAF', 'LRR', 'LDEV', 'BDEV', 'BDEV_Phase')
cmd <- paste('bcftools query --format', fmt, args$vcf)
if (!is.null(args$samples)) cmd <- paste(cmd, '--samples', args$samples)

chroms <- unlist(lapply(regions[regions!='all'], function(x) unlist(strsplit(x, ':'))[1]))
begs <- as.numeric(unlist(lapply(regions[regions!='all'], function(x) unlist(strsplit(unlist(strsplit(x, ':'))[2], '-'))[1])))
ends <- as.numeric(unlist(lapply(regions[regions!='all'], function(x) unlist(strsplit(unlist(strsplit(x, ':'))[2], '-'))[2])))
lefts <- round(pmax(1.5 * begs - .5 * ends, 0))
rights <- round(pmin(1.5 * ends - .5 * begs, chrlen[chroms]))
if (!('all' %in% regions)) cmd <- paste(cmd, '--regions', paste0(chroms,':',lefts,'-',rights, collapse=','))
chroms <- as.factor(gsub('^chr', '', gsub('^chrM', 'MT', chroms)))

write(paste('Command:', cmd), stderr())
df <- setNames(fread(cmd, sep = '\t', header = FALSE, na.strings = '.', data.table = FALSE), names)

# fix non-reference SNPs
nonref <- grepl(',', df$ALT)
df$GT[nonref & df$GT=='1/1'] <- '0/0'
df$GT[nonref & df$GT=='1/2'] <- '0/1'
df$GT[nonref & df$GT=='2/2'] <- '1/1'

# fills in variables of interest
df$CHROM <- as.factor(gsub('^chr', '', gsub('^chrM', 'MT', df$CHROM)))
ord <- order(as.numeric(gsub('MT', '26', gsub('Y', '24', gsub('X', '23', levels(df$CHROM))))))
df$CHROM <- factor(df$CHROM, levels(df$CHROM)[ord])
df$UNPHASED_GT <- df$GT
df$UNPHASED_GT[df$GT=='0/0' | df$GT=='0|0'] <- '0/0'
df$UNPHASED_GT[df$GT=='1/1' | df$GT=='1|1'] <- '1/1'
df$UNPHASED_GT[df$GT=='0/1' | df$GT=='1/0' | df$GT=='0|1' | df$GT=='1|0'] <- '0/1'
df$PHASE = (df$GT == '0|1') - (df$GT == '1|0')
df$pBAF <- NaN
df$pBAF[df$GT=='0|1' | df$GT=='1|0'] <- (df$BAF[df$GT=='0|1' | df$GT=='1|0'] - 0.5) * df$PHASE[df$GT=='0|1' | df$GT=='1|0'] + 0.5
df$pBAF[df$GT=='1/1' | df$GT=='1|1'] <- df$BAF[df$GT=='1/1' | df$GT=='1|1']
df$pBAF[df$GT=='0/0' | df$GT=='0|0'] <- 1 - df$BAF[df$GT=='0/0' | df$GT=='0|0']

df$LDEV <- as.factor(df$LDEV)
df$LDEV[df$LDEV==0] <- NaN
df$BDEV <- as.factor(df$BDEV)
df$BDEV[df$BDEV==0] <- NaN

if (!is.null(args$pdf)) {
  pdf(args$pdf)
} else {
  png(args$png, width = 7, height = 7, units = 'in', res = 150)
}

if ('all' %in% regions) {
  write('Plotting region: all', stderr())
  p <- ggplot(df[!is.na(df$BAF),], aes(x = POS/1e6, y = BAF, color = LDEV)) +
    scale_x_continuous('Mbp position', expand = c(0, 0)) +
    scale_y_continuous('BAF', breaks = c()) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_color_discrete(guide = FALSE) +
    ggtitle(args$samples) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    facet_grid(CHROM ~ .) +
    theme_bw(base_size = args$fontsize) +
    theme(strip.background = element_rect(color=NA, fill=NA), plot.title = element_text(hjust = 0.5))
  if (all(is.na(df$LDEV))) {
    p <- p + geom_bin2d(binwidth = c(1, .05), alpha = 1/4, show.legend = FALSE, color = 'gray50')
  } else {
    p <- p + geom_bin2d(binwidth = c(1, .05), alpha = 1/4, show.legend = FALSE)
  }
  print(p)
  regions <- regions[regions!='all']
}

if (length(regions)>0) {
  df$CHROM <- as.character(df$CHROM)
  for (i in 1:length(regions)) {
    write(paste('Plotting region:', regions[i]), stderr())

    p_grid <- list()

    title <- paste0('chr', chroms[i], ':', format(begs[i], scientific = FALSE, big.mark=","), '-', format(ends[i], scientific = FALSE, big.mark=","))

    idx <- df$CHROM == chroms[i] & df$POS >= lefts[i] & df$POS <= rights[i]
    
    p <- ggplot(df[idx & !is.na(df$LRR),], aes(x = POS/1e6, y = exp(LRR), color = LDEV, shape = UNPHASED_GT)) +
      geom_vline(xintercept = c(begs[i]/1e6, ends[i]/1e6), alpha=1/2, linetype = 'dashed', color = 'gray') +
      scale_x_continuous('', limits = c(lefts[i], rights[i])/1e6, expand = c(0, 0)) +
      scale_y_continuous(expression(e^LRR), expand = c(0, 0)) +
      coord_cartesian(ylim = c(0, 2)) +
      ggtitle(title) +
      scale_color_discrete(guide = FALSE) +
      scale_shape_manual(guide = FALSE, values = c('0/0' = 3, '0/1' = 8, '1/1' = 4, './.' = 1)) +
      theme_bw(base_size = args$fontsize) +
      theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5)) +
      facet_grid(. ~ SAMPLE)
    if (all(is.na(df$LDEV[idx]))) {
      p_grid[['lrr']] <- p + geom_point(size = 1, color = 'gray50')
    } else {
      p_grid[['lrr']] <- p + geom_point(size = 1)
    }
  
    p <- ggplot(df[idx & !is.na(df$BAF),], aes(x = POS/1e6, y = BAF, color = LDEV, shape = UNPHASED_GT)) +
      geom_hline(yintercept = c(1/3, 1/2, 2/3), alpha=1/2, linetype = 'dashed', color = 'gray') +
      geom_vline(xintercept = c(begs[i]/1e6, ends[i]/1e6), alpha=1/2, linetype = 'dashed', color = 'gray') +
      scale_x_continuous('', limits = c(lefts[i], rights[i])/1e6, expand = c(0, 0)) +
      scale_y_continuous('BAF', expand = c(0, 0)) +
      coord_cartesian(ylim = c(0, 1)) +
      scale_color_discrete(guide = FALSE) +
      scale_shape_manual(guide = FALSE, values = c('0/0' = 3, '0/1' = 8, '1/1' = 4, './.' = 1)) +
      theme_bw(base_size = args$fontsize) +
      theme(strip.background = element_blank(), strip.text.x = element_blank()) +
      facet_grid(. ~ SAMPLE)
    if (all(is.na(df$LDEV[idx]))) {
      p_grid[['baf']] <- p + geom_point(size = 1, color = 'gray50')
    } else {
      p_grid[['baf']] <- p + geom_point(size = 1)
    }
    
    p <- ggplot(df[idx & !is.na(df$pBAF),], aes(x = POS/1e6, y = pBAF, color = LDEV, shape = UNPHASED_GT)) +
      geom_hline(yintercept = c(1/3, 1/2, 2/3), alpha=1/2, linetype = 'dashed', color = 'gray') +
      geom_vline(xintercept = c(begs[i]/1e6, ends[i]/1e6), alpha=1/2, linetype = 'dashed', color = 'gray') +
      scale_x_continuous('', limits = c(lefts[i], rights[i])/1e6, expand = c(0, 0)) +
      scale_y_continuous('phased BAF', expand = c(0, 0)) +
      coord_cartesian(ylim = c(.25, 1)) +
      scale_color_discrete(guide = FALSE) +
      scale_shape_manual('', breaks = c('0/0', '0/1', '1/1', './.'), values = c('0/0' = 3, '0/1' = 8, '1/1' = 4, './.' = 1), labels = c('0/0' = 'Homozygous\nreference', '0/1' = 'Heterozygous', '1/1' = 'Homozygous\nalternate', './.' = 'Missing')) +
      theme_bw(base_size = args$fontsize) +
      theme(strip.background = element_blank(), strip.text.x = element_blank(), legend.position = 'bottom', legend.box = 'horizontal') +
      facet_grid(. ~ SAMPLE)
    if (all(is.na(df$LDEV[idx]))) {
      p_grid[['pbaf']] <- p + geom_point(size = 1, color = 'gray50')
    } else {
      p_grid[['pbaf']] <- p + geom_point(size = 1)
    }
    
    do.call(grid.arrange, c(p_grid, nrow = 3, ncol = 1))
  }
}

invisible(dev.off())
