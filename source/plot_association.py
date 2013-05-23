#!/usr/bin/env python
#
# $File: plot_association.py $
# $LastChangedDate: 2012-06-05 12:31:19 -0500 (Tue, 05 Jun 2012) $
# $Rev: 1179 $
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://varianttools.sourceforge.net for details.
#
# Copyright (C) 2012 Zongxiao He, Bo Peng and Gao Wang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import argparse, sys, os, re
try:
    # python 2 has pickle and cPickle
    import cPickle as pickle
except:
    # python 3 has pickle
    import pickle
from collections import OrderedDict
from .utils import env, runCommand, mkdir_p, downloadFile
from .rtester import Str4R

def plotAssociation(args):
    env.logger.info("Reading from standard input ...")
    data = []
    data_size = 0
    for x in sys.stdin.readlines():
        try:
            data_size += len(x.encode("utf-8"))
        except:
            data_size += len(x)
        data.append(x.rstrip().split())
    env.logger.info("Processing {} of input data ...".format(size(len(data) + data_size)))
    p = Plot(args, data)
    pinput = eval(args.method.upper() + '_FOO') + eval('p.{}'.format(args.method if args.method == 'qq' else 'manhattan'))() + eval(args.method.upper() + '_MAIN')
    env.logger.info("Generating graph(s) ...")
    cmd = "R --slave --no-save --no-restore"
    out = runCommand(cmd, pinput)
    env.logger.info("Complete!")
    #print(pinput)
    return

def qq(args):
	plot(args, 'qq')
	return

def manhattan(args):
	plot(args, 'manhattan')
	return

def manhattan_plain(args):
	# have to adjust font for this plot
	args.font_size = args.font_size / 2.5
	plot(args, 'manhattan')
	return

def size(bytes):
    system = [
    (1024 ** 5, 'P'),
    (1024 ** 4, 'T'),
    (1024 ** 3, 'G'),
    (1024 ** 2, 'M'),
    (1024 ** 1, 'K'),
    (1024 ** 0, 'B'),
    ]

    for factor, suffix in system:
        if bytes >= factor:
            break
    amount = int(bytes/factor)
    if isinstance(suffix, tuple):
        singular, multiple = suffix
        if amount == 1:
            suffix = singular
        else:
            suffix = multiple
    return str(amount) + suffix


class PlotAssociationOpt:
    def __init__(self, master_parser):
        self.master_parser = master_parser
        subparsers = self.master_parser.add_subparsers()
        # subparser 1
        parserQQ = subparsers.add_parser('qq', help='QQ plot via ggplot2')
        self.qqArguments(parserQQ)
        self.commonArguments(parserQQ)
        parserQQ.set_defaults(func=qq)
        # subparser 2
        parserMan = subparsers.add_parser('manhattan', help='Manhattan plot via ggplot2')
        self.manArguments(parserMan)
        self.commonArguments(parserMan)
        parserMan.set_defaults(func=manhattan)
        # subparser 3
        parserManPlain = subparsers.add_parser('manhattan_plain',
                                               help='Manhattan plot implementation not using ggplot2')
        self.manArguments(parserManPlain)
        self.commonArguments(parserManPlain)
        parserManPlain.set_defaults(func=manhattan_plain)

    def get(self):
        return self.master_parser

    def qqArguments(self, parser):
        parser.add_argument('--shape',
                metavar='INTEGER',
                            type=int,
                default=1,
                help='''Choose a shape theme
                (integer 1 to 16) for dots on QQ plot.
                Default set to 1.''')
        parser.add_argument('--fixed_shape',
                            action='store_true',
                help='''Use the same dot-shape theme for all plots''')
        parser.add_argument('--no_slope',
                            action='store_true',
                help='''Do not plot the diagonal line''')

    def manArguments(self, parser):
        parser.add_argument('--chrom',
                metavar = 'CHROMOSOME',
                nargs = '+',
                default=list(map(str, range(1,23))) + ['X','Y','Un'],
                help='''Specify the particular chromosome(s) to display. Can be
                one or multiple in this list: "{}". Slicing syntax "?:?" is 
                supported. For example "1:22" is equivalent to displaying 
                all autosomes; "1:Y" is equivalent to displaying 
                all mapped chromosomes. Default set to all including unmapped 
                chromosomes.'''.format(' '.join(list(map(str, range(1,23))) + ['X','Y','Un', '?:?'])))
        parser.add_argument('--chrom_prefix',
                metavar = 'PREFIX',
                type = str,
                default = 'chr',
                help='''Prefix chromosome ID with a string.
                Default is set to "chr" (X-axis will be displayed
                as "chr1", "chr2", etc). Use "None" for no prefix.
                ''')
        parser.add_argument('--gene_map',
                metavar = 'FILE',
                type = str,
                help='''If the plot units are genes and the program fails to map certain genes to 
                chromosomes, you can fix it by providing a text file of genomic coordinate 
                information of these genes. Each gene in the file is a line of 3 columns
                specifying "GENENAME CHROM MIDPOINT_POSITION", e.g., "IKBKB 8 42128820".
                ''')

    def commonArguments(self, parser):
        parser.add_argument("--method", default = sys.argv[2] if len(sys.argv) > 2 else '', help=argparse.SUPPRESS)
        settings = parser.add_argument_group('graph properties')
        settings.add_argument('-t', '--title',
                            type=str,
                default='',
                            help='''Title of plot.''')
        settings.add_argument('--color',
                            type=str,
                choices=CTHEME,
                            help='''Choose a color theme from the list above to apply
                to the plot. (via the 'RColorBrewer' package:
                cran.r-project.org/web/packages/RColorBrewer)''')
        settings.add_argument('--width_height',
                metavar = 'INCHES',
                nargs = 2,
                help='''The width and height of the graphics region in inches''')
        settings.add_argument('-s', '--same_page',
                            action='store_true',
                            help='''Plot multiple groups of p-values on the same graph''')
        settings.add_argument('-o', '--output',
                metavar = 'FILE',
                type = str,
                help='''Specify output graph filename. 
                Output is in pdf format. It can be converted to jpg format
                via the 'convert' command in Linux (e.g., convert -density 180 p.pdf p.jpg)''')
        labelling = parser.add_argument_group('variants/genes highlighting')
        labelling.add_argument('-b', '--bonferroni',
                            action='store_true',
                            help='''Plot the horizontal line at 0.05/N on Y-axis
                (significance level after Bonferroni correction)''')
        labelling.add_argument('-l', '--hlines',
                metavar = 'POSITION',
                nargs = '+',
                type=float,
                help='''Additional horizontal line(s) to
                be drawn on the Y-axis.''')
        labelling.add_argument('--label_top',
                metavar='INTEGER',
                            type=int,
                default=1,
                help='''Specify how many top hits (smallest p-values by rank)
                you want to highlight with their identifiers in text.''')
        labelling.add_argument('--label_these',
                metavar='NAME',
                            type=str,
                nargs = '+',
                help='''Specify the names of variants (chr:pos, e.g., 1:87463) 
                or genes (genename, e.g., IKBKB) you want to
                highlight with their identifiers in text.''')
        labelling.add_argument('-f','--font_size',
                metavar='SIZE',
                            type=float,
                default=2.5,
                help='''Font size of text labels. Default set to '2.5'.''')


CTHEME = ["Dark2", "grayscale", "default", "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", 
"RdYlBu", "RdYlGn", "Spectral","Accent", "Paired", 
"Pastel1", "Pastel2", "Set1", "Set2", "Set3", "Blues", 
"BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges",
"OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", 
"Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"]

class Plot:
    def __init__(self, args, inlines):
        '''Prepare R script to generate qq and manhattan plots'''
        if len(inlines) < 2:
            sys.exit("ERROR: Input file is empty.")
        self.a = args
        # processing chroms
        if hasattr(self.a, 'chrom'):
            allchroms = list(map(str, range(1,23))) + ['X','Y','Un']
            userange = [idx for idx, item in enumerate(self.a.chrom) if ':' in item and len(item.split(":")) == 2]
            for item in userange:
                start, end = self.a.chrom[item].split(':')
                start = allchroms.index(start)
                end = allchroms.index(end)
                if start < end:
                    self.a.chrom += allchroms[start:end+1]
                else:
                    self.a.chrom += allchroms[end:start+1]
            chrom = []
            for item in self.a.chrom:
                if item in allchroms:
                    if item not in chrom:
                        chrom.append(item)
            if len(chrom) == 0:
                sys.exit("ERROR: Invalid --chrom input")
            self.a.chrom = chrom
        # processing gene map
        if hasattr(self.a, 'gene_map'):
            self.gene_map = self.a.gene_map
        else:
            self.gene_map = None
        if self.gene_map and not os.path.exists(self.gene_map):
            sys.exit("ERROR: {0} does not exist".format(self.gene_map))
        self.fname = self.a.output if self.a.output else args.method
        self.m_read(inlines)

    def m_read(self, inlines):
        # fields headers
        headline = [x.replace('#','') for x in inlines[0]]
        self.fields = ['_'.join(x.split('_')[1:]) for x in headline if x.startswith('pvalue')]
        if len(self.fields) == 0:
            sys.exit("ERROR: Invalid headers (see header convention of 'vtools associate' command output).")
        idxes = []
        for field in self.fields:
            idxes.append(headline.index('pvalue_' + field))
        self.fields.insert(0,'GENE')
        # data
        self.is_snv = False
        self.dat = []
        badlines = []
        if 'chr' in headline and 'pos' in headline:
            if headline.index('chr') == 0 and headline.index('pos') == 1:
                # merge the 1st and 2nd columns if input data is for single SNV
                for idx, x in enumerate(inlines[1:]):
                    try:
                        self.dat.append([':'.join(x[0:2])] + [x[i] for i in idxes])
                    except IndexError:
                        badlines.append(str(idx+2))
                self.is_snv = True
        else:
                for idx, x in enumerate(inlines[1:]):
                    try:
                        self.dat.append([x[0]] + [x[i] for i in idxes])
                    except IndexError:
                        badlines.append(str(idx+2))
        if len(badlines):
            env.logger.warning("The following lines are ignored due to having empty fields: {}. You may want to manually fill them up with placeholder 'NaN'.".format(', '.join(badlines)))
        return

    def m_qqdata(self):
        rdat = OrderedDict()
        for x, y in zip(self.fields, list(zip(*self.dat))):
            if x == 'GENE':
                rdat[x] = list(y)
            else:
                try:
                    rdat[x] = [i if i == i else None for i in list(map(float, list(y)))]
                except ValueError as e:
                    sys.exit("ERROR: {}".format(e))
        self.rdat = Str4R(rdat)
        return

    def m_manhattandata(self):
        self.fields = [self.fields[0]] + ['CHR', 'BP'] + self.fields[1:]
        dat = []
        if self.is_snv:
            # for single variants can get their positions easily
            dat = [[x[0]] + x[0].split(":") + x[1:] for x in self.dat]
        else:
            # for genes have to annotate their positions
            # refgene is previously written to pickle
            #>>> output = open('refgene.pkl', 'wb')
            #>>> pickle.dump(gendict, output)
            #>>> output.close()
            pfname = os.path.join(env.local_resource, 'resource/refgene.pkl')
            if not os.path.exists(pfname):
                mkdir_p(os.path.join(env.local_resource, 'resource'))
                downloadFile('http://vtools.houstonbioinformatics.org/resource/refgene.pkl',
                             os.path.join(env.local_resource, 'resource')) 
            with open(pfname, 'rb') as f:
                gdict = pickle.load(f)
            if self.gene_map:
                # include user provided annotation list
                gd = {}
                with open(self.gene_map, 'r') as f:
                    try:
                        for item in f.readlines():
                            item = item.rstrip().split()
                            gd[item[0]] = [item[1],item[2]]
                    except:
                        sys.exit('ERROR: {0} format is not valid'.format(self.gene_map))
                gdict.update(gd)
                #
            failed = []
            multi_chrom = []
            i = 0
            for x in self.dat:
                try:
                    if len(gdict[x[0]]) > 2:
                        multi_chrom.append(x[0])
                        dat.append([x[0]] + gdict[x[0]][0:2] + x[1:])
                    else:
                        dat.append([x[0]] + gdict[x[0]] + x[1:])
                except KeyError:
                    failed.append(x[0])
                    dat.append([x[0]] + ['Un', '{0}'.format(i)] + x[1:])
                    i += 1000
            if len(multi_chrom):
                env.logger.warning('The following genes belong to more than one chromosomes. Using whatever first entry (alphanumeric order) as the coordinate: {0}'.format(', '.join(multi_chrom)))
            if len(failed):
                env.logger.warning('The following genes are not found in aviewer database. You may want to provide your own list of genomic coordinates of genes: {0}'.format(', '.join(failed)))
        #
        rdat = OrderedDict()
        for x, y in zip(self.fields, list(zip(*dat))):
            if x == 'BP':
                rdat[x] = [i if i == i else None for i in list(map(int, list(y)))]
            elif x in ['GENE', 'CHR']:
                rdat[x] = list(y)
            else:
                try:
                    rdat[x] = [i if i == i else None for i in list(map(float, list(y)))]
                except ValueError as e:
                    sys.exit("ERROR: {}".format(e))
        self.rdat = Str4R(rdat)
        return

    def qq(self):
        self.m_qqdata()
        astr = '''
        options(warn=1)
        pdfname <- paste("{0}", ".pdf", sep="")
        onePlot <- as.logical("{1}")
        color <- "{2}"
        shapeFixed <- as.logical("{3}")
        shapeValue <- as.numeric("{4}")
        title <- "{5}"
        gwLine <- as.logical("{6}")
        slopeLine <- as.logical("{7}")
        optLinestr <- "{8}"
        optLines <- -log10(as.numeric(strsplit(optLinestr, " ")[[1]]))
        topFont <- as.numeric("{9}")

        labelTopGenes <- as.numeric("{10}")
        annotatestr <- "{11}"
        annotate <- strsplit(annotatestr, " ")[[1]]
        tryCatch({{dat <- as.data.frame({12})}}, warning = function(ex) {{ paste("WARNING from Str4R:", str(ex)) }})
        qq_width <- {13}
        qq_height <- {14}
        '''.format(self.fname,
                'true' if self.a.same_page else 'false',
                'grayscale' if not self.a.color else self.a.color,
                'true' if self.a.fixed_shape or not self.a.same_page else 'false',
                self.a.shape,
                self.a.title,
                'true' if self.a.bonferroni else 'false',
                'false' if self.a.no_slope else 'true',
                ' ' if not self.a.hlines else ' '.join(list(map(str,self.a.hlines))),
                self.a.font_size,
                self.a.label_top,
                ' ' if not self.a.label_these else ' '.join(self.a.label_these),
                self.rdat,
                self.a.width_height[0] if self.a.width_height else 9,
                self.a.width_height[1] if self.a.width_height else 8,
                )
        return astr

    def manhattan(self):
        self.m_manhattandata()
        astr = '''
        options(warn=1)
        pdfname <- paste("{0}", ".pdf", sep="")
        facet <- as.logical("{1}")
        color <- "{2}"
        title <- "{3}"
        chromstr <- "{4}"
        chroms <- strsplit(chromstr, " ")[[1]]
        gwLine <- as.logical("{5}")
        optLinestr <- "{6}"
        optLines <- -log10(as.numeric(strsplit(optLinestr, " ")[[1]]))
        topFont <- as.numeric("{7}")
        labelTopGenes <- as.numeric("{8}")
        annotatestr <- "{9}"
        annotate <- strsplit(annotatestr, " ")[[1]]
        chrprefix <- {10}
        tryCatch({{dat <- as.data.frame({11})}}, warning = function(ex) {{ paste("WARNING from Str4R:", str(ex)) }})
        man_width <- {12}
        man_height <- {13}
        '''.format(self.fname,
                'true' if self.a.same_page else 'false',
                'grayscale' if not self.a.color else self.a.color,
                self.a.title,
                ' '.join(self.a.chrom),
                'true' if self.a.bonferroni else 'false',
                ' ' if not self.a.hlines else ' '.join(list(map(str,self.a.hlines))),
                self.a.font_size,
                self.a.label_top,
                ' ' if not self.a.label_these else ' '.join(self.a.label_these),
                'NULL' if self.a.chrom_prefix == "None" else repr(self.a.chrom_prefix),
                self.rdat,
                self.a.width_height[0] if self.a.width_height else 12,
                self.a.width_height[1] if self.a.width_height else 6,
                )
        return astr

QQ_FOO = '''
#! QQplot function
QQplot <- function(dat, multiple=T, color='default', shapeFixed=F, shapeValue=1, title="", gwLine=T, slopeLine=T, optLines=c(), topFont=3, labelTopGenes=0, annotate=NULL, index=2)
{
	suppressMessages(library(ggplot2))
	suppressMessages(library(RColorBrewer))
	if (labelTopGenes>0 || length(annotate)>0) label=T
	else label=F
	qqDat <- data.frame()
	qqTopHit <- data.frame()
	for( i in 2:ncol(dat) ) {
		one <- dat[, c(1,i)]
		names(one) <- c("gene", "obs")
		one$obs <- as.numeric(one$obs)
		one <- subset(na.omit(one[order(one$obs), ]), (obs>0 & obs<=1))
		one$logObs <- -log10(as.numeric(one$obs))
		one$ept <- -log10(c(1:length(one$obs))/(1+length(one$obs)))
		one$method <- colnames(dat)[i]
		if (label) {
			topGenes <- one[which(one$logObs>=sort(one$logObs, decreasing=TRUE, na.last=NA)[labelTopGenes]), ]$gene
			qqTopHit<- rbind(qqTopHit, one[which(one$gene %in% union(topGenes, annotate)), ])
		}
		qqDat <- rbind(qqDat, one)
	}
	shapeValue <- as.numeric(shapeValue)
	if(is.na(shapeValue)) shapeValue <- 1
	else shapeValue <- shapeValue%%20
	if (multiple) {
		if(shapeFixed) {
			pobj <- ggplot(qqDat, aes(x=ept, y=logObs, colour=factor(method))) + geom_point(binwidth = 0.8, size=1.5, shape=c(shapeValue))
		} else {
			pobj <- ggplot(qqDat, aes(x=ept, y=logObs, colour=factor(method), shape=factor(method))) + geom_point(binwidth = 0.8,size=1.5) +
				scale_shape_manual(values = seq(shapeValue, shapeValue+length(unique(qqDat$method))-1, 1)%%20)
		}
		if (color!='default') {
			mycols <- rep(c("gray10", "gray50"),20)[1:length(unique(qqDat$method))]
			if(color!='grayscale') mycols <- rep(brewer.pal(brewer.pal.info[color,]$maxcolors, name = color), 4)[1:length(unique(qqDat$method))]
			pobj <- pobj + scale_colour_manual(values=mycols)
		}
		pobj <- pobj +
		opts(legend.key = theme_rect(colour = 'white', fill = 'white', size = 0.5, linetype='blank')) +
		opts(legend.justification=c(1,0), legend.position=c(1,0)) +
		opts(legend.background = theme_rect(fill="transparent", size=.5, linetype="blank")) +
		opts(legend.title=theme_blank()) +
		scale_x_continuous(limits=c(0,max(qqDat$ept)), breaks=0:max(qqDat$ept), labels=0:max(qqDat$ept))
		if(label && nrow(qqTopHit)>0) pobj <- pobj + geom_text(data = qqTopHit, aes(ept, logObs, label = gene, colour=factor(method)), vjust = 1.5, size = as.numeric(topFont), show_guide = FALSE, legend = FALSE)
	} else {
		if(!shapeFixed) shapeValue <- (shapeValue+as.numeric(index))%%20
        mycols <- NA
        if (color!='default') {
			mycols <- rep(c("gray10", "gray50"),20)[1:length(unique(qqDat$method))]
			if(color!='grayscale') mycols <- rep(brewer.pal(brewer.pal.info[color,]$maxcolors, name = color), 4)[1:length(unique(qqDat$method))]
		    pobj <- ggplot(qqDat, aes(x=ept, y=logObs)) + geom_point(binwidth = 0.8, size=2.0, colour=mycols, shape=shapeValue)
		} else {
		    pobj <- ggplot(qqDat, aes(x=ept, y=logObs)) + geom_point(binwidth = 0.8, size=2.0, shape=shapeValue)
        }
		pobj <- pobj + opts(legend.position="none") + scale_x_continuous(limits=c(0,max(qqDat$ept)), breaks=0:max(qqDat$ept), labels=0:max(qqDat$ept))
		if(label && nrow(qqTopHit)>0) {
            if (!is.na(mycols)) {
                pobj <- pobj + geom_text(data = qqTopHit, aes(ept, logObs, label = gene), colour=mycols, vjust = 1.5, size = as.numeric(topFont), show_guide = FALSE, legend = FALSE)
            } else {
                pobj <- pobj + geom_text(data = qqTopHit, aes(ept, logObs, label = gene), vjust = 1.5, size = as.numeric(topFont), show_guide = FALSE, legend = FALSE)
            }
        }
	}
	pobj <- pobj +
	opts(title=paste(title,'\\n')) +
	opts(axis.text.x = theme_text(angle=0, size=12)) +
	opts(axis.text.y = theme_text(angle=0, size=12)) +
	opts(axis.title.y = theme_text(face="bold", colour= "black", angle=90, size=13)) +
	opts(axis.title.x = theme_text(face="bold", colour= "black", size=13)) +
	#opts(plot.title = theme_text(size=16, face="bold")) +
	opts(plot.title = theme_text(size=14)) +
    opts(panel.background = theme_rect(fill = "white",colour = "black")) +
	opts(panel.grid.minor=theme_blank()) +
	opts(plot.background = theme_rect(fill = "transparent",colour = NA)) +
	#opts(strip.text.x = theme_text(face = 'bold')) +
	#xlab(expression(-log[10](italic(p[expected])))) + ylab(expression(-log[10](italic(p[observed]))))
	xlab(expression(paste(-log[10], " " ,italic(p)[expected]))) + ylab(expression(paste(-log[10], " " ,italic(p)[observed])))

	if (slopeLine) pobj <- pobj + geom_abline(intercept = 0, slope = 1, col = "black", linetype="dotted")
	if (gwLine) {
		gwPval = -log10(0.05/length(unique(qqDat$gene)))
		pobj <- pobj + geom_hline(yintercept=gwPval,colour="red", linetype="dashed")
	}
	if (length(optLines) > 0){
		for( i in optLines ) {
			pobj <- pobj + geom_hline(yintercept=i,colour="blue", linetype="dashed")
		}
	}
	#label significant
	return(pobj)
}
'''
MANHATTAN_FOO = '''
# Stephen Turner
# http://StephenTurner.us/
# http://GettingGeneticsDone.blogspot.com/
# See license at http://gettinggeneticsdone.blogspot.com/p/copyright.html

# Modified by Zongxiao He, Thu Jun  7 16:34:07 CDT 2012
# http://gettinggeneticsdone.blogspot.com/2011/04/annotated-manhattan-plots-and-qq-plots.html
# This is input data format
#    GENE CHR    BP       P
#    ZAN   1 235800006 0.62220
#  DCP1B   1  46100028 0.06195
#  ABCA7   1 143700035 0.10700
#   SINA   1 202300047 0.47280
#   SOHU   1  66400050      NA
#   HAHA   1  64900051 0.53430
#       Accent    8 
#       Dark2     8 
#       Paired   12 
#       Pastel1   9 
#       Pastel2   8 
#       Set1      9 
#       Set2      8 
#       Set3     12 

#! manhattan plot using ggplot2
manhattanplot <- function(dataframe, facet=F, color='default', title='', chroms=c(1:22,'X','Y','Un'), chrprefix=chrprefix, gwLine=T, optLines=c(), topFont=2, labelTopGenes=0, annotate=NULL) {
    d=dataframe
    if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
    if ((labelTopGenes>0 || length(annotate)>0) && !("GENE" %in% names(d))) stop("You requested annotation but your data doesn't contains column GENE")
    if (facet && !("FACET" %in% names(d))) stop("You requested facet the plot but but your data doesn't contains column FACET")
    suppressMessages(library(ggplot2))
    suppressMessages(library(plyr))
    suppressMessages(library(RColorBrewer))
    d$GENE <- as.character(d$GENE)
    d$CHR <- as.character(d$CHR)
    d$CHR[-which(d$CHR %in% as.character(c(1:22,'X','Y')))] <- 'Un'
    d=d[d$CHR %in% chroms, ]
    d$CHR[d$CHR=="X"] <- 23
    d$CHR[d$CHR=="Y"] <- 24
    d$CHR[d$CHR=="Un"] <- 25
    d$CHR <- as.numeric(d$CHR)
    d$P <- as.numeric(d$P)
    d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
    d$logp = -log10(d$P)
    d$pos=NA
    ticks=NULL
    lastbase=0
    ymax<-ceiling(max(d$logp))
    #if (ymax<6) ymax<-6.5 # exome rare variant association, -log10(0.05/20000) = 5.60206
    numchroms=length(unique(d$CHR))
    if (numchroms==1) {
        d$pos=d$BP
        ticks=floor(length(d$pos))/2+1
    } else {
        index <- sort(unique(d$CHR))
        for (i in unique(d$CHR)) {
            if (which(index==i)==1) {
                d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
            } else {
                lastbase=lastbase+tail(subset(d,CHR==index[which(index==i)-1])$BP, 1)
                d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
            }
            ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
        }
    }
    if (numchroms==1) {
        plotObj <- ggplot(d, aes(x=pos, y=logp)) + geom_point(binwidth = 0.8,size=1.5) +
        ylab(expression(paste(-log[10], " ", italic(p)))) + xlab(paste("Chromosome",unique(d$CHR),"position")) + scale_x_continuous(labels=NULL)
    }else {
        plotObj <- ggplot(d, aes(x=pos, y=logp, group=factor(CHR), colour=factor(CHR)))  + geom_point(binwidth = 0.8,size=1.5)
        label <- unique(d$CHR)
        label[which(label==23)] <- 'X'
        label[which(label==24)] <- 'Y'
        label[which(label==25)] <- 'Un'
        label <- lapply(label, function(x) paste(chrprefix,x,sep=''))
        plotObj <- plotObj +
        ylab(expression(paste(-log[10], " ", italic(p)))) +
        scale_x_continuous(name="Chromosome", breaks=ticks, labels=label) +
        scale_y_continuous(limits=c(0,ymax), breaks=seq(1, ymax, by=1), labels=seq(1, ymax, by=1))
        if(color!='default') {
            mycols <- rep(c("gray10", "gray50"),20)[1:length(unique(d$CHR))]
            if(color!='grayscale') mycols <- rep(brewer.pal(brewer.pal.info[color,]$maxcolors, name = color), 4)[1:length(unique(d$CHR))]
            plotObj <- plotObj + scale_colour_manual(values=mycols)
        }
    }
    #print(d[c(1:20),])

    if (facet) plotObj <- plotObj + facet_grid(FACET ~ ., scales = "free", space="free")
    if(length(annotate)>0){
        d.annotate <- d[which(d$GENE %in% annotate), ]
        if (nrow(d.annotate)) plotObj <- plotObj + geom_text(data = d.annotate, aes(pos, logp, label = GENE), face='bold', colour='black', vjust = -1.0, size = as.numeric(topFont)+1, show_guide = FALSE, legend = FALSE)
    }
    if(labelTopGenes > 0){
        d.tops=data.frame()
        if(facet && "FACET" %in% names(d)) {
            d.tops <- ddply(d, "FACET", function(x) {cutoff = sort(x$logp, decreasing=TRUE, na.last=NA)[labelTopGenes];  x[which(x$logp>=cutoff), ]})
            if(length(annotate)>0) d.tops <- ddply(d.tops, "FACET", function(x) {x[which(!(x$GENE %in% annotate)), ]})
        } else {
            d.tops <- d[which(d$logp>=sort(d$logp, decreasing=TRUE, na.last=NA)[labelTopGenes]), ]
            if(length(annotate)>0) d.tops <- d.tops[which(!(d.tops$GENE %in% annotate)), ]
        }
        if (nrow(d.tops)) plotObj <- plotObj + geom_text(data = d.tops, aes(pos, logp, label = GENE, colour=factor(CHR)), vjust = -1.0, size = as.numeric(topFont), show_guide = FALSE, legend = FALSE)
    }
    agl = 90
    if (is.null(chrprefix)) agl = 0
    plotObj <- plotObj +
    opts(legend.position = "none") +
    opts(title=paste(title, '\\n')) +
    opts(panel.background=theme_rect(fill="white", colour="black")) +
    opts(panel.grid.minor=theme_blank()) +
    opts(axis.text.x = theme_text(angle=agl, size=11)) +
    opts(axis.text.y = theme_text(angle=0, size=12)) +
    opts(axis.title.y = theme_text(colour= "black", angle=90, size=13)) +
    opts(axis.title.x = theme_text(colour= "black", size=13))
    #opts(axis.ticks=theme_segment(colour=NA))

    gwPval = -log10(0.05/length(unique(d$GENE)))
    if (gwLine) plotObj <- plotObj + geom_hline(yintercept=gwPval,colour="red", linetype="dashed")
    if (length(optLines) > 0){
        for( i in optLines ) {
            plotObj <- plotObj + geom_hline(yintercept=i,colour="blue", linetype="dashed")
        }
    }
    return(plotObj)
}
'''
MANHATTAN_PLAIN_FOO = '''
# Stephen Turner
# http://StephenTurner.us/
# http://GettingGeneticsDone.blogspot.com/
# See license at http://gettinggeneticsdone.blogspot.com/p/copyright.html

# Last updated: Tuesday, April 19, 2011
# R code for making manhattan plots and QQ plots from plink output files. 
# manhattan() with GWAS data this can take a lot of memory, recommended for use on 64bit machines only, for now. 
# Altnernatively, use bmanhattan() , i.e., base manhattan. uses base graphics. way faster.

# Modified by Zongxiao He, Thu July 19 16:34:07 CDT 2012
# http://gettinggeneticsdone.blogspot.com/2011/04/annotated-manhattan-plots-and-qq-plots.html

# manhattan plot using base graphics
manhattanplainplot <- function(dataframe, facetLastOne=T, color='default', title='', subtitle=subtitle, chroms=c(1:22,'X','Y','Un'), chrprefix=chrprefix, gwLine=T, optLines=c(), topFont=1, labelTopGenes=0, annotate=NULL) {
    d=dataframe
    if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
    if ((labelTopGenes>0 || length(annotate)>0) && !("GENE" %in% names(d))) stop("You requested annotation but your data doesn't contains column GENE")
    suppressMessages(library(plyr))
    suppressMessages(library(RColorBrewer))
    d$GENE <- as.character(d$GENE)
    d$CHR <- as.character(d$CHR)
    d$CHR[-which(d$CHR %in% as.character(c(1:22,'X','Y')))] <- 'Un'
    d=d[d$CHR %in% chroms, ]
    d$CHR[d$CHR=="X"] <- 23
    d$CHR[d$CHR=="Y"] <- 24
    d$CHR[d$CHR=="Un"] <- 25
    d$CHR <- as.numeric(d$CHR)
    d$P <- as.numeric(d$P)
    d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
    d$logp = -log10(d$P)
    d$pos=NA
    ticks=NULL
    lastbase=0
    ymax<-ceiling(max(d$logp))
    #if (ymax<6) ymax<-6.5 # exome rare variant association, -log10(0.05/20000) = 5.60206
    numchroms=length(unique(d$CHR))

    mycols <- rep(c("gray10", "gray50"),20)[1:length(unique(d$CHR))]
    if(color!='grayscale') {
        tryCatch( {mycols <- rep(brewer.pal(brewer.pal.info[color,]$maxcolors, name = color), 4)[1:length(unique(d$CHR))]}, error = function(e) { mycols <- rep(c("gray10", "gray50"),20)[1:length(unique(d$CHR))] } )
}

    index <- sort(unique(d$CHR))
    if (numchroms==1) {
        d$pos=d$BP
        ticks=floor(length(d$pos))/2+1
    } else {
        for (i in unique(d$CHR)) {
            if (which(index==i)==1) {
                d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
            } else {
                lastbase=lastbase+tail(subset(d,CHR==index[which(index==i)-1])$BP, 1)
                d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
            }
            ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
        }
    }

    label <- sort(unique(d$CHR))
    label[which(label==23)] <- 'X'
    label[which(label==24)] <- 'Y'
    label[which(label==25)] <- 'Un'
    label <- lapply(label, function(x) paste(chrprefix,x,sep=''))

    if (numchroms == 1) {
        with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(paste(-log[10], " ", italic(p))), cex.lab=1, xlab=paste("Chromosome",unique(d$CHR),"position"), xaxt="n", main=title))
        with(d, points(pos, logp, col=mycols[1]))
    } else {
        if(facetLastOne) {
            agl = 3
            if (is.null(chrprefix)) agl = 1
            with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(paste(-log[10], " ", italic(p))), cex.lab=1, xlab="Chromosome", xaxt="n", type="n", main=title, font.main= 4))
            axis(1, at=ticks, lab=label, las=agl)
        } else {
            with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(paste(-log[10], " ", italic(p))), cex.lab=1, xlab="", xaxt="n", type="n", main=NULL))
            axis(1, at=ticks, labels=FALSE, ann=FALSE)
        }
        icol=1
        for (i in sort(unique(d$CHR))) {
            with(d[d$CHR==i, ],points(pos, logp, col=mycols[icol]))
            icol=icol+1
        }
    }

    mtext(subtitle,outer=F,side=3,line=0.2) #subtitle
    chr2col=data.frame(chr=index, col=mycols)
    if(length(annotate)>0){
        d.annotate <- d[which(d$GENE %in% annotate), ]
        if (nrow(d.annotate)) {
            cols <- d.annotate$CHR
            for(i in 1:length(cols)) cols[i] <- as.character(chr2col$col[chr2col$chr==cols[i]])
            with(d.annotate, text(pos, logp, GENE, adj=c(0.5,1.5), cex=as.numeric(topFont),  col=cols))
            with(d.annotate, points(pos, logp, col=cols, pch=4))
        }
    }

    if(labelTopGenes > 0){
        d.tops <- d[which(d$logp>=sort(d$logp, decreasing=TRUE, na.last=NA)[labelTopGenes]), ]
        if(length(annotate)>0) d.tops <- d.tops[which(!(d.tops$GENE %in% annotate)), ] # in case double highlight
        if (nrow(d.tops)){
            #print(chr2col$col[chr2col$chr %in% d.tops$CHR])
            cols <- d.tops$CHR
            for(i in 1:length(cols)) cols[i] <- as.character(chr2col$col[chr2col$chr==cols[i]])
            with(d.tops, text(pos, logp, GENE, adj=c(0.5,1.5), cex=as.numeric(topFont),  col=cols))
            with(d.tops, points(pos, logp, col=cols, pch=19))
        }
    }

    gwPval = -log10(0.05/length(d$GENE))
    if (gwLine) abline(h=gwPval, col="red", lty=2) 
    if (length(optLines) > 0){
        for( i in optLines ) {
            abline(h=i,,col="blue", lty=3)
        }
    }
}
'''
QQ_MAIN = '''
if(topFont<1.5) topFont <- 1.5
if(topFont>8) topFont <- 8
#tryCatch(  {dat <- read.delim(pipe('cat /dev/stdin'), header=T, stringsAsFactors=F)}, error = function(e) { quit("no") } )
if (ncol(dat)<2) stop("data doesn't have enought columns")
pdf(pdfname, qq_width, qq_height)
if (onePlot) {
    pObj <- QQplot(dat, multiple=T, color=color, shapeFixed=shapeFixed, shapeValue=shapeValue, title=title, gwLine=gwLine, slopeLine=slopeLine, optLines=optLines, topFont=topFont, labelTopGenes=labelTopGenes, annotate=annotate)
    print(pObj)
} else {
    for( i in 2:ncol(dat) ) {
        pval <- dat[, c(1,i)]
        subtitle <- title
        if (ncol(dat)>2) subtitle <- paste(title, " ", colnames(dat)[i], " method", sep="")
        pObj <- QQplot(pval, multiple=F, color=color, shapeFixed=shapeFixed, shapeValue=shapeValue, title=subtitle, gwLine=gwLine, slopeLine=slopeLine, optLines=optLines, topFont=topFont, labelTopGenes=labelTopGenes, annotate=annotate, index=i)
        print(pObj)
    }
}
dev.off()
'''
MANHATTAN_MAIN = '''
if(topFont<1.5) topFont <- 1.5
if(topFont>10) topFont <- 10
#tryCatch(  {dat <- read.delim(pipe('cat /dev/stdin'), header=T, stringsAsFactors=F)}, error = function(e) { quit("no") } )
if (ncol(dat)<4) stop("data doesn't have enought columns")
if (ncol(dat)==4) {
    names(dat) <- c("GENE", "CHR", "BP", "P")
    pObj <- manhattanplot(dat, facet=F, color=color, title=title, chroms=chroms, chrprefix=chrprefix, gwLine=gwLine, optLines=optLines, topFont=topFont, labelTopGenes=labelTopGenes, annotate=annotate)
    pdf(pdfname, man_width, man_height)
    print(pObj)
    dev.off()
} else {
    if (facet) {
        manDat <- data.frame()
        for( i in 4:ncol(dat) ) {
            pval <- dat[, c(1,2,3,i)]
            names(pval) <- c("GENE", "CHR", "BP", "P")
            pval$FACET <- colnames(dat)[i]
            manDat <- rbind(manDat, pval)
        }
        names(manDat) <- c("GENE", "CHR", "BP", "P", "FACET")
        pdf(pdfname, man_width, man_height*length(levels(factor(manDat$FACET))))
        pObj <- manhattanplot(manDat,facet=T, color=color, title=title, chroms=chroms, chrprefix=chrprefix, gwLine=gwLine, optLines=optLines, topFont=topFont, labelTopGenes=labelTopGenes, annotate=annotate)
        print(pObj)
        dev.off()
    } else {
        pdf(pdfname, man_width, man_height)
        for( i in 4:ncol(dat) ) {
            pval <- dat[, c(1,2,3,i)]
            names(pval) <- c("GENE", "CHR", "BP", "P")
            subtitle <- title
            if (ncol(dat)>4) subtitle <- paste(title, " ", colnames(dat)[i], " method", sep="")
            pObj <- manhattanplot(pval,facet=F, color=color, title=subtitle, chroms=chroms, chrprefix=chrprefix, gwLine=gwLine, optLines=optLines, topFont=topFont, labelTopGenes=labelTopGenes, annotate=annotate)
            print(pObj)
        }
        dev.off()
    }
}
'''
MANHATTAN_PLAIN_MAIN = '''
if(topFont<0.5 | topFont>2) topFont <- 1
#tryCatch(  {dat <- read.delim(pipe('cat /dev/stdin'), header=T, stringsAsFactors=F)}, error = function(e) { quit("no") } )

if (ncol(dat)<4) stop("data doesn't have enought columns")
if (ncol(dat)==4) facet=F

if (facet) {
    pdf(pdfname, man_width, man_height*(ncol(dat)-3))
    par(mfrow=c(ncol(dat)-3, 1))
    for( i in 4:ncol(dat) ) {
        pval <- dat[, c(1,2,3,i)]
        names(pval) <- c("GENE", "CHR", "BP", "P")
        subtitle <- paste(colnames(dat)[i], " method", sep="")
        if(i==ncol(dat)) {
            lastOne=T
            par(mar=c(5,5,2,2)+0.1)
        } else {
            lastOne=F  #no x lab
            if(i==4) par(mar=c(2,5,5,2)+0.1) 
            else par(mar=c(2,5,2,2)+0.1)
        }
        manhattanplainplot(pval, facetLastOne=lastOne, color=color, title="", subtitle=subtitle, chroms=chroms, chrprefix=chrprefix, gwLine=gwLine, optLines=optLines, topFont=topFont, labelTopGenes=labelTopGenes, annotate=annotate)
    }
    title(title, outer=T, line=-2, font.main= 4, cex.main=2)
    dev.off()
} else {
	pdf(pdfname, man_width, man_height)
	par(mar=c(5,5,4,2)+0.1)
    for( i in 4:ncol(dat) ) {
        pval <- dat[, c(1,2,3,i)]
        names(pval) <- c("GENE", "CHR", "BP", "P")
        subtitle <- ""
        if (ncol(dat) > 4) subtitle <- paste(colnames(dat)[i], " method", sep="")
        manhattanplainplot(pval,facetLastOne=T, color=color, title=title, subtitle=subtitle, chroms=chroms, chrprefix=chrprefix, gwLine=gwLine, optLines=optLines, topFont=topFont, labelTopGenes=labelTopGenes, annotate=annotate)
    }
	dev.off()
}
'''
