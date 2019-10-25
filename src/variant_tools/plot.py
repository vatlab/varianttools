#!/usr/bin/env python
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit https://github.com/vatlab/varianttools for details.
#
# Copyright (C) 2012 Bo Peng and Gao Wang
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

import os
import pickle
import subprocess
import sys
from collections import OrderedDict

from .rtester import Str4R
from .utils import downloadFile, env, mkdir_p, runCommand, whereisRPackage


def size(bytes):
    system = [
        (1024**5, 'P'),
        (1024**4, 'T'),
        (1024**3, 'G'),
        (1024**2, 'M'),
        (1024**1, 'K'),
        (1024**0, 'B'),
    ]

    for factor, suffix in system:
        if bytes >= factor:
            break
    amount = int(bytes / factor)
    if isinstance(suffix, tuple):
        singular, multiple = suffix
        if amount == 1:
            suffix = singular
        else:
            suffix = multiple
    return str(amount) + suffix


def ismissing(x):
    return (x != x or x.lower() in [
        'none', 'null', 'nan', 'na', '.', '-', '', '\t', '\n', ' '
    ])


CTHEME = [
    "Dark2", "grayscale", "default", "BrBG", "PiYG", "PRGn", "PuOr", "RdBu",
    "RdGy", "RdYlBu", "RdYlGn", "Spectral", "Accent", "Paired", "Pastel1",
    "Pastel2", "Set1", "Set2", "Set3", "Blues", "BuGn", "BuPu", "GnBu",
    "Greens", "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples",
    "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"
]


def RdeviceFromFilename(filename, width=800, height=600):
    # guess the device used to plot R plot from filename
    basename, ext = os.path.splitext(filename)
    # default
    device = 'postscript'
    params = ''
    try:
        # functions are not available in, for example, R 2.6.2
        if ext.lower() == '.pdf':
            device = 'pdf'
            params = ', width={}/90, height={}/90'.format(width, height)
        elif ext.lower() == '.png':
            device = 'png'
            params = ', width={}, height={}'.format(width, height)
        elif ext.lower() == '.bmp':
            device = 'bmp'
            params = ', width={}, height={}'.format(width, height)
        elif ext.lower() in ['.jpg', '.jpeg']:
            device = 'jpeg'
            params = ', width={}, height={}'.format(width, height)
        elif ext.lower() in ['.tif', '.tiff']:
            device = 'tiff'
            params = ', width={}, height={}'.format(width, height)
        elif ext.lower() == '.eps':
            device = 'postscript'
            params = ', width={}/90, height={}/90'.format(width, height)
    except Exception as e:
        env.logger.warning(
            'Can not determine which device to use to save file {}. A postscript driver is used: {}'
            .format(filename, e))
        device = 'postscript'
    return '{}("{}" {})'.format(device, filename, params)


def resolvePlotFilename(name, fields):
    if len(fields) == 1:
        return [name]
    else:
        if name is None:
            return [None for x in fields]
        else:
            return [
                os.path.splitext(name)[0] + '_{}'.format(x) +
                os.path.splitext(name)[1] for x in fields
            ]


def executeRScript(script):
    # write script to log file for debugging and customization purposes
    env.logger.info(
        'Running R script (complete script available in vtools_report.log)'
        .format(script))
    env.logger.debug(script)
    # start R
    process = subprocess.Popen(['R', '--slave', '--no-save', '--no-restore'],
                               stdin=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    # send script and get standard output and error
    out, err = process.communicate(script)
    if err.strip():
        env.logger.warning(err)
    return process.wait()


def loadGgplot(script):
    # load R libraries
    for l in ['ggplot2', 'scales', 'RColorBrewer', 'plyr']:
        rlib = whereisRPackage('{}'.format(l))
        script = ('\nsuppressMessages(library("{}", lib.loc="{}"))'.format(
            l, rlib) if rlib else
                  '\nsuppressMessages(library("{}"))'.format(l)) + script
    return script


def rhist(data,
          output,
          width,
          height,
          vlines=None,
          normcurve=True,
          save_data=None,
          save_script=None):
    '''draw histogram using ggplot2'''
    if vlines:
        vlines = Str4R(vlines)
    else:
        vlines = 'NULL'
    rstr = '''
        tryCatch( {{dat <- read.delim(pipe('cat /dev/stdin'), header=T, stringsAsFactors=F)
        }}, error = function(e) {{ quit("no") }} )
        stat_foo <- NULL
        vlines <- {0}
        fns <- c({2})
        for (i in 1:ncol(dat)) {{
        if ({1}) stat_foo <- function(x) dnorm(x, mean(dat[,i], na.rm=T), sd(dat[,i], na.rm=T))
        d <- subset(dat, select=c(colnames(dat)[i]))
        p <- gghist(d, stat_foo, vlines=vlines, xname=colnames(dat)[i])
        eval(parse(text=fns[i]))
        print(p)
        graphics.off()
        }}'''.format(
        vlines, int(normcurve),
        ','.join([repr(RdeviceFromFilename(x, width, height)) for x in output]))
    #
    # Here we pipe data from standard input (which can be big),
    # create a script dynamically, and pump it to a R process.
    #
    if save_data:
        with open(save_data, 'w') as f:
            f.write(data)
    rstr = HIST_FOO + loadGgplot(rstr)
    rscript = save_script if save_script else os.path.join(
        env.cache_dir, 'hist.R')
    with open(rscript, 'w') as f:
        f.write(rstr)
    env.logger.info('Generating histogram {} ...'.format(repr(output)))
    cmd = "Rscript {} --slave --no-save --no-restore".format(rscript)
    out = runCommand(cmd, data)
    if out:
        sys.stdout.write(out)
    env.logger.info("Complete!")
    # clean up
    if save_script is None:
        os.remove(rscript)
    return


def rdot(data,
         output,
         width,
         height,
         psize=2.5,
         color=None,
         save_data=None,
         save_script=None):
    '''draw dotplot using ggplot2
    input can be 1, 2 or 3 columns for x, y, and z (i.e., dot colors) axis'''
    rstr = '''
        tryCatch( {{dat <- read.delim(pipe('cat /dev/stdin'), header=T, stringsAsFactors=F)
        }}, error = function(e) {{ quit("no") }} )
        p <- ggdot(dat, psize = {0}, color = {1})
        eval(parse(text={2}))
        print(p)
        graphics.off()
        '''.format(psize,
                   repr(color) if color else "NULL",
                   repr(RdeviceFromFilename(output, width, height)))
    #
    # Here we pipe data from standard input (which can be big),
    # create a script dynamically, and pump it to a R process.
    #
    if save_data:
        with open(save_data, 'w') as f:
            f.write(data)
    rstr = DOT_FOO + loadGgplot(rstr)
    rscript = save_script if save_script else os.path.join(
        env.cache_dir, 'hist.R')
    with open(rscript, 'w') as f:
        f.write(rstr)
    env.logger.info('Generating dot plot {} ...'.format(repr(output)))
    cmd = "Rscript {} --slave --no-save --no-restore".format(rscript)
    out = runCommand(cmd, data)
    if out:
        sys.stdout.write(out)
    env.logger.info("Complete!")
    # clean up
    if save_script is None:
        os.remove(rscript)
    return


def rbox(raw_data,
         fields,
         stratify,
         output,
         width,
         height,
         psize=2,
         color=None,
         save_data=None,
         save_script=None):
    '''draw box using ggplot2
    input should be 3 columns: 1st column is values, 2nd column is bin ID's and 3rd column is bin order'''
    if stratify is not None:
        if len(fields) > 2:
            raise ValueError(
                'Less than 2 input fields are allowed with --stratify option')
        # create strata for given field
        strata = sorted(stratify)
        data = '{}\t{}\torder\n'.\
          format(
              fields[0].replace('_','.'),
              (fields[-1] if fields[-1]!=fields[0] else fields[-1] + '.strata').replace('_','.')
              )
        for item in raw_data.split('\n')[1:]:
            item = item.split('\t')
            if item[0] == 'NA':
                continue
            else:
                found = False
                for i, s in enumerate(strata):
                    if float(item[-1]) < s:
                        found = True
                        if i == 0:
                            xaxis = ('Below {}'.format(s))
                        else:
                            xaxis = ('{}-{}'.format(strata[i - 1], s))
                        order = i + 1
                        break
                if found is not True:
                    xaxis = '{} or more'.format(strata[-1])
                    order = len(strata) + 1
            data += '{}\t{}\t{}\n'.format(item[0], xaxis, order)
    else:
        # stack data to 3 columns retaining the input order of fields
        data = 'values\tvariables\torder\n'
        for item in [x.split() for x in raw_data.split('\n')[1:]]:
            for i, x in enumerate(item):
                data += '{}\t{}\t{}\n'.format(x, fields[i], i + 1)
    rstr = '''
        tryCatch( {{dat <- read.delim(pipe('cat /dev/stdin'), header=T, stringsAsFactors=F)
        }}, error = function(e) {{ quit("no") }} )
        p <- ggbox(dat, psize = {0}, color = {1})
        eval(parse(text={2}))
        print(p)
        graphics.off()
        '''.format(psize,
                   repr(color) if color else "NULL",
                   repr(RdeviceFromFilename(output, width, height)))
    #
    # Here we pipe data from standard input (which can be big),
    # create a script dynamically, and pump it to a R process.
    #
    if save_data:
        with open(save_data, 'w') as f:
            f.write(data)
    rstr = BOX_FOO + loadGgplot(rstr)
    rscript = save_script if save_script else os.path.join(
        env.cache_dir, 'hist.R')
    with open(rscript, 'w') as f:
        f.write(rstr)
    env.logger.info('Generating dot plot {} ...'.format(repr(output)))
    cmd = "Rscript {} --slave --no-save --no-restore".format(rscript)
    out = runCommand(cmd, data)
    if out:
        sys.stdout.write(out)
    env.logger.info("Complete!")
    # clean up
    if save_script is None:
        os.remove(rscript)
    return


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
    env.logger.info("Processing {} of input data ...".format(
        size(len(data) + data_size)))
    p = PlotAssociation(args, data)

    pinput = eval(args.method.upper() + '_FOO') + eval(
        'p.{}'.format(args.method if args.method == 'qq' else 'manhattan'))(
        ) + eval(args.method.upper() + '_MAIN')

    pinput = loadGgplot(pinput)
    env.logger.info("Generating graph(s) ...")
    cmd = "R --slave --no-save --no-restore"
    out = runCommand(cmd, pinput)
    if out:
        sys.stdout.write(out)
    env.logger.info("Complete!")
    return


class PlotAssociation:

    def __init__(self, args, inlines):
        '''Prepare R script to generate qq and manhattan plots'''
        if len(inlines) < 2:
            raise ValueError("Input file is empty.")
        self.a = args
        if args.method == 'manhattan_plain':
            self.a.font_size = self.a.font_size / 2.5
        # processing chroms
        if hasattr(self.a, 'chrom'):
            allchroms = list(map(str, list(range(1, 23)))) + ['X', 'Y', 'Un']
            userange = [
                idx for idx, item in enumerate(self.a.chrom)
                if ':' in item and len(item.split(":")) == 2
            ]
            for item in userange:
                start, end = self.a.chrom[item].split(':')
                start = allchroms.index(start)
                end = allchroms.index(end)
                if start < end:
                    self.a.chrom += allchroms[start:end + 1]
                else:
                    self.a.chrom += allchroms[end:start + 1]
            chrom = []
            for item in self.a.chrom:
                if item in allchroms:
                    if item not in chrom:
                        chrom.append(item)
            if len(chrom) == 0:
                raise ValueError("Invalid --chrom input")
            self.a.chrom = chrom
        else:
            self.a.chrom = ''
        if not hasattr(self.a, 'chrom_prefix'):
            self.a.chrom_prefix = "None"
        # processing gene map
        if hasattr(self.a, 'gene_map'):
            self.gene_map = self.a.gene_map
        else:
            self.gene_map = None
        if self.gene_map and not os.path.exists(self.gene_map):
            raise ValueError("{0} does not exist".format(self.gene_map))
        self.fname = self.a.output if self.a.output else args.method
        self.m_read(inlines)

    def m_read(self, inlines):
        # fields headers
        headline = [x.replace('#', '') for x in inlines[0]]
        self.fields = [
            '_'.join(x.split('_')[1:])
            for x in headline
            if x.startswith('pvalue')
        ]
        if len(self.fields) == 0:
            raise ValueError(
                "Invalid headers (see header convention of 'vtools associate' command output)."
            )
        idxes = []
        for field in self.fields:
            idxes.append(headline.index('pvalue_' + field))
        self.fields.insert(0, 'GENE')
        # data
        self.is_snv = False
        self.dat = []
        badlines = []
        if 'chr' in headline[0] and 'pos' in headline[1]:
            # merge the 1st and 2nd columns if input data is for single SNV
            for idx, x in enumerate(inlines[1:]):
                try:
                    self.dat.append([':'.join(x[0:2])] + [x[i] for i in idxes])
                except IndexError:
                    badlines.append(str(idx + 2))
            self.is_snv = True
        else:
            for idx, x in enumerate(inlines[1:]):
                try:
                    self.dat.append([x[0]] + [x[i] for i in idxes])
                except IndexError:
                    badlines.append(str(idx + 2))
        if len(badlines):
            env.logger.warning(
                "The following lines are ignored due to having empty fields: {}. You may want to manually fill them up with placeholder 'NaN'."
                .format(', '.join(badlines)))
        return

    def m_qqdata(self):
        rdat = OrderedDict()
        for x, y in zip(self.fields, list(zip(*self.dat))):
            if x == 'GENE':
                rdat[x] = list(y)
            else:
                try:
                    rdat[x] = [
                        i if i == i else None
                        for i in list(map(float, list(y)))
                    ]
                except ValueError as e:
                    raise ValueError("{}".format(e))
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
                downloadFile('resource/refgene.pkl',
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
                            gd[item[0]] = [item[1], item[2]]
                    except:
                        sys.exit('ERROR: {0} format is not valid'.format(
                            self.gene_map))
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
                env.logger.warning(
                    'There are {1} genes belonging to more than one chromosomes. This might due to discrepancies between refgene database versions. Using whatever first entry (alphanumeric order) as the coordinate: {0}'
                    .format(', '.join(multi_chrom), len(multi_chrom)))
            if len(failed):
                env.logger.warning(
                    'There are {1} genes not found in local gene name database. You may want to provide your own list of genomic coordinates of genes: {0}'
                    .format(', '.join(failed), len(failed)))
        #
        rdat = OrderedDict()

        for x, y in zip(self.fields, list(zip(*dat))):
            if x == 'BP':
                rdat[x] = [
                    i if i == i else None for i in list(map(int, list(y)))
                ]
            elif x in ['GENE', 'CHR']:
                rdat[x] = list(y)
            else:
                try:
                    rdat[x] = [
                        i if i == i else None
                        for i in list(map(float, list(y)))
                    ]
                except ValueError as e:
                    raise ValueError("{}".format(e))
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
        '''.format(
            self.fname,
            'true' if self.a.same_page else 'false',
            'grayscale' if not self.a.color else self.a.color,
            'true' if self.a.fixed_shape or not self.a.same_page else 'false',
            self.a.shape,
            self.a.title,
            'true' if self.a.bonferroni else 'false',
            'false' if self.a.no_slope else 'true',
            ' ' if not self.a.hlines else ' '.join(
                list(map(str, self.a.hlines))),
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
        '''.format(
            self.fname,
            'true' if self.a.same_page else 'false',
            'grayscale' if not self.a.color else self.a.color,
            self.a.title,
            ' '.join(self.a.chrom),
            'true' if self.a.bonferroni else 'false',
            ' ' if not self.a.hlines else ' '.join(
                list(map(str, self.a.hlines))),
            self.a.font_size,
            self.a.label_top,
            ' ' if not self.a.label_these else ' '.join(self.a.label_these),
            'NULL'
            if self.a.chrom_prefix == "None" else repr(self.a.chrom_prefix),
            self.rdat,
            self.a.width_height[0] if self.a.width_height else 12,
            self.a.width_height[1] if self.a.width_height else 6,
        )
        return astr


QQ_FOO = '''
#! QQplot function
QQplot <- function(dat, multiple=T, color='default', shapeFixed=F, shapeValue=1, title="", gwLine=T, slopeLine=T, optLines=c(), topFont=3, labelTopGenes=0, annotate=NULL, index=2)
{
	if (labelTopGenes>0 || length(annotate)>0) label=T
	else label=F
	qqDat <- data.frame()
	qqTopHit <- data.frame()
    med <- vector()
	for( i in 2:ncol(dat) ) {
		one <- dat[, c(1,i)]
		names(one) <- c("gene", "obs")
		one$obs <- as.numeric(one$obs)
		one <- subset(na.omit(one[order(one$obs), ]), (obs>0 & obs<=1))
		one$logObs <- -log10(as.numeric(one$obs))
		one$ept <- -log10(c(1:length(one$obs))/(1+length(one$obs)))
		one$method <- colnames(dat)[i]
        lambda <- median(one$logObs)/median(one$ept)
        med[i-1] = paste("Genomic inflation factor for method '", one$method, "' is: ", lambda, sep='')
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
		theme(legend.key = element_rect(colour = 'white', fill = 'white', size = 0.5, linetype='blank')) +
		theme(legend.justification=c(1,0), legend.position=c(1,0)) +
		theme(legend.background = element_rect(fill="transparent", size=.5, linetype="blank")) +
		theme(legend.title=element_blank()) +
		scale_x_continuous(limits=c(0,max(qqDat$ept)), breaks=0:max(qqDat$ept), labels=0:max(qqDat$ept))
		if(label && nrow(qqTopHit)>0) pobj <- pobj + geom_text(data = qqTopHit, aes(ept, logObs, label = gene, colour=factor(method)), vjust = 1.5, size = as.numeric(topFont), show_guide = FALSE)
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
		pobj <- pobj + theme(legend.position="none") + scale_x_continuous(limits=c(0,max(qqDat$ept)), breaks=0:max(qqDat$ept), labels=0:max(qqDat$ept))
		if(label && nrow(qqTopHit)>0) {
            if (!is.na(mycols)) {
                pobj <- pobj + geom_text(data = qqTopHit, aes(ept, logObs, label = gene), colour=mycols, vjust = 1.5, size = as.numeric(topFont), show_guide = FALSE)
            } else {
                pobj <- pobj + geom_text(data = qqTopHit, aes(ept, logObs, label = gene), vjust = 1.5, size = as.numeric(topFont), show_guide = FALSE)
            }
        }
	}
	pobj <- pobj +
	ggtitle(paste(title,'\\n')) +
	theme(axis.text.x = element_text(angle=0, size=12)) +
	theme(axis.text.y = element_text(angle=0, size=12)) +
	theme(axis.title.y = element_text(face="bold", colour= "black", angle=90, size=13)) +
	theme(axis.title.x = element_text(face="bold", colour= "black", size=13)) +
	#theme(plot.title = element_text(size=16, face="bold")) +
	theme(plot.title = element_text(size=14)) +
    theme(panel.background = element_rect(fill = "white",colour = "black")) +
	theme(panel.grid.minor=element_blank()) +
	theme(plot.background = element_rect(fill = "transparent",colour = NA)) +
	#theme(strip.text.x = element_text(face = 'bold')) +
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
	return(list(plot=pobj, median=med))
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
        if (nrow(d.annotate)) plotObj <- plotObj + geom_text(data = d.annotate, aes(pos, logp, label = GENE), face='bold', colour='black', vjust = -1.0, size = as.numeric(topFont)+1, show_guide = FALSE)
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
        if (nrow(d.tops)) plotObj <- plotObj + geom_text(data = d.tops, aes(pos, logp, label = GENE, colour=factor(CHR)), vjust = -1.0, size = as.numeric(topFont), show_guide = FALSE)
    }
    agl = 90
    if (is.null(chrprefix)) agl = 0
    plotObj <- plotObj +
    theme(legend.position = "none") +
    ggtitle(paste(title, '\\n')) +
    theme(panel.background=element_rect(fill="white", colour="black")) +
    theme(panel.grid.minor=element_blank()) +
    theme(axis.text.x = element_text(angle=agl, size=11)) +
    theme(axis.text.y = element_text(angle=0, size=12)) +
    theme(axis.title.y = element_text(colour= "black", angle=90, size=13)) +
    theme(axis.title.x = element_text(colour= "black", size=13))
    #theme(axis.ticks=element_segment(colour=NA))

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
    print(pObj$plot)
    for (i in 1:length(pObj$median)) write(pObj$median[i], stdout())
} else {
    for( i in 2:ncol(dat) ) {
        pval <- dat[, c(1,i)]
        subtitle <- title
        if (ncol(dat)>2) subtitle <- paste(title, " ", colnames(dat)[i], " method", sep="")
        pObj <- QQplot(pval, multiple=F, color=color, shapeFixed=shapeFixed, shapeValue=shapeValue, title=subtitle, gwLine=gwLine, slopeLine=slopeLine, optLines=optLines, topFont=topFont, labelTopGenes=labelTopGenes, annotate=annotate, index=i)
        print(pObj$plot)
        for (i in 1:length(pObj$median)) write(pObj$median[i], stdout())
    }
}
graphics.off()
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
    graphics.off()
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
        graphics.off()
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
        graphics.off()
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
    graphics.off()
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
	graphics.off()
}
'''
HIST_FOO = '''
gghues <- function(n) {
        hues = seq(15, 375, length=n+1)
        hcl(h=hues, l=65, c=100)[1:n]
}
gghue <- function(n) { gghues(n)[n] }
gghist <- function(dat, stat_foo = NULL, vlines = NULL, xname='x') {
        average <- round(mean(dat[,xname], na.rm=T),4)
        stdev <- round(sd(dat[,xname], na.rm=T),4)
        #med <- round(median(dat[,xname], na.rm=T), 4)
        r1 <- round(min(dat[,xname], na.rm=T),4)
        r2 <- round(max(dat[,xname], na.rm=T),4)
        # convert dat obj from numeric to data.frame
        myplot <- ggplot(dat) +
                        aes_string(x = xname) +
                        geom_histogram(aes_string(y = '..density..', fill = '..density..'), color = 'white', binwidth = (r2-r1)/30) +
                        scale_fill_gradient('bin mass', low = 'darkolivegreen3', high = colors()[552]) +
                        geom_line(stat = 'density', size = 0.5, linetype = 2, color = 'grey50') +
                        geom_rug(color = 'grey80') +
                        scale_x_continuous(name = paste('\\n', xname)) +
                        scale_y_continuous(name = 'Density\\n') +
                        theme_bw()
        if (!is.null(vlines)) myplot <- myplot + geom_vline(xintercept = vlines, color = '#9ACD32', linetype = 2)
        if (!is.null(stat_foo)) {
        myplot <- myplot + stat_function(fun = stat_foo, color = 'red')
        plottitle <- 'Histogram & fitted vs. normal distribution density curve for "'
        } else {
        plottitle <- 'Histogram & fitted density curve for "'
        }
        myplot <- myplot + labs(title = paste(plottitle, xname, '"\\n', 'mean = ', toString(average), '; ', 'stdev = ', toString(stdev), '; ', 'range = ', '[', toString(r1), ',', toString(r2), ']', '\\n', sep=''))
        return(myplot)
}
'''

DOT_FOO = '''
ggdot <- function(dat, psize=2.5, color=NULL) {
    xyz = colnames(dat)
	#dat = as.data.frame(dat[order(dat[,1]),])
    if (length(xyz) >= 3) {
        if (!is.null(color)) dat[,3] = as.factor(dat[,3])
        myplot = ggplot(dat, aes_string(x=xyz[1], y=xyz[2], colour=xyz[3])) +
            xlab(paste('\\n', xyz[1], sep='')) + ylab(paste(xyz[2], '\\n', sep=''))
        if (!is.null(color)) {
            ncolors = length(unique(dat[,3]))
            mycolors = rep(brewer.pal(brewer.pal.info[color,]$maxcolors, name = color), ncolors)[1:ncolors]
            myplot = myplot + scale_colour_manual(values = mycolors)
        }
    } else if (length(xyz) == 2) {
        myplot = ggplot(dat, aes_string(x=xyz[1], y=xyz[2])) +
            xlab(paste('\\n', xyz[1], sep='')) + ylab(paste(xyz[2], '\\n', sep=''))
    } else {
        dat = cbind(seq(1:nrow(dat)), dat)
        colnames(dat) = c('index', xyz[1])
        myplot = ggplot(dat, aes_string(x='index', y=xyz[1])) +
            xlab(paste('\\n', 'index', sep='')) + ylab(paste(xyz[1], '\\n', sep=''))
    }
    myplot = myplot + geom_point(size = psize, binwidth = range(dat[,1])/30) +
	    theme(legend.text=element_text(size=8)) +
	    theme(legend.title=element_text(size=10)) +
	    theme(axis.title.x=element_text(size=10)) +
	    theme(axis.title.y=element_text(size=10)) +
	    theme_bw()
    return(myplot)
}
'''
BOX_FOO = '''
ggbox <- function(dat, psize=2, color=NULL) {
    tmp <- matrix(apply(dat[,c(2,3)], 2, unique), ncol=2)
    tmp <- tmp[order(tmp[,2]), ]
    levnames <- as.character(tmp[,1])
    dat[,3] <- as.factor(dat[,3])
    levels(dat[,3]) <- levnames
    names = colnames(dat)
    myplot = ggplot(dat, aes_string(x = names[3], y = names[1], fill = names[2])) +
            xlab(paste('\\n', names[2], sep='')) + ylab(paste(names[1], '\\n', sep=''))
    if (!is.null(color)) {
        ncolors = length(unique(dat[,2]))
        mycolors = rep(brewer.pal(brewer.pal.info[color,]$maxcolors, name = color), ncolors)[1:ncolors]
        myplot = myplot + scale_fill_manual(values = mycolors)
    }
    myplot = myplot + geom_boxplot(outlier.size = psize, varwidth = TRUE) +
	    theme(legend.text=element_text(size=8)) +
	    theme(legend.title=element_text(size=10)) +
	    theme(axis.title.x=element_text(size=10)) +
	    theme(axis.title.y=element_text(size=10)) +
	    theme_bw()
    return(myplot)
}
'''
