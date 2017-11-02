#!/usr/bin/env python

""" MultiQC module to parse output from Salmon """

from __future__ import print_function
from collections import OrderedDict
import json
import logging
import os
from GCModel import GCModel

from multiqc import config
from multiqc.plots import linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Salmon', anchor='salmon',
        href='http://combine-lab.github.io/salmon/',
        info="is a tool for quantifying the expression of transcripts using RNA-seq data.")

        # Parse meta information. JSON win!
        self.salmon_meta = dict()
        self.salmon_bias_FirstSample = dict()
        self.salmon_bias_MiddleSample = dict()
        self.salmon_bias_LastSample = dict()
        for f in self.find_log_files('salmon/meta'):
            # Get the s_name from the parent directory
            s_name = os.path.basename( os.path.dirname(f['root']) )
            firstRatio= OrderedDict()
            middleRatio = OrderedDict()
            lastRatio = OrderedDict()
            s_name = self.clean_s_name(s_name, f['root'])
            if os.path.basename(os.path.dirname(f['root'])) == 'bias' :
                gc = GCModel()
                gc.from_file(os.path.dirname(f['root']))
                if gc.valid_ :
                    first_Row = gc.obs_[0] / gc.exp_[0]
                    middle_Row = gc.obs_[1] / gc.exp_[1]
                    last_Row = gc.obs_[2] / gc.exp_[2]
                    for i in range(len(first_Row)):
                        firstRatio[i]=first_Row[i]
                        middleRatio[i]=middle_Row[i]
                        lastRatio[i]=last_Row[i]
                    self.salmon_bias_FirstSample[s_name]=firstRatio
                    self.salmon_bias_MiddleSample[s_name]=middleRatio
                    self.salmon_bias_LastSample[s_name]=lastRatio
                else :
                    log.error('Failed to open {}'.format(os.path.dirname(f['root'])))

            #s_name = self.clean_s_name(s_name, f['root'])
            self.salmon_meta[s_name] = json.loads(f['f'])
        # Parse Fragment Length Distribution logs
        self.salmon_fld = dict()
        for f in self.find_log_files('salmon/fld'):
            # Get the s_name from the parent directory
            if os.path.basename(f['root']) == 'libParams':
                s_name = os.path.basename( os.path.dirname(f['root']) )
                s_name = self.clean_s_name(s_name, f['root'])
                parsed = OrderedDict()
                for i, v in enumerate(f['f'].split()):
                    parsed[i] = float(v)
                if len(parsed) > 0:
                    if s_name in self.salmon_fld:
                        log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                    self.add_data_source(f, s_name)
                    self.salmon_fld[s_name] = parsed

        # Parse Fragment Length Distribution logs

        # Filter to strip out ignored sample names
        self.salmon_meta = self.ignore_samples(self.salmon_meta)
        self.salmon_fld = self.ignore_samples(self.salmon_fld)

        if len(self.salmon_meta) == 0 and len(self.salmon_fld) == 0:
            raise UserWarning

        if len(self.salmon_meta) > 0:
            log.info("Found {} meta reports".format(len(self.salmon_meta)))
            log.info("hfffgfdfda")
            log.error("errrror")
            log.debug("debug")
            self.write_data_file(self.salmon_meta, 'multiqc_salmon')
        if len(self.salmon_fld) > 0:
            log.info("Found {} fragment length distributions".format(len(self.salmon_fld)))

        # Add alignment rate to the general stats table
        headers = OrderedDict()
        headers['percent_mapped'] = {
            'title': '% Aligned',
            'description': '% Mapped reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn'
        }
        headers['num_mapped'] = {
            'title': 'M Aligned',
            'description': 'Mapped reads (millions)',
            'min': 0,
            'scale': 'PuRd',
            'modify': lambda x: float(x) / 1000000,
            'shared_key': 'read_count'
        }
        self.general_stats_addcols(self.salmon_meta, headers)

        # Fragment length distribution plot
        pconfig = {
            'smooth_points': 500,
            'id': 'salmon_plot',
            'title': 'Salmon: Fragment Length Distribution',
            'ylab': 'Fraction',
            'xlab': 'Fragment Length (bp)',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section( plot = linegraph.plot(self.salmon_fld, pconfig) )

        pconfig2 = {
            'smooth_points': 500,
            'id': 'salmon_plot1',
            'title': 'Salmon : GC Bias Ratio in Beginning of Read',
            'ylab': 'Ratio',
            'xlab': 'GC Biases',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section(plot=linegraph.plot(self.salmon_bias_FirstSample, pconfig2))

        pconfig3 = {
            'smooth_points': 500,
            'id': 'salmon_plot2',
            'title': 'Salmon : GC Bias Ratio in Middle of Read',
            'ylab': 'Ratio',
            'xlab': 'GC Biases',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section(plot=linegraph.plot(self.salmon_bias_MiddleSample, pconfig3))

        pconfig4 = {
            'smooth_points': 500,
            'id': 'salmon_plot2',
            'title': 'Salmon : GC Bias Ratio in End of Read',
            'ylab': 'Ratio',
            'xlab': 'GC Biases',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section(plot=linegraph.plot(self.salmon_bias_LastSample, pconfig4))
