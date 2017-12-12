#!/usr/bin/env python

""" MultiQC module to parse output from Salmon """

from __future__ import print_function
from collections import OrderedDict
import json
import logging
import os
from GCModel import GCModel
from SeqModel import SeqModel
from QuantReader import QuantModel

from multiqc import config
from multiqc.plots import linegraph
from multiqc.plots import heatmap,bargraph
from multiqc.modules.base_module import BaseMultiqcModule
import numpy as np

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
        # Declaring dicts to hold sample names and ratios for first,midddle and last rows without weights and with weights
        self.salmon_bias_FirstSample = dict()
        self.salmon_bias_MiddleSample = dict()
        self.salmon_bias_LastSample = dict()
        self.salmon_bias_FirstSampleWeights = dict()
        self.salmon_bias_MiddleSampleWeights = dict()
        self.salmon_bias_LastSampleWights = dict()
        self.salmon_bias_Average=dict()
        self.salmon_seq3 = dict()
        self.salmon_seq5 = dict()
        self.salmon_seqAverage=dict()
        self.salmon_quant =dict()
        # Dictionary for sample wise averages for all bias ratios
        self.salmon_bias_TotalAverage = dict()
        # Declaring lists to hold avergaes of first , middle and last row
        self.heatmapFirstrow=[]
        self.heatMapMiddleRow=[]
        self.heatMapLastRow=[]
        self.averageBiasHeatMap=[]
        self.seq3HeatMap=[]
        self.seq5HeatMap=[]
        self.sample_names=[]
        heatMapListFirst=[]
        heatMapListMiddle=[]
        heatMapListLast=[]
        count = 0
        for f in self.find_log_files('salmon/meta'):
            # Get the s_name from the parent directory
            s_name = os.path.basename( os.path.dirname(f['root']) )
            #Dicts for every sample for all the bucket(25) ratios
            firstRatio= OrderedDict()
            middleRatio = OrderedDict()
            lastRatio = OrderedDict()
            firstRatioWeight = OrderedDict()
            middleRatioWeight = OrderedDict()
            lastRatioWeight = OrderedDict()
            average = OrderedDict()
            sampleAverage = OrderedDict()
            seq3 = OrderedDict()
            seq5 = OrderedDict()
            seq_Average = OrderedDict()
            quant_Dict = OrderedDict()
            s_name = self.clean_s_name(s_name, f['root'])
            # Check only for bias folder
            if os.path.basename(os.path.dirname(f['root'])) == 'bias':
                gc = GCModel() # Instantiate GCModel class
                s_name_trimmed = s_name.partition('|')[0].split()
                self.sample_names.append(s_name_trimmed)
                # Call the GCModel method to get all observed and expected values
                gc.from_file(os.path.dirname(f['root']))
                # Calculate the ratio of all rows for observed by expected
                seq = SeqModel()
                seq.from_file(os.path.dirname(f['root']))
                seq3_prob = seq.obs3_/seq.exp3_
                seq5_prob = seq.obs5_/seq.exp5_

                quant = QuantModel()
                quant.from_file(os.path.dirname(f['root']))
                quant_ratio = quant.ratio

                first_Row = gc.obs_[0] / gc.exp_[0]
                firstRowWeightsRatio=gc.obs_weights_[0]/gc.exp_weights_[0]
                middle_Row = gc.obs_[1] / gc.exp_[1]
                middleRowWeightsRatio = gc.obs_weights_[1]/gc.exp_weights_[1]
                last_Row = gc.obs_[2] / gc.exp_[2]
                lastRowWeightsRatio = gc.obs_weights_[2]/gc.exp_weights_[2]
                # Avergaing all the ratios for the entire sample
                totalSampleAverage = ( (sum(first_Row)+sum(middle_Row)+sum(last_Row))/(len(first_Row)+len(middle_Row)+len(last_Row)))
                sampleAverage[count]=totalSampleAverage
                count = count +1
                self.salmon_bias_TotalAverage[s_name_trimmed[0]] = sampleAverage
                #Avergaing ratios for each row
                heatMapListFirst.append(sum(first_Row)/len(first_Row))
                heatMapListMiddle.append(sum(middle_Row) / len(middle_Row))
                heatMapListLast.append(sum(last_Row) / len(last_Row))
                self.heatmapFirstrow.append(first_Row.tolist())
                self.heatMapMiddleRow.append(middle_Row.tolist())
                self.heatMapLastRow.append(last_Row.tolist())
                self.seq3HeatMap.append(seq3_prob.tolist())
                self.seq5HeatMap.append(seq5_prob.tolist())
                heatmapAverage=[]
                # Iterating over all the buckets to create Ordered Dicts
                for i in range(len(first_Row)):
                    index = i*(100/len(first_Row))
                    firstRatio[i]=first_Row[i]
                    middleRatio[i]=middle_Row[i]
                    lastRatio[i]=last_Row[i]
                    firstRowRatios = firstRowWeightsRatio*first_Row[i]
                    middleRowRatios = middleRowWeightsRatio * middle_Row[i]
                    lastRowRatios = lastRowWeightsRatio * last_Row[i]
                    firstRatioWeight[index] = firstRowRatios
                    middleRatioWeight[index] = middleRowRatios
                    lastRatioWeight[index] = lastRowRatios
                    average[index]=np.mean([firstRowRatios,middleRowRatios,lastRowRatios])
                    heatmapAverage.append(average[index])

                self.averageBiasHeatMap.append(heatmapAverage)
                for i in range(len(seq3_prob)):
                    index=i*(100/len(seq3_prob))
                    seq3[index]=seq3_prob[i]
                    seq5[index]=seq5_prob[i]
                    seq_Average[i]=np.mean([seq3_prob[i],seq5_prob[i]])


                for i in range(len(quant_ratio)):
                    quant_Dict[i] = quant_ratio[i]

                # Setting all the ordered dicts to the outermost Dictionaries
                self.salmon_bias_FirstSample[s_name]=firstRatio
                self.salmon_bias_MiddleSample[s_name]=middleRatio
                self.salmon_bias_LastSample[s_name]=lastRatio
                self.salmon_bias_FirstSampleWeights[s_name]=firstRatioWeight
                self.salmon_bias_MiddleSampleWeights[s_name] = middleRatioWeight
                self.salmon_bias_LastSampleWights[s_name] = lastRatioWeight
                self.salmon_bias_Average[s_name] = average
                self.salmon_seq3[s_name]=seq3
                self.salmon_seq5[s_name]=seq5
                self.salmon_seqAverage[s_name]=seq_Average
                self.salmon_quant[s_name] = quant_Dict

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


        pconfig1 = {
            'smooth_points': 500,
            'id': 'salmon_plot1',
            'title': 'Salmon : GC Bias Ratio in Beginning of Read',
            'ylab': 'Ratio',
            'xlab': 'GC Biases',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        #self.add_section(plot=linegraph.plot(self.salmon_bias_FirstSample, pconfig1))

        pconfig2 = {
            'smooth_points': 500,
            'id': 'salmon_plot2',
            'title': 'Salmon : GC Bias Ratio in Middle of Read',
            'ylab': 'Ratio',
            'xlab': 'GC Biases',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        #self.add_section(plot=linegraph.plot(self.salmon_bias_MiddleSample, pconfig2))

        pconfig3 = {
            'smooth_points': 500,
            'id': 'salmon_plot3',
            'title': 'Salmon : GC Bias Ratio in Last of Read',
            'ylab': 'Ratio',
            'xlab': 'GC Biases',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        #self.add_section(plot=linegraph.plot(self.salmon_bias_LastSample, pconfig3))

        pconfig4 = {
            'smooth_points': 500,
            'id': 'salmon_plot4',
            'title': 'Salmon : GC Bias Ratio in Beginning of Read With Weights',
            'ylab': 'Ratio',
            'xlab': 'GC Biases',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section(plot=linegraph.plot(self.salmon_bias_FirstSampleWeights, pconfig4))

        pconfig5 = {
            'smooth_points': 500,
            'id': 'salmon_plot5',
            'title': 'Salmon : GC Bias Ratio in Middle of Read With Weights',
            'ylab': 'Ratio',
            'xlab': 'GC Biases',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section(plot=linegraph.plot(self.salmon_bias_MiddleSampleWeights, pconfig5))

        pconfig6 = {
            'smooth_points': 500,
            'id': 'salmon_plot6',
            'title': 'Salmon : GC Bias Ratio in Last of Read With Weights',
            'ylab': 'Ratio',
            'xlab': 'GC Biases',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section(plot=linegraph.plot(self.salmon_bias_LastSampleWights, pconfig6))

        pconfig7 = {
            'smooth_points': 500,
            'id': 'salmon_plot6',
            'title': 'Salmon : Average GC Bias of all samples',
            'ylab': 'Ratio',
            'xlab': 'Bias',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section(plot=linegraph.plot(self.salmon_bias_Average, pconfig7))

        pconfig7 = {
            'smooth_points': 500,
            'id': 'salmon_plot6',
            'title': 'Salmon : Seq 3',
            'ylab': 'Ratio',
            'xlab': 'Seq ',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section(plot=linegraph.plot(self.salmon_seq3, pconfig7))
        pconfig8 = {
            'smooth_points': 500,
            'id': 'salmon_plot6',
            'title': 'Salmon : Seq 5',
            'ylab': 'Ratio',
            'xlab': 'Seq ',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section(plot=linegraph.plot(self.salmon_seq5, pconfig8))

        pconfig8 = {
            'smooth_points': 500,
            'id': 'salmon_plot6',
            'title': 'Salmon : Seq Average',
            'ylab': 'Ratio',
            'xlab': 'Seq ',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section(plot=linegraph.plot(self.salmon_seqAverage, pconfig8))

        pconfig7 = {
            'smooth_points': 500,
            'id': 'salmon_plot7',
            'title': 'Salmon : GC Bias Average for every samples',
            'ylab': 'Ratios',
            'xlab': 'Samples',
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section(plot=bargraph.plot(self.salmon_bias_TotalAverage, pconfig=pconfig7))

        pconfig7 = {
            'smooth_points': 500,
            'id': 'salmon_plot7',
            'title': 'Salmon : Quant plot',
            'ylab': 'Ratios',
            'xlab': 'Samples',
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section(plot=linegraph.plot(self.salmon_quant, pconfig7))
        """
        # Calculating Correlatioon Coefficients for all the rows across all samples
        FirstRowCoff = np.corrcoef(self.heatmapFirstrow)
        self.add_section(name='Sample Similarity', description='Heatmap to display variance between first row ratios of all the samples',
                         plot=heatmap.plot(FirstRowCoff, self.sample_names, self.sample_names))

        MiddleRowCoff = np.corrcoef(self.heatMapMiddleRow)
        self.add_section(name='Sample Similarity', description='Heatmap to display variance between middle row ratios of all the samples',
                         plot=heatmap.plot(MiddleRowCoff, self.sample_names, self.sample_names))

        LastRowCoff = np.corrcoef(self.heatMapLastRow)
        self.add_section(name='Sample Similarity', description='Heatmap to display variance between last row ratios of all the samples',
                         plot=heatmap.plot(LastRowCoff, self.sample_names, self.sample_names))

        # Calculating correlation coefficients between the different rows
        heatMapData = []
        heatMapData.append(heatMapListFirst)
        heatMapData.append(heatMapListMiddle)
        heatMapData.append(heatMapListLast)
        corr = np.corrcoef(heatMapData)
        x=['FirstRow','MiddleRow','LastRow']
        self.add_section(name='Sample Similarity',description='Heatmap to display variance between the first , middle and last row of all the samples',plot=heatmap.plot(corr,x,pconfig3))

        AverageCoff = np.corrcoef(self.averageBiasHeatMap)
        self.add_section(name='Sample Similarity', description='Heatmap to display average bias across all samples',
                         plot=heatmap.plot(AverageCoff, self.sample_names, self.sample_names))


        Seq3HeatMap = np.corrcoef(self.seq3HeatMap)
        self.add_section(name='Sample Similarity', description='Heatmap to display Sequence 3 prime across all samples',
                         plot=heatmap.plot(Seq3HeatMap, self.sample_names, self.sample_names))

        Seq5HeatMap = np.corrcoef(self.seq5HeatMap)
        self.add_section(name='Sample Similarity', description='Heatmap to display Sequence 5 prime across all samples',
                         plot=heatmap.plot(Seq5HeatMap, self.sample_names, self.sample_names))"""
        
