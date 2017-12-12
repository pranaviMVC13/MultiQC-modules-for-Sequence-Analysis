import pandas as pd
import numpy as np
import matplotlib.pyplot as py

class QuantModel:
    def __init__(self):
        self.ratio = None

    def populate_model_(self, data_):
        import struct
        from numpy.linalg import norm

        weights = None
        model = None
        offset = 0
        int_struct = struct.Struct('@i')
        long_struct = struct.Struct('@q')

        mspace = int_struct.unpack_from(data_[offset:])[0]
        offset += int_struct.size

        nrow = long_struct.unpack_from(data_[offset:])[0]
        offset += long_struct.size

        ncol = long_struct.unpack_from(data_[offset:])[0]
        offset += long_struct.size

        weight_struct = struct.Struct('@' + nrow * 'd')
        weights = weight_struct.unpack_from(data_[offset:])
        offset += weight_struct.size

        model_struct = struct.Struct('@' + nrow * ncol * 'd')
        model = model_struct.unpack_from(data_[offset:])
        model = np.array(model)
        model = model.reshape(ncol, nrow).T
        model = (model.T / model.sum(axis=1)).T
        return weights, model

    # dname is the root directory of salmon output
    def from_file(self, dname):
        import os
        import gzip
        quant_name = os.path.sep.join([dname, 'quant.sf'])
        df = pd.read_csv(quant_name, sep='\t', header=(0))
        self.ratio = (df['Length']/df['EffectiveLength']).tolist()
        return True
        """
        # Observed and Expected for GC Bias ############################
        obs_dat, exp_dat = None, None
        try:
            with gzip.open(obs_name) as obs_file:
                obs_dat = obs_file.read()
            self.obs_weights_, self.obs_ = self.populate_model_(obs_dat)
        except IOError:
            print("Could not open file {}".format(obs_name))
            return False

        try:
            with gzip.open(exp_name) as exp_file:
                exp_dat = exp_file.read()
            self.exp_weights_, self.exp_ = self.populate_model_(exp_dat)
        except IOError:
            print("Could not open file {}".format(exp_name))
            return False

        self.valid_ = True
        return True
        """
