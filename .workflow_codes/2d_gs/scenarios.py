#!/usr/bin/env python

import numpy as np

base_dic              = {}
train_points          = {}
test_points           = {}

sample_mesh_fractions = {}
leverage_scores_betas = {}

rom_energies          = {}
rom_basis_sets        = {}

'''
GalerkinFull
GalerkinCollocation
GalerkinGappy
GalerkinMaskedCollocation
GalerkinMaskedGappy
'''
rom_algos             = {}

'''
SplitVars or DontSplitVars
'''
rom_basis_styles      = {}

interpolation_sets    = {}
polyfit_sets          = {}

################################
################################

'''
------------------------------------------------------
SCENARIO 1
------------------------------------------------------
'''
base_dic[1] = {
  'general' : {
    'meshDir': "tbd",
  },

  'fom' : {
    'finalTime': 1000.,
    'odeScheme': "RungeKutta4",
    'dt' : 0.25,
    'stateSamplingFreq' : 5,
    'velocitySamplingFreq' : 5
  },

  'rom' : {
    'finalTime': 1000.,
    'odeScheme': "RungeKutta4",
    'dt' : 0.25,
    'stateSamplingFreq' : 5,
    # the following params are changed by driver scripts
    'numModes'            : 'tbd',
    'podFile'             : 'tbd',
    'algoName'            : 'tbd',
    'romInitialStateFile' : 'tbd',
    'referenceStateFile'  : 'tbd'
  },

  'physicalCoefficients' : {
    'diffusionA': 0.0002,
    'diffusionB': 0.00005,
    'feedRate'   : "tbd",
    'killRate'   : 0.062
  }
}

train_points[1] = {
  0: 0.042
}

test_points[1]  = {
  0: 0.042
}

sample_mesh_fractions[1] = [0.05]
leverage_scores_betas[1] = [0.9]

rom_algos[1]        = ["GalerkinFull", "GalerkinGappy", "GalerkinMaskedGappy"]
rom_basis_styles[1] = ["DontSplitVars"]
rom_energies[1]     = [99.999999]
rom_basis_sets[1]   = {
  0: [0]
}


'''
male list of valid scenarios, figure out from keys in dic
'''
valid_scenarios_ids = list(base_dic.keys())
