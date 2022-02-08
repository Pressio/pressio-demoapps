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
    'finalTime': 5.0,
    'inviscidFluxReconstruction' : "Weno3",
    'odeScheme': "RungeKutta4",
    'dt' : 0.002,
    'stateSamplingFreq' : 5,
    'velocitySamplingFreq' : 5
  },

  'rom' : {
    'finalTime': 5.0,
    'inviscidFluxReconstruction' : "Weno3",
    'odeScheme': "RungeKutta4",
    'dt' : 0.002,
    'stateSamplingFreq' : 5,
    # the following params are here but left tbd
    # because these are changed by driver scripts
    'numModes'            : 'tbd',
    'podFile'             : 'tbd',
    'algoName'            : 'tbd',
    'romInitialStateFile' : 'tbd',
    'referenceStateFile'  : 'tbd'
  },

  'physicalCoefficients' : {
    'amplitude'   : "tbd"
  }
}

train_points[1] = {
  0: 0.0337,
  1: 0.05,
  2: 0.1
}

test_points[1]  = {
  0: 0.05,
  1: 0.08
}

sample_mesh_fractions[1] = [0.01, 0.05]
# leverage_scores_betas[1] = [0.9]

rom_algos[1]        = ["GalerkinGappy"]
rom_basis_styles[1] = ["DontSplitVars"]
rom_energies[1]     = [99.99]
rom_basis_sets[1]   = {
  0: [0, 1, 2]
}


'''
male list of valid scenarios, figure out from keys in dic
'''
valid_scenarios_ids = list(base_dic.keys())
