#!/usr/bin/env python

#==============================================================
# imports
#==============================================================

import copy

# LOCAL imports
from .scenarios import \
  base_dic, test_points, train_points
from .constants import \
  problemName, numDofsPerCell


#==============================================================
# functions
#==============================================================

def _replace_dic_field(scenario, dic, val):
  if scenario == 1:
    dic['physicalCoefficients']['feedRate'] = float(val)

  # !!!!!
  # if you add another scenario, add specifics here
  # !!!!!

  else:
    sys.exit("2d_gs: scenario = {} not supported yet".format(scenario))


def _create_fom_input_dics_impl(scenario, meshPath, kind):
  assert(meshPath != None)

  values = train_points[scenario] if kind=="train" else test_points[scenario]

  result = []
  for runNumber, v in values.items():
    # make a copy of the base dic
    thisDic = copy.deepcopy(base_dic[scenario])

    thisDic['runId'] = runNumber

    # we are making dics for the FOM, so remove the rom key
    if 'rom' in thisDic: del thisDic['rom']

    # always need to set mesh and problem
    thisDic['general']['meshDir'] = meshPath
    thisDic['general']['problem'] = problemName
    thisDic['general']['numDofsPerCell'] = numDofsPerCell

    _replace_dic_field(scenario, thisDic, v)
    result.append(thisDic)

  return result

#============================================================
def _create_rom_input_dics_impl(scenario, meshPath, algoName, \
                                numModes, lsvFile, refStateFile, \
                                projectorFile = None,\
                                # we need an optional sampleMeshPath because
                                # for masked galerkin the mesh to use is the FULL mesh
                                # but we still need to use gids from the sample mesh
                                sampleMeshPath = None):

  # we only care about test points here
  values = test_points[scenario]

  result = []
  for runNumber, v in values.items():
    # make a copy of the base dic
    thisDic = copy.deepcopy(base_dic[scenario])

    # the run ID MUST be the runNumber here because it has to
    # match the enumeration of the FOM
    thisDic['runId'] = runNumber

    # we are making dics for ROM, so remove the fom key
    if 'fom' in thisDic: del thisDic['fom']

    # always need to set mesh and problem
    thisDic['general']['meshDir'] = meshPath
    thisDic['general']['problem'] = problemName
    thisDic['general']['numDofsPerCell'] = numDofsPerCell

    thisDic['rom']['numModes'] = numModes
    thisDic['rom']['podFile']  = lsvFile
    thisDic['rom']['algoName'] = algoName
    thisDic['rom']['referenceStateFile'] = refStateFile

    if projectorFile != None:
      # if the projector is nontrivial
      thisDic['rom']['projectorFile'] = projectorFile

      if "Masked" not in algoName:
        # for real hyper-reduction, the mesh IS a sample mesh that contains gids
        thisDic['rom']['sampleMeshGidsFile']  = meshPath+"/sample_mesh_gids.dat"
        thisDic['rom']['stencilMeshGidsFile'] = meshPath+"/stencil_mesh_gids.dat"

      elif "Masked" in algoName:
        # if we are doing MASKED, then the mesh to use is a FULL mesh
        # so the GIDs must be taken from the sampleMeshPath argument
        assert( sampleMeshPath != None )
        # gappy Galerkin needs to know about gids files
        thisDic['rom']['sampleMeshGidsFile']  = sampleMeshPath+"/sample_mesh_gids.dat"
        thisDic['rom']['stencilMeshGidsFile'] = sampleMeshPath+"/stencil_mesh_gids.dat"

    _replace_dic_field(scenario, thisDic, v)
    result.append(thisDic)

  return result
