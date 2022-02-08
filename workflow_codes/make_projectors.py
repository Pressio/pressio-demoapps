
import numpy as np
from scipy import linalg

from pod_functions import \
  _load_basis_from_binary_file

def _make_gappy_projector_dont_split_vars_impl(outputFile, sampleMeshPath, \
                                               statePodDir, rhsPodDir, \
                                               K, numDofsPerCell):

  smIndices = np.loadtxt(sampleMeshPath+"/sample_mesh_gids.dat", dtype=np.int64)
  smCount   = len(smIndices)

  # load phi
  fullPhi = _load_basis_from_binary_file(statePodDir +"/lsv")[:, 0:K]

  # load theta
  fullTheta = _load_basis_from_binary_file(rhsPodDir +"/lsv")[:, 0:K]
  theta     = np.zeros((smCount*numDofsPerCell, K), order='F')
  for k in range(numDofsPerCell):
    theta[k::numDofsPerCell, :] = fullTheta[numDofsPerCell*smIndices + k, :]

  A = fullPhi.transpose() @ fullTheta
  projector = A @ linalg.pinv(theta)
  print("SHAPE = ", projector.shape)

  # write to file
  # here when writing w need to consider that project above is computed
  # such that it is short and wide, so to write it to file we need to
  # "view" it as flipped. the actual num rows is the cols and vice versa.
  numRows = np.int64(projector.shape[1])
  numCols = np.int64(projector.shape[0])
  fileo = open(outputFile, "wb")
  np.array([numRows]).tofile(fileo)
  np.array([numCols]).tofile(fileo)
  '''
  NOTE: tofile write an array in rowwise fashion, REGARDLESS of the layout of the matrix.
  so here because projector is already transposed, we are good
  '''
  projector.tofile(fileo)
  fileo.close()

# #=============================================================================
# def make_collocation_projector_dont_split_vars(outputFile, sampleMeshPath, \
#                                                statePodDir, K, numDofsPerCell):

#   smIndices = np.loadtxt(sampleMeshPath+"/sample_mesh_gids.dat", dtype=np.int64)
#   smCount   = len(smIndices)

#   # load phi
#   fullPhi = load_basis_from_binary_file(statePodDir +"/lsv")[:, 0:K]

#   myProjector = None
#   if numDofsPerCell==1:
#     myProjector = fullPhi[smIndices, :]
#   else:
#     N = len(smIndices)*numDofsPerCell
#     myProjector = np.zeros((N, K), order='F')
#     for k in range(numDofsPerCell):
#       myProjector[k::numDofsPerCell, :] = fullPhi[numDofsPerCell*smIndices + k, :]

#   # write to file
#   numRows = np.int64(myProjector.shape[0])
#   numCols = np.int64(myProjector.shape[1])
#   fileo = open(outputFile, "wb")
#   np.array([numRows]).tofile(fileo)
#   np.array([numCols]).tofile(fileo)
#   '''
#   NOTE: tofile write an array in rowwise fashion, REGARDLESS of the layout of the matrix.
#   So here we need to pass .T to tofile so that tofile writes to file
#   the project "column by column" because this is required format for reading in cpp code
#   '''
#   np.transpose(myProjector).tofile(fileo)
#   fileo.close()
