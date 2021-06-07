
import numpy as np
from pressiodemoapps import mesh as mesh
import pathlib, sys
file_path = pathlib.Path(__file__).parent.absolute()

def test1():
  meshPath = str(file_path)
  #meshO = mesh.CellCenteredUniformMesh(meshPath)
  meshO = mesh.loadCellCenterUniformMesh(meshPath)

  assert( meshO.dimensionality() == 2 )
  assert( meshO.stencilMeshSize() == 25 )
  assert( meshO.sampleMeshSize() == 25 )
  assert( meshO.stencilSize() == 3 )
  assert( meshO.graph().shape[0] == 25 )
  assert( meshO.graph().shape[1] == 5 )

  assert( meshO.dx() == 0.2 )
  assert( meshO.dxInv() == 5. )
  assert( meshO.dy() == 0.2 )
  assert( meshO.dyInv() == 5. )

  goldX = [0.1, 0.3, 0.5, 0.7, 0.9, 0.1, 0.3, 0.5,\
           0.7, 0.9, 0.1, 0.3, 0.5, 0.7, 0.9, 0.1,\
           0.3, 0.5, 0.7, 0.9, 0.1, 0.3, 0.5, 0.7, 0.9]
  goldY = [0.1, 0.1, 0.1, 0.1, 0.1,
           0.3, 0.3, 0.3, 0.3, 0.3,
           0.5, 0.5, 0.5, 0.5, 0.5,
           0.7, 0.7, 0.7, 0.7, 0.7,
           0.9, 0.9, 0.9, 0.9, 0.9]
  assert( np.allclose(meshO.viewX(), goldX) )
  assert( np.allclose(meshO.viewY(), goldY) )

def testGrahp():
  meshPath = str(file_path)
  meshO = mesh.CellCenteredUniformMesh(meshPath)
  G = meshO.graph()
  print(G)

  goldG = [
    [0,        4,        5,        1,       20],
    [1,        0,        6,        2,       21],
    [2,        1,        7,        3,       22],
    [3,        2,        8,        4,       23],
    [4,        3,        9,        0,       24],
    [5,        9,       10,        6,        0],
    [6,        5,       11,        7,        1],
    [7,        6,       12,        8,        2],
    [8,        7,       13,        9,        3],
    [9,        8,       14,        5,        4],
    [10,       14,       15,       11,        5],
    [11,       10,       16,       12,        6],
    [12,       11,       17,       13,        7],
    [13,       12,       18,       14,        8],
    [14,       13,       19,       10,        9],
    [15,       19,       20,       16,       10],
    [16,       15,       21,       17,       11],
    [17,       16,       22,       18,       12],
    [18,       17,       23,       19,       13],
    [19,       18,       24,       15,       14],
    [20,       24,        0,       21,       15],
    [21,       20,        1,       22,       16],
    [22,       21,        2,       23,       17],
    [23,       22,        3,       24,       18],
    [24,       23,        4,       20,       19]
  ]
  assert( np.allclose(G, goldG) )

if __name__ == '__main__':
  test1()
  testGraph()
