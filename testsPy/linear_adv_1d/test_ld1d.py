
import pathlib, sys
file_path = pathlib.Path(__file__).parent.absolute()

import numpy as np
from pressiodemoapps import mesh as mesh
from pressiodemoapps.advection import LinearAdvection1d as LA1d

def test_create_vel():
  meshPath = str(file_path)
  meshObj = mesh.CellCenteredUniformMesh(meshPath)
  appObj = LA1d(meshObj)
  v = appObj.createVelocity()
  print(v.shape)

if __name__ == '__main__':
  test_create_vel()
