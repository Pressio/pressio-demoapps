
import numpy as np

def removeStepCells(meshObj, xBd, yBd):

  stepXb = [xBd[2], xBd[3]]
  stepYb = [yBd[1], yBd[2]]
  print("stepXb ", stepXb)
  print("stepYb ", stepYb)

  totNumCells = meshObj.getTotCells()
  [x, y] = meshObj.getCoordinates()
  gids = meshObj.getGIDs()
  G = meshObj.getGraph()

  a = np.where(x>stepXb[0])
  b = np.where(x<stepXb[1])
  c = np.where(y>stepYb[0])
  d = np.where(y<stepYb[1])
  r1 = np.intersect1d(a,b)
  r2 = np.intersect1d(r1,c)
  stepCellsGids = np.intersect1d(r2,d)
  domainGids = list(set(np.arange(totNumCells)).difference(set(stepCellsGids)))
  #print("stepCellsGids = ", stepCellsGids)
  #print("domainGids    = ", domainGids)

  G2 = {}
  for k,v in G.items():
    if k not in stepCellsGids:
      v2 = [i if i not in stepCellsGids else -1 for i in v]
      G2[k] = v2
  #printDicPretty(G2)

  # reindex
  newIds = {}
  count = 0
  for k in G2.keys():
    newIds[k] = count
    count+=1

  # fix graph
  G3 = {}
  for k,v in G2.items():
    v2 = [newIds[i] if i != -1 else -1 for i in v]
    G3[newIds[k]] = v2

  gids = np.array(list(G3.keys()))
  return x[domainGids], y[domainGids], G3, gids
