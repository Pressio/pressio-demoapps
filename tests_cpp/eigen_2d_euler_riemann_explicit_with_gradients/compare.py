
import numpy as np
import sys, os

if __name__== "__main__":
  goldInit  = np.loadtxt("grad_gold_init.txt")
  goldFinal = np.loadtxt("grad_gold_final.txt")

  Dinit  = np.loadtxt("grad_result_init.txt")
  Dfinal = np.loadtxt("grad_result_final.txt")

  assert(np.allclose(goldInit.shape, Dinit.shape))
  assert(np.isnan(Dinit).all() == False)
  assert(np.allclose(goldInit, Dinit, rtol=1e-8, atol=1e-8))

  assert(np.allclose(goldFinal.shape, Dfinal.shape))
  assert(np.isnan(Dfinal).all() == False)
  assert(np.allclose(goldFinal, Dfinal, rtol=1e-8, atol=1e-8))
