

#==============================================================
# functions
#==============================================================

def reconstruct_final_state_on_full_mesh(dryRun, workDir, romRunDir):
  from rom_reconstruction_and_error import _reconstruct_final_state_on_full_mesh_impl
  _reconstruct_final_state_on_full_mesh_impl(dryRun, workDir, romRunDir)

def compute_final_state_error(dryRun, fomTestDir, romDir):
  from rom_reconstruction_and_error import _compute_final_state_error_impl
  _compute_final_state_error_impl(dryRun, fomTestDir, romDir)
