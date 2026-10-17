import numpy as np

def compute_basis_changing_overlap(old_esds, new_esds, old_psi, new_psi):
    """
    Prototype for Phase 2: Handle overlaps when the ESD basis changes.
    If a determinant exists in both steps, use its overlap. If it vanishes, overlap is 0.
    """
    n_old = len(old_esds)
    n_new = len(new_esds)
    
    old_dict = {tuple(esd): i for i, esd in enumerate(old_esds)}
    new_dict = {tuple(esd): i for i, esd in enumerate(new_esds)}
    
    O_det = np.zeros((n_old, n_new), dtype=complex)
    for esd_str, i in old_dict.items():
        if esd_str in new_dict:
            j = new_dict[esd_str]
            O_det[i, j] = 1.0 # Dummy overlap for basis-change prototype
            
    S_ij = old_psi.T @ O_det @ new_psi
    return S_ij
