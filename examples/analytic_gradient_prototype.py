import numpy as np

def hellmann_feynman_gradient(ci_coeffs, dH_dR):
    # FIX: Initialize with dH_dR[0, 0] which is (3,3) instead of dH_dR[0] which is (2,3,3)
    E_grad = np.zeros_like(dH_dR[0, 0])
    for state in range(len(ci_coeffs[0])):
        c = ci_coeffs[:, state]
        grad = np.einsum('i,ijkl,j->kl', c, dH_dR, c)
        E_grad += grad
    return E_grad

if __name__ == "__main__":
    print("Hellmann-Feynman gradient prototype loaded.")
