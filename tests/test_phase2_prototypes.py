import numpy as np
import pytest
from src.basis_change import compute_basis_changing_overlap
from examples.analytic_gradient_prototype import hellmann_feynman_gradient
from src.soc import compute_fragment_soc

def test_basis_change_overlap():
    """Test that basis-changing overlaps correctly zero out unmatched ESDs."""
    old_esds = [(1, 2), (3, 4)]
    new_esds = [(1, 2), (5, 6)] # (3,4) is lost, (5,6) is new
    old_psi = np.array([[1.0, 0.0], [0.0, 1.0]]) 
    new_psi = np.array([[1.0, 0.0], [0.0, 1.0]])
    
    S_ij = compute_basis_changing_overlap(old_esds, new_esds, old_psi, new_psi)
    
    assert S_ij[0, 0] == pytest.approx(1.0)
    assert S_ij[0, 1] == pytest.approx(0.0)
    assert S_ij[1, 0] == pytest.approx(0.0)
    assert S_ij[1, 1] == pytest.approx(0.0)

def test_hellmann_feynman_gradient():
    """Test the analytic gradient contraction."""
    ci_coeffs = np.array([[1.0, 0.0], [0.0, 1.0]])
    dH_dR = np.random.rand(2, 2, 3, 3) 
    
    grad = hellmann_feynman_gradient(ci_coeffs, dH_dR)
    assert grad.shape == (3, 3)
    assert np.allclose(grad, dH_dR[0, 0] + dH_dR[1, 1])

def test_soc_integrals():
    """Test SOC contraction logic using a mock molecule to avoid PySCF build issues."""
    # Mock the molecule and MO coefficients to bypass missing int1e_psoc in PySCF build
    class MockMol:
        nao = 2
        def intor(self, name):
            if name == 'int1e_psoc':
                # Return a fake 3-component SOC matrix
                return np.random.rand(3, 2, 2) + 1j * np.random.rand(3, 2, 2)
            return np.zeros((2, 2))
            
    mo_coeff = np.eye(2)
    
    soc_mo = compute_fragment_soc(MockMol(), mo_coeff)
    
    # Check that it properly contracted to (3, nao, nao)
    assert soc_mo.shape == (3, 2, 2)
