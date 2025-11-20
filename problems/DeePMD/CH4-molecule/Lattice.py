from ase import Atoms
import numpy as np

def create_fcc(kind : str, a : float, n1 : int, n2 : int, n3 : int):
    # Define the basis for an FCC structure: 4 atoms per unit cell
    fcc_basis = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
    
    # Generate positions for all atoms in the supercell
    positions = []
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                for basis_atom in fcc_basis:
                    position = np.array(basis_atom) + np.array([i, j, k])
                    positions.append(position)
    
    # Create the Atoms object with the specified positions and symbols (assuming a single element, e.g., Cu)
    atoms = Atoms(kind * len(positions), 
                  positions=np.array(positions) * a,  # Scale by lattice parameter 'a'
                  cell=[n1*a, n2*a, n3*a], 
                  pbc=True)
    
    return atoms

def create_bcc(kind : str, a : float, n1 : int, n2 : int, n3 : int):
    # Define the basis for a BCC structure: 2 atoms per unit cell
    bcc_basis = [(0, 0, 0), (0.5, 0.5, 0.5)]
    
    # Generate positions for all atoms in the supercell
    positions = []
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                for basis_atom in bcc_basis:
                    position = np.array(basis_atom) + np.array([i, j, k])
                    positions.append(position)
    
    # Create the Atoms object with the specified positions and symbols (assuming a single element, e.g., Fe)
    atoms = Atoms(kind * len(positions), 
                  positions=np.array(positions) * a,  # Scale by lattice parameter 'a'
                  cell=[n1*a, n2*a, n3*a], 
                  pbc=True)
    
    return atoms

def create_zincblende(kinds : list, a : float, n1 : int, n2 : int, n3 : int):
    """
    Creates a Zinc Blende (ZnS) structure.

    Parameters:
    a (float): Lattice parameter in Angstroms.
    n1, n2, n3 (int): Supercell dimensions. Default is 1x1x1 unit cell.

    Returns:
    Atoms object representing the ZnS structure.
    """
    fcc_basis = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
    # Define the basis for a Zinc Blende structure: 2 atoms per unit cell
    zincblende_basis = [(0, 0, 0), (0.25, 0.25, 0.25)]

    positions = []
    symbols   = []
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                for fcc_el in fcc_basis:
                    for zb_atom in zincblende_basis:
                        if zb_atom == (0, 0, 0):
                            symbols.append(kinds[0])
                        else:
                            symbols.append(kinds[1])
                        position = np.array(zb_atom) + np.array([i, j, k]) + np.array(fcc_el)
                        positions.append(position)
    
    atoms = Atoms(symbols, 
                  positions=np.array(positions) * a,
                  cell=[n1*a, n2*a, n3*a], 
                  pbc=True)
    
    return atoms

def create_cscl(kinds: list, a: float, n1: int, n2: int, n3: int):
    """
    Creates a CsCl (BCC) structure.

    Parameters:
    kinds (list): List of two elements to be used in the structure.
    a (float): Lattice parameter in Angstroms.
    n1, n2, n3 (int): Supercell dimensions. Default is 1x1x1 unit cell.

    Returns:
    Atoms object representing the CsCl structure.
    """
    cscl_basis = [(0, 0, 0), (0.5, 0.5, 0.5)]

    positions = []
    symbols = []
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                for basis_atom in cscl_basis:
                    if basis_atom == (0, 0, 0):
                        symbols.append(kinds[0])
                    else:
                        symbols.append(kinds[1])
                    position = np.array(basis_atom) + np.array([i, j, k])
                    positions.append(position)

    atoms = Atoms(symbols,
                  positions=np.array(positions) * a,
                  cell=[n1*a, n2*a, n3*a],
                  pbc=True)

    return atoms

if __name__ == "__main__":
    from ase.visualize import view

    a = 4.0  # Lattice parameter in Angstroms (e.g., for Fe which has a BCC structure with a=2.87 A)
    bcc_system = create_bcc('Fe', a, 3, 3, 3)
    #view(bcc_system)
    assert len(bcc_system) == 54, "Expected 54 atoms"
    assert np.allclose(bcc_system.get_cell(), [[3*a, 0, 0], [0, 3*a, 0], [0, 0, 3*a]]), "Cell dimensions do not match expectations"

    a = 4.0
    fcc_system = create_fcc('Cu', a, 3, 3, 3)
    assert len(fcc_system) == 108, "Expected 108 atoms"
    assert np.allclose(fcc_system.get_cell(), [[3*a, 0, 0], [0, 3*a, 0], [0, 0, 3*a]]), "Cell dimensions do not match expectations"
    #view(fcc_system)

    a = 5.42
    n1, n2, n3 = 3, 3, 3
    zincblende_system = create_zincblende(['S', 'Zn'], a, n1, n2, n3)
    # print(zincblende_system)
    # print(len(zincblende_system.positions))
    expected_atoms = 2 * 4 * (n1 * n2 * n3)
    assert len(zincblende_system) == expected_atoms, f"Expected {expected_atoms} atoms"
    assert np.allclose(zincblende_system.get_cell(), [[n1*a, 0, 0], [0, n2*a, 0], [0, 0, n3*a]]), "Cell dimensions do not match expectations"
    #view(zincblende_system)

    a = 4.12
    n1, n2, n3 = 3, 3, 3
    cscl_system = create_cscl(['Cs', 'Cl'], a, n1, n2, n3)
    print(cscl_system)
    expected_atoms = 2 * (n1 * n2 * n3)
    assert len(cscl_system) == expected_atoms, f"Expected {expected_atoms} atoms"
    assert np.allclose(cscl_system.get_cell(), [[n1*a, 0, 0], [0, n2*a, 0], [0, 0, n3*a]]), "Cell dimensions do not match expectations"
    view(cscl_system)