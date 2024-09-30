import numpy as np
from ase.io import read, write
from ase import Atoms

def read_xyz_file(filename):
    """Function to read an XYZ file and return an ASE Atoms object."""
    return read(filename)

def generate_random_orientations(molecule, num_orientations):
    """Function to generate different random orientations of the molecule."""
    molecules_oriented = []
    for _ in range(num_orientations):
        # Copy the molecule to avoid modifying the original object
        mol_copy = molecule.copy()
        
        # Generate random rotation angles for the x, y, z axes
        angle_x, angle_y, angle_z = np.random.uniform(0, 360, 3)

        # Apply rotations to the molecule
        mol_copy.rotate(angle_x, 'x', center='COM')
        mol_copy.rotate(angle_y, 'y', center='COM')
        mol_copy.rotate(angle_z, 'z', center='COM')
        
        molecules_oriented.append(mol_copy)
    
    return molecules_oriented

def position_molecule_near_surface(surface, molecule, distance, num_structures):
    """Positions the molecule at a specified distance from the substrate in different orientations."""
    oriented_molecules = generate_random_orientations(molecule, num_structures)
    structures = []
    
    for mol in oriented_molecules:
        # Calculate the center of mass (COM) of the surface and the molecule
        surface_com = surface.get_center_of_mass()
        mol_com = mol.get_center_of_mass()
        
        # Calculate the maximum height of the molecule relative to its COM
        mol_height = np.max(mol.get_positions()[:, 2]) - mol_com[2]

        # Define the new initial position for the molecule
        new_position = surface_com + np.array([0, 0, distance + mol_height])
        mol.translate(new_position - mol_com)
        
        # Check if the molecule is overlapping with the surface
        if is_overlapping(surface, mol):
            # Adjust the position to avoid overlap
            adjustment = 0.1  # Smaller adjustment
            while is_overlapping(surface, mol):
                new_position[2] += adjustment
                mol.translate(new_position - mol_com)

        # Create a new system containing the surface and the positioned molecule
        combined_system = surface + mol
        structures.append(combined_system)
    
    return structures

def is_overlapping(surface, molecule):
    """Checks if there is overlap between the surface and the molecule."""
    for atom in molecule:
        for surface_atom in surface:
            distance = np.linalg.norm(atom.position - surface_atom.position)
            if distance < 1.0:  # Overlap tolerance (1.0 Å)
                return True
    return False

def main():
    # Prompt the user for input files
    surface_file = input("Enter the path to the XYZ file of the substrate: ")
    molecule_file = input("Enter the path to the XYZ file of the molecule to be adsorbed: ")
    distance = float(input("Enter the desired adsorption distance (in Å): "))
    num_structures = int(input("Enter the number of structures you want to generate: "))

    # Read the XYZ files
    surface = read_xyz_file(surface_file)
    molecule = read_xyz_file(molecule_file)

    # Generate the structures
    structures = position_molecule_near_surface(surface, molecule, distance, num_structures)

    # Save the generated structures in XYZ files
    for i, structure in enumerate(structures):
        output_file = f'structure_{i + 1}.xyz'
        write(output_file, structure)
        print(f"Structure {i + 1} saved as {output_file}")

if __name__ == "__main__":
    main()

