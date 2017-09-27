#!/usr/bin/env python3
import distmesh as dm
import numpy as np
import fem
import shapes
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.sparse.linalg import eigsh


# Step 1: Get the nodes and elements from a given shape
# Which you just do by calling a shape from shapes.py
circle = shapes.circle(0.1)


# Step 2: Get the mesh data - nodes, elements, boundary and internal points
def get_meshdata(shape, savefile='meshdata.npz'):
    """An element is a triangular area with the three vertices called nodes.
    The shape is meshed and has nodes at coordinates (x, y). We need:
    Node Coordinates (coordinates, nx2 array); Elements (nodes, nx3 array);
    Boundary Nodes (nodes, 1d array); and Internal Nodes (nodes, 1d array)
    """
    node_coords, elements = shape
    boundary_nodes = np.unique(dm.boundedges(node_coords, elements))
    internal_nodes = np.setdiff1d(elements, boundary_nodes, assume_unique=True)

    if savefile:    # They could always set it to None!
        np.savez(savefile, node_coords=node_coords, elements=elements,
                boundary_nodes=boundary_nodes, internal_nodes=internal_nodes)

    return node_coords, elements, boundary_nodes, internal_nodes


# Step 3: Analyze the data from step 2
def analyze_data(shape, states=12, draw=True, nth_state=0,
                 meshfile=None, saveMeshFile=None, eigenfile=None, saveEigenFile=None
                 ):
    """Analyzes meshdata from a shape using the function get_meshdata(). Can use a provided
    meshfile or will run get_meshdata() to gather info, and can use/or create an eigenfile with
    with stored solutions.

    Gets the allowed energies for 'states' number of states, and if 'draw' is True will plot the
    mesh (twice, I think one comes from distmesh) and the wave function of the 'nth_state.'
    """
    if meshfile:
        print("Using meshdata from a numpy-packaged (.npz) mesh-file...")
        data = np.load(meshfile)
        nodes_coords = data['node_coords']
        elements = data['elements']
        boundary_nodes = data['boundary_nodes']
        internal_nodes = data['internal_nodes']     # Unused?
    else:
        print("Gathering mesh data anew...")
        nodes_coords, elements, boundary_nodes, internal_nodes = get_meshdata(shape, saveMeshFile)

    if eigenfile:
        solutions = np.load(eigenfile)
        energies = solutions['energies']
        u = solutions['u']  # I really don't know what u is about so for now it's staying as u
    else:
        # What are Tm and B? I don't know! Energy matrices, probably
        Tm = 0.5 * fem.matrix(nodes_coords, elements)
        Tm = np.delete(Tm, boundary_nodes, axis=0)  # delete boundary rows
        Tm = np.delete(Tm, boundary_nodes, axis=1)  # delete boundary columns

        B = fem.overlap_matrix(nodes_coords, elements)
        B = np.delete(B, boundary_nodes, axis=0)
        B = np.delete(B, boundary_nodes, axis=1)

        energies, u = eigsh(Tm, states, B, which='SA')  # Solve the eigen-equation
        if saveEigenFile:
            np.savez(saveEigenFile, energies=energies, u=u)

    print("Energies:\n", energies)

    if draw:
        nodes = np.asarray(nodes_coords)
        wave_function = u[:, nth_state]  # len(wave_function) = len(internal_points)
        for i in boundary_nodes:
            wave_function = np.insert(wave_function, i, 0)

        plt.figure()
        plt.subplot(111, aspect='equal')
        plt.triplot(nodes[:, 0], nodes[:, 1], elements, 'o-', linewidth=1)
        plt.xlabel('$x$', size=20)
        plt.ylabel('$y$', size=20)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_trisurf(nodes[:, 0], nodes[:, 1], wave_function, cmap=plt.cm.jet, linewidth=0.2)
        plt.axis('off')
        plt.title('n={}'.format(nth_state))

        plt.show()


if __name__ == '__main__':
    analyze_data(shape=circle, meshfile=None)
