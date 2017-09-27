"""J Wang's finite element method code
Edited by Adam Whitney (primarily increasing readability)
"""

from numpy import zeros


def abg(p1, p2, p3):
    """Returns alpha, beta, gamma, and area of a triangular element, which are related to the
    vertex points p1, p2, p3.
    """
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    alpha = [
        x2*y3 - x3*y2,
        x3*y1 - x1*y3,
        x1*y2 - x2*y1
    ]
    beta = [
        y2 - y3,
        y3 - y1,
        y1 - y2
    ]
    gamma = [
        x3 - x2,
        x1 - x3,
        x2 - x1
    ]
    area = 0.5 * sum(alpha)
    return alpha, beta, gamma, area


def overlap(i, j, p1, p2, p3):
    """Returns the integral of (phi_i * phi_j dx dy) over element triangle including endpoints"""
    a, b, g, area = abg(p1, p2, p3)     # alpha, beta, gamma, area
    X, Y, XY, X2, Y2 = 0, 0, 0, 0, 0
    for [x, y] in [p1, p2, p3]:
        X += x
        Y += y
        XY += x * y
        X2 += x * x
        Y2 += y * y
    var1 = (a[i]*b[j] + b[i]*a[j]) * X
    var2 = (a[i]*g[j] + g[i]*a[j]) * Y
    var3 = (b[i]*b[j]*(X2+(X*X)) + (b[i]*g[j] + g[i]*b[j])*(XY+(X*Y)) + (g[i]*g[j]*(Y2+(Y*Y)))) / 4
    var4 = 3 * a[i] * a[j]
    return (var1 + var2 + var3 + var4) / (12 * area)


def matrix(nodes, elements):
    """Fills the matrix from: \int \nabla \phi_i \cdot \nabla \phi_j dx dy
    The matrix value = 2 * KE of the Quantum Dot?
    """
    matrix = zeros((len(nodes), len(nodes)))    # Initialize square matrix the size of the # of nodes
    for elem in elements:
        x1, y1 = nodes[elem[0]]
        x2, y2 = nodes[elem[1]]
        x3, y3 = nodes[elem[2]]
        area = 2 * (x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))  # Element area?
        if area <= 0:
            # TODO I think if this issue comes up it will crash anyway, but here it is for now
            print("Zero or negative element area")
        beta = [
            y2 - y3,
            y3 - y1,
            y1 - y2
        ]
        gamma = [
            x3 - x2,
            x1 - x3,
            x2 - x1
        ]
        for i in range(3):
            for j in range(i, 3):   # Account for symmetry around diagonal
                matrix[elem[i], elem[j]] += (beta[i]*beta[j] + gamma[i]*gamma[j]) / area
                if i != j:  # Symmetry around diagonal
                    matrix[elem[j], elem[i]] = matrix[elem[i], elem[j]]
    return matrix


def overlap_matrix(nodes, elements):
    """Overlap matrix, integral of phi_i phi_j dx dy"""
    matrix = zeros((len(nodes), len(nodes)))
    for elem in elements:
        p1 = nodes[elem[0]]
        p2 = nodes[elem[1]]
        p3 = nodes[elem[2]]
        for i in range(3):
            for j in range(i, 3):
                matrix[elem[i], elem[j]] += overlap(i, j, p1, p2, p3)
                if i != j:
                    matrix[elem[j], elem[i]] = matrix[elem[i], elem[j]]
    return matrix
