#!/usr/bin/python3
"""A library of shapes to call on in main.py, each with
    a set variable for h0.
    distmesh2d vars: fd, fh, h0, bbox, pfix
"""
import distmesh as dm
import numpy as np


def show(shape):
    """Runs fstats function, presents the data arrays, shows the mesh"""
    p, t = shape
    print('{} nodes, {} elements, min quality {:.2}'.format(len(p), len(t), dm.simpqual(p, t).min()))
    plt.show(shape)


def circle(h0=0.05):
    """Uniform Mesh on Unit Circle"""
    fd = lambda p: np.sqrt((p ** 2).sum(1)) - 1.0
    return dm.distmesh2d(fd, dm.huniform, h0, (-1, -1, 1, 1))


def rectangleHole(h0=0.05):
    """Rectangle with circular hole, refined at circle boundary"""
    fd = lambda p: dm.ddiff(dm.drectangle0(p, -1, 1, -1, 1), dm.dcircle(p, 0, 0, 0.5))
    fh = lambda p: 0.05 + 0.3 * dm.dcircle(p, 0, 0, 0.5)
    return dm.distmesh2d(fd, fh, h0, (-1, -1, 1, 1))


def doughnut(h0=0.01):
    """Circle with circular hole, refined around it"""
    fd = lambda p: dm.ddiff(dm.dcircle(p, 0, 0, 1), dm.dcircle(p, 0, 0, 0.2))
    fh = lambda p: 0.01 + 0.3 * dm.dcircle(p, 0, 0, 0.2)
    return dm.distmesh2d(fd, fh, h0, (-1, -1, 1, 1))


def doughnut2(h0=0.01):
    """Circle with circular hole, not? refined around it"""
    fd = lambda p: dm.ddiff(dm.dcircle(p, 0, 0, 1), dm.dcircle(p, 0, 0, 0.2))
    # fh = lambda p: 0.01+0.3*dm.dcircle(p, 0,0, 0.2)
    return dm.distmesh2d(fd, dm.huniform, h0, (-1, -1, 1, 1))


def hexagon(h0=0.05):
    """A hexagon with a rotated hexagonal hole. d1 is the larger,
    d2 the smaller, and fd the difference."""
    n = 6
    outer = []
    for i in range(0, 2 * n + 1, 2):
        phi = (2 * np.pi * i) / (2 * n)
        outer.append((np.cos(phi), np.sin(phi)))
    inner = []
    for i in range(1, 2 * n + 2, 2):
        phi = (2 * np.pi * i) / (2 * n)
        inner.append((0.5 * np.cos(phi), 0.5 * np.sin(phi)))
    d1 = lambda p: dm.dpoly(p, outer)  # Larger Hexagon
    d2 = lambda p: dm.dpoly(p, inner)  # Smaller, Rotated Hexagon
    fd = lambda p: dm.ddiff(d1(p), d2(p))  # Difference Between them
    return dm.distmesh2d(d1, dm.huniform, h0, (-1, -1, 1, 1))  # Use d1, d2, or fd


def hexagonHole(h0=0.05):
    """A hexagon with a rotated hexagonal hole. d1 is the larger,
    d2 the smaller, and fd the difference."""
    n = 6
    outer = []
    for i in range(0, 2 * n + 1, 2):
        phi = (2 * np.pi * i) / (2 * n)
        outer.append((np.cos(phi), np.sin(phi)))
    inner = []
    for i in range(1, 2 * n + 2, 2):
        phi = (2 * np.pi * i) / (2 * n)
        inner.append((0.5 * np.cos(phi), 0.5 * np.sin(phi)))
    d1 = lambda p: dm.dpoly(p, outer)  # Larger Hexagon
    d2 = lambda p: dm.dpoly(p, inner)  # Smaller, Rotated Hexagon
    fd = lambda p: dm.ddiff(d1(p), d2(p))  # Difference Between them
    return dm.distmesh2d(fd, dm.huniform, h0, (-1, -1, 1, 1))  # Use d1, d2, or fd


def stadium(h0=0.05):
    """Two Circles + a Rectangle between them!"""
    d1 = lambda p: np.sqrt((p[:, 0] + 1) ** 2 + (p[:, 1]) ** 2) - 1.0
    d2 = lambda p: np.sqrt(((p[:, 0] - 1) ** 2) + (p[:, 1]) ** 2) - 1.0
    d3 = lambda p: dm.drectangle0(p, -1, 1, -1, 1)
    fd = lambda p: dm.dunion(d1(p), dm.dunion(d2(p), d3(p)))

    return dm.distmesh2d(fd, dm.huniform, h0, (-2, -2, 2, 2))


def halfstadium(h0=0.05):
    """Very similar to the full stadium"""
    d1 = lambda p: np.sqrt((p[:, 0] + 1) ** 2 + (p[:, 1]) ** 2) - 1.0
    d2 = lambda p: dm.drectangle0(p, -1, 1, -1, 1)
    fd = lambda p: dm.dunion(d1(p), d2(p))

    return dm.distmesh2d(fd, dm.huniform, h0, (-2, -2, 2, 2))


if __name__ == '__main__':
    """Test that a function returns correctly by plotting it with show()
    Imports matplotlib here rather than global, because it isn't needed otherwise
    """
    import matplotlib.pyplot as plt

    #show(circle(0.1))
