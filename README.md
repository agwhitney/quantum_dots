# Finite Element Analysis

This is a snapshot of a project that I did during my senior year at UMass Dartmouth, with Dr. Jianyi Jay Wang as my advisor. The goal was to find or create a method in Python for finite element analysis, which would be then applied to solve equations for bound quantum systems called quantum dots. Essentially: find a mesh package, learn how to use it, and apply it to gather and analyze data.

I've recently come back to clean the code up a little bit. It used to be much, *much* rougher. Though it's low priority, there are still plans to improve the project - particularly in terms of flexibility and user-friendliness, with the overall intent of reducing having to dig through the code.

As is, running `main.py` will produce results for a unit circle - allowed energies up to 12 states, an image of the wave function, and the mesh used. The shapes in `shapes.py` should all produce results as is, as well.


## PyDistMesh

The mesh package used is called **pydistmesh** and is hosted [here.](https://github.com/bfroehle/pydistmesh) It can be installed via pip: `pip install pydistmesh` and requires numpy and matplotlib.

As far as I am aware the package does not compile on Windows. I worked in Ubuntu 14.04 and encountered no errors, but 17.04 is missing some packages, which are installed from terminal as below.

`sudo apt-get install libblas-dev liblapack-dev python3-tk`
