#!/usr/bin/env python3

# Plot nearest and second nearest neighbors found by test_delaunay_c.x and
# test_delaunay_f.x
#
# Felipe H. da Jornada (2019)

import numpy as np
import matplotlib.pyplot as plt

def plot_points(fname):
    with open(fname, 'r') as f:
        ndims, npts, npts_orig, nnei = map(int, f.readline().split())
        data = np.loadtxt(f)
    data = data[:,1:]
    pts = data[:npts]
    pts_orig = pts[:npts_orig]
    pts_nei = data[npts:npts+nnei]
    pt_test = data[-1]

    plt.scatter(pts[:,0], pts[:,1])
    plt.scatter(pts_orig[:,0], pts_orig[:,1])
    plt.scatter(pts_nei[:,0], pts_nei[:,1])
    plt.scatter(pt_test[0], pt_test[1], s=5)
    plt.gca().set_aspect('equal')
    plt.show()


if __name__ == '__main__':
    import sys
    fname = sys.argv[1]
    plot_points(fname)
