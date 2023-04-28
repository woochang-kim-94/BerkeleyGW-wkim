#!/usr/bin/env python3

# Script to generate file points.dat with test coordinates for
# test_delaunay_c.x and test_delaunay_f.x
#
# Felipe H. da Jornada (2019)

import numpy as np


def gen_grid(ndims, grid, bvec=None, plot=False, save=True, periodic=True):
    n1, n2, n3 = grid

    if bvec is None:
        bvec = np.eye(3)
    bvec = bvec[:ndims, :ndims]

    # Define point in "crystal" coordinates.
    # We will only convert to cartesian at the very end.
    pts_orig = (np.mgrid[:n1,:n2,:n3].reshape((3,-1))).T / np.array([n1,n2,n3])[None,:]
    pts_orig = pts_orig[:,:ndims]
    pts_orig = np.ascontiguousarray(pts_orig)
    npts_orig = pts_orig.shape[0]

    # Note: n_ghost=1 means only the original set of points
    n_ghost = 1
    if periodic:
        n_ghost = 3**ndims

    npts = npts_orig * n_ghost
    pts = np.empty((npts, ndims))

    ghost_G = np.zeros((ndims,), dtype=int)
    for ighost in range(n_ghost):
      offset = npts_orig*ighost
      # FHJ: the outermost mod()`s map (0, 1, 2) to (0, 1, -1)
      ghost_G[0] = (((ighost%3) + 1) % 3) - 1
      if ndims > 1:
        ghost_G[1] = ((((ighost//3) % 3) + 1) % 3) - 1
      if ndims > 2:
        ghost_G[2] = ((((ighost//9) % 3) + 1) % 3) - 1

      pts[offset:offset+npts_orig,:] = pts_orig + ghost_G[None,:]

    pts = np.dot(pts, bvec)
    pts_orig = np.dot(pts_orig, bvec)

    center = np.mean(pts_orig, axis=0)
    pts -= center[None,:]
    pts_orig -= center[None,:]

    if plot:
        import matplotlib.pyplot as plt
        plt.scatter(pts[:,0], pts[:,1])
        plt.scatter(pts_orig[:,0], pts_orig[:,1])
        plt.gca().set_aspect('equal')
        plt.show()

    np.savetxt('points.dat', pts, comments='', fmt='%.9f '*ndims,
               header='{} {} {}'.format(ndims, npts, npts_orig))


def gen_grid_2d(plot=False):
    args = dict(
        plot = plot,
        periodic = True,
        ndims = 2,
        grid = (4, 4, 1),
        bvec = np.array([
            [np.sqrt(3)/2, -1/2, 0],
            [np.sqrt(3)/2, 1/2, 0],
            [0, 0, 1]
            ], dtype = np.float64) * 2.
        #bvec = np.array([
        #    [1, 0, 0],
        #    [1/2, np.sqrt(3)/2, 0],
        #    [0, 0, 1],
        #    ], dtype = np.float64) * 2.
        )
    gen_grid(**args)


def gen_grid_3d(plot=False):

    # Simple cubic
    bvec = None
    # bcc
    a0 = np.array([1/2, 1/2, 1/2]) #cener
    a1 = np.array([1, 0, 0])
    a2 = np.array([0, 1, 0])
    a3 = np.array([0, 0, 1])
    bvec = np.row_stack((a1-a0, a2-a0, a3-a0))
    args = dict(
        plot = plot,
        periodic = True,
        ndims = 3,
        grid = (3, 3, 3),
        bvec = bvec,
        )
    gen_grid(**args)


if __name__ == '__main__':
    plot = False
    plot = True
    gen_grid_2d(plot)
    #gen_grid_3d(plot)
