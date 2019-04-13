# -*- coding: utf-8 -*-

import numpy as np
import pytest
from snf import cv

rs = np.random.RandomState(1234)


def test_compute_SNF(simdata):
    # don't define cluster number (find using compute.get_n_clusters)
    zaff, labels = cv.compute_SNF(simdata.data, metric='euclidean',
                                  n_perms=100)
    # define cluster number
    zaff, labels = cv.compute_SNF(simdata.data, metric='euclidean',
                                  n_clusters=3, n_perms=100)
    assert np.unique(labels).size == 3
    # provide list of cluster numbers
    zaff, labels = cv.compute_SNF(simdata.data, metric='euclidean',
                                  n_clusters=[3, 4], n_perms=100)
    assert isinstance(labels, list)
    for n, f in enumerate(labels, 3):
        assert np.unique(f).size == n


def test_snf_gridsearch(simdata):
    # only a few parameters to test
    zaff, labels = cv.snf_gridsearch(*simdata.data, metric='euclidean',
                                     mu=[0.35, 0.85], K=[10, 20],
                                     n_clusters=[2, 3], n_perms=100, seed=1234)
    # get optimal parameters based on diff corners
    for neighbors in ['edges', 'corners']:
        mu, K = cv.get_optimal_params(zaff, labels, neighbors=neighbors)


@pytest.mark.parametrize('x, y, edges, corners', [
    (0, 0, [0, 4, 1], [0, 4, 1, 5]),
    (1, 2, [5, 2, 6, 10, 7], [1, 5, 9, 2, 6, 10, 3, 7, 11]),
    (2, 1, [8, 5, 9, 13, 10], [4, 8, 12, 5, 9, 13, 6, 10, 14]),
    (3, 3, [14, 11, 15], [10, 14, 11, 15])
])
def test_neighbors(x, y, edges, corners):
    X = np.arange(4**2).reshape(4, 4)

    assert np.allclose(edges, X[cv.get_neighbors(x, y, 'edges', X.shape)])
    assert np.allclose(corners, X[cv.get_neighbors(x, y, 'corners', X.shape)])

    with pytest.raises(ValueError):
        cv.get_neighbors(x, y, 'badneighbors')


def test_extract_max_inds():
    X2d, X3d = np.arange(4**2).reshape(4, 4), np.arange(4**3).reshape(4, 4, 4)

    assert np.allclose([np.array([0, 1, 2, 3]), np.array([3, 3, 3, 3])],
                       cv.extract_max_inds(X2d))
    assert np.allclose(cv.extract_max_inds(X2d, axis=0),
                       cv.extract_max_inds(X2d, axis=-2))
    assert np.allclose([np.array([[0, 0, 0, 0], [1, 1, 1, 1], [2, 2, 2, 2],
                                  [3, 3, 3, 3]]),
                        np.array([[0, 1, 2, 3], [0, 1, 2, 3], [0, 1, 2, 3],
                                  [0, 1, 2, 3]]),
                        np.array([[3, 3, 3, 3], [3, 3, 3, 3], [3, 3, 3, 3],
                                  [3, 3, 3, 3]])],
                       cv.extract_max_inds(X3d))


def test_zrand():
    # make the same two-group community assignments (with different labels)
    label = np.ones((100, 1))
    X, Y = np.vstack((label, label * 2)), np.vstack((label * 2, label))
    # compare
    assert cv.zrand(X, Y) == cv.zrand(X, Y[::-1])
    assert cv.zrand(X, Y) > cv.zrand(X, rs.choice([0, 1], size=X.shape))
    assert cv.zrand(X, Y) == cv.zrand(X[:, 0], Y[:, 0])


def test_zrand_partitions():
    # make random communities
    comm = rs.choice(range(6), size=(10, 100))
    all_diff = cv.zrand_partitions(comm)
    all_same = cv.zrand_partitions(np.repeat(comm[:, [0]], 10, axis=1))

    # partition of labels that are all the same should have higher average
    # zrand and lower stdev zrand
    assert all_same[0] > all_diff[0]
    assert all_same[1] < all_diff[1]


def test_zrand_convolve():
    # random grid of community labels, testing against diff neighbor cases
    grid = rs.choice(range(6), size=(10, 10, 100))
    for neighbors in ['edges', 'corners']:
        cv.zrand_convolve(grid, neighbors=neighbors)
