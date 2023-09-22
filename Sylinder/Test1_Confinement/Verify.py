import numpy as np
import scipy as sp
import glob

files = glob.glob('./result/result*-*/SylinderAscii_*.dat')
radSy = 0.25
eps = 4e-4


def checkWall1(point):
    center = np.array([5, 5, 5])
    norm = np.array([1, 1, 0])
    norm_hat = norm / np.linalg.norm(norm)
    rvec = point - center
    error = rvec.dot(norm_hat) - radSy
    if error < -eps:
        print("Failed wall1", point, error)

    return


def checkWall2(point):
    center = np.array([0, 0, 20])
    norm = np.array([0, 0, -1])
    norm_hat = norm / np.linalg.norm(norm)
    rvec = point - center
    error = rvec.dot(norm_hat) - radSy
    if error < -eps:
        print("Failed wall2", point, error)

    return


def checkTube(point):
    center = np.array([12, 12, 12])
    radius = 8
    axis = np.array([1, 1, 1])
    axis_hat = axis / np.linalg.norm(axis)
    rvec = point - center
    radvec = rvec - rvec.dot(axis_hat) * axis_hat
    error = np.linalg.norm(radvec) - (radius - radSy)
    if error > eps:
        print("Failed tube", point, error)

    return


def checkSphere(point):
    center = np.array([12, 12, 12])
    radius = 10
    error = np.linalg.norm(point - center) - (radius - radSy)
    if error > eps:
        print("Failed sphere", point, error)
    return


def checkPoint(point):
    checkWall1(point)
    checkWall2(point)
    checkTube(point)
    checkSphere(point)


# for f in files:
#     print(f)
#     rods = np.loadtxt(f, skiprows=2, delimiter=' ',
#                       usecols=(3, 4, 5, 6, 7, 8))
#     for sy in rods:
#         checkPoint(sy[:3])
#         checkPoint(sy[3:])
