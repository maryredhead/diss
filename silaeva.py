import decimal
import numpy as np
import csv
import seaborn as sns
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# import plotly.graph_objects as go
from itertools import product
from scipy.integrate import odeint, ode, RK45

x_00, y_00, z_00 = 76291, 101997, 5365


# закон изменения для величин по системе модели города
# def lorenzsys(XYZ, t, a1, a2, c1, c2, d1, d2):
#
#    x, y, z = XYZ
#    x_dt = a1*y*z - a2 * x
#    y_dt = c1 * (y + z) - c2 * x * z
#    z_dt = d1 * y- d2 * z
#    return x_dt, y_dt, z_dt

# закон изменения для величин по системе модели города
def lorenzsys_0(XYZ, t, a1, a2, a3, c1, c2, c3, c4, d1, d2):
    # cust_koef.append((0.19, 1.44, 1.07, 1.66, 1.15, 0.79, 1.7, 2.19, 0.66))
    x, y, z = XYZ
    x_dt = a1 * (a2 * y - a3 * x)
    y_dt = c1 * (c2 * x - c3 * y) - c4 * x * z
    z_dt = d1 * x * y - d2 * z
    return x_dt, y_dt, z_dt


# управление по ренте постановка 1
# макроперемнная Y-ro*Z
def lorenzsys_upr_z1(XYZ, t, ro, T):
    x, y, z = XYZ
    a1, a2, a3, c1, c2, c3, c4, d1, d2 = 5.4, 0.8, 0.6, 1.6, 0.6, 0.3, 0.9, 1.2, 0.3
    psi = y - ro * z
    u = (c1 * (c2 * x - c3 * y) - c4 * x * z) / ro - d1 * x * y - d2 * z + psi / ro * T
    x_dt = a1 * (a2 * y - a3 * x)
    y_dt = c1 * (c2 * x - c3 * y) - c4 * x * z
    z_dt = d1 * x * y - d2 * z + u
    return x_dt, y_dt, z_dt


def lorenz_sys(koef):
    # a1, a2,  c1, c2,d1, d2 = koef[0], koef[1], koef[2], koef[3], koef[4], koef[5]
    a1, a2, a3, c1, c2, c3, c4, d1, d2 = koef[0], koef[1], koef[2], koef[3], koef[4], koef[5], koef[6], koef[7], koef[8]

    tmax, n = 10, 500
    t = np.linspace(0, tmax, 10)

    f = odeint(lorenzsys_0, (0.186, 0.268, 0.305), t, args=(a1, a2, a3, c1, c2, c3, c4, d1, d2))
    # f = ode.(lorenzsys, (x_00, y_00, z_00), t, f_args=(a1, a2, c1, c2, d1, d2))
    X, Y, Z = f.T
    return X, Y, Z

def lorenz_sys_upr(koef):
    ro, T = koef[0], koef[1]

    tmax, n = 10, 500
    t = np.linspace(0, tmax, 10)

    f = odeint(lorenzsys_upr_z1, (0.186, 0.268, 0.305), t, args=(ro,T))
    X, Y, Z = f.T
    return X, Y, Z

# статистика по годам - реальные данные
def get_statistic():
    # продукция
    # arr_x = [17895, 19658, 18958, 23468, 24955, 35999, 32980, 37075, 38096, 41876]
    arr_x = [0.186541495, 0.204919402, 0.197622445, 0.244635697, 0.260136519, 0.375261653, 0.343790919, 0.386478118,
             0.397121251, 0.436524819]
    # Население
    # arr_y = [78967, 83489, 83456, 88956, 87678, 95987, 94678, 103457, 103678, 104987]
    arr_y = [0.268614068, 0.283996098, 0.283883845, 0.30259264, 0.298245396, 0.326509282, 0.322056589, 0.351919227,
             0.35267098, 0.357123673]
    # рента
    arr_z = [0.305396529, 0.301089401, 0.304095839, 0.301552968, 0.310095214, 0.2749025, 0.325599975, 0.320599746,
             0.330528193, 0.333485124]
    # arr_z = [67856, 66899, 67567, 67002, 68900, 72765, 72345, 71234, 73440, 74097]
    return arr_x, arr_y, arr_z


# 5,  = 8,  = 1.9,  = 2.1,  = 3.16,  = 0.9
def main():
    x_y_z = (get_statistic())
    # lineplot(range(1,6), x_y_z[i], 't', 'ось значений', 'X')
    fig, axes = plt.subplots(3, 3)
    len_arr = len(x_y_z[0])
    for i in range(0, 3):
        axes[0, i].plot(range(0, len_arr), x_y_z[i])
    cust_koef = []

    # cust_koef.append((13.54, 10.92, 0.28, 0.72, 1.78, 1.58))

    cust_koef.append((1.18, 2.68, 1.28, 1.33, 0.7, 0.9, 0.77, 15.54, 8.92))
    cust_koef.append((1.08, 2.78, 1.28, 1.33, 0.7, 0.9, 0.77, 15.54, 6.92))
    cust_koef.append((1.08, 3.18, 1.28, 1.33, 0.6, 0.7, 0.57, 17.54, 9.92))
    cust_koef.append((1.11, 2.98, 1.11, 1.33, 0.6, 0.7, 0.97, 18.54, 8.92))
    cust_koef.append((1.11, 2.98, 1.11, 1.33, 0.6, 0.7, 0.97, 1.54, 0.92))
    cust_koef.append((1.11, 2.98, 1.11, 1.83, 0.6, 0.7, 0.97, 1.54, 0.92))
    cust_koef.append((1.07, 2.88, 1.33, 1.53, 0.8, 0.9, 0.79, 1.74, 0.97))
    # cust_koef.append((1.07, 2.88, 1.43, 1.13, 0.9, 0.9, 0.59, 11.94, 3.77))
    # cust_koef.append((1.07, 2.58, 1.73, 1.23, 0.9, 0.8, 0.69, 1.94, 0.83))
    # cust_koef.append((1.07, 2.58, 1.73, 1.23, 0.9, 0.8, 0.69, 22.94, 10.83))
    # cust_koef.append((1.05, 2.67, 1.13, 1.43, 0.77, 0.8, 0.69, 6.94, 2.83))

    # cust_koef.append((1.19, 2.01, 1.07, 1.38, 1.8, 0.7, 0.57, 17.34, 9.72))
    # cust_koef.append((1.19, 3.31, 1.27, 1.08, 1.8, 1.7, 0.97, 16.34, 9.72))
    # cust_koef.append((1.19, 2.31, 1.27, 1.77, 0.91, 0.67, 1.57, 1.94, 0.92))
    cust_koef.append((1.07, 2.88, 1.33, 1.53, 0.8, 0.9, 0.79, 1.74, 0.97))
    cust_koef.append((0.19, 1.44, 1.07, 1.66, 1.15, 0.79, 1.7, 2.19, 0.66))  # !!!!
    # cust_koef.append((1.12, 1.4, 1.35, 1.77, 1.1, 0.71, 1.55, 2.11, 0.71))
    arr_xyz_lorenz_cust = []
    arr_xyz_lorenz_upr = []
    for i in range(0, len(cust_koef)):
        arr_xyz_lorenz_cust.append(lorenz_sys(cust_koef[i]))
    arr_xyz_lorenz_upr.append(lorenz_sys([5.4, 0.8, 0.6, 1.6, 0.6, 0.3, 0.9, 1.2, 0.3]))
    # for j in range(0, len(arr_xyz_lorenz_cust)):
    #    for n in range(0, 3):
    #        axes[j+1, n].plot(range(0, 7), arr_xyz_lorenz_cust[j][n])

    # for n in range(0, 3):
    #    axes[2, n].plot(range(0, 7), arr_xyz_lorenz_cust[1][n])
    # for j in range(1, len(arr_xyz_lorenz_cust)):
    for i in range(0, 3):
        axes[0, i].plot(range(0, 10), arr_xyz_lorenz_upr[0][i])
        axes[1, i].plot(range(0, 10), arr_xyz_lorenz_upr[0][i])

    X, Y, Z = lorenz_sys_upr([0.6,5])
    axes[1, 0].plot(range(0, 10), X)

    axes[1, 1].plot(range(0, 10), Y)
    axes[1, 2].plot(range(0, 10), Z)
    #psi
    axes[2, 0].plot(range(0, 10), Y-0.6*Z)
    plt.show()


if __name__ == '__main__':
    main()
