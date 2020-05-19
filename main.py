import decimal
import numpy as np
import csv
import seaborn as sns
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# import plotly.graph_objects as go
from itertools import product
from scipy.integrate import odeint,ode,RK45


def float_range(start, stop, step):
    while start < stop:
        yield float(start)
        start += float(decimal.Decimal(step))


# создаем массив всех возможных наборов коэффициентов в заданном диапазоне для каждого
# с заданным шагом изменения внутри диапазона
#  = 5,   = 8,   = 1.9,   = 2.1,   = 3.16,   = 0.9
#   = 13.54,  = 0.72,  = 10.92,   = 0.28,  = 1.58,  = 1.78
def make_array():
    a1 = float_range(0.9, 1.3, 0.1)
    a2 = float_range(0.6, 3, 0.2)
    a3 = float_range(2.4, 3.9, 0.1)

    c1 = float_range(2.1, 2.5, 0.4)
    c2 = float_range(0.9, 1.1, 0.1)
    c3 = float_range(0.9, 1.1, 0.1)
    c4 = float_range(7, 8, 0.5)

    d1 = float_range(5, 6, 0.5)
    d2 = float_range(0.8, 2, 0.4)
    arr = []

    lists = [a1, a2, a3, c1, c2, c3, c4, d1, d2]
    for items in product(*lists):
        yield items
        #arr.append(items)
    #return arr


# статистика по годам - реальные данные
def get_statistic():
    # продукция
    #arr_x = [17895, 19658, 18958, 23468, 24955, 35999, 32980, 37075, 38096, 41876]
    arr_x = [0.615032994, 0.675625516, 0.651567226, 0.80657135, 0.857678031, 1.237249106, 1.133489139, 1.274230135, 1.309320869, 1.439235634]
    # Население
    #arr_y = [78967, 83489, 83456, 88956, 87678, 95987, 94678, 103457, 103678, 104987]
    arr_y = [0.853390077, 0.902258971, 0.901902342, 0.961340404, 0.94752916, 1.037323861, 1.023177602, 1.118051555, 1.120439885, 1.134586144]
    # рента
    arr_z = [0.96646513, 0.95283469, 0.962348936, 0.954301707, 0.981334701, 1.036383447, 1.030401436, 1.014577592, 1.045997394, 1.055354968]
    #arr_z = [67856, 66899, 67567, 67002, 68900, 72765, 72345, 71234, 73440, 74097]
    return arr_x, arr_y, arr_z


#x_00 = 17895
#y_00 = 78967
#z_00 = 67856
x_00 = 0.615032994
y_00 = 0.853390077
z_00 =0.96646513
t_len = 10


# закон изменения для величин по системе модели города
def lorenzsys(XYZ, t, a1, a2, a3, c1, c2, c3, c4, d1, d2):
    x, y, z = XYZ
    x_dt = a1 * (a2 * y - a3 * x)
    y_dt = c1 * (c2 * x - c3 * y) - c4 * x * z
    z_dt = d1 * x * y - d2 * z
    return x_dt, y_dt, z_dt


# управление по ренте постановка 1
# макроперемнная Y-ro*Z
def lorenzsys_upr_z1(XYZ, t, a1, a2, a3, c1, c2, c3, c4, d1, d2, ro, T):
    x, y, z = XYZ
    psi = y- ro*z
    #u =
    x_dt = a1 * (a2 * y - a3 * x)
    y_dt = c1 * (c2 * x - c3 * y) - c4 * x * z
    z_dt = d1 * x * y - d2 * z+u
    return x_dt, y_dt, z_dt


# отображение графика
def lineplot(x_data, y_data, x_label="", y_label="", title=""):
    # Create the plot object
    _, ax = plt.subplots()
    ax.plot(x_data, y_data, lw=2, color='#539caf', alpha=1)
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)


def ode_array():
    x_r, y_r, z_r = (get_statistic())
    arr = make_array()
    dict = {}
    count = 0
    for j in arr:
        count+=1
    # начальные значения модельных данных = начальным реальных
    # рассчитываем систему для каждого набора коэффициентов
    # считаем сумму квадратов отклонения для каждого набора относительно реальных данных
    #for r in range(0, len(arr)):
        #row = arr[r]
        q_sum = 0
        len_t = len(x_r)

        X, Y, Z = lorenz_sys(j, [x_00, y_00, z_00])
        for i in range(0, len_t):
            q_sum += (X[i] - x_r[i]) ** 2 + (Y[i] - y_r[i]) ** 2 + (Z[i] - z_r[i]) ** 2
        dict[j] = q_sum / len_t
        with open('dict.csv', "a", newline="") as csvfile:
            elem = [j, q_sum]
            writer = csv.writer(csvfile, dialect='excel')
            writer.writerow([j,q_sum])
        print(count)
    return dict


def lorenz_sys(koef, arr_0):
    a1, a2, a3, c1, c2, c3, c4, d1, d2 = koef[0], koef[1], koef[2], koef[3], koef[4], koef[5], koef[6], koef[7], koef[8]
    x_0, y_0, z_0 = arr_0[0], arr_0[1], arr_0[2]

    tmax, n = t_len, t_len
    t = np.linspace(0, tmax, n)

    f = odeint(lorenzsys, (x_0, y_0, z_0), t, args=(a1, a2, a3, c1, c2, c3, c4, d1, d2))
    X, Y, Z = f.T
    return X, Y, Z

def lorenz_sys_rk(koef):
    a1, a2, a3, c1, c2, c3, c4, d1, d2 = koef[0], koef[1], koef[2], koef[3], koef[4], koef[5], koef[6], koef[7], koef[8]

    tmax, n = t_len, t_len
    t = np.linspace(0, tmax, n)

    f = RK45(lorenzsys, 0.0, (x_00, y_00, z_00), t_len, f_args =(a1, a2, a3, c1, c2, c3, c4, d1, d2))
    X, Y, Z = f.T
    return X, Y, Z

def main():
    x_y_z = (get_statistic())
    # lineplot(range(1,6), x_y_z[i], 't', 'ось значений', 'X')
    fig, axes = plt.subplots(6, 3)
    len_arr = len(x_y_z[0])
    for i in range(0, 3):
        axes[0, i].plot(range(0, len_arr), x_y_z[i])

    cust_koef=[]

    cust_koef.append((1.0, 2.8, 3.4, 2.1, 1.0, 0.9, 7.5, 5.0, 0.8))
    cust_koef.append((1.2, 0.6, 2.7, 2.1, 0.9, 0.9, 7.5, 5.0, 1.6))
    cust_koef.append((0.9, 1.0, 3.4, 2.1, 0.9, 1.0, 7.5, 5.0, 0.8))
    cust_koef.append((1.2, 2.19, 3.7, 2.1, 0.9, 1.0, 7.5, 5.0, 1.6))
    cust_koef.append((1, 1.78, 1.58, 0.28, 1, 1, 0.72, 13.54, 10.92))
    arr_xyz_lorenz_cust=[]
    for i in range(0, len(cust_koef)):
        arr_xyz_lorenz_cust.append(lorenz_sys_rk(cust_koef[i]))
    for j in range(1, len(arr_xyz_lorenz_cust)):
        for i in range(0, 3):
            axes[j, i].plot(range(0, len_arr), arr_xyz_lorenz_cust[j-1][i])
   # поиск лучших коэф. возьмет первые 5 с мин суммой квадратов отклонений
   #sum_dict = ode_array()
   #min_sum = sorted(sum_dict.items(), key=lambda kv: kv[1])[:5]
   #best_coef_arr = []
   #arr_xyz_lorenz = []

   ## выводим 5 наилучших сочетаний коэф с минимальной суммой квадратов ошибки
   #for i in range(0, len(min_sum)):
   #    best_coef_arr.append(min_sum[i][0])
   #    arr_xyz_lorenz.append(lorenz_sys(min_sum[i][0], [x_00, y_00, z_00]))
   #    print(min_sum[i][0])
   #for j in range(1, len(arr_xyz_lorenz)):
   #    for i in range(0, 3):
   #        axes[j, i].plot(range(0, len_arr), arr_xyz_lorenz[j][i])

    plt.show()


if __name__ == '__main__':
    main()
