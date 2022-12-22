import numpy as np
import matplotlib.pyplot as plt


# Количество точек, удовлетворяющих условию для круга
def count_for_cells(x, y, r=0.25):
    counter_for_cells = 0
    for i in x:
        for j in y:
            if (i ** 2 + j ** 2 <= r ** 2):
                counter_for_cells += 1
    return counter_for_cells


def q(x, h_x, y, h_y, r, counter_for_cells):
    if (x ** 2 + y ** 2 <= r ** 2):
        return 1 / (counter_for_cells * h_x * h_y)
    return 0


# Начальное условие
def phi(x, y):
    return 0


# теплоизолированные стенки
# ГУ для оси x
def alpha_x(y, t):
    alpha_x1 = 1
    alpha_x2 = 1

    return [alpha_x1, alpha_x2]


def beta_x(y, t):
    beta_x1 = 0
    beta_x2 = 0

    return [beta_x1, beta_x2]


def gamma_x(y, t):
    gamma_x1 = 0
    gamma_x2 = 0

    return [gamma_x1, gamma_x2]


# ГУ для оси y
def alpha_y(x, t):
    alpha_y1 = 1
    alpha_y2 = 1

    return [alpha_y1, alpha_y2]


def beta_y(x, t):
    beta_y1 = 0
    beta_y2 = 0

    return [beta_y1, beta_y2]


def gamma_y(x, t):
    gamma_y1 = 0
    gamma_y2 = 0
    return [gamma_y1, gamma_y2]


# условие устойчивости Куранта
def Kurant_condition(h_x, h_y, a):
    t_x = h_x ** 2 / 2 / a ** 2 / 2
    t_y = h_y ** 2 / 2 / a ** 2 / 2

    t = min(t_x, t_y)

    return t

def solution(x, y, h_x, h_y, Nt, a = 1, r = 0.25):
    t = Kurant_condition(h_x, h_y, a)
    counter_for_cells = count_for_cells(x, y, r)

    c_x = a ** 2 * t / h_x ** 2
    c_y = a ** 2 * t / h_y ** 2

    u = []

    y_x_layer = np.zeros((len(y), len(x)))#phi == 0

    t_prev_layer = y_x_layer
    u.append(y_x_layer)
    q_ = np.zeros((len(y), len(x)))
    for j in range(1, len(y) - 1):
        for i in range(1, len(x) - 1):
            q_[j][i] = q(x[i], h_x, y[j], h_y, r, counter_for_cells)
    for t_i in range(Nt):
        y_x_layer = np.zeros((len(y), len(x)))
        for j in range(1,len(y)-1):
            for i in range(1,len(x)-1):
                y_x_layer[j][i] = (c_x * (t_prev_layer[j][i + 1] - 2 * t_prev_layer[j][i] + t_prev_layer[j][i - 1]) +
                          c_y * (t_prev_layer[j + 1][i] - 2 * t_prev_layer[j][i] + t_prev_layer[j - 1][i]) +
                          t_prev_layer[j][i] + q(x[i], h_x, y[j], h_y, r, counter_for_cells))
        for j in range(1,len(y)-1):
            y_x_layer[j][0] = (gamma_x(y[j], t_i)[0] * h_x- alpha_x(y[j], t_i)[0] * y_x_layer[j][1] ) / \
                              (beta_x(y[j], t_i)[0] * h_x- alpha_x(y[j], t_i)[0])
            y_x_layer[j][len(x)-1] = (gamma_x(y[j], t_i)[-1] * h_x - alpha_x(y[j], t_i)[-1] * y_x_layer[j][len(x)-2]) /\
                              (beta_x(y[j], t_i)[-1] * h_x - alpha_x(y[j], t_i)[-1])

        for i in range(len(x)):
            y_x_layer[0][i] = (gamma_y(x[i], t_i)[0] * h_y - alpha_y(x[i], t_i)[0] * y_x_layer[1][i]) /\
                      (beta_y(x[i], t_i)[0] * h_y - alpha_y(x[i], t_i)[0])
            y_x_layer[len(y)-1][i] = (gamma_y(x[i], t_i)[-1] * h_y - alpha_y(x[i], t_i)[-1] * y_x_layer[len(y)-2][i]) /\
                      (beta_y(x[i], t_i)[-1] * h_y - alpha_y(x[i], t_i)[-1])
        t_prev_layer = np.copy(y_x_layer)
        u.append(t_prev_layer)

    return u


xmin = -1
xmax = 1
ymin = -2
ymax = 2
h_x = 0.1
h_y = 0.2
t = 0.05
Nt = 10
x = np.arange(xmin, xmax+ h_x, h_x)
y = np.arange(ymin, ymax+ h_y, h_y)
counter_for_cells = count_for_cells(x, y, 0.25)
#print(q(x[36], h_x, y[71], h_y, 0.25, counter_for_cells))
my = solution(x, y, h_x, h_y, Nt, r = 0.25)
print(1)

fig, ax = plt.subplots(7)
fig.set_size_inches((20,20))

for i in range(7):
    ax[i].imshow(my[i*10], extent=[xmin, xmax, ymin, ymax], interpolation='bilinear', origin='lower',
                     cmap='jet')
    ax[i].set_title(f'Шаг {i*100}')

plt.show()