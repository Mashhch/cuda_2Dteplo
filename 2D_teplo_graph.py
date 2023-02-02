import numpy as np
import matplotlib.pyplot as plt




def q(x, y,t):
    return x+y-4


def resh(x,y,t):
    return x**2+x*t+y**2+y*t


# Начальное условие
def phi(x, y):
    return x**2+y**2


# теплоизолированные стенки
# ГУ для оси x
def alpha_x():
    alpha_x1 = 1
    alpha_x2 = 1

    return [alpha_x1, alpha_x2]


def beta_x():
    beta_x1 = 0
    beta_x2 = 0

    return [beta_x1, beta_x2]


def gamma_x(y, t):
    gamma_x1 = -2+t
    gamma_x2 = 2+t

    return [gamma_x1, gamma_x2]


# ГУ для оси y
def alpha_y():
    alpha_y1 = 1
    alpha_y2 = 1

    return [alpha_y1, alpha_y2]


def beta_y():
    beta_y1 = 0
    beta_y2 = 0

    return [beta_y1, beta_y2]


def gamma_y(x, t):
    gamma_y1 = t-4
    gamma_y2 = t+4
    return [gamma_y1, gamma_y2]


# условие устойчивости Куранта
def Kurant_condition(h_x, h_y, a):
    t_x = h_x ** 2 / 2 / a ** 2 / 2
    t_y = h_y ** 2 / 2 / a ** 2 / 2

    t = min(t_x, t_y)

    return t

def solution_1ord(x, y, h_x, h_y, Nt, a = 1):
    t = Kurant_condition(h_x, h_y, a)
    c_x = a ** 2 * t / h_x ** 2
    c_y = a ** 2 * t / h_y ** 2
    u = []
    y_x_layer = np.zeros((len(y), len(x)))
    for j in range(len(y)):
        for i in range(len(x)):
            y_x_layer[j][i] = phi(x[i], y[j])
    t_prev_layer = y_x_layer
    #u.append(y_x_layer)
    for t_i in range(1, Nt):
        y_x_layer = np.zeros((len(y), len(x)))
        for j in range(1,len(y)-1):
            for i in range(1,len(x)-1):
                y_x_layer[j][i] = (c_x * (t_prev_layer[j][i + 1] - 2 * t_prev_layer[j][i] + t_prev_layer[j][i - 1]) +
                                   c_y * (t_prev_layer[j + 1][i] - 2 * t_prev_layer[j][i] + t_prev_layer[j - 1][i])
                                   + t_prev_layer[j][i] + t*q(x[i],y[j],t_i*t))
        for j in range(1,len(y)-1):
            y_x_layer[j][0] = (gamma_x(y[j], t_i*t)[0] * h_x- alpha_x()[0] * y_x_layer[j][1] ) / \
                              (beta_x()[0] * h_x- alpha_x()[0])
            y_x_layer[j][len(x)-1] = (gamma_x(y[j], t_i*t)[-1] * h_x + alpha_x()[-1] * y_x_layer[j][len(x)-2]) /\
                              (beta_x()[-1] * h_x + alpha_x()[-1])

        for i in range(len(x)):
            y_x_layer[0][i] = (gamma_y(x[i], t_i*t)[0] * h_y - alpha_y()[0] * y_x_layer[1][i]) /\
                      (beta_y()[0] * h_y - alpha_y()[0])
            y_x_layer[len(y)-1][i] = (gamma_y(x[i], t_i*t)[-1] * h_y + alpha_y()[-1] * y_x_layer[len(y)-2][i]) /\
                      (beta_y()[-1] * h_y + alpha_y()[-1])
        t_prev_layer = np.copy(y_x_layer)
        u= np.copy(t_prev_layer)

    return u


def solution(x, y, h_x, h_y, Nt, a = 1):
    t = Kurant_condition(h_x, h_y, a)
    c_x = a ** 2 * t / h_x ** 2
    c_y = a ** 2 * t / h_y ** 2
    u = []
    y_x_layer = np.zeros((len(y), len(x)))
    t_prev_layer = np.zeros((len(y), len(x)))
    for j in range(len(y)):
        for i in range(len(x)):
            t_prev_layer[j][i] = phi(x[i], y[j])
    #u.append(y_x_layer)
    for t_i in range(1, Nt+1):
        y_x_layer = np.zeros((len(y), len(x)))
        for j in range(1,len(y)-1):
            for i in range(1,len(x)-1):
                y_x_layer[j][i] = (c_x * (t_prev_layer[j][i + 1] - 2 * t_prev_layer[j][i] + t_prev_layer[j][i - 1]) +
                                   c_y * (t_prev_layer[j + 1][i] - 2 * t_prev_layer[j][i] + t_prev_layer[j - 1][i])
                                   + t_prev_layer[j][i] + t*q(x[i],y[j],t_i*t))
        for j in range(1,len(y)-1):
            y_x_layer[j][0] = c_x * (2*t_prev_layer[j][1] - 2 * t_prev_layer[j][0] +
                                     2*h_x*(beta_x()[0]*t_prev_layer[j][0]-gamma_x(y[j], t_i*t)[0])/alpha_x()[0]) + t_prev_layer[j][0]
            #y_x_layer[j][0] = (gamma_x(y[j], t_i*t)[0] * h_x- alpha_x()[0] * y_x_layer[j][1] ) / \
            #                  (beta_x()[0] * h_x- alpha_x()[0])
            y_x_layer[j][len(x)-1] = c_x * (2 * t_prev_layer[j][len(x)-2] - 2 * t_prev_layer[j][len(x)-1] +
                                     2 * h_x * ( - beta_x()[-1] * t_prev_layer[j][len(x)-1] + gamma_x(y[j], t_i * t)[-1]) /
                                     alpha_x()[-1]) + t_prev_layer[j][len(x)-1]
            #y_x_layer[j][len(x)-1] = (gamma_x(y[j], t_i*t)[-1] * h_x + alpha_x()[-1] * y_x_layer[j][len(x)-2]) /\
            #                  (beta_x()[-1] * h_x + alpha_x()[-1])

        for i in range(len(x)):
            y_x_layer[0][i] = c_y * (2 * t_prev_layer[1][i] - 2 * t_prev_layer[0][i] +
                                     2 * h_y * (beta_y()[0] * t_prev_layer[0][i] - gamma_y(x[i], t_i * t)[0]) /
                                     alpha_y()[0]) + t_prev_layer[0][i]
            # y_x_layer[0][i] = (gamma_y(x[i], t_i*t)[0] * h_y - alpha_y()[0] * y_x_layer[1][i]) /\
            #           (beta_y()[0] * h_y - alpha_y()[0])
            y_x_layer[len(y) - 1][i] = c_y * (2 * t_prev_layer[len(y) - 2][i] - 2 * t_prev_layer[len(y) - 1][i] +
                                       2 * h_y * (- beta_y()[-1] * t_prev_layer[len(y)-1][i] + gamma_y(x[i], t_i * t)[-1]) /
                                       alpha_y()[-1]) + t_prev_layer[len(y) - 1][i]
            # y_x_layer[len(y)-1][i] = (gamma_y(x[i], t_i*t)[-1] * h_y + alpha_y()[-1] * y_x_layer[len(y)-2][i]) /\
            #           (beta_y()[-1] * h_y + alpha_y()[-1])
        t_prev_layer = np.copy(y_x_layer)
        u= np.copy(t_prev_layer)

    return u


xmin = -1
xmax = 1
ymin = -2
ymax = 2
h_x = 0.01
h_y = 0.01
Nt = 5
x = np.arange(xmin, xmax+ h_x, h_x)
y = np.arange(ymin, ymax+ h_y, h_y)

t = Kurant_condition(h_x, h_y, a=1)
my = solution(x, y, h_x, h_y, Nt)
er= np.zeros((len(y), len(x)))
maxx = -100
for j in range(len(y)):
    for i in range(len(x)):
        er[j][i]= abs(my[j][i] - resh(x[i], y[j], Nt*t))
        if maxx < er[j][i]:
            maxx = er[j][i]

for j in range(len(y)):
    for i in range(len(x)):
        if er[j][i] == maxx:
            print(j, i)

hx = [0.02, 0.01, 0.005, 0.0025]
hx2 = [hx[i]**2 for i in range(len(hx))]
hy = [0.02, 0.01, 0.005, 0.0025]
errs = []
for i in range(len(hx)):
    h_x = hx[i]
    h_y = hy[i]
    x_ = np.arange(xmin, xmax + h_x, h_x)
    y_ = np.arange(ymin, ymax + h_y, h_y)
    my1 = solution(x_, y_, h_x, h_y, Nt)
    mysol = np.array(my1)
    t = Kurant_condition(h_x, h_y, a=1)
    er= np.zeros((len(y_), len(x_)))
    maxx = -100
    for j in range(len(y_)):
        for i in range(len(x_)):
            er[j][i]= abs(mysol[j][i] - resh(x_[i], y_[j], Nt*t))
            if maxx < er[j][i]:
                maxx=er[j][i]
            if (er[j][i] >= 0.01):
                print(er[j][i], i , j)
    errs.append(maxx)
print(errs)
h_log = [np.log(hx[i]) for i in range(len(hx))]
h_log2 = [np.log(hx[i]**2) for i in range(len(hx))]
errs_log = [np.log(errs[i]) for i in range(len(hx))]
plt.figure()
plt.plot(h_log, errs_log, color='yellow', label = 'maxerr')
plt.plot(h_log, h_log2, color='green', label = 'y=x^2')
plt.title(' график ошибки')
plt.legend()
plt.grid(True)
plt.show()