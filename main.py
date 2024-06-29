import matplotlib.pyplot as plt
import numpy as np
import math
class cubic_spline:
    splines = []
    n = None        # количество узлов сетки
    # Построение сплайна
    # x - узлы сетки, должны быть упорядочены по возрастанию, кратные узлы запрещены
    # y - значения функции в узлах сетки
    # n - количество узлов сетки
    def build_spline(self, x, y, n):
        self.n = n
        self.splines = [{'a': y[i], 'b': 0, 'c': 0, 'd': 0, 'x': x[i]} for i in range(self.n)]
        self.splines[0]['c'] = 0
        # Решение СЛАУ относительно коэффициентов сплайнов c[i] методом прогонки для трехдиагональных матриц
        # Вычисление прогоночных коэффициентов - прямой ход метода прогонки
        alpha = []
        beta = []
        A, B, C, F, h_i, z = None,None,None,None,None,None
        alpha.append(0)
        beta.append(0)
        for i in range(1, n-1):
            h_i = x[i] - x[i - 1]
            h_i1 = x[i + 1] - x[i]
            A = h_i
            C = 2 * (h_i + h_i1)
            B = h_i1
            F = 6 * ((y[i+1]-y[i])/h_i1 - (y[i] - y[i-1]) / h_i)
            z = (A * alpha[i-1] + C)
            alpha.append(-B/z)
            beta.append((F-A * beta[i-1]) / z)
        self.splines[n - 1]['c'] = (F - A * beta[n - 2]) / (C + A * alpha[n - 2])
        # Нахождение решения - обратный ход метода прогонки
        for i in range(n-2,0,-1):
            self.splines[i]['c'] = alpha[i] * self.splines[i + 1]['c'] + beta[i]
        # По известным коэффициентам c[i] находим значения b[i] и d[i]
        for i in range(n-1,0,-1):
            h_i = x[i] - x[i - 1];
            self.splines[i]['d'] = (self.splines[i]['c'] - self.splines[i - 1]['c']) / h_i
            self.splines[i]['b'] = h_i * (2 * self.splines[i]['c'] + self.splines[i - 1]['c']) / 6 + (y[i] - y[i - 1]) / h_i
    
    # Вычисление значения интерполированной функции в произвольной точке
    def f(self, x):
        if (not self.splines):
            return None; # Если сплайны ещё не построены - возвращаем None
        if (x <= self.splines[0]['x']):             # Если x меньше точки сетки x[0] - пользуемся первым эл-тов массива
            s = self.splines[1];
        elif (x >= self.splines[self.n - 1]['x']):  # Если x больше точки сетки x[n - 1] - пользуемся последним эл-том массива
            s = self.splines[self.n - 1]
        else:                                       # Иначе x лежит между граничными точками сетки - производим бинарный поиск нужного эл-та массива
            i = 0
            j = self.n - 1
            while (i + 1 < j):
                k = i + (j - i) // 2
                if (x <= self.splines[k]['x']):
                    j = k
                else:
                    i = k
                s = self.splines[j]
        dx = (x - s['x'])
        return s['a'] + (s['b'] + (s['c'] / 2 + s['d'] * dx / 6) * dx) * dx; # Вычисляем значение сплайна в заданной точке по схеме Горнера

def read_file():
    f = open('input.txt', 'r')
    # чтение строки с x-ами
    temp = [f.readline()]
    temp[0] = [float(x) for x in temp[0].split()]
    # чтение строки с y-ами
    temp.append(f.readline())
    temp[1] = [float(x) for x in temp[1].split()]
    return temp
# Функция
def function(x):
    #return math.fabs(x)    # y = |x|
    #return math.exp(-x)    # y = e^(-x)
    return math.sin(x)      # y = sin x
print("Выберете режим (1 или 2): ", end='')
mode = float(input())
a = cubic_spline()
info = read_file()
x = info[0]
if(mode == 1):
    y = info[1]
    a.build_spline(x,y,len(x))
    print("Введите x: ", end='')
    inp = float(input())
    print("F(x) = ",a.f(inp))
else:
    y = [function(i) for i in x]
    a.build_spline(x,y,len(x))
    print("Выберете режим (1 или 2): ", end='')
    mode = float(input())
    # Вычисление точек для отрисовки графиков
    step = 0.01
    pointsX = [x[0] + i * step for i in range(int((x[len(x)-1] - x[0]) / step) + 1)]
    pointsY = [function(i) for i in pointsX]        # Поиск y 
    interpolarPointsY = [a.f(i) for i in pointsX]   # y-ки для интерполированного графика
    # изменение значение y для одной из точек
    if(mode == 2):
        y[len(y)//2] = y[len(y)//2]*2   # центральная точка поднимается в 2 раза
        a.build_spline(x,y,len(x))
        incorrectInterpolarPointsY = [a.f(i) for i in pointsX]  # расчёт y для интерполированной функции с неправильной точкой
    # Поиск отклонения
    maxDeviation = 0
    for i in range(len(pointsX)):
        if(abs(abs(pointsY[i]) - abs(interpolarPointsY[i])) > maxDeviation): 
            maxDeviation = abs(abs(pointsY[i]) - abs(interpolarPointsY[i]))
    print("Максимальная погрешность = ",maxDeviation)
    # Отрисовка графиков
    plt.plot(pointsX, pointsY, label="y = sin x")
    plt.plot(pointsX, interpolarPointsY, label="Интерполированный y = sin x")
    if(mode == 2):
        plt.plot(pointsX, incorrectInterpolarPointsY, label="С изменённой таблицей")
    plt.legend()
    plt.grid(True)
    plt.show()
