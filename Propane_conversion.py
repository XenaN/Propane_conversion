import numpy as np

'''Входящие данные:
        T0 = int              - начальная температура (K) 
        q = np.ndarray(5,)    - мольные потоки i компонетов(моль/с)
        с = np.ndarray(2,)    - конверсия пропана (%)
        P = 1                 - давление (атм)
   
   Другие переменные:
        po = int              - плотность смеси (г/см3)
        Cp = []               - теплоемкость i компонетов (кДж/моль*К)
        Cp_mix = int          - теплоемкость смеси (Дж/моль*К)
        T = int               - выходная температура (К)
        Q = np.ndarray(2,)    - тепловой эффект реакции (кДж/моль)
        Q0 = np.ndarray(2,)   - тепловой эффект реакции для стандартных условий(кДж/моль)
        с0 = np.ndarray(2,)   - начальная концентрация пропана (моль/л)
        Mr = np.ndarray(5,)   - молярная масса i компонетов (г/моль)

   Считаем, что массивы np.ndarray(5,) передаются в виде:
   ([propane, propene, hydrogen, ethene, methane])

   Основываясь на законе Гессе, напишем уравнение теплового баланса:
        po * Cp * (T - T0) = sum(Q * C0 * c)
        
    Отсюда находим температуру для экзотермической реакции
        T = T0 - sum(Q * c0 * c)/(po * Cp)
'''


def calculate_temperature(T0, q, c, P, Q0, Cp, Mr):
    #Calculate enthalpy
    T_standart = 298.15
    Cp_initial = np.array([Cp[0], Cp[0]])                        # теплоемкость пропана
    Cp_product = np.array([Cp[1] + Cp[2], Cp[3] + Cp[4]])        # теплоемкость продуктов
    Q = Q0 + (Cp_initial - Cp_product) * (T0 - T_standart)

    #calculate propane concentrate
    R = 8.31                             # Дж / (моль*К)
    P = P * 101325                       # Па
    c0 = P / (R * T0)

    #Calculate heat capacity
    n_i = q / sum(q)                     # мольная доля i компонента
    Cp_mix = sum(Cp * n_i)

    #Calculate density
    Mr_mix = sum(Mr * n_i)
    po = P * Mr_mix / (R * T0)           # г/м3 или кг/л

    # Cp_mix из кДж/(моль*К) в кДж/(кг*К)
    Cp_mix = Cp_mix * 1000 / Mr_mix

    #Calculate temperature
    T = T0 - sum(Q * c0 * c) / (po * Cp_mix)

    return T



# Q0 считаем как разницу энтальпий образования исходных веществ и продуктов
# 1 реакции Q0 = -103.85 - (20.41 + 0) = -124.26
# 2 реакции Q0 = -103.85 - (-74.85 + 52.30) = -81.3
Q0 = np.array([-124.26, -81.3])

Mr = np.array([44.1, 42.08, 2.02, 28.05, 16.04])
Cp = [0.0736, 0.0643, 0.0288, 0.0429, 0.03569]
P = 1

T0 = 350
c = np.array([19, 25])
q = np.array([35.3, 8.0, 11.5, 6.2, 7.1])

result = calculate_temperature(T0, q, c, P, Q0, Cp, Mr)
print('Temperature: %.2f' % result)
