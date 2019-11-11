import numpy as np

'''Входящие данные:
        T0 = int              - начальная температура (K) 
        q = np.ndarray(5,)    - мольные потоки i компонетов(моль/с)
        с = np.ndarray(2,)    - конверсия пропана (%)
        P = 1                 - давление (атм)
   
   Другие переменные:
        n_i = np.ndarray(5,)  - мольные доли компонентов смеси
        Cp = np.ndarray(5,)   - теплоемкость i компонетов (кДж/моль*К)
        T = int               - выходная температура (К)
        Q0_input = int        - энтальпия входящей смеси (кДж/моль)
        Q0_output = int       - энтальпия выходящей смеси(кДж/моль)


   Считаем, что массивы np.ndarray(5,) передаются в виде:
   ([propane, propene, hydrogen, ethene, methane])

'''


def calculate_temperature(T0, q, c, Q0, a):
    n_i = q / sum(q)                                               # Мольная доля i компонента по входящему потоку
    Q0_input = Q0 * n_i
    Q0_input = np.split(Q0_input, [1, 5])
    Q0_input = sum(Q0_input[1]) - sum(2*Q0_input[0])               # Стандартная энтальпия входящей смеси

    a_input = np.zeros(a.shape)
    Cp = np.array([0, 0, 0, 0, 0])
    for i in range(0, len(a_input)):
        for j in range(0, len(a_input[i])):
            a_input[i][j] = a[i][j] * (T0 - 298.15)**(j+1)

        Cp[i] = sum(a_input[i])
    integrated_Cp_input = sum(Cp * n_i)                             # Теплоемкость входящей смеси

    c = c / 100
    q_1 = np.array([q[0] * (1 - c[0]), q[1] + q[0]*c[0],
                    q[2] + q[0]*c[0], q[3], q[4]])                  # Мольный поток подуктов по первой реакции
    q_1 = np.array([q_1[0] * (1 - c[1]), q_1[1], q_1[2],
                    q_1[3] + q_1[0]*c[1], q_1[4] + q_1[0]*c[1]])    # Мольный поток подуктов по второй реакции

    n_i = q_1 / sum(q_1)                                            # Мольная доля i компонента по выходящему потоку

    Q0_output = Q0 * n_i
    Q0_output = np.split(Q0_output, [1, 5])
    Q0_output = sum(Q0_output[1]) - sum(2*Q0_output[0])             # Стандартная энтальпия выходящей смеси

    a_poly_constant = a * n_i
    a_poly_constant = np.sum(a_poly_constant, axis=0)
    a_poly_constant = np.flip(a_poly_constant, axis=0)
    a_poly_constant = np.append(a_poly_constant,
                        (Q0_input+integrated_Cp_input-Q0_output)*(-1))
    delta_T_poly = np.roots(a_poly_constant)                        # Находим корни полинома
    T = delta_T_poly[4].real + 298.15                               # и выбираем не вещественное число
    return T





Q0 = np.array([-103.85, 20.41, 0, 52.30, -74.85])
a = np.array([[3.847, 5.131*10**(-3), 6.011*10**(-5), -7.793*10**(-8), 3.079*10**(-11)],
             [3.834, 3.893*10**(-3), 4.688*10**(-5), -6.013*10**(-8), 2.283*10**(-11)],
             [2.883, 3.681*10**(-3), -0.772*10**(-5), 0.692*10**(-8), -0.213*10**(-11)],
             [4.221, -8.782*10**(-3), 5.795*10**(-5), -6.729*10**(-8), 2.511*10**(-11)],
             [4.568, -8.975*10**(-3), 3.631*10**(-5), -3.407*10**(-8), 1.091*10**(-11)]])

T0 = int(input())
c = np.array([10, 10])
q = np.array([35.3, 8.0, 11.5, 6.2, 7.1])

result = calculate_temperature(T0, q, c, Q0, a)
print('Temperature: %.2f' % result)
