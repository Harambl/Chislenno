#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <clocale>

long double func(long double x) {
    return x - sin(x);
}

long double func_pr(long double x) {
    return 1 - cos(x);
}

long double func_2pr(long double x) {
    return sin(x);
}

long double func_3pr(long double x) {
    return cos(x);
}

long double LagranF1(long double x_, long double xi,long  double xi_1) {
    return func(xi) * (x_ - xi_1) / (xi - xi_1) + func(xi_1) * (x_ - xi) / (xi_1 - xi);
}

long double LagranF2(long double x_, long double xi, long double xi_1, long double x1) {
    return func(x1) * ((x_ - xi) * (x_ - xi_1)) / ((x1 - xi) * (x1 - xi_1)) 
        + func(xi) * ((x_ - x1) * (x_ - xi_1)) / ((xi - x1) * (xi - xi_1)) 
        + func(xi_1) * ((x_ - x1) * (x_ - xi)) / ((xi_1 - x1) * (xi_1 - xi));
}



using namespace std;

int main() {
    setlocale(LC_ALL, "Russian");

    long double a = 0.6;
    long double b = 1.1;

    vector<long double> x_ = { 0.88, 0.63, 1.08, 0.83 };

    long double h = (b - a) / 10;

    long double x__ = x_[0];
    int i_ = 0;
    long double xi = 0;

    vector<long double> x(11, 0);
    for (int i = 0; i < 11; ++i) {
        x[i] = a + i * h;
        if (x[i] < x__) {
            if (x[i] > xi) {
                xi = x[i];
                i_ = i;
            }
        }
    }

    cout << "Таблица x_i" << endl;
    int count = 1;
    for (long double xi : x) {
        cout << "x[" << count << "]= " << xi << endl;
        count++;
    }

    cout << endl;

    {
        long double exactf = func(x__);
        long double L1 = LagranF1(x__, xi, x[i_ + 1]);
        long double true_error1 = L1 - exactf;


        long double min_f = min(func_2pr(xi), func_2pr(x[i_ + 1]));
        long double max_f = max(func_2pr(xi), func_2pr(x[i_ + 1]));

        long double omega1 = (x__ - xi) / (x__ - x[i_ + 1]);
        long double minR1 = min_f * omega1 / 2;
        long double maxR1 = max_f * omega1 / 2;


        long double L2 = LagranF2(x__, xi, x[i_ + 1], x[i_ - 1]);
        long double true_error2 = L2 - exactf;

        long double min_f2 = min(func_3pr(xi), func_3pr(x[i_ + 1]));
        long double max_f2 = max(func_3pr(xi), func_3pr(x[i_ + 1]));

        long double omega2 = (x__ - x[i_ - 1]) * (x__ - xi) * (x__ - x[i_ + 1]);
        long double minR2 = min_f * omega2 / 6;
        long double maxR2 = max_f * omega2 / 6;

        cout << "При x_ = " << x__ << endl;
        cout << "L1= " << L1 << endl;
        cout << "min f''(x) = " << min_f << endl;
        cout << "max f''(x) = " << max_f << endl;
        cout << "min R1 = " << minR1 << endl;
        cout << "max R1 = " << maxR1 << endl;
        cout << "L1(x*) - f(x*) = " << true_error1 << " = R1(x*) " << endl;
        if (minR1 < true_error1 < maxR1) {
            cout << "Неравенство min R1 < R1 < max R1 выполнено " << endl;
        }
        else {
            cout << "Неравенство min R1 < R1 < max R1 не выполнено " << endl;
        }
        if (abs(true_error1) < 0.0001) {
            cout << "Погрешность = " << abs(true_error1) << " подходит " << endl;
        }
        else {
            cout << "Погрешность = " << abs(true_error1) << "  не подходит " << endl;
        }

        cout << endl;

        cout << "L2= " << L2 << endl;
        cout << "min f'''(x) = " << min_f2 << endl;
        cout << "max f'''(x) = " << max_f2 << endl;
        cout << "min R2 = " << minR2 << endl;
        cout << "max R2 = " << maxR2 << endl;
        cout << "L2(x*) - f(x*) = " << true_error2 << " = R2(x*) " << endl;
        if (minR2 < true_error2 < maxR2) {
            cout << "Неравенство min R2 < R2 < max R2 выполнено " << endl;
        }
        else {
            cout << "Неравенство min R2 < R2 < max R2 не выполнено " << endl;
        }
        if (abs(true_error2) < 0.00001) {
            cout << "Погрешность = " << abs(true_error2) << " подходит " << endl;
        }
        else {
            cout << "Погрешность = " << abs(true_error2) << "  не подходит " << endl;
        }
        cout << endl;
    }

    x__ = x_[1];
    i_ = 0;
    xi = 0;

    for (int i = 0; i < 11; ++i) {
        if (x[i] < x__) {
            if (x[i] > xi) {
                xi = x[i];
                i_ = i;
            }
        }
    }

    {
        long double exactf = func(x__);
        long double L1 = LagranF1(x__, xi, x[i_ + 1]);
        long double true_error1 = L1 - exactf;


        long double min_f = min(func_2pr(xi), func_2pr(x[i_ + 1]));
        long double max_f = max(func_2pr(xi), func_2pr(x[i_ + 1]));

        long double omega1 = (x__ - xi) / (x__ - x[i_ + 1]);
        long double minR1 = min_f * omega1 / 2;
        long double maxR1 = max_f * omega1 / 2;


        long double L2 = LagranF2(x__, xi, x[i_ + 1], x[i_ + 2]);
        long double true_error2 = L2 - exactf;

        long double min_f2 = min(func_3pr(xi), func_3pr(x[i_ + 1]));
        long double max_f2 = max(func_3pr(xi), func_3pr(x[i_ + 1]));

        long double omega2 = (x__ - x[i_ + 2]) * (x__ - xi) * (x__ - x[i_ + 1]);
        long double minR2 = min_f * omega2 / 6;
        long double maxR2 = max_f * omega2 / 6;

        cout << "При x_ = " << x__ << endl;
        cout << "L1= " << L1 << endl;
        cout << "min f''(x) = " << min_f << endl;
        cout << "max f''(x) = " << max_f << endl;
        cout << "min R1 = " << minR1 << endl;
        cout << "max R1 = " << maxR1 << endl;
        cout << "L1(x*) - f(x*) = " << true_error1 << " = R1(x*) " << endl;
        if (minR1 < true_error1 < maxR1) {
            cout << "Неравенство min R1 < R1 < max R1 выполнено " << endl;
        }
        else {
            cout << "Неравенство min R1 < R1 < max R1 не выполнено " << endl;
        }
        if (abs(true_error1) < 0.0001) {
            cout << "Погрешность = " << abs(true_error1) << " подходит " << endl;
        }
        else {
            cout << "Погрешность = " << abs(true_error1) << "  не подходит " << endl;
        }

        cout << endl;

        cout << "L2= " << L2 << endl;
        cout << "min f'''(x) = " << min_f2 << endl;
        cout << "max f'''(x) = " << max_f2 << endl;
        cout << "min R2 = " << minR2 << endl;
        cout << "max R2 = " << maxR2 << endl;
        cout << "L2(x*) - f(x*) = " << true_error2 << " = R2(x*) " << endl;
        if (minR2 < true_error2 < maxR2) {
            cout << "Неравенство min R2 < R2 < max R2 выполнено " << endl;
        }
        else {
            cout << "Неравенство min R2 < R2 < max R2 не выполнено " << endl;
        }
        if (abs(true_error2) < 0.00001) {
            cout << "Погрешность = " << abs(true_error2) << " подходит " << endl;
        }
        else {
            cout << "Погрешность = " << abs(true_error2) << "  не подходит " << endl;
        }
        cout << endl;
    }

    x__ = x_[2];
    i_ = 0;
    xi = 0;

    for (int i = 0; i < 11; ++i) {
        if (x[i] < x__) {
            if (x[i] > xi) {
                xi = x[i];
                i_ = i;
            }
        }
    }

    {
        long double exactf = func(x__);
        long double L1 = LagranF1(x__, xi, x[i_ + 1]);
        long double true_error1 = L1 - exactf;


        long double min_f = min(func_2pr(xi), func_2pr(x[i_ + 1]));
        long double max_f = max(func_2pr(xi), func_2pr(x[i_ + 1]));

        long double omega1 = (x__ - xi) / (x__ - x[i_ + 1]);
        long double minR1 = min_f * omega1 / 2;
        long double maxR1 = max_f * omega1 / 2;


        long double L2 = LagranF2(x__, xi, x[i_ + 1], x[i_ - 1]);
        long double true_error2 = L2 - exactf;

        long double min_f2 = min(func_3pr(xi), func_3pr(x[i_ + 1]));
        long double max_f2 = max(func_3pr(xi), func_3pr(x[i_ + 1]));

        long double omega2 = (x__ - x[i_ - 1]) * (x__ - xi) * (x__ - x[i_ + 1]);
        long double minR2 = min_f * omega2 / 6;
        long double maxR2 = max_f * omega2 / 6;

        cout << "При x_ = " << x__ << endl;
        cout << "L1= " << L1 << endl;
        cout << "min f''(x) = " << min_f << endl;
        cout << "max f''(x) = " << max_f << endl;
        cout << "min R1 = " << minR1 << endl;
        cout << "max R1 = " << maxR1 << endl;
        cout << "L1(x*) - f(x*) = " << true_error1 << " = R1(x*) " << endl;
        if (minR1 < true_error1 < maxR1) {
            cout << "Неравенство min R1 < R1 < max R1 выполнено " << endl;
        }
        else {
            cout << "Неравенство min R1 < R1 < max R1 не выполнено " << endl;
        }
        if (abs(true_error1) < 0.0001) {
            cout << "Погрешность = " << abs(true_error1) << " подходит " << endl;
        }
        else {
            cout << "Погрешность = " << abs(true_error1) << "  не подходит " << endl;
        }

        cout << endl;

        cout << "L2= " << L2 << endl;
        cout << "min f'''(x) = " << min_f2 << endl;
        cout << "max f'''(x) = " << max_f2 << endl;
        cout << "min R2 = " << minR2 << endl;
        cout << "max R2 = " << maxR2 << endl;
        cout << "L2(x*) - f(x*) = " << true_error2 << " = R2(x*) " << endl;
        if (minR2 < true_error2 < maxR2) {
            cout << "Неравенство min R2 < R2 < max R2 выполнено " << endl;
        }
        else {
            cout << "Неравенство min R2 < R2 < max R2 не выполнено " << endl;
        }
        if (abs(true_error2) < 0.00001) {
            cout << "Погрешность = " << abs(true_error2) << " подходит " << endl;
        }
        else {
            cout << "Погрешность = " << abs(true_error2) << "  не подходит " << endl;
        }
        cout << endl;
    }

    x__ = x_[3];
    i_ = 0;
    xi = 0;

    for (int i = 0; i < 11; ++i) {
        if (x[i] < x__) {
            if (x[i] > xi) {
                xi = x[i];
                i_ = i;
            }
        }
    }

    {
        long double exactf = func(x__);
        long double L1 = LagranF1(x__, xi, x[i_ + 1]);
        long double true_error1 = L1 - exactf;


        long double min_f = min(func_2pr(xi), func_2pr(x[i_ + 1]));
        long double max_f = max(func_2pr(xi), func_2pr(x[i_ + 1]));

        long double omega1 = (x__ - xi) / (x__ - x[i_ + 1]);
        long double minR1 = min_f * omega1 / 2;
        long double maxR1 = max_f * omega1 / 2;


        long double L2 = LagranF2(x__, xi, x[i_ + 1], x[i_ - 1]);
        long double true_error2 = L2 - exactf;

        long double min_f2 = min(func_3pr(xi), func_3pr(x[i_ + 1]));
        long double max_f2 = max(func_3pr(xi), func_3pr(x[i_ + 1]));

        long double omega2 = (x__ - x[i_ - 1]) * (x__ - xi) * (x__ - x[i_ + 1]);
        long double minR2 = min_f * omega2 / 6;
        long double maxR2 = max_f * omega2 / 6;

        cout << "При x_ = " << x__ << endl;
        cout << "L1= " << L1 << endl;
        cout << "min f''(x) = " << min_f << endl;
        cout << "max f''(x) = " << max_f << endl;
        cout << "min R1 = " << minR1 << endl;
        cout << "max R1 = " << maxR1 << endl;
        cout << "L1(x*) - f(x*) = " << true_error1 << " = R1(x*) " << endl;
        if (minR1 < true_error1 < maxR1) {
            cout << "Неравенство min R1 < R1 < max R1 выполнено " << endl;
        }
        else {
            cout << "Неравенство min R1 < R1 < max R1 не выполнено " << endl;
        }
        if (abs(true_error1) < 0.0001) {
            cout << "Погрешность = " << abs(true_error1) << " подходит " << endl;
        }
        else {
            cout << "Погрешность = " << abs(true_error1) << "  не подходит " << endl;
        }

        cout << endl;

        cout << "L2= " << L2 << endl;
        cout << "min f'''(x) = " << min_f2 << endl;
        cout << "max f'''(x) = " << max_f2 << endl;
        cout << "min R2 = " << minR2 << endl;
        cout << "max R2 = " << maxR2 << endl;
        cout << "L2(x*) - f(x*) = " << true_error2 << " = R2(x*) " << endl;
        if (minR2 < true_error2 < maxR2) {
            cout << "Неравенство min R2 < R2 < max R2 выполнено " << endl;
        }
        else {
            cout << "Неравенство min R2 < R2 < max R2 не выполнено " << endl;
        }
        if (abs(true_error2) < 0.00001) {
            cout << "Погрешность = " << abs(true_error2) << " подходит " << endl;
        }
        else {
            cout << "Погрешность = " << abs(true_error2) << "  не подходит " << endl;
        }
        cout << endl;
    }

    return 0;
}
