#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <clocale>

long double func(long double x) {
    return x - sin(x);
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
    cout << "При x_ = " << x__ << endl;
    cout << "L1= " << LagranF1(x__, xi, x[i_ + 1]) << endl;
    cout << "L2= " << LagranF2(x__, xi, x[i_ + 1], x[i_ - 1]) << endl;
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
        cout << "При x_ = " << x__ << endl;
        cout << "L1= " << LagranF1(x__, xi, x[i_ + 1]) << endl;
        cout << "L2= " << LagranF2(x__, xi, x[i_ + 1], x[i_ + 2]) << endl;
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
        cout << "При x_ = " << x__ << endl;
        cout << "L1= " << LagranF1(x__, xi, x[i_ + 1]) << endl;
        cout << "L2= " << LagranF2(x__, xi, x[i_ + 1], x[i_ - 1]) << endl;
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
        cout << "При x_ = " << x__ << endl;
        cout << "L1= " << LagranF1(x__, xi, x[i_ + 1]) << endl;
        cout << "L2= " << LagranF2(x__, xi, x[i_ + 1], x[i_ - 1]) << endl;
        cout << endl;
    }

    return 0;
}
