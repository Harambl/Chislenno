#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <clocale>

using namespace std;


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

long double center(long double a, long double b) {
    return (a + b) / 2;
}

long double Gauss(long double a, long double b, int n) {
    vector<long double> t(n, 0);
    for (int i = 0; i < n; ++i) {
        if (i % 2 == 1) {
            t[i] = -0.577350;
        }
        if (i % 2 == 0) {
            t[i] = 0.577350;
        }
    }
    vector<long double> x(n, 0);
    for (int i = 0; i < n; ++i) {
        x[i] = (b + a) / 2 + (b - a) * t[i] / 2;
    }
    long double sum = 0;
    for (int i = 0; i < n; ++i) {
        sum += func(x[i]);
    }
    long double INTE = ((b - a) / 2) * sum;
    return INTE;
}


int main() {
    setlocale(LC_ALL, "Russian");

    long double a = 0.6;
    long double b = 1.1;

    long double epsilon = 0.000001;

    vector<long double> x_ = { 0.88, 0.63, 1.08, 0.83 };

    long double h = (b - a) / 10;

    int n = 2;

    vector<long double> x(11, 0);
    for (int i = 0; i < 11; ++i) {
        x[i] = a + i * h;
    }
    long double INTE1 = Gauss(a, b, n);
    cout << "Вариант 14(по формуле Гаусса для n = 2)" << endl;
    cout << "Интеграл = " << INTE1 << endl;
    cout << "Берем эпсилон = " << epsilon << endl;
    long double cent1 = center(a, b);
    long double cent2 = center(a, b);
    long double Gauss1 = Gauss(a, cent1, n);
    long double Gauss2 = Gauss(cent2, b, n);
    vector<long double> inte2 = {Gauss1, Gauss2};
    long double INTE2 = accumulate(inte2.begin(), inte2.end(), 0.0);
    cout << "Интеграл 2n = " << INTE2 << endl;
    cout << endl;
    if (abs(INTE1 - INTE2) <= epsilon) {
        cout << "Uspex" << endl;
    }
    else {
        long double c = cent1;
        while (abs(INTE1 - INTE2) >= epsilon) {
            cent1 = center(a, cent1);
            cent2 = center(cent1, b);
            INTE1 = INTE2;
            INTE2 = 0.0;
            inte2.clear();
            Gauss1 = Gauss(a, cent1, n);
            Gauss2 = Gauss(cent1, c, n);
            inte2.push_back(Gauss1);
            inte2.push_back(Gauss2);
            Gauss1 = Gauss(c, cent2, n);
            Gauss2 = Gauss(cent2, b, n);
            inte2.push_back(Gauss1);
            inte2.push_back(Gauss2);
            INTE2 = accumulate(inte2.begin(), inte2.end(), 0.0);
        }
        cout << "Интеграл 2n = " << INTE1 << endl;
        cout << "Интеграл 4n = " << INTE2 << endl;
        cout << endl;
        cout << "Uspex 2" << endl;
    }


    cout << endl;

    return 0;
}
