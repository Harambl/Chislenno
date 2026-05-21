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



int main() {
    setlocale(LC_ALL, "Russian");

    long double a = 0.6;
    long double b = 1.1;

    int n = 2;

     vector<long double> t(n, 0);
    for (int i = 0; i < n; ++i) {
        if (i % 2 == 1) {
            t[i] = 0.577350;
        }
        if (i % 2 == 0) {
            t[i] = -0.577350;
        }
    }

    vector<long double> x(n, 0);
    for (int i = 0; i < n; ++i) {
        x[i] = (b + a) / 2 + (b - a) * t[i] / 2;
    }
    vector<long double> h = { x[1] - x[0] };

    long double inte = h[0] * (func(x[0]) + func(x[1])) - (pow(h[0], 3) / 24) * (a - b);

    long double inte2 = (pow(x[0], 2) / 2) - (pow(x[1], 2) / 2) + (cos(x[0]) - cos(x[1]));


    cout << "Вариант 14 (сплайн через моменты, задача интегрирования, n = 2, 2 краевое условие)" << endl;
    cout << "Шаг равен " << h[0] << endl;
    cout << "Интеграл 1 равен " << inte << endl;
    cout << "Интеграл 2 равен " << inte2 << endl;

    cout << endl;

    return 0;
}
