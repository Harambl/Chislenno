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
        if (i % 2 == 0) {
            t[i] = 0.577350;
        }
        if (i % 2 == 1) {
            t[i] = -0.577350;
        }
    }

    vector<long double> x(n + 2, 0);
    x.push_back(a);
    for (int i = 0; i < n; ++i) {
        x[i] = (b + a) / 2 + (b - a) * t[i] / 2;
        cout << x[i] << endl;
    }
    x.push_back(b);
    long double h = 1;

    vector<long double> M{func_2pr(a), func_2pr(x[0]), func_2pr(x[1]), func_2pr(b)};
    for (int i = 0; i < M.size(); ++i) {
        cout << M[i] << 'j' << endl;
    }

    long double sum = 0;
    for (int i = 0; i < M.size() - 1; ++i) {
         sum += (func(x[i]) + func(x[i+1]))/2 - (pow(h, 3) / 24) * (M[i] + M[i+1]);
         cout << sum << 's' << endl;
    }
    
    long double inte1 = (func(a) + func(b)) / 2 - (pow(h, 3) / 24) * (func_2pr(a) + func_2pr(b));

    long double inte2 = (pow(b, 2) / 2) - (pow(a, 2) / 2) + (cos(b) - cos(a));

    cout << "Вариант 14 (сплайн через моменты, задача интегрирования, n = 2, 2 краевое условие)" << endl;
    cout << "Шаг равен " << h << endl;
    cout << "Интеграл 1 равен " << inte1 << endl;
    cout << "Интеграл 2 равен " << inte2 << endl;

    long double inte3 = inte1 - inte2;

    cout << "Разность между теоретическим и сплайновым интегралом = " << inte3 << " погрешность" << endl;

    cout << endl;

    return 0;
}
