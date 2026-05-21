#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <clocale>

using namespace std; // polozhitelnye korni


long double func(long double x) {
    return 2*x - log10(x) - 7;
}

long double func_pr(long double x) {
    return 2 - 1/(x * log(10));
}

long double func_2pr(long double x) {
    return 1/(pow(x, 2) * log(10));
}

bool verify(long double x_, long double xk, long double epsilon) {
    return abs(x_ - xk) < epsilon;
}

long double chord(long double x_, long double x) {
    long double xn = x - (func(x)*(x_ - x)) / (func(x_) - func(x));
    return xn;
}

long double casat(long double xn_1) {
    long double xn = xn_1 - func(xn_1) / func_pr(xn_1);
    return xn;
}


int main() {
    setlocale(LC_ALL, "Russian");

    long double epsilon = 0.000000000001;

    long double a = 1.0;
    long double b = 5.0;

    long double x_n = b - func(b) / func_pr(b);
    long double xn = a - func(a)*(x_n - a) / func(x_n) - func(a);

    long double x1 = xn;
    long double x1_ = x_n;

    long double k = xn - x_n;

    cout << x_n << ' ' << xn << ' ' << k << endl;

    cout << func(a) << ' ' << func_2pr(a) << endl;
    cout << func(b) << ' ' << func_2pr(b) << endl;

    while (abs(k) > epsilon) {
        x_n = casat(x1_);
        xn = chord(x1_, x1);
        x1_ = x_n;
        x1 = xn;
        k = xn - x_n;
    }

    cout << "Yspex  " << xn << ' ' << x_n << endl;
    
    return 0;
}
