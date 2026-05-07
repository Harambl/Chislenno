#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <clocale>

using namespace std;

long double lagrange_poly(long double x, vector<long double>& xs, vector<long double>& ys);
long double lagrange_deriv1(long double x, vector<long double>& xs, vector<long double>& ys);
long double lagrange_deriv2(long double x, vector<long double>& xs, vector<long double>& ys);

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

long double lagrange_poly(long double x, vector<long double>& xs,  vector<long double>& ys) {
    long double res = 0.0;
    int n = xs.size();
    for (int j = 0; j < n; ++j) {
        long double lj = 1.0;
        for (int m = 0; m < n; ++m) {
            if (m == j) continue;
            lj *= (x - xs[m]) / (xs[j] - xs[m]);
        }
        res += ys[j] * lj;
    }
    return res;
}


long double lagrange_deriv1(long double x, vector<long double>& xs, vector<long double>& ys) {
    long double res = 0.0;
    int n = xs.size();
    for (int j = 0; j < n; ++j) {
        long double lj = 1.0;
        long double sum_inv = 0.0;
        for (int m = 0; m < n; ++m) {
            if (m == j) continue;
            long double diff = x - xs[m];
            if (fabs(diff) < 1e-12) continue;
            lj *= (x - xs[m]) / (xs[j] - xs[m]);
            sum_inv += 1.0 / diff;
        }
        res += ys[j] * lj * sum_inv;
    }
    return res;
}


long double lagrange_deriv2(long double x, vector<long double>& xs, vector<long double>& ys) {
    long double res = 0.0;
    int n = xs.size();
    for (int j = 0; j < n; ++j) {
        long double lj = 1.0;
        long double sum_inv = 0.0;
        long double sum_inv2 = 0.0;
        for (int m = 0; m < n; ++m) {
            if (m == j) continue;
            long double diff = x - xs[m];
            if (fabs(diff) < 1e-12) continue;
            lj *= (x - xs[m]) / (xs[j] - xs[m]);
            sum_inv += 1.0 / diff;
            sum_inv2 += 1.0 / (diff * diff);
        }
        res += ys[j] * lj * (sum_inv * sum_inv - sum_inv2);
    }
    return res;
}



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
        if (x[i] < x__ && x[i+1] < x__ && x[i + 2] > x__) {
                xi = x[i];
                i_ = i;
        }
    }

    vector<long double> xs = { xi, x[i_ + 1], x[i_ + 2], x[i_ + 3] };
    vector<long double> ys = { func(xs[0]), func(xs[1]), func(xs[2]), func(xs[3]) };

    cout << "Таблица x_i" << endl;
    int count = 1;
    for (long double xi : x) {
        cout << "x[" << count << "]= " << xi << endl;
        count++;
    }

    cout << endl;

    {
        long double exactf = func(x__);
        long double L3 = lagrange_poly(x__, xs, ys);

        

        long double L3d = lagrange_deriv1(x__, xs, ys);
        long double r1 = abs(func_pr(x__) - L3d);

        long double L3dd = lagrange_deriv2(x__, xs, ys);
        long double r2 = abs(func_2pr(x__) - L3dd);

        cout << "При x_ = " << x__ << endl;
        cout << "L3= " << L3 << endl;
        cout << "L3' = " << L3d << endl;
        cout << "L3'' = " << L3dd << endl;
        cout << "R1 = " << r1 << endl;
        cout << "R2 = " << r2 << endl;

        cout << endl;

    }

    x__ = x_[1];
    i_ = 0;
    xi = 0;

    for (int i = 0; i < 11; ++i) {
        x[i] = a + i * h;
        if (x[i] < x__ && x[i + 1] < x__ && x[i + 2] > x__) {
            xi = x[i];
            i_ = i;
        }
    }

    xs = { xi, x[i_ + 1], x[i_ + 2], x[i_ + 3] };
    ys = { func(xs[0]), func(xs[1]), func(xs[2]), func(xs[3]) };

    {
        long double exactf = func(x__);
        long double L3 = lagrange_poly(x__, xs, ys);



        long double L3d = lagrange_deriv1(x__, xs, ys);
        long double r1 = abs(func_pr(x__) - L3d);

        long double L3dd = lagrange_deriv2(x__, xs, ys);
        long double r2 = abs(func_2pr(x__) - L3dd);

        cout << "При x_ = " << x__ << endl;
        cout << "L3= " << L3 << endl;
        cout << "L3' = " << L3d << endl;
        cout << "L3'' = " << L3dd << endl;
        cout << "R1 = " << r1 << endl;
        cout << "R2 = " << r2 << endl;

        cout << endl;

    }

    x__ = x_[2];
    i_ = 0;
    xi = 0;



    xs = { x[7], x[8], x[9], x[10]};
    ys = { func(xs[0]), func(xs[1]), func(xs[2]), func(xs[3]) };

    {
        long double exactf = func(x__);
        long double L3 = lagrange_poly(x__, xs, ys);



        long double L3d = lagrange_deriv1(x__, xs, ys);
        long double r1 = abs(func_pr(x__) - L3d);

        long double L3dd = lagrange_deriv2(x__, xs, ys);
        long double r2 = abs(func_2pr(x__) - L3dd);

        cout << "При x_ = " << x__ << endl;
        cout << "L3= " << L3 << endl;
        cout << "L3' = " << L3d << endl;
        cout << "L3'' = " << L3dd << endl;
        cout << "R1 = " << r1 << endl;
        cout << "R2 = " << r2 << endl;

        cout << endl;

    }

    x__ = x_[3];
    i_ = 0;
    xi = 0;

    for (int i = 0; i < 11; ++i) {
        x[i] = a + i * h;
        if (x[i] < x__ && x[i + 1] < x__ && x[i + 2] > x__) {
            xi = x[i];
            i_ = i;
        }
    }

    xs = { xi, x[i_ + 1], x[i_ + 2], x[i_ + 3] };
    ys = { func(xs[0]), func(xs[1]), func(xs[2]), func(xs[3]) };

    {
        long double exactf = func(x__);
        long double L3 = lagrange_poly(x__, xs, ys);



        long double L3d = lagrange_deriv1(x__, xs, ys);
        long double r1 = abs(func_pr(x__) - L3d);

        long double L3dd = lagrange_deriv2(x__, xs, ys);
        long double r2 = abs(func_2pr(x__) - L3dd);

        cout << "При x_ = " << x__ << endl;
        cout << "L3= " << L3 << endl;
        cout << "L3' = " << L3d << endl;
        cout << "L3'' = " << L3dd << endl;
        cout << "R1 = " << r1 << endl;
        cout << "R2 = " << r2 << endl;

        cout << endl;

    }

    return 0;
}
