#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <clocale>

double func(double x) {
    return x - sin(x);
}

double LagranF1(double x_, double xi, double xi_1) {
    return func(xi) * (x_ - xi_1) / (xi - xi_1) + func(xi_1) * (x_ - xi) / (xi_1 - xi);
}

double LagranF2(double x_, double xi, double xi_1, double x1) {
    return func(x1) * ((x_ - xi) * (x_ - xi_1)) / ((x1 - xi) * (x1 - xi_1)) 
        + func(xi) * ((x_ - x1) * (x_ - xi - 1)) / ((xi - x1) * (xi - xi_1)) 
        + func(xi_1) * ((x_ - xi_1) * (x_ - xi)) / ((xi_1 - x1) * (xi_1 - xi));
}



using namespace std;

int main() {
    setlocale(LC_ALL, "Russian");

    double a = 0.6;
    double b = 1.1;

    vector<double> x_ = {0.88, 0.63, 1.08, 0.83};

    double h = (b - a) / 10;

    double x__ = x_[0];
    int i_ = 0;
    double xi = 0;

    vector<double> x (11, 0);
    for (int i = 0; i < 11; ++i) {
        x[i] = a + i * h;
        if(x[i] < x__){
            if (x[i] > xi) {
                xi = x[i];
                i_ = i;
            }
        }
    }

    cout << "Таблица x_i" << endl;
    int count = 1;
    for (double xi : x) {
        cout << "x[" << count << "]= " << xi << endl;
        count++;
    }

    for (double x : x_) {

    }


    return 0;
}
