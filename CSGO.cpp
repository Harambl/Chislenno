#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <clocale>

using namespace std;


double norm1(const vector<vector<double>>& A);


double vectnorm(const vector<double>& x);

vector<vector<double>> buildB(vector<vector<double>>& A, int l, int m, double c, double s);

double calculate_c(vector<vector<double>>& A, int l, int m);

double calculate_s(vector<vector<double>>& A, int l, int m, double c);


double norm1(const vector<vector<double>>& A) {
    int n = A.size();
    double maxColSum = 0;
    for (int j = 0; j < n; j++) {
        double colSum = 0;
        for (int i = 0; i < n; i++) {
            colSum += abs(A[i][j]);
        }
        maxColSum = max(maxColSum, colSum);
    }
    return maxColSum;
}



double vectnorm(const vector<double>& x) {
    double sum = 0;
    for (int i = 0; i < x.size(); ++i) {
        sum += pow(x[i], 2);
    }
    return sqrt(sum);
}



vector<int> find_pos_of_max(vector<vector<double>>& A) {
    int curr_i = 0;
    int curr_j = 0;
    double max_elem = 0;
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A.size(); ++j) {
            if (A[i][j] > max_elem && i != j) {
                max_elem = A[i][j];
                curr_i = i;
                curr_j = j;
            }
        }
    }
    vector<int> lm= { curr_i, curr_j };
    return lm;
}

vector<vector<double>> buildB(vector<vector<double>>& A, int l, int m, double c, double s) {
    vector<vector<double>> B = A;
    for (int i = 0; i < A.size(); ++i) {
        if (i != l && m != i) {
            B[l][i] = c * A[l][i] + s * A[m][i];
            B[m][i] = -1 * s * A[l][i] + c * A[m][i];
        }
    }
    B[l][l] = pow(c, 2) * A[l][l] + 2 * c * s * A[l][m] + pow(s, 2) * A[m][m];
    B[m][m] = pow(s, 2) * A[l][l] - 2 * c * s * A[l][m] + pow(c, 2) * A[m][m];
    B[l][m] = c * s * (A[m][m] - A[l][l]) + (pow(c, 2) - pow(s, 2)) * A[l][m];
    return B;
}


double calculate_c(vector<vector<double>>& A, int l, int m) {
    double c = sqrt((1 + ((abs(A[l][l] - A[m][m])) / sqrt(pow(2 * A[l][m], 2) + pow(A[l][l] - A[m][m], 2)))) / 2;
    return c;
}

double calculate_s(vector<vector<double>>& A, int l, int m, double c) {
    double sign = 2 * A[l][m] * (A[l][l] - A[m][m]);
    if (sign > 0) {
        sign = 1;
    }
    else {
        sign = -1;
    }
    double s = (sign * abs(2 * A[l][m])) / 2 * c * sqrt(pow(2 * A[l][m], 2) + pow(A[l][l] - A[m][m], 2));
    return s;
}

bool verify(vector<vector<double>>& A, double epsilon) {
    double sum = 0;
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A.size(); ++j) {
            if (i != j) {
                sum += pow(A[i][j], 2);
            }
        }
    }
    if (sqrt(sum) < epsilon) {
        return true;
    }
    else {
        return false;
    }
}

int main() {
    setlocale(LC_ALL, "Russian");

    double alpha = 0.75;
    int delta = 14;

    cout << "=== ВАРИАНТ 14 ===" << endl;
    cout << "alpha = " << alpha <<  ", delta = " << delta << endl;

    int n1 = 5;
    vector<vector<double>> A1 = {
        {5.18 + alpha, 1.12, 0.95, 1.32, 0.83},
        {1.12, 4.28 - alpha, 2.12, 0.57, 0.91},
        {0.95, 2.12, 6.13 + alpha, 1.29, 1.57},
        {1.32, 0.57, 1.29, 4.57 - alpha, 1.25},
        {0.83, 0.91, 1.57, 1.25, 5.21 + alpha}
    };

    vector<vector<double>> A2(11, vector<double>(11));
    for (int i = 0; i < 11; i++) {
        for (int j = 0; j < 11; j++) {
            A2[i][j] = 1.0 / (i + j + 1) + 0.1 * delta;
        }
    }

    vector<double> epsi{ 0.001, 0.00001, 0.0000001 };
    
    cout << "=== Первый Эксперимент" << endl;

    for (double epsilon : epsi) {

        int count = 1;
         
        vector<int> lm = find_pos_of_max(A1);
        int l = lm[0];
        int m = lm[1];

        double c = calculate_c(A1, l, m);
        double s = calculate_s(A1, l, m, c);
        
        vector<vector<double>> B = buildB(A1, l, m, c, s);

        cout << "Итерация номер " << count << endl;
        for (int i = 0; i < xk1.size(); ++i) {
            cout << "xk1[" << i << "] = " << xk1[i] << endl;
        }
        cout << endl;

        while (vectnorm(minusAB(xk, xk1)) >= epsilon) {
            count += 1;
            xk = xk1;
            xk1 = SimpleIterationMethod(xk, A1, I, tau, b1);
            cout << "Итерация номер " << count << endl;
            for (int i = 0; i < xk1.size(); ++i) {
                cout << "xk1[" << i << "] = " << xk1[i] << endl;
            }
            cout << endl;
        }
        cout << "\nТочность (epsilon): " << epsilon << endl;
        cout << "Финальное решение под номером " << count << endl;
        for (int i = 0; i < xk1.size(); ++i) {
            cout << "xk1[" << i << "] = " << xk1[i] << endl;
        }
        cout << endl;
    }

    cout << "\n=== ЗАВЕРШЕНО ===" << endl;
    return 0;
}
