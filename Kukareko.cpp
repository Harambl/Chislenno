#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;

vector<double> IterationMethod(vector<double>& xk, vector<vector<double>>& A, double tau, vector<double> f) {

}


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

vector<double> multiplication(vector<vector<double>>& A, vector<double>& B) {
    vector<double> AB(B.size());
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A.size(); ++j) {
            AB[i] += A[i][j] * B[j];
        }
    }
    return AB;
}

void multiplicationConst(vector<vector<double>>& A, double con) {
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A.size(); ++j) {
            A[i][j] *= con;
        }
    }
}

vector<double> minusAB(vector<double>& A, vector<double>& B) {
    vector<double> MAB(A.size());
    for (int i = 0; i < A.size(); ++i) {
        MAB[i] = A[i] - B[i];
    }
    return MAB;
}

int main() {
    setlocale(LC_ALL, "Russian");

    int n1 = 4;
    vector<vector<double>> A1 = {
        {0.85, -0.05, 0.08, -0.14},
        {-0.32, 1.13, 0.12, -0.11},
        {-0.17, -0.06, 1.08, -0.12},
        {-0.21, 0.16, -0.36, 1.00}
    };

    vector<double> b1 = { -0.4800, 1.2400, 1.1500, -0.8800 };

    vector<double> x0(4, 0);

    vector<vector<double>> I(4, vector<double>(4, 0));

    for (int i = 0; i < I.size(); ++i) {
        I[i][i] = 1;
    }

    cout << "=== Первый Эксперимент" << endl;

    vector<double> epsi{ 0.001, 0.00001, 0.0000001 };

    double tau = 2 / norm1(A1);

    for (double epsilon : epsi) {

    }


    cout << "\n=== ЗАВЕРШЕНО ===" << endl;
    return 0;
}
