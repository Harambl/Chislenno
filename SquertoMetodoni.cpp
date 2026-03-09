#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;

void choleskyDecomposition(const vector<vector<double>>& A, vector<vector<double>>& W, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0;

            if (j == i) {
                for (int k = 0; k < j; k++) {
                    sum += pow(W[j][k], 2);
                }
                W[j][j] = sqrt(A[j][j] - sum);
            }
            else {       
                for (int k = 0; k < j; k++) {
                    sum += W[i][k] * W[j][k];
                }
                W[i][j] = (A[i][j] - sum) / W[j][j];
            }
        }
    }
}

void forwardSubstitution(const vector<vector<double>>& W, const vector<double>& b, vector<double>& y, int n) {
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += W[i][j] * y[j];
        }
        y[i] = (b[i] - sum) / W[i][i];
    }
}

void backwardSubstitution(const vector<vector<double>>& W, const vector<double>& y, vector<double>& x, int n) {
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += W[j][i] * x[j];  // W^T[i][j] = W[j][i]
        }
        x[i] = (y[i] - sum) / W[i][i];
    }
}


double calculateConditionNumber(const vector<vector<double>>& A, int n) {
    double maxVal = abs(A[0][0]);
    double minVal = abs(A[0][0]);

    for (int i = 0; i < n; i++) {
        double rowSum = 0;
        for (int j = 0; j < n; j++) {
            rowSum += abs(A[i][j]);
        }
        if (rowSum > maxVal) maxVal = rowSum;
        if (rowSum < minVal) minVal = rowSum;
    }

    return maxVal / minVal;
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

vector<double> minusAB(vector<double>& A, vector<double>& B) {
    vector<double> MAB(A.size());
    for (int i = 0; i < A.size(); ++i) {
        MAB[i] = A[i] - B[i];
    }
    return MAB;
}

int main() {
    setlocale(LC_ALL, "Russian");

    double alpha = 0.75;   
    double beta = 0.7;     
    int delta = 14;        

    cout << "=== ВАРИАНТ 14 ===" << endl;
    cout << "alpha = " << alpha << ", beta = " << beta << ", delta = " << delta << endl;

    int n1 = 5;
    vector<vector<double>> A1 = {
        {5.18 + alpha, 1.12, 0.95, 1.32, 0.83},
        {1.12, 4.28 - alpha, 2.12, 0.57, 0.91},
        {0.95, 2.12, 6.13 + alpha, 1.29, 1.57},
        {1.32, 0.57, 1.29, 4.57 - alpha, 1.25},
        {0.83, 0.91, 1.57, 1.25, 5.21 + alpha}
    };

    vector<double> b1 = { 6.19 + beta, 3.21, 4.28 - beta, 6.25, 4.95 + beta };

    cout << "\n=== ПЕРВАЯ СИСТЕМА ===" << endl;

    vector<vector<double>> W1(n1, vector<double>(n1, 0));
    choleskyDecomposition(A1, W1, n1);

    cout << "\nМатрица Холецкого W1:" << endl;

    vector<double> y1(n1);
    forwardSubstitution(W1, b1, y1, n1);

    vector<double> x1(n1);
    backwardSubstitution(W1, y1, x1, n1);

    cout << "\nРешение первой системы:" << endl;
    for (int i = 0; i < n1; i++) {
        cout << "x1[" << i << "] = " << fixed << setprecision(6) << x1[i] << endl;
    }

    double cond1 = calculateConditionNumber(A1, n1);
    cout << "\nЧисло обусловленности A1: " << fixed << setprecision(4) << cond1 << endl;

    
    cout << "\n\n=== ВТОРАЯ СИСТЕМА ===" << endl;

    vector<int> sizes = { 5, 7, 9, 11 };

    for (int n : sizes) {
        cout << "\n--- Размер n = " << n << " ---" << endl;

        vector<vector<double>> A2(n, vector<double>(n));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A2[i][j] = 1.0 / (i + j + 1) + 0.1 * delta;
            }
        }

        vector<double> b2(n, 0);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                b2[i] += A2[i][j];
            }
        }

        vector<vector<double>> W2(n, vector<double>(n, 0));
        choleskyDecomposition(A2, W2, n);

        vector<double> y2(n);
        forwardSubstitution(W2, b2, y2, n);

        vector<double> x2(n);
        backwardSubstitution(W2, y2, x2, n);

        double cond = calculateConditionNumber(A2, n);

        cout << "Число обусловленности: " << fixed << setprecision(4) << cond << endl;
        cout << "Первые значения решения:" << endl;
        for (int i = 0; i < n; i++) {
            cout << "x[" << i << "] = " << fixed << setprecision(6) << x2[i] << endl;
        }

        vector<double> Prov1 = multiplication(A2, x2);
        vector<double> Prov2 = minusAB(Prov1, b2);
        cout << "Проверка Ax - B";
        for (int i = 0; i < n; i++) {
            cout << "Prov[" << i << "] = " << fixed << setprecision(6) << Prov2[i] << endl;
        }
        
    }

    cout << "\n=== ЗАВЕРШЕНО ===" << endl;
    return 0;
}
