#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <clocale>

using namespace std;

vector<double> multiplicationConst( vector<double>& A, double con);

vector<vector<double>> multiplicationConst( vector<vector<double>>& A, double con);

vector<double> multiplication( vector<vector<double>>& A,  vector<double>& B);

vector<vector<double>> minusABe( vector<vector<double>>& A,  vector<vector<double>>& B);

vector<double> minusAB( vector<double>& A,  vector<double>& B);

vector<double> additive( vector<double>& A,  vector<double>& B);

double norm1( const vector<vector<double>>& A);

double vectnorm( const vector<double>& x);

vector<double> IterationMethod( vector<double>& xk, vector<vector<double>>& A, vector<vector<double>>& I, double tau,  vector<double>& f);



vector<double> IterationMethod(vector<double>& xk, vector<vector<double>>& A, vector<vector<double>>& I,  double tau, vector<double>& f) {
    vector<double> ftau = multiplicationConst(f, tau);
    vector<vector<double>> Atau = multiplicationConst(A, tau);
    vector<vector<double>> AItau = minusABe(I, Atau);
    vector<double> Axk = multiplication(AItau, xk);
    vector<double> xk1 = additive(Axk, ftau);

    return xk1;
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



double vectnorm(const vector<double>& x) {
    double sum = 0;
    for (int i = 0; i < x.size(); ++i) {
        sum += pow(x[i], 2);
    }
    return sqrt(sum);
}



vector<double> multiplication(vector<vector<double>>& A, vector<double>& B) {
    vector<double> AB(B.size(), 0);
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A.size(); ++j) {
            AB[i] += A[i][j] * B[j];
        }
    }
    return AB;
}



vector<vector<double>> multiplicationConst(vector<vector<double>>& A, double con) {
    vector<vector<double>> A_new = A;
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A.size(); ++j) {
            A_new[i][j] *= con;
        }
    }
    return A_new;
}


vector<double> multiplicationConst(vector<double>& A, double con) {
    vector<double> A_new = A;
    for (int i = 0; i < A.size(); ++i) {
            A_new[i] *= con;
    }
    return A_new;
}



vector<double> minusAB(vector<double>& A, vector<double>& B) {
    vector<double> MAB(A.size());
    for (int i = 0; i < A.size(); ++i) {
        MAB[i] = A[i] - B[i];
    }
    return MAB;
}




vector<vector<double>> minusABe(vector<vector<double>>& A, vector<vector<double>>& B) {
    vector<vector<double>> MABe(A.size(), vector<double>(A.size(), 0));
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A.size(); ++j) {
            MABe[i][j] = A[i][j] - B[i][j];
        }
    }
    return MABe;
}



vector<double> additive(vector<double>& A, vector<double>& B) {
    vector<double> A_new = A;
    for (int i = 0; i < A.size(); ++i) {
        A_new[i] += B[i];
    }
    return A_new;
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

        int count = 1;
        vector<double> xk = x0;
        vector<double> xk1 = IterationMethod(xk, A1, I, tau, b1);

        cout << "Итерация номер " << count << endl;
        for (int i = 0; i < xk1.size(); ++i) {
            cout << "xk1[" << i << "] = " << xk1[i] << endl;
        }
        cout << endl;

        while (vectnorm(minusAB(xk, xk1)) >= epsilon) {
            count += 1;
            xk = xk1;
            xk1 = IterationMethod(xk, A1, I, tau, b1);
            cout << "Итерация номер " << count << endl;
            for (int i = 0; i < xk1.size(); ++i) {
                cout << "xk1[" << i << "] = " << xk1[i] << endl;
            }
            cout << endl;
        }
        cout << "\nТочность (epsilon): " << epsilon << endl;
        cout << "Финальное решение под номером " << count << endl;
        for (int i = 0; i < xk1.size(); ++i) {
            cout << "xk1[" << i << "] = "  << xk1[i] << endl;
        }
        cout << endl;
    }


    cout << "=== Второй Эксперимент" << endl;

    tau = 0.5;
    double epsilon = epsi[1];
    int count = 1;

    cout << "\nТочность (epsilon): " << epsilon << endl;

    for (int i = 0; i < 9; ++i) {
        cout << "Попытка номер " << i << endl;
        vector<double> xk = x0;
        vector<double> xk1 = IterationMethod(xk, A1, I, tau, b1);

        while (vectnorm(minusAB(xk, xk1)) >= epsilon) {
            count += 1;
            xk = xk1;
            xk1 = IterationMethod(xk, A1, I, tau, b1);
        }

        cout << "Финальное решение под номером " << count << endl;
        for (int i = 0; i < xk1.size(); ++i) {
            cout << "xk1[" << i << "] = " << xk1[i] << endl;
        }

        tau -= 0.05;
    }

    epsilon = epsi[2];
    tau = 0.5;

    cout << "\nТочность (epsilon): " << epsilon << endl;

    for (int i = 0; i < 9; ++i) {
        cout << "Попытка номер " << i << endl;
        vector<double> xk = x0;
        vector<double> xk1 = IterationMethod(xk, A1, I, tau, b1);

        while (vectnorm(minusAB(xk, xk1)) >= epsilon) {
            count += 1;
            xk = xk1;
            xk1 = IterationMethod(xk, A1, I, tau, b1);
        }

        cout << "Финальное решение под номером " << count << endl;
        for (int i = 0; i < xk1.size(); ++i) {
            cout << "xk1[" << i << "] = " << xk1[i] << endl;
        }

        tau -= 0.05;
    }

    cout << "\n=== ЗАВЕРШЕНО ===" << endl;
    return 0;
}
