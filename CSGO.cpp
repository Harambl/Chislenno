#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <clocale>

using namespace std;


long double norm1(const vector<vector<long double>>& A);


long double vectnorm(const vector<long double>& x);

vector<vector<long double>> buildB(vector<vector<long double>>& A, int l, int m, long double c, long double s);

long double calculate_c(vector<vector<long double>>& A, int l, int m);

long double calculate_s(vector<vector<long double>>& A, int l, int m, long double c);


long double norm1(const vector<vector<long double>>& A) {
    int n = A.size();
    long double maxColSum = 0;
    for (int j = 0; j < n; j++) {
        long double colSum = 0;
        for (int i = 0; i < n; i++) {
            colSum += abs(A[i][j]);
        }
        maxColSum = max(maxColSum, colSum);
    }
    return maxColSum;
}



long double vectnorm(const vector<long double>& x) {
    long double sum = 0;
    for (int i = 0; i < x.size(); ++i) {
        sum += pow(x[i], 2);
    }
    return sqrt(sum);
}



vector<int> find_pos_of_max(vector<vector<long double>>& A) {
    int curr_i = 0;
    int curr_j = 1;
    long double max_elem = 0;
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A.size(); ++j) {
            if (std::abs(A[i][j]) > std::abs(max_elem) && i != j) {
                max_elem = A[i][j];
                curr_i = i;
                curr_j = j;
            }
        }
    }
    vector<int> lm= { curr_i, curr_j };
    return lm;
}

vector<vector<long double>> buildB(vector<vector<long double>>& A, int l, int m, long double c, long double s) {
    vector<vector<long double>> B = A;
    for (int i = 0; i < A.size(); ++i) {
        if (i != l && m != i) {
            B[l][i] = c * A[l][i] + s * A[m][i];
            B[i][l] = B[l][i];
            B[m][i] = -1 * s * A[l][i] + c * A[m][i];
            B[i][m] = B[m][i];
        }
    }
    B[l][l] = pow(c, 2) * A[l][l] + 2 * c * s * A[l][m] + pow(s, 2) * A[m][m];
    B[m][m] = pow(s, 2) * A[l][l] - 2 * c * s * A[l][m] + pow(c, 2) * A[m][m];
    B[l][m] = c * s * (A[m][m] - A[l][l]) + (pow(c, 2) - pow(s, 2)) * A[l][m];
    B[m][l] = c * s * (A[m][m] - A[l][l]) + (pow(c, 2) - pow(s, 2)) * A[l][m];
    return B;
}


long double calculate_c(vector<vector<long double>>& A, int l, int m) {
    long double c = sqrt((1 + ((std::abs(A[l][l] - A[m][m])) / sqrt(pow(2 * A[l][m], 2) + pow(A[l][l] - A[m][m], 2)))) / 2);
    //cout << c << endl;
    return c;
}

long double calculate_s(vector<vector<long double>>& A, int l, int m, long double c) {
    long double sign = 2 * A[l][m] * (A[l][l] - A[m][m]);
    if (sign > 0) {
        sign = 1.0;
    }
    else {
        sign = -1.0;
    }
    long double s = (sign * std::abs(2 * A[l][m])) / (2 * c * sqrt(pow(2 * A[l][m], 2) + pow(A[l][l] - A[m][m], 2)));
   // cout << s << endl;
    return s;
}

bool verify(vector<vector<long double>>& A, long double epsilon) {
    long double sum = 0;
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A.size(); ++j) {
            if (i != j) {
                sum += pow(A[i][j], 2);
            }
        }
    }
    return sqrt(sum) < epsilon;

}

int main() {
    setlocale(LC_ALL, "Russian");

    long double alpha = 0.75;
    int delta = 14;

    cout << "=== ВАРИАНТ 14 ===" << endl;
    cout << "alpha = " << alpha <<  ", delta = " << delta << endl;

    int n1 = 5;
    vector<vector<long double>> A1 = {
        {5.18 + alpha, 1.12, 0.95, 1.32, 0.83},
        {1.12, 4.28 - alpha, 2.12, 0.57, 0.91},
        {0.95, 2.12, 6.13 + alpha, 1.29, 1.57},
        {1.32, 0.57, 1.29, 4.57 - alpha, 1.25},
        {0.83, 0.91, 1.57, 1.25, 5.21 + alpha}
    };


    vector<vector<long double>> A2(11, vector<long double>(11));
    for (int i = 0; i < 11; i++) {
        for (int j = 0; j < 11; j++) {
            A2[i][j] = 1.0 / (i + j + 1) + 0.1 * delta;
        }
    }

    vector<long double> epsi{ 0.001, 0.00001, 0.0000001 };
    
    cout << "=== Первый Эксперимент" << endl;

    for (long double epsilon : epsi) {

        int count = 1;
         
        vector<int> lm = find_pos_of_max(A1);
        int l = lm[0];
        int m = lm[1];

        long double c = calculate_c(A1, l, m);
        long double s = calculate_s(A1, l, m, c);
        
        vector<vector<long double>> B = buildB(A1, l, m, c, s);

        cout << "Итерация номер " << count << endl;
        for (int i = 0; i < B.size(); ++i) {

            cout << "Собственные значения B[" << i << "][" << i << "] = " << B[i][i] << endl;
        }
        cout << endl;


        while (!verify(B, epsilon)  ) {
            count += 1;
            lm = find_pos_of_max(B);
            int l = lm[0];
            int m = lm[1];

            long double c = calculate_c(B, l, m);
            long double s = calculate_s(B, l, m, c);

            B = buildB(B, l, m, c, s);

            cout << "Итерация номер " << count << endl;
            for (int i = 0; i < B.size(); ++i) {
                cout << "Собственные значения B[" << i << "][" << i << "] = " << B[i][i] << endl;
            }
            cout << endl;
        }
        cout << "\nТочность (epsilon): " << epsilon << endl;
        cout << "Финальное решение под номером " << count << endl;
        for (int i = 0; i < B.size(); ++i) {
            cout << "Собственные значения B[" << i << "][" << i << "] = " << B[i][i] << endl;
        }
        cout << endl;
    }


    long double epsilon = epsi[2];

    cout << "=== Второй Эксперимент" << endl;

    {
        int count = 1;

        vector<int> lm = find_pos_of_max(A1);
        int l = lm[0];
        int m = lm[1];

        long double c = calculate_c(A1, l, m);
        long double s = calculate_s(A1, l, m, c);

        vector<vector<long double>> B = buildB(A2, l, m, c, s);

        cout << "Итерация номер " << count << endl;
        for (int i = 0; i < B.size(); ++i) {
            cout << "Собственные значения B[" << i << "][" << i << "] = " << B[i][i] << endl;
        }
        cout << endl;

        while (!verify(B, epsilon)) {
            count += 1;
            lm = find_pos_of_max(B);
            int l = lm[0];
            int m = lm[1];

            long double c = calculate_c(B, l, m);
            long double s = calculate_s(B, l, m, c);

            B = buildB(B, l, m, c, s);
            cout << "Итерация номер " << count << endl;
            for (int i = 0; i < B.size(); ++i) {
                cout << "Собственные значения B[" << i << "][" << i << "] = " << B[i][i] << endl;
            }
            cout << endl;
        }
        cout << "\nТочность (epsilon): " << epsilon << endl;
        cout << "Финальное решение под номером " << count << endl;
        for (int i = 0; i < B.size(); ++i) {
            cout << "Собственные значения B[" << i << "][" << i << "] = " << B[i][i] << endl;
        }
        cout << endl;
    }
    cout << "\n=== ЗАВЕРШЕНО ===" << endl;
    return 0;
}
