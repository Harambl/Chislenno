#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

double d, s;

void transtA(vector<vector<double>> &A, int n) {
	vector<vector<double>> At(n + 1, vector<double>(n + 1));
	for (int i = 0; i <= n; ++i) {
		for (int j = 0; j <= n; ++j) {
			At[i][j] = A[j][i];
		}
	}
	A = At;
}

void formW(vector<vector<double>> &N, vector<vector<double>> &A) {
	double prevw = 0;
	for (int i = 0; i <= N.size(); ++i) {
		if (i == 0) {
			N[i][i] = sqrt(A[i][i]);
		}
		else {
			prevw += N[i - 1][i - 1];
		}
		N[i][i] = sqrt(A[i][i] - pow(prevw, 2));
	}
}

vector<double> multiplication(vector<vector<double>> &A, vector<double> &B) {
	vector<double> AB;
	for (int i = 0; i <= A.size(); ++i) {
		for (int j = 0; j <= A.size(); ++j) {
			AB[i] += A[i][j] * B[j];
		}
	}
	return AB;
}



int main()
{
	double alpha = 0.75;
	double beta = 0.7;

	int n;
	cin >> n;
	n -= 1;


	vector<vector<double>> A{ {5.18 + alpha, 1.12, 0.95, 1.32, 0.83},
							  {1.12, 4.28 - alpha, 2.12, 0.57, 0.91},
							  {0.95, 2.12, 6.13 + alpha, 1.29, 1.57},
							  {1.32, 0.57, 1.29, 4.57 - alpha, 1.25},
							  {0.83, 0.91, 1.57, 1.25, 5.21 + alpha} };

	vector<vector<double>> W1(6, vector<double>(6, 0));

	formW(W1, A);

	vector<vector<double>> W2 = W1;

	transtA(W2, W2.size());

	vector<double> B{ 6.19 + beta, 3.21, 4.28 - beta, 6.25, 4.95 + beta};

	vector<vector<double>> A2(n + 1, vector<double>(n + 1));
	for (int i = 0; i <= n; ++i) {
		for (int j = 0; j <= n; ++j) {
			A2[i][j] = 1.0 / (i + j + 1) + 0.1 * 14;
		}
	}

	vector<double> B2(n + 1);

	vector<double> X(n + 1);

	vector<double> Y(n + 1);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			B2[i] += A2[i][j];
		}
	}

	cin >> n;
	n -= 1;


	cout << endl;
}

