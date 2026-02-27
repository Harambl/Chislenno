#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

double d, s;

void gauss(vector<vector<double>> A, vector<double> B, vector<double> X, int n) {
	for (int k = 0; k <= n; k++) // прямой ход
	{
		for (int j = k + 1; j <= n; j++)
		{
			d = A[j][k] / A[k][k]; // формула (1)

			for (int i = k; i <= n; i++)
			{

				A[j][i] = A[j][i] - d * A[k][i]; // формула (2)

			}

			B[j] = B[j] - d * B[k]; // формула (3)

		}

	}

	for (int k = n; k >= 0; k--) // обратный ход
	{
		d = 0;

		for (int j = k + 1; j <= n; j++)
		{

			s = A[k][j] * X[j]; // формула (4)

			d = d + s; // формула (4)

		}

		X[k] = (B[k] - d) / A[k][k]; // формула (4)

	}

	cout << "Korni sistemy: " << endl;

	for (int i = 0; i <= n; i++)

		cout << "x[" << i << "]=" << X[i] << " " << endl;
}

int main()
{
	int n;
	cin >> n;
	n -= 1;
	vector<vector<double>> A{ {0.411, 0.421, -0.333, 0.313, -0.141, -0.381, 0.245},
							{0.241, 0.705, 0.139, -0.409, 0.321, 0.0625, 0.101},
							{0.123, -0.239, 0.502, 0.901, 0.243, 0.819, 0.321},
							{0.413, 0.309, 0.801, 0.865, 0.423, 0.118, 0.183},
							{0.241, -0.221, -0.243, 0.134, 1.274, 0.712, 0.423},
							{0.281, 0.525, 0.719, 0.118, -0.974, 0.808, 0.923},
							{0.246, -0.301, 0.231, 0.813, -0.702, 1.223, 1.105} };
	vector<vector<double>> A1{ { 5, 4, 7, 5, 6, 7, 5},
								{4, 12, 8, 7, 8, 8, 6},
								{7, 8, 10, 9, 8, 7, 7},
								{5, 7, 9, 11, 9, 7, 5},
								{6, 8, 8, 9, 10, 8, 9},
								{7, 8, 7, 7, 8, 10, 10},
								{5, 6, 7, 5, 9, 10, 10} };

	vector<vector<double>> A2(n + 1, vector<double>(n + 1));
	for (int i = 0; i <= n; ++i) {
		for (int j = 0; j <= n; ++j) {
			A2[i][j] = 1.0 / (i + j + 1); 
		}
	}

	
	vector<vector<double>> A3(n + 1, vector<double>(n + 1));
	for (int i = 0; i <= n; ++i) {
		for (int j = 0; j <= n; ++j) {
			int ii = i + 1, jj = j + 1;
			A3[i][j] = static_cast<double>(min(ii, jj)) / max(ii, jj); 
		}
	}

	vector<double> B{ 0.096, 1.252, 1.024, 1.023, 1.155, 1.937, 1.673 };
	vector<double> B1(n + 1);
	vector<double> B2(n + 1);
	vector<double> B3(n + 1);

	vector<double> X(n + 1);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			B1[i] += A1[i][j];
		}
	}

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			B2[i] += A2[i][j];
		}
	}

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			B3[i] += A3[i][j];
		}
	}

	gauss(A, B, X, n);
	cout << endl;
	gauss(A1, B1, X, n);
	cout << endl;
	gauss(A2, B2, X, n);
	cout << endl;
	gauss(A3, B3, X, n);
}
