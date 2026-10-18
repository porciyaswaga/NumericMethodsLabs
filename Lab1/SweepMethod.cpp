#include <iostream>
#include <vector>


int n = 5;

std::vector<double> SweetMethod(std::vector<std::vector<double>> A, std::vector<double> d) {
    std::vector<double> x(n);
    std::vector<double> P(n);
    std::vector<double> Q(n);
    
    P[0] = -A[0][1]/A[0][0];
    Q[0] = d[0]/A[0][0];

    for (int i = 1; i < n; i++) {
        P[i] = -A[i][i + 1] / (A[i][i] + A[i][i - 1] * P[i - 1]);
        Q[i] = (d[i] - A[i][i - 1] * Q[i - 1]) / (A[i][i] + A[i][i - 1] * P[i - 1]);
    }

    x[n - 1] = Q[n - 1];

    for (int i = n - 2; i >= 0; i--) {
        x[i] = P[i] * x[i + 1] + Q[i];
    }

    return x;
}

int main() {
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0));
    std::vector<double> d(n);
    std::cout << "Введите матрицу A: \n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cin >> A[i][j];
        }
    }
    std::cout << "Введите вектор d: \n";
    for (int i = 0; i < n; i++) {
        std::cin >> d[i];
    }

    std::vector<double> x = SweetMethod(A, d);

    std::cout << "\nИсходная матрица A:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << A[i][j] << "\t";
        }
        std::cout << "\n";
    }

    std::cout << "\nВектор правой части d:\n";
    for (int i = 0; i < n; i++) {
        std::cout << d[i] << "\n";
    }

    std::cout << "\nРешение СЛАУ:\n";
    for (int i = 0; i < n; i++) {
        std::cout << "x[" << i + 1 << "] = " << x[i] << "\n";
    }
}
