#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include <stdexcept>
#include <tuple>


int n = 4;

double Determinator(std::vector<std::vector<double>> U, std::vector<std::vector<int>> P, int countOfSwaps) {
    double determinator = 1.0;

    if (countOfSwaps % 2 != 0) {
        determinator = -1.0;
    }

    for (int i = 0; i < n; i++) {
        determinator *= U[i][i];
    }
    return determinator;
}

std::vector<double> Solve(
    std::vector<std::vector<double>> L,
    std::vector<std::vector<double>> U,
    std::vector<std::vector<int>> P,
    std::vector<double> b) {

    std::vector<double> Pb(n, 0.0);
    std::vector<double> y(n, 0.0);
    std::vector<double> x(n, 0.0);

    for (int i = 0; i < n; i++) {
        Pb[i] = 0.0;
        for (int j = 0; j < n; j++) {
            Pb[i] += P[i][j] * b[j];
        }
    }

    for (int i = 0; i < n; i++) {
        y[i] = Pb[i];
        for (int j = 0; j < i; j++) {
            y[i] -= L[i][j] * y[j];
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];

        for (int j = i + 1; j < n; j++) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }

    return x;
}

std::vector<std::vector<double>> InverseMatrix(
    std::vector<std::vector<double>> L,
    std::vector<std::vector<double>> U,
    std::vector<std::vector<int>> P) {

    std::vector<std::vector<double>> E(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        E[i][i] = 1;
    }

    std::vector<double> z;
    std::vector<std::vector<double>> inverse(n, std::vector<double>(n, 0.0));

    for (int k = 0; k < n; k++) {
        std::vector<double> b(n);

        for (int i = 0; i < n; i++) {
            b[i] = E[i][k];
        }

        z = Solve(L, U, P, b);

        for (int i = 0; i < n; i++) {
            inverse[i][k] = z[i];
        }
    }

    return inverse;
}

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<int>>, int> LUDecomposition(std::vector<std::vector<double>> A) {
    
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> U = A;
    std::vector<std::vector<int>> P(n, std::vector<int>(n, 0));

    int countOfSwaps = 0;

    for (int i = 0; i < n; i++) {
        L[i][i] = 1.0;
        P[i][i] = 1;
    }

    for (int k = 0; k < n; k++) {
        int maxRow = k;
        double maxChoose = std::abs(U[k][k]);

        for (int i = k + 1; i < n; i++) {
            if (std::abs(U[i][k]) > maxChoose) {
                maxChoose = std::abs(U[i][k]);
                maxRow = i;
            }
        }
        
        if (maxRow != k) {
            std::swap(U[k], U[maxRow]);
            std::swap(P[k], P[maxRow]);

            for (int j = 0; j < k; j++) {
                std::swap(L[k][j], L[maxRow][j]);
            }
            countOfSwaps += 1;
        }

        if (U[k][k] == 0.0) {
            throw std::runtime_error("LU decomposition failed: zero pivot");
        }

        for (int i = k + 1; i < n; i++) {
            L[i][k] = U[i][k] / U[k][k];
            for (int j = k; j < n; j++) {
                U[i][j] = U[i][j] - L[i][k] * U[k][j];
            }
        }
    }
    return {L, U, P, countOfSwaps};
}

int main() {
    std::vector<std::vector<double>> A(n, std::vector<double>(n));
    std::vector<double> b(n);
    std::cout << "Введите матрицу A: \n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cin >> A[i][j];
        }
    }
    std::cout << "Введите вектор b: \n";
    for (int i = 0; i < n; i++) {
        std::cin >> b[i];
    }

    auto [L, U, P, countOfSwaps] = LUDecomposition(A);

    double determinantor = Determinator(U, P, countOfSwaps);
    std::vector<std::vector<double>> Inverse = InverseMatrix(L, U, P);
    std::vector<double> x = Solve(L, U, P, b);

    std::cout << "\nИсходная матрица A:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << A[i][j] << "\t";
        }
        std::cout << "\n";
    }

    std::cout << "\nВектор правой части b:\n";
    for (int i = 0; i < n; i++) {
        std::cout << b[i] << "\n";
    }

    std::cout << "\nРешение СЛАУ:\n";
    for (int i = 0; i < n; i++) {
        std::cout << "x[" << i + 1 << "] = " << x[i] << "\n";
    }

    std::cout << "\nОбратная матрица A^(-1):\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << Inverse[i][j] << "\t";
        }
        std::cout << "\n";
    }

    std::cout << "\nОпределитель det(A) = " << determinantor << "\n";

    return 0;
}
