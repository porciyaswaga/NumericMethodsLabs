
#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>


int n = 3;

std::vector<std::vector<double>> JacobiRotationMethod(std::vector<std::vector<double>> A, double parameter, std::vector<std::vector<double>>& V, int& iterations) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (std::abs(A[i][j] - A[j][i]) > 1e-10) {
                throw std::runtime_error("RotationMethod failed: Matrix A is not symmetric");
            }
        }
    }

    V = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; i++) {
        V[i][i] = 1.0;
    }

    for (int iter = 0; iter < 1000; iter++) {
        double maxEl = 0.0;
        int p = 0;
        int q = 1;

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (std::abs(A[i][j]) > maxEl) {
                    maxEl = std::abs(A[i][j]);
                    p = i;
                    q = j;
                }
            }
        }

        if (maxEl <= parameter) {
            iterations = iter;
            return A;
        }

        double phi = 0.5 * std::atan2(2.0 * A[p][q], A[q][q] - A[p][p]);
        double c = std::cos(phi);
        double s = std::sin(phi);

        std::vector<std::vector<double>> U(n, std::vector<double>(n, 0.0));

        for (int i = 0; i < n; i++) {
            U[i][i] = 1.0;
        }

        U[p][p] = c;
        U[q][q] = c;
        U[p][q] = s;
        U[q][p] = -s;

        std::vector<std::vector<double>> temp(n, std::vector<double>(n, 0.0));
        std::vector<std::vector<double>> ANew(n, std::vector<double>(n, 0.0));
        std::vector<std::vector<double>> VNew(n, std::vector<double>(n, 0.0));

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    temp[i][j] += A[i][k] * U[k][j];
                }
            }
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    ANew[i][j] += U[k][i] * temp[k][j];
                }
            }
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    VNew[i][j] += V[i][k] * U[k][j];
                }
            }
        }

        A = ANew;
        V = VNew;
    }

    throw std::runtime_error("RotationMethod failed: iteration limit exceeded");
}

int main() {
    std::vector<std::vector<double>> A(n, std::vector<double>(n));

    std::cout << "Введите матрицу A:\n";

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cin >> A[i][j];
        }
    }

    double parameter;

    std::cout << "Введите точность вычислений:\n";
    std::cin >> parameter;

    if (parameter <= 0.0) {
        throw std::runtime_error("Precision must be positive");
    }

    std::vector<std::vector<double>> V;
    int iterations = 0;

    std::vector<std::vector<double>> ANew = JacobiRotationMethod(A, parameter, V, iterations);

    std::cout << "\nСобственные значения:\n";

    for (int i = 0; i < n; i++) {
        std::cout << ANew[i][i] << "\n";
    }

    std::cout << "\nСобственные векторы (по столбцам):\n";

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << V[i][j] << " ";
        }
        std::cout << "\n";
    }

    std::cout << "\nКоличество итераций: " << iterations << "\n";

    double maxEl = 0.0;

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (std::abs(ANew[i][j]) > maxEl) {
                maxEl = std::abs(ANew[i][j]);
            }
        }
    }

    std::cout << "Максимальный внедиагональный элемент: " << maxEl << "\n";

    return 0;
}
