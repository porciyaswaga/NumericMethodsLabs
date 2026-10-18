#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <utility>


int n = 4;

std::pair<std::vector<double>, int> SimpleIteration(std::vector<std::vector<double>> A, std::vector<double> b, double parameter) {
    std::vector<std::vector<double>> B(n, std::vector<double>(n, 0.0));
    std::vector<double> c(n, 0.0);
    std::vector<double> x(n, 0.0);

    for (int i = 0; i < n; i++) if (A[i][i] == 0) throw std::runtime_error("Simple Iterarions failed: zero pivot");

    for (int i = 0; i < n; i++) {
        c[i] = b[i]/A[i][i];
        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            B[i][j] = -A[i][j]/A[i][i];
        }
    }

    double alpha = 0.0;

    for (int i = 0; i < n; i++) {
        double newAlpha = 0.0;

        for (int j = 0; j < n; j++) {
            newAlpha += std::abs(B[i][j]);
        }

        if (newAlpha > alpha) {
            alpha = newAlpha;
        }
    }

    if (alpha >= 1.0) {
        throw std::runtime_error("Alpha must be less than 1");
    }

    std::vector<double> xNew(n, 0.0);
    int iter = 0;
    while(iter <= 1000) {
        iter += 1;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                xNew[i] += B[i][j] * x[j];
            }
            xNew[i] += c[i];
        }

        double maxDiff = 0.0;

        for (int i = 0; i < n; i++) {
            double diff = std::abs(xNew[i] - x[i]);

            if (diff > maxDiff) {
                maxDiff = diff;
            }
        }

        double errorEstimate = alpha / (1.0 - alpha) * maxDiff;

        if (errorEstimate <= parameter) {
            break;
        }

        x = xNew;

        for (int i = 0; i < n; i++) {
            xNew[i] = 0.0;
        }
    }
    return {xNew, iter};
}

std::pair<std::vector<double>, int> Zeidel(std::vector<std::vector<double>> A, std::vector<double> b, double parameter) {
    std::vector<std::vector<double>> B(n, std::vector<double>(n, 0.0));
    std::vector<double> x(n, 0.0);
    std::vector<double> c(n, 0.0);

    for (int i = 0; i < n; i++) {
        if (A[i][i] == 0) throw std::runtime_error("Zeidel failed: zero pivot");
    }

    for (int i = 0; i < n; i++) {
        c[i] = b[i] / A[i][i];

        for (int j = 0; j < n; j++) {
            if (i != j) B[i][j] = -A[i][j] / A[i][i];
        }
    }

    int iter = 0;

    while (iter < 1000) {
        iter++;

        std::vector<double> xOld = x;

        for (int i = 0; i < n; i++) {
            x[i] = c[i];

            for (int j = 0; j < n; j++) {
                if (j < i) {
                    x[i] += B[i][j] * x[j];
                }
                else if (j > i) {
                    x[i] += B[i][j] * xOld[j];
                }
            }
        }

        double maxDiff = 0.0;

        for (int i = 0; i < n; i++) {
            double diff = std::abs(x[i] - xOld[i]);

            if (diff > maxDiff) {
                maxDiff = diff;
            }
        }

        if (maxDiff <= parameter) {
            return {x, iter};
        }
    }

    throw std::runtime_error("Zeidel failed: iteration limit exceeded");
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

    double parameter;
    std::cout << "Введите точность вычислений:\n";
    std::cin >> parameter;

    auto [x1, iter1] = SimpleIteration(A, b, parameter);

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

    std::cout << "\nРешение СЛАУ (методом простых итераций):\n";
    for (int i = 0; i < n; i++) {
        std::cout << "x[" << i + 1 << "] = " << x1[i] << "\n";
    }

    std::cout << "\nКол-во итераций в методе простых итераций:\n";
    std::cout << iter1 << std::endl;

    auto [x2, iter2] = Zeidel(A, b, parameter);

    std::cout << "\nРешение СЛАУ (методом Зейделя):\n";
    for (int i = 0; i < n; i++) {
        std::cout << "x[" << i + 1 << "] = " << x2[i] << "\n";
    }

    std::cout << "\nКол-во итераций в методе Зейделя:\n";
    std::cout << iter2 << std::endl;
}
