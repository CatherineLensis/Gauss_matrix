#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

// Функция для вычисления элементов матрицы A
void calculateMatrixA(vector<vector<double>>& A, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A[i][j] = atan(i + pow(j, 2));
        }
    }
}

// Функция для вывода матрицы в консоль и файл
void printMatrix(const vector<vector<double>>& matrix, const string& filename) {
    ofstream file(filename);
    if (file.is_open()) {
        for (const auto& row : matrix) {
            for (const auto& elem : row) {
                cout << elem << "\t";
                file << elem << "\t";
            }
            cout << endl;
            file << endl;
        }
        file.close();
    }
    else {
        cerr << "Unable to open file for writing.";
    }
}

// Функция для вычисления определителя матрицы
double determinant(const vector<vector<double>>& A) {
    int n = A.size();
    vector<vector<double>> temp = A;
    double det = 1;

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            while (temp[j][i] != 0) {
                double ratio = temp[i][i] / temp[j][i];
                for (int k = 0; k < n; ++k) {
                    temp[i][k] -= ratio * temp[j][k];
                }
                if (temp[i][i] == 0) break;
                swap(temp[i], temp[j]);
                det *= -1;
            }
        }
        det *= temp[i][i];
    }
    return det;
}

// Функция для вычисления обратной матрицы с помощью метода Гаусса-Жордана
vector<vector<double>> inverseMatrix(const vector<vector<double>>& A) {
    int n = A.size();
    vector<vector<double>> R(n, vector<double>(2 * n, 0));
    vector<vector<double>> RR = A;

    // Создание матрицы R
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            R[i][j] = A[i][j];
        }
        R[i][n + i] = 1;
    }

    // Прямой ход метода Гаусса-Жордана
    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < 2 * n; ++j) {
            if (j == k) continue;
            RR[k][j] /= RR[k][k];
        }
        RR[k][k] = 1;

        for (int i = 0; i < n; ++i) {
            if (i == k) continue;
            for (int j = 0; j < 2 * n; ++j) {
                if (j == k) continue;
                RR[i][j] -= RR[i][k] * RR[k][j];
            }
            RR[i][k] = 0;
        }
    }

    // Извлечение обратной матрицы из RR
    vector<vector<double>> A_inverse(n, vector<double>(n, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A_inverse[i][j] = RR[i][n + j];
        }
    }

    return A_inverse;
}

// Функция для вычисления матрицы H и вывода ее в консоль и файл
void calculateMatrixH(const vector<vector<double>>& A, const vector<vector<double>>& A_inverse, const string& filename) {
    int n = A.size();
    vector<vector<double>> H(n, vector<double>(n, 0));

    // Вычисление матрицы H = A * A^-1
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                H[i][j] += A[i][k] * A_inverse[k][j];
            }
        }
    }

    // Вывод матрицы H
    printMatrix(H, filename);
}

// Функция для вычисления числа обусловленности матрицы A
double conditionNumber(const vector<vector<double>>& A, const vector<vector<double>>& A_inverse) {
    double normA = 0, normA_inverse = 0;
    int n = A.size();

    // Вычисление норм матрицы A и A^-1
    for (int j = 0; j < n; ++j) {
        double sumA = 0, sumA_inverse = 0;
        for (int i = 0; i < n; ++i) {
            sumA += abs(A[i][j]);
            sumA_inverse += abs(A_inverse[i][j]);
        }
        normA = max(normA, sumA);
        normA_inverse = max(normA_inverse, sumA_inverse);
    }

    return normA * normA_inverse;
}

int main() {
    int n;
    cout << "Enter the size of the matrix: ";
    cin >> n;

    // Создание матрицы A
    vector<vector<double>> A(n, vector<double>(n, 0));
    calculateMatrixA(A, n);

    // Вывод матрицы A
    cout << "Matrix A:" << endl;
    printMatrix(A, "matrix_A.txt");

    // Вычисление определителя матрицы A
    double det = determinant(A);
    cout << "Determinant of A: " << det << endl;

    // Вычисление обратной матрицы A^-1
    vector<vector<double>> A_inverse = inverseMatrix(A);

    // Вывод матрицы A^-1
    cout << "Matrix A^-1:" << endl;
    printMatrix(A_inverse, "matrix_A_inverse.txt");

    // Вычисление матрицы H = A * A^-1
    cout << "Matrix H = A * A^-1:" << endl;
    calculateMatrixH(A, A_inverse, "matrix_H.txt");

    // Вычисление числа обусловленности матрицы A
    double cond = conditionNumber(A, A_inverse);
    cout << "Condition number of A: " << cond << endl;

    return 0;
}
