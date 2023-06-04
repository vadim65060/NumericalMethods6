//
// Created by vadim on 03.06.2023.
//

#include "Cauchy.h"
#include <iomanip>

Cauchy::Cauchy(const std::function<double(double, double)> &func, double h, int n, double x0, double y0,
               const std::function<double(double)> &trueY) {
    if (trueY == nullptr) {
        this->trueY = [](double x) -> double {
            return 0;
        };
    } else {
        this->trueY = trueY;
    }

    this->func = func;
    this->h = h;
    this->N = n;
    this->x0 = x0;
    this->y0 = y0;

    points = std::vector<double>(N + 1);
    points[0] = x0;
    for (int i = 1; i <= N; i++) {
        points[i] = x0 + i * h;
    }

    values = std::vector<double>(N + 1);
    values[0] = y0;
    for (int i = 1; i <= N; i++) {
        values[i] = this->trueY(points[i]);
    }

    pointsn = std::vector<double>(N + 3);
    for (int i = -2; i <= N; i++) {
        pointsn[i + 2] = x0 + i * h;
    }

    valuesn = std::vector<double>(N + 3);
    for (int i = -2; i <= N; i++) {
        valuesn[i + 2] = this->trueY(pointsn[i + 2]);
    }

}

void Cauchy::Solve() {
    GetTaylor();
    GetEulerOne();
    GetEulerTwo();
    GetEulerThree();
    GetRunge();
    GetAdams();

    std::cout.setf(std::ios::left);
    std::cout << "true y = " << values[N]<<'\n';
    std::cout << "Table for y_N for all methods:\n";
    std::cout << "method      y" << N << "             ABSOLUTE ERROR\n";
    std::cout << "Talor       " << std::setw(wrap) << taylor[N] << "       " << std::setw(wrap)
              << std::abs(taylor[N] - values[N]) << '\n';
    std::cout << "Euler 1     " << std::setw(wrap) << eulerOne[N] << "       " << std::setw(wrap)
              << std::abs(eulerOne[N] - values[N]) << '\n';
    std::cout << "Euler 2     " << std::setw(wrap) << eulerTwo[N] << "       " << std::setw(wrap)
              << std::abs(eulerTwo[N] - values[N]) << '\n';
    std::cout << "Euler 3     " << std::setw(wrap) << eulerThree[N] << "       " << std::setw(wrap)
              << std::abs(eulerThree[N] - values[N]) << '\n';
    std::cout << "Runge       " << std::setw(wrap) << runge[N] << "       " << std::setw(wrap)
              << std::abs(runge[N] - values[N]) << '\n';
    std::cout << "Adams       " << std::setw(wrap) << adams[N] << "       " << std::setw(wrap)
              << std::abs(adams[N] - values[N]) << '\n';
    std::cout << "------------------------------------------------------------------------\n";
    std::cout.unsetf(std::ios::left);
}

void Cauchy::GetTaylor() {
    taylor = std::vector<double>(N + 3);
    taylor[0] = 0;
    der = std::vector<double>(7);
    der[0] = 1;
    der[1] = -1;
    der[2] = 2;
    der[3] = -2;
    der[4] = 1;
    der[5] = -1;
    der[6] = 2;
    for (int k = -2; k <= N; k++) {
        for (int j = 0; j <= 6; j++) {
            taylor[k + 2] += (der[j] * pow((pointsn[k + 2] - points[0]), j)) / fact(j);
        }
    }
}

void Cauchy::GetEulerOne() {
    eulerOne = std::vector<double>(N + 1);
    eulerOne[0] = y0;
    for (int i = 1; i <= N; i++) {
        eulerOne[i] = eulerOne[i - 1] + h * func(points[i - 1], eulerOne[i - 1]);
    }
}

void Cauchy::GetEulerTwo() {
    eulerTwo = std::vector<double>(N + 1);
    eulerTwo[0] = y0;
    for (int i = 1; i <= N; i++) {
        eulerTwo[i] = eulerTwo[i - 1] +
                      h * func(points[i - 1] + h / 2, eulerTwo[i - 1] + h / 2 * func(points[i - 1], eulerTwo[i - 1]));
    }
}

void Cauchy::GetEulerThree() {
    eulerThree = std::vector<double>(N + 1);
    eulerThree[0] = y0;
    for (int i = 1; i <= N; i++) {
        eulerThree[i] = eulerThree[i - 1] +
                        h / 2 * (func(points[i - 1], eulerThree[i - 1]) +
                                 func(points[i], eulerThree[i - 1] + h * func(points[i - 1], eulerThree[i - 1])));
    }
}

void Cauchy::GetRunge() {
    runge = std::vector<double>(N + 1);
    runge[0] = y0;
    double k1, k2, k3, k4;
    for (int i = 1; i <= N; i++) {
        k1 = h * func(points[i - 1], runge[i - 1]);
        k2 = h * func(points[i - 1] + h / 2, runge[i - 1] + k1 / 2);
        k3 = h * func(points[i - 1] + h / 2, runge[i - 1] + k2 / 2);
        k4 = h * func(points[i - 1] + h, runge[i - 1] + k3);
        runge[i] = runge[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6.;
    }
}

void Cauchy::GetAdams() {
    std::vector<std::vector<double> > matrix(N + 100, std::vector<double>(100 + N));
    adams = std::vector<double>(N + 1);
    matrix[0] = taylor;
    for (int i = 0; i < 5; i++) {
        matrix[1][i] = h * func(pointsn[i], matrix[0][i]);
    }
    for (int i = 2; i < 6; i++)
        for (int j = 0; j < 6 - i; j++)
            matrix[i][j] = matrix[i - 1][j + 1] - matrix[i - 1][j];

    for (int i = 1; i < N - 1; i++) {
        matrix[0][4 + i] = matrix[0][4 + i - 1] + matrix[1][4 + i - 1] + matrix[2][4 + i - 2] / 2 +
                           5 * matrix[3][4 + i - 3] / 12 + 3 * matrix[4][4 + i - 4] / 8 +
                           251 * matrix[5][4 + i - 5] / 720;

        matrix[1][4 + i] = h * func(pointsn[4 + i], matrix[0][4 + i]);

        matrix[2][4 + i - 1] = matrix[1][4 + i] - matrix[1][4 + i - 1];

        matrix[3][4 + i - 2] = matrix[2][4 + i - 1] - matrix[2][4 + i - 2];

        matrix[4][4 + i - 3] = matrix[3][4 + i - 2] - matrix[3][4 + i - 3];

        matrix[5][4 + i - 4] = matrix[4][4 + i - 3] - matrix[4][4 + i - 4];
    }
    if (N - 2 >= 0)
        for (int i = 0; i < N - 2; i++)
            adams[i + 3] = matrix[0][5 + i];
}

void Cauchy::PrintTaylor() {
    std::cout << "Taylor\n";
    for (int i = 0; i < N + 3; i++) {
        PrintMethodValue(pointsn[i], taylor[i], valuesn[i]);
    }
    std::cout << "------------------------------------------------------------------------\n";
}

void Cauchy::PrintEulerOne() {
    std::cout << "Euler 1:\n";
    for (int i = 0; i < N + 1; i++) {
        PrintMethodValue(points[i], eulerOne[i], values[i]);
    }
    std::cout << "------------------------------------------------------------------------\n";
}

void Cauchy::PrintEulerTwo() {
    std::cout << "Euler 2:\n";
    for (int i = 0; i < N + 1; i++) {
        PrintMethodValue(points[i], eulerTwo[i], values[i]);
    }
    std::cout << "------------------------------------------------------------------------\n";
}

void Cauchy::PrintEulerThree() {
    std::cout << "Euler 3:\n";
    for (int i = 0; i < N + 1; i++) {
        PrintMethodValue(points[i], eulerThree[i], values[i]);
    }
    std::cout << "------------------------------------------------------------------------\n";
}

void Cauchy::PrintRunge() {
    std::cout << "Runge:\n";
    for (int i = 0; i < N + 1; i++) {
        PrintMethodValue(points[i], runge[i], values[i]);
    }
    std::cout << "------------------------------------------------------------------------\n";
}

void Cauchy::PrintAdams() {
    std::cout << "Adams:\n";
    for (int i = 3; i < N + 1; i++) {
        PrintMethodValue(points[i], adams[i], values[i]);
    }
    std::cout << "------------------------------------------------------------------------\n";
}

void Cauchy::PrintValues() {
    std::cout << "True value\n";
    for (int i = 0; i < N + 3; i++) {
        PrintMethodValue(pointsn[i], valuesn[i]);
    }
    std::cout << "------------------------------------------------------------------------\n";
}

void Cauchy::PrintMethodValue(double x, double y, double trueValue) const {
    std::ostringstream sstream;
    sstream << x;
    std::string strX = sstream.str();

    std::string y_x = "y(" + strX + ")";
    std::cout << std::setw(8) << std::right << y_x << " = " << std::left << std::setw(wrap) << y << "   ";
    if (trueValue != 3e-20)
        std::cout << std::setw(24) << "|" + y_x + " - trueY(" + strX + ")|" << " = " << std::abs(y - trueValue);
    std::cout << '\n';
}
