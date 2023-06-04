//
// Created by vadim on 03.06.2023.
//

#ifndef NM6_CAUCHY_H
#define NM6_CAUCHY_H


#include <utility>
#include <vector>
#include <iostream>
#include <cmath>
#include <functional>

class Cauchy {
private:
    std::function<double(double, double)> func;
    std::function<double(double)> trueY;
    double x0 = 0;
    double y0 = 1;
    double h;
    int N;
    int wrap = 10;
    std::vector<double> points;
    std::vector<double> der;
    std::vector<double> values;
    std::vector<double> taylor;
    std::vector<double> eulerOne;
    std::vector<double> eulerTwo;
    std::vector<double> eulerThree;
    std::vector<double> runge;
    std::vector<double> adams;
    std::vector<double> pointsn;
    std::vector<double> valuesn;

public:
    Cauchy(const std::function<double(double, double)> &func, double h, int n, double x0, double y0,
           const std::function<double(double)> &trueY = nullptr);

    void Solve();

    void GetTaylor();

    void GetEulerOne();

    void GetEulerTwo();

    void GetEulerThree();

    void GetRunge();

    void GetAdams();

    void PrintTaylor();

    void PrintEulerOne();

    void PrintEulerTwo();

    void PrintEulerThree();

    void PrintRunge();

    void PrintAdams();

    void PrintValues();

private:

    void PrintMethodValue(double x, double y, double trueValue = 3e-20) const;

    double fact(int j) {
        int res = 1;
        for (int i = 2; i <= j; i++) {
            res *= i;
        }
        return res;
    }
};


#endif //NM6_CAUCHY_H
