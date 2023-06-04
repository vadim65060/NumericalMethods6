#include <iostream>
#include "Cauchy.h"

double f(double x, double y) {
    return -y * y;
}

double TrueY(double x) {
    return 1 / (x + 1);
}

int main() {
    int n = 10;
    double h = 0.1, x0 = 0, y0 = 1;
    std::cout << n << ' ' << h << '\n';
    Cauchy cauchy = Cauchy(f, h, n, x0, y0, TrueY);
    cauchy.Solve();
    while (true) {
        std::cout
                << "0 - all\n1 - Taylor\n2 - Euler 1\n3 - Euler 2\n4 - Euler 3\n5 - Runge\n6 - Adams\n7 - True value\n-1 - exit\n";
        std::string s;
        std::cin >> s;
        system("cls");
        int key;
        try {
            key = std::stoi(s);
        } catch (const std::exception &e) {
            continue;
        }
        switch (key) {
            case 0:
                cauchy.Solve();
                break;
            case 1:
                cauchy.PrintTaylor();
                break;
            case 2:
                cauchy.PrintEulerOne();
                break;
            case 3:
                cauchy.PrintEulerTwo();
                break;
            case 4:
                cauchy.PrintEulerThree();
                break;
            case 5:
                cauchy.PrintRunge();
                break;
            case 6:
                cauchy.PrintAdams();
                break;
            case 7:
                cauchy.PrintValues();
                break;
            default:
                break;
        }
        if (key == -1) return 0;
    }
}
