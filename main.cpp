#include <iostream>
#include <ctime>
#include <random>
#include <stdlib.h>

std::vector<std::vector<double>> res(2, std::vector<double>(21));

double Function(double c, double d, int N) {
    double E = 0;
    for (auto i = 0; i < N; i++) {
        auto x = res[0][i];
        auto t = res[1][i];
        double y = c * x + d;
        E += (y - t) * (y - t);
    }
    return E;
}

double Dichotomy(double a, double b, double d,int N){
    const double e = 0.1;
    const double delta = 0.01;
    double Xleft, Xright;
    double Ylower, Yupper;
    do{
        Xleft = 0.5 * (b + a) - delta;
        Xright = 0.5 * (b + a) + delta;
        Ylower = Function(Xleft, d, N);
        Yupper = Function(Xright, d, N);
        if (Ylower > Yupper)
        {
            a = Xleft;
        }
        else {
            b = Xright;
        }

    } while ((b - a) > e);

    return (a + b) / 2;
}


double passive(double a,double b, double c, int N) {
    double e = 0.1; 
    int n = 1;
    double d = a;
    double sr = (b - a) / 2;
    while (sr > e) { 
        double x, x_min = a;
        double f;
        double f_min = Function(c, a, N);
        for (auto i = 1; i <= n; i++) {
            x = (b - a) * i / (n + 1) + a;
            f = Function(c, x, N);
            if (f < f_min) {
                f_min = f;
                x_min = x;
            }
        }
        n++;
        sr = (b - a) / (n + 1);
        d = x_min;
    }
    return d;
}



int main() {
    srand(time(nullptr));
    double c = 0.5, d = 0, a = -2, b = 1;
    int N = 20;
    int A = 1;
    double step = (b - a) / N;
    int i = 0;

    for (double x = a; x < b; x += step) {
        res[0][i] = x;
        res[1][i] = c * x + d;
        i++;
    }
    double C = Dichotomy(-5, 5, d, N);
    double D = passive(-5, 5, c, N);
    std::cout << "C, D without noise:" << C << ", " << D<<std::endl;
    for (auto i = 0; i < N; i++) {
        std::cout << "(" << res[0][i] << ";" << res[1][i] << ") " << std::endl;
    }
    std::cout << std::endl;

    i = 0;
    for (double x = a; x < b; x += step) {
        res[0][i] = x;
        res[1][i] = c * x + d + A * (-0.5 + ((double)(rand() % 1000)/1000));
        i++;
    }

    C = Dichotomy(-5, 5, d, N);
    D = passive(-5, 5, c, N);
    std::cout <<"C, D with noise:"<< C << ", " << D << std::endl;
    for (auto i = 0; i < N; i++) {
        std::cout << "(" << res[0][i] << ";" << res[1][i] << ") "<<std::endl;
    }
	return 0;
}
