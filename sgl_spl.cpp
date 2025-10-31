#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <limits>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;
class RandGen {
private:
    double mean;
    double variance;
    default_random_engine generator;
    normal_distribution<double> distribution;
public:
    RandGen(double m, double d) : mean(m), variance(d), generator(time(0)), distribution(m, d) {}

    double generate() {
        return distribution(generator);
    }
};

class Point {
private:
    double x_val, y_val;

public:
    Point(double x = 0, double y = 0) : x_val(x), y_val(y) {}

    double x() const { return x_val; }
    double y() const { return y_val; }
    void setX(double x) { x_val = x; }
    void setY(double y) { y_val = y; }
};

class SmoothingSpline {
private:
    vector<Point> points;
    double p; 
    vector<double> alpha; 
    double basis_function(int i, double psi) const {
        if (i == 1) return 1 - psi;
        if (i == 2) return psi;
        return 0;
    }
    double der_basis_function(int i, double psi) const {
        if (i == 1) return -1;
        if (i == 2) return 1;
        return 0;
    }

    void transition_to_master_element(int segment, double x, double& psi) const {
        double h = points[segment + 1].x() - points[segment].x();
        psi = (x - points[segment].x()) / h;
    }

    vector<double> thomas_algorithm(const vector<double>& a, const vector<double>& b,
        const vector<double>& c, const vector<double>& d) {
        int n = d.size();
        vector<double> alpha(n);
        vector<double> beta(n);
        vector<double> x(n);

        beta[0] = c[0] / b[0];
        alpha[0] = d[0] / b[0];

        for (int i = 1; i < n; i++) {
            double denominator = b[i] - a[i] * beta[i - 1];
            beta[i] = c[i] / denominator;
            alpha[i] = (d[i] - a[i] * alpha[i - 1]) / denominator;
        }
        x[n - 1] = alpha[n - 1];
        for (int i = n - 2; i >= 0; i--) {
            x[i] = alpha[i] - beta[i] * x[i + 1];
        }

        return x;
    }

public:
    SmoothingSpline(const vector<Point>& pts, double smoothing_param)
        : points(pts), p(smoothing_param) {
        build_spline();
    }

    void build_spline() {
        int n = points.size();
        int num_segments = n - 1;

        if (n < 2) {
            alpha = vector<double>(n, 0);
            return;
        }

        vector<double> a(n, 0), b(n, 0), c(n, 0), d(n, 0);

        for (int i = 0; i < num_segments; i++) {
            double h = points[i + 1].x() - points[i].x();

            double m11 = h / 3.0, m12 = h / 6.0;
            double m21 = h / 6.0, m22 = h / 3.0;

            double k11 = 1.0 / h, k12 = -1.0 / h;
            double k21 = -1.0 / h, k22 = 1.0 / h;

            b[i] += (1 - p) * m11 + p * k11;
            c[i] += (1 - p) * m12 + p * k12;
            a[i + 1] += (1 - p) * m21 + p * k21;
            b[i + 1] += (1 - p) * m22 + p * k22;

            d[i] += (1 - p) * (m11 * points[i].y() + m12 * points[i + 1].y());
            d[i + 1] += (1 - p) * (m21 * points[i].y() + m22 * points[i + 1].y());
        }

        b[0] += p / (points[1].x() - points[0].x());
        b[n - 1] += p / (points[n - 1].x() - points[n - 2].x());
        alpha = thomas_algorithm(a, b, c, d);
    }
    double evaluate(double x) const {
        double eps = 1e-7;
        int num_segments = points.size() - 1;
        for (int i = 0; i < num_segments; i++) {
            if ((x > points[i].x() && x < points[i + 1].x()) ||
                fabs(x - points[i].x()) < eps ||
                fabs(x - points[i + 1].x()) < eps) {

                double psi;
                transition_to_master_element(i, x, psi);

                return alpha[i] * basis_function(1, psi) +
                    alpha[i + 1] * basis_function(2, psi);
            }
        }
        if (x <= points[0].x()) return alpha[0];
        if (x >= points.back().x()) return alpha.back();

        return 0;
    }
    double derivative(double x) const {
        double eps = 1e-7;
        int num_segments = points.size() - 1;

        for (int i = 0; i < num_segments; i++) {
            if ((x > points[i].x() && x < points[i + 1].x()) ||
                fabs(x - points[i].x()) < eps ||
                fabs(x - points[i + 1].x()) < eps) {

                double h = points[i + 1].x() - points[i].x();
                double psi;
                transition_to_master_element(i, x, psi);

                return (alpha[i] * der_basis_function(1, psi) +
                    alpha[i + 1] * der_basis_function(2, psi)) * 2.0 / h;
            }
        }

        return 0;
    }
};

int main() {
    setlocale(LC_ALL, "Russian");

    int n;
    double mean, variance;
    cout << "Введите количество наблюдений N = ";
    cin >> n;
    cout << "Введите мат. ожидание M = ";
    cin >> mean;
    cout << "Введите среднее квадратичное отклонение = ";
    cin >> variance;
    RandGen generator(mean, variance);
    vector<double> y(n);
    vector<Point> points(n);
    for (int i = 0; i < n; ++i) {
        y[i] = generator.generate();
        points[i] = Point(i, y[i]);
        cout << "Точка " << i + 1 << ": x = " << i << ", y = " << y[i] << endl;
    }
    vector<double> p_values = { 0.0, 0.1, 0.5, 0.7, 0.99 };
    vector<vector<double>> results(p_values.size(), vector<double>(n));

    for (int i = 0; i < p_values.size(); ++i) {
        SmoothingSpline spline(points, p_values[i]);

        for (int j = 0; j < n; ++j) {
            results[i][j] = spline.evaluate(j);
        }
    }

    ofstream csvFile("spline.txt");
    csvFile << fixed << setprecision(5);
    for (int i = 0; i < n; ++i) {
        csvFile << i + 1 << "," << y[i];
        for (const auto& result : results)
            csvFile << "," << result[i];
        csvFile << "\n";
    }

    csvFile.close();
    cout << "\nДанные в файл записаны\n";
}
