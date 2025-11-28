#include <iostream>
#include <cmath>
#include <functional>
#include <sstream>
#include <map>
#include <cctype>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <clocale>
using namespace std;
class RealFunctionParser {
private:
    string expression;
    size_t pos;

    void skipWhitespace() {
        while (pos < expression.length() && isspace(expression[pos])) {
            pos++;
        }
    }

    double parseNumber() {
        skipWhitespace();
        size_t start = pos;

        while (pos < expression.length() &&
            (isdigit(expression[pos]) || expression[pos] == '.' ||
                expression[pos] == '-' || expression[pos] == 'e' || expression[pos] == 'E')) {
            pos++;
        }

        string numStr = expression.substr(start, pos - start);
        return stod(numStr);
    }

    double parseFunction() {
        skipWhitespace();
        size_t start = pos;

        while (pos < expression.length() && isalpha(expression[pos])) {
            pos++;
        }

        string funcName = expression.substr(start, pos - start);
        skipWhitespace();

        if (expression[pos] != '(') {
            throw runtime_error("Expected '(' after function " + funcName);
        }

        pos++;
        double arg = parseExpression();
        skipWhitespace();

        if (expression[pos] != ')') {
            throw runtime_error("Expected ')' after function argument " + funcName);
        }
        pos++;

        // Math functions
        if (funcName == "sin") return sin(arg);
        if (funcName == "cos") return cos(arg);
        if (funcName == "exp") return exp(arg);
        if (funcName == "log") return log(arg);
        if (funcName == "sqrt") return sqrt(arg);
        if (funcName == "tan") return tan(arg);
        if (funcName == "abs") return fabs(arg);

        throw runtime_error("Unknown function: " + funcName);
    }

    double parseFactor() {
        skipWhitespace();

        if (expression[pos] == '(') {
            pos++;
            double result = parseExpression();
            skipWhitespace();
            if (expression[pos] == ')') {
                pos++;
            }
            return result;
        }
        else if (expression[pos] == 'x' || expression[pos] == 'X') {
            pos++;
            return 0;
        }
        else if (isdigit(expression[pos]) || expression[pos] == '-' || expression[pos] == '.') {
            return parseNumber();
        }
        else if (isalpha(expression[pos])) {
            return parseFunction();
        }

        throw runtime_error("Unexpected character: " + string(1, expression[pos]));
    }

    double parsePower() {
        double left = parseFactor();
        skipWhitespace();

        while (pos < expression.length() && expression[pos] == '^') {
            pos++;
            double right = parseFactor();
            left = pow(left, right);
            skipWhitespace();
        }

        return left;
    }

    double parseTerm() {
        double left = parsePower();
        skipWhitespace();

        while (pos < expression.length() && (expression[pos] == '*' || expression[pos] == '/')) {
            char op = expression[pos];
            pos++;
            double right = parsePower();

            if (op == '*') left *= right;
            else if (op == '/') {
                if (right == 0) throw runtime_error("Division by zero");
                left /= right;
            }
            skipWhitespace();
        }

        return left;
    }

    double parseExpression() {
        double left = parseTerm();
        skipWhitespace();

        while (pos < expression.length() && (expression[pos] == '+' || expression[pos] == '-')) {
            char op = expression[pos];
            pos++;
            double right = parseTerm();

            if (op == '+') left += right;
            else if (op == '-') left -= right;
            skipWhitespace();
        }

        return left;
    }

public:
    function<double(double)> parse(const string& expr) {
        return [expr](double x) {
            RealFunctionParser parser;
            parser.expression = expr;
            parser.pos = 0;

            string processed;
            for (size_t i = 0; i < expr.length(); i++) {
                if (expr[i] == 'x' || expr[i] == 'X') {
                    string x_str = to_string(x);
                    processed += x_str;
                }
                else {
                    processed += expr[i];
                }
            }

            parser.expression = processed;
            return parser.parseExpression();
            };
    }
};

class NumericalIntegrator {
private:
    static const int GAUSS_NODES = 4;
    static const double gauss_points[GAUSS_NODES];
    static const double gauss_weights[GAUSS_NODES];

public:
    double simpson(const function<double(double)>& f, double a, double b, int n_segments = 4) {
        if (n_segments < 1) n_segments = 1;

        vector<double> nodes(n_segments + 1);
        vector<double> h(n_segments);

        for (int i = 0; i <= n_segments; i++) {
            nodes[i] = a + i * (b - a) / n_segments;
        }

        for (int i = 0; i < n_segments; i++) {
            h[i] = nodes[i + 1] - nodes[i];
        }
        double sum = 0.0;

        sum += h[0] * f(nodes[0]);

        sum += h[n_segments - 1] * f(nodes[n_segments]);

        for (int i = 0; i < n_segments; i++) {
            double mid_point = (nodes[i] + nodes[i + 1]) / 2.0;
            sum += 4.0 * h[i] * f(mid_point);
        }

        for (int i = 1; i < n_segments; i++) {
            sum += (h[i] + h[i - 1]) * f(nodes[i]);
        }

        return sum / 6.0;
    }

    double gauss(const function<double(double)>& f, double a, double b, int n_segments = 4) {
        if (n_segments < 1) n_segments = 1;

        vector<double> nodes(n_segments + 1);
        vector<double> h(n_segments);

        for (int i = 0; i <= n_segments; i++) {
            nodes[i] = a + i * (b - a) / n_segments;
        }

        for (int i = 0; i < n_segments; i++) {
            h[i] = nodes[i + 1] - nodes[i];
        }

        double sum = 0.0;
        for (int k = 0; k < n_segments; k++) {
            double segment_sum = 0.0;

            for (int i = 0; i < GAUSS_NODES; i++) {
                double x_ik = h[k] * (gauss_points[i] + 1.0) / 2.0 + nodes[k];
                segment_sum += gauss_weights[i] * f(x_ik);
            }

            sum += h[k] * segment_sum;
        }

        return sum / 2.0;
    }

    void compareMethods(const function<double(double)>& f, double a, double b) {
        double result_simpson = simpson(f, a, b, 4);
        double result_gauss = gauss(f, a, b, 4);

        cout << fixed << setprecision(6);
        const char* old_locale = setlocale(LC_ALL, "Russian");
        cout << "Метод Симпсона (4 узла): " << result_simpson << endl;
        cout << "Метод Гаусса (4 узла): " << result_gauss << endl;
        setlocale(LC_ALL, old_locale);

        cout.unsetf(ios_base::floatfield);
    }
};

const double NumericalIntegrator::gauss_points[4] = {
    -0.8611363115940526,
    -0.3399810435848563,
     0.3399810435848563,
     0.8611363115940526
};

const double NumericalIntegrator::gauss_weights[4] = {
    0.3478548451374538,
    0.6521451548625461,
    0.6521451548625461,
    0.3478548451374538
};

int main() {
    NumericalIntegrator integrator;
    RealFunctionParser parser;

        string user_expr;
        cout << "f(x) = ";
        getline(cin, user_expr);

        if (user_expr == "exit" || user_expr == "quit") {
            cout << "Exiting program.\n";

        }

        function<double(double)> func;
        try {
            func = parser.parse(user_expr);
        }
        catch (const exception& e) {
            cout << "Error parsing function: " << e.what() << ". Try again.\n";
        }

        double a, b;
        cout << "a, b: ";
        cin >> a >> b;
        cin.ignore();

        if (a > b) {
            swap(a, b);
            cout << "a=" << a << ", b=" << b << endl;
        }
        cout << "f(x)= " << user_expr << endl;
        cout << "[a,b]: [" << a << ", " << b << "]\n";
        integrator.compareMethods(func, a, b);
    }
