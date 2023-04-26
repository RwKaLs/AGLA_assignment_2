#include <iostream>
#include <cstdio>
#include <utility>
#include <vector>
#include <cmath>
#include <iomanip>

#include <math.h>

#define GNUPLOT_LSA "C:\\gnuplot\\bin\\gnuplot -perist"

using namespace std;

template<typename T>
class Matrix
{
protected:
    int rows;
    int cols;
    vector<vector<T>> data;
public:
    Matrix(int rows, int cols) : rows(rows), cols(cols)
    {}

    Matrix(int rows, int cols, vector<vector<T>> data) : rows(rows), cols(cols), data(std::move(data))
    {}

    friend std::istream &operator>>(std::istream &is, Matrix &mat)
    {
        mat.data = vector<vector<T>>(mat.rows);
        for (int i = 0; i < mat.rows; i++) {
            mat.data[i] = vector<T>(mat.cols);
            for (int j = 0; j < mat.cols; j++) {
                is >> mat.data[i][j];
            }
        }
        return is;
    }

    Matrix<T> operator+(const Matrix<T> &m)
    {
        if (this->rows != m.rows || this->cols != m.cols) {
            return {0, 0, vector<vector<T>>()};
        }
        vector<vector<T>> res = vector<vector<T>>(this->rows);
        for (int i = 0; i < this->rows; ++i) {
            res[i] = vector<T>(this->cols);
            for (int j = 0; j < this->cols; ++j) {
                res[i][j] = this->data[i][j] + m.data[i][j];
            }
        }
        Matrix<T> result(this->rows, this->cols, res);
        return result;
    }

    Matrix<T> operator-(const Matrix<T> &m)
    {
        if (this->rows != m.rows || this->cols != m.cols) {
            return {0, 0, vector<vector<T>>()};
        }
        vector<vector<T>> res = vector<vector<T>>(this->rows);
        for (int i = 0; i < m.rows; ++i) {
            res[i] = vector<T>(m.cols);
            for (int j = 0; j < m.cols; ++j) {
                res[i][j] = this->data[i][j] - m.data[i][j];
            }
        }
        Matrix<T> result(this->rows, this->cols, res);
        return result;
    }

    Matrix<T> operator*(const Matrix<T> &m)
    {
        if (this->cols != m.rows) {
            return {0, 0, vector<vector<T>>()};
        }
        vector<vector<T>> res = vector<vector<T>>(this->rows);
        for (int i = 0; i < this->rows; ++i) {
            res[i] = vector<T>(m.cols);
            for (int j = 0; j < m.cols; ++j) {
                res[i][j] = 0;
                for (int k = 0; k < this->cols; ++k) {
                    res[i][j] += this->data[i][k] * m.data[k][j];
                }
                if (round(res[i][j] * 10000) / 10000 == 0) {
                    res[i][j] = 0.0000;
                }
            }
        }
        Matrix<T> result(this->rows, m.cols, res);
        return result;
    }

    Matrix<T> &operator=(const Matrix<T> &m)
    = default;

    Matrix<T> transpose()
    {
        vector<vector<T>> res = vector<vector<T>>(this->cols);
        for (int i = 0; i < this->cols; ++i) {
            res[i] = vector<T>(this->rows);
            for (int j = 0; j < this->rows; ++j) {
                res[i][j] = this->data[j][i];
            }
        }
        Matrix<T> result(this->cols, this->rows, res);
        return result;
    }

    T getElement(int i, int j)
    {
        return data[i][j];
    }

    void setElement(int i, int j, T x)
    {
        data[i][j] = x;
    }

    vector<vector<T>> getData() {
        return this->data;
    }

    friend std::ostream &operator<<(std::ostream &os, const Matrix &m)
    {
        if (m.rows == 0 || m.cols == 0) {
            os << "Error: the dimensional problem occurred" << endl;
        } else {
            for (int i = 0; i < m.rows; ++i) {
                for (int j = 0; j < m.cols; ++j) {
                    if (j == m.cols - 1) {
                        if (round(m.data[i][j] * 10000) / 10000 == 0) {
                            os << "0.0000";
                        } else {
                            os << fixed << setprecision(4) << round(m.data[i][j] * 10000) / 10000;
                        }
                    } else {
                        if (round(m.data[i][j] * 10000) / 10000 == 0) {
                            os << "0.0000 ";
                        } else {
                            os << fixed << setprecision(4) << round(m.data[i][j] * 10000) / 10000 << " ";
                        }
                    }
                }
                os << endl;
            }
        }
        return os;
    }
};

template<class T>
class SquareMatrix : public Matrix<T>
{
protected:
    int n{};
public:
    explicit SquareMatrix(int n) : Matrix<T>(n, n)
    {
        this->n = n;
    }

    SquareMatrix(int n, Matrix<T> m) : Matrix<T>(n, n, m.getData()) {
        this->n = n;
    }

    virtual ~SquareMatrix() = default;
};

template<class T>
class IdentityMatrix : public SquareMatrix<T>
{
public:
    explicit IdentityMatrix(int n) : SquareMatrix<T>(n)
    {
        identify();
    }

    void identify()
    {
        this->data = vector<vector<T>>(this->n);
        for (int i = 0; i < this->n; ++i) {
            this->data[i] = vector<T>(this->n);
            for (int j = 0; j < this->n; ++j) {
                if (i == j) {
                    this->data[i][j] = 1;
                } else {
                    this->data[i][j] = 0;
                }
            }
        }
    }
};

template<typename T>
class EliminationMatrix : public IdentityMatrix<T>
{
public:
    explicit EliminationMatrix(int size) : IdentityMatrix<T>(size)
    {}

    void setValues(Matrix<double> m, int i, int j)
    {
        if (i > this->n || j > this->n) {
            cout << "Size error" << endl;
        }
        for (int k = 0; k < this->n; ++k) {
            this->data[k][k] = 1;
        }
        this->data[i - 1][j - 1] = -(m.getElement(i - 1, j - 1) / m.getElement(j - 1, j - 1));
    }
};

template<typename T>
class PermutationMatrix : public IdentityMatrix<T>
{
public:
    explicit PermutationMatrix(int n) : IdentityMatrix<T>(n)
    {}

    void doPermutation(int i, int j)
    {
        vector<T> temp = this->data[i - 1];
        this->data[i - 1] = this->data[j - 1];
        this->data[j - 1] = temp;
    }
};

template<typename T>
class ColumnVector : public Matrix<T>
{
public:
    explicit ColumnVector(int n, vector<vector<double>> content) : Matrix<T>(n, 1, content)
    {
        this->rows = n;
        this->cols = 1;
    }
    void setElement(int i, T x)
    {
        this->data[i][0] = x;
    }
    T getElement(int i)
    {
        return this->data[i][0];
    }
    T norm()
    {
        T sum = 0;
        for (int i = 0; i < this->rows; ++i) {
            sum += pow(this->data[i][0], 2);
        }
        return sqrt(sum);
    }
};

template<typename T>
class AugmentedMatrix : public Matrix<T>
{
private:
    Matrix<T> a{1, 1};
    Matrix<T> b{1, 1};
    int step;
    int n{};
public:
    AugmentedMatrix(SquareMatrix<T> a, SquareMatrix<T> b, int n) : Matrix<T>(n, 2 * n)
    {
        this->a = a;
        this->b = b;
        this->n = n;
        this->rows = n;
        this->cols = 2 * n;
        this->data = vector<vector<T>>(n);
        this->initializeAugmented();
    }

    Matrix<T> getA() {
        return this->a;
    }

    Matrix<T> getB() {
        return this->b;
    }

    void initializeAugmented()
    {
        this->data = vector<vector<T>>(this->n);
        for (int i = 0; i < this->n; ++i) {
            this->data[i] = vector<T>(2 * this->n);
            for (int j = 0; j < this->n; ++j) {
                this->data[i][j] = a.getElement(i, j);
            }
            for (int j = 0; j < this->n; ++j) {
                this->data[i][this->n + j] = b.getElement(i, j);
            }
        }
    }

    void updateAugmented()
    {
        for (int i = 0; i < this->n; ++i) {
            for (int j = 0; j < this->n; ++j) {
                this->data[i][j] = a.getElement(i, j);
            }
            this->data[i][this->n] = b.getElement(i, 0);
        }
    }

    void directWay()
    {
        PermutationMatrix<double> p(n);
        EliminationMatrix<double> e(n);
        T max;
        int maxRow;
        this->step = 1;
        for (int c = 0; c < this->n; ++c) {
            max = this->a.getElement(c, c);
            maxRow = c;
            for (int i = c + 1; i < this->n; ++i) {
                if (abs(round(this->a.getElement(i, c) * 10000) / 10000) > abs(round(max * 10000) / 10000)) {
                    max = abs(this->a.getElement(i, c));
                    maxRow = i;
                }
            }
            if (c != maxRow && round(max * 10000) / 10000 != 0) {
                PermutationMatrix<double> temp_p(n);
                temp_p.doPermutation(c + 1, maxRow + 1);
                this->a = temp_p * this->a;
                this->b = temp_p * this->b;
                updateAugmented();
            }
            for (int i = c + 1; i < this->n; ++i) {
                if (round(a.getElement(c, c) * 10000) / 10000 == 0 || round(a.getElement(i, c) * 10000) / 10000 == 0
                    || a.getElement(i, c) == -(a.getElement(i, c) / a.getElement(c, c))) {
                    continue;
                }
                EliminationMatrix<double> temp_e(n);
                temp_e.setValues(this->a, i + 1, c + 1);
                this->a = temp_e * this->a;
                this->b = temp_e * this->b;
                updateAugmented();
            }
        }
    }

    void wayBack()
    {
        for (int i = n - 1; i > 0; --i) {
            for (int j = i - 1; j >= 0; --j) {
                if (round(a.getElement(i, i) * 10000) / 10000 == 0 || round(a.getElement(j, i) * 10000) / 10000 == 0
                    || a.getElement(j, i) == -(a.getElement(j, i) / a.getElement(i, i))) {
                    continue;
                }
                EliminationMatrix<double> e(this->n);
                e.setValues(this->a, j + 1, i + 1);
                this->a = e * this->a;
                this->b = e * this->b;
                this->updateAugmented();
//                cout << "step #" << step++ << ": elimination" << endl;
//                cout << *this;
            }
        }
    }

    void diagNorm()
    {
        for (int i = 0; i < this->n; ++i) {
            T d = this->a.getElement(i, i);
            this->a.setElement(i, i, 1);
            for (int j = 0; j < this->n; ++j) {
                this->b.setElement(i, j, this->b.getElement(i, j) / d);
            }
        }
        this->updateAugmented();
//        cout << "Diagonal normalization:" << endl;
//        cout << *this;
    }
};

int main()
{
    FILE* pipe = _popen(GNUPLOT_LSA, "w");
    int m, n;
    cin >> m;
    vector<vector<double>> bv(m);
    vector<double> in(m);
    for (int i = 0; i < m; ++i) {
        cin >> in[i];
        bv[i] = *new vector<double>(1);
        cin >> bv[i][0];
    }
    cin >> n;
    vector<vector<double>> mContent(m);
    for (int i = 0; i < m; ++i) {
        mContent[i] = *new vector<double>(n+1);
        for (int j = 0; j <= n; ++j) {
            mContent[i][j] = pow(in[i], j);
        }
    }
    ColumnVector<double> b(m, bv);
    Matrix<double> A(m, n+1, mContent);
    cout << "A:" << endl << A;
    Matrix<double> A_T = A.transpose();
    SquareMatrix<double> m_AT_A(n+1, A_T * A);
    cout << "A_T*A:" << endl << m_AT_A;
    SquareMatrix<double> *I = new IdentityMatrix<double>(n+1);
    AugmentedMatrix<double> inv_mAT_A(m_AT_A, *I, n+1);
    inv_mAT_A.directWay();
    inv_mAT_A.wayBack();
    inv_mAT_A.diagNorm();
    cout << "(A_T*A)^-1:" << endl << inv_mAT_A.getB();
    Matrix<double> m_AT_b = A_T * b;
    cout << "A_T*b:" << endl << m_AT_b;
    Matrix<double> mulX = inv_mAT_A.getB() * m_AT_b;
    cout << "x~:" << endl << mulX;
    fprintf(pipe, "%s\n", "set multiplot");
    fprintf(pipe, "%s\n", "set xrange [-10:10]");
    fprintf(pipe, "%s\n", "set yrange [-10:10]");
    string ans = "";
    for (int i = n; i >= 0; --i) {
        if (i == n) {
            ans += to_string(mulX.getElement(i, 0)) + "*x**" + to_string(i);
        } else if (i == 0) {
            if (mulX.getElement(0, 0) < 0) {
                ans += " - " + to_string(abs(mulX.getElement(0, 0)));
            } else {
                ans += " + " + to_string(mulX.getElement(0, 0));
            }
        } else {
            if (mulX.getElement(i, 0) < 0) {
                ans += " - " + to_string(abs(mulX.getElement(i, 0))) + "*x**" + to_string(i);
            } else {
                ans += " + " + to_string(mulX.getElement(i, 0)) + "*x**" + to_string(i);
            }
        }
    }
    fprintf(pipe, "%s %s %s\n", "plot ", ans.c_str(), "with lines linecolor 'blue'");
    fprintf(pipe, "%s\n", "plot '-' using 1:2 title '' with points pointtype 6 linecolor 'red'");
    for (int i = 0; i < m; ++i) {
        fprintf(pipe, "%lf %lf\n", in[i], bv[i][0]);
    }
    fprintf(pipe, "%s\n", "e");
    fprintf(pipe, "%s\n", "unset multiplot");
    fflush(pipe);
    _pclose(pipe);
    return 0;
}
