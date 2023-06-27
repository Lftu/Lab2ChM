#include <bits/stdc++.h>
#include <locale>;
#include "windows.h";

using namespace std;

const int n = 7;
const int N = 7;
const double e = 0.00001;

void inversion(double **A)
{
    double temp;

    double **E = new double *[N];

    for (int i = 0; i < N; i++)
        E[i] = new double [N];

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            E[i][j] = 0.0;

            if (i == j)
                E[i][j] = 1.0;
        }

    for (int k = 0; k < N; k++)
    {
        temp = A[k][k];

        for (int j = 0; j < N; j++)
        {
            A[k][j] /= temp;
            E[k][j] /= temp;
        }

        for (int i = k + 1; i < N; i++)
        {
            temp = A[i][k];

            for (int j = 0; j < N; j++)
            {
                A[i][j] -= A[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }

    for (int k = N - 1; k > 0; k--)
    {
        for (int i = k - 1; i >= 0; i--)
        {
            temp = A[i][k];

            for (int j = 0; j < N; j++)
            {
                A[i][j] -= A[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            A[i][j] = E[i][j];

    for (int i = 0; i < N; i++)
        delete [] E[i];

    delete [] E;
}

void PrintMatrix(double **A, double *F)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            cout << A[i][j] << " ";
        cout << F[i];
        cout << endl;
    }
}
void PrintMatrix(double **A)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            cout << A[i][j] << " ";
        cout << endl;
    }
}
void PrintVec(double *A)
{
    for (int i = 0; i < n; i++)
    {
        cout << A[i] << " ";
    }
    cout << endl;
}
void Print(double *A)
{
    for (int i = 0; i < n; i++)
    {
        cout <<"x"<<i+1<<"="<< fixed << setprecision(8) << A[i] << " ";
    }
    cout << endl;
}
bool converge(double *xk, double *xkp, int n)
{
    double norm = 0;
    for (int i = 0; i < n; i++)
        norm += (xk[i] - xkp[i])*(xk[i] - xkp[i]);
    return (sqrt(norm)<e);
}


double *Zeydel(double **a,double *b, double*x)
{
    double *p;
    p = new double[n];
    double *TempX = new double[n];
    for (int k = 0; k < n; k++)
        TempX[k] = x[k];
    int cnt = 0;
    do
    {
        for (int i = 0; i < n; i++)
            p[i] = x[i];
        for (int i = 0; i < n; i++)
        {
            double var = 0;
            for (int j = 0; j < i; j++)
                var += (a[i][j] * x[j]);
            for (int j = i + 1; j < n; j++)
                var += (a[i][j] * p[j]);
            x[i] = (b[i] - var) / a[i][i];
        }
        double norm = 0;
        for (int i = 0; i < n; i++)
            norm += (x[i] - p[i])*(x[i] - p[i]);
        cout << "Норма = " << sqrt(norm) << endl;
        cnt++;
    }
    while (!converge(x, p, n));
    double *nev = new double[n];
    cout << "Вектор нев'язки: (";
    for (int i = 0; i < n; i++)
    {
        nev[i] = -b[i];
        for (int k = 0; k < n; k++)
            nev[i] += a[i][k] * x[k];
        cout << nev[i];
        if (i != n - 1)
            cout << ", ";
        else
            cout << ")" << endl;
    }
    cout << "Кількість ітерацій методу Зейделя = " << cnt << endl;
    delete[] TempX;
    return x;
}

void input(double **&A, double *&F) /// A - Matrix
{
    F = new double[n];
    A = new double *[n];
    for (int i = 0; i < n; i++)
        A[i] = new double[n];
    double n1 = n;
    for (int i = 1; i <= n; i++)
    {
        double i1 = i;
        for(int j = 1; j <= n; j++)
        {
            double j1 = j;
            if(i == j)
            {
                A[i - 1][j - 1] = n1 / 2 + 10 + i1 / n1 + j1 / 10;
            }
            else
            {
                A[i - 1][j - 1] = (i1 * j1) / (10 + n1 * n1);
            }
        }
        F[i - 1] = 2 * i1 * i1 - n1 * n1 / 2;
    }
}
bool checkMatrix(double **A, double *F) /// Достатня умова збіжності
{
    for(int i = 0; i < n; i++)
    {
        int sum = 0;
        for(int j = 0; j < n; j++)
        {
            if(i != j)
            {
                sum += abs(A[i][j]);
            }
        }
        if(abs(A[i][i]) < sum)
        {
            return false;
        }
    }
    return true;
}

void print(int n, int t, double **z)// вывести матрицу
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < t; j++)
        {
            cout << setw(15) << z[i][j];
        }
        cout << endl;
    }
}

double** transp(double** a, int n, int m)// транспонировать матрицу
{
    int i, j;
    double **b;
    b = new double *[n];
    for (i = 0; i< n; i++)
    {
        b[i] = new double[m];
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
            b[j][i] = a[i][j];
    }
    return b;
}

int sign(double z)
{
    if (z > 0)
        return 1;
    else
        return -1;
}

double Norm(double **a)
{
    double temp = 0, norm_m = 0;
    for(int i = 0; i < n; i++)
    {
        temp = 0;
        for(int j = 0; j < n; j++)
            temp += fabs(a[i][j]);
        if(temp > norm_m)
            norm_m = temp;
    }
    return norm_m;
}

void SquareRoot(double **a, double *b)
{
    double **at = new double*[n], **rab = new double*[n], **s = new double*[n], **d = new double*[n];
    double opred = 1, *nev = new double[n], *y = new double[n], *x = new double[n], *b1 = new double[n];
    int i, j, k;
    bool flag = true;
    double kst;
    for (i = 0; i < n; i++)
    {
        d[i] = new double[n];
        s[i] = new double[n];
        rab[i] = new double[n];
        at[i] = new double[n];
        for (j = 0; j < n; j++)
        {
            if (a[i][j] != a[j][i])
                flag = false;
            d[i][j] = 0;
            s[i][j] = 0;
        }
    }
    if (!flag)
    {
        cout << "Матриця не симетрична, необхідно домножити на транспоновану до неї" << endl;
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
                at[j][i] = a[i][j];
        }
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                rab[i][j] = 0;
                b1[i] = 0;
                for (k = 0; k < n; k++)
                {
                    rab[i][j] += at[i][k] * a[k][j];
                    b1[i] += at[i][k] * b[k];
                }
            }
        }
    }
    else
    {
        cout << "Матриця симетрична" << endl;
        for (i = 0; i < n; i++)
        {
            b1[i] = b[i];
            for (j = 0; j < n; j++)
                rab[i][j] = a[i][j];
        }
    }
    cout << "Робоча матриця:" << endl;
    print(n, n, rab);
    for (int i = 0; i< n; i++)
    {
        for (int k = 0; k < (i + 1); k++)
        {
            double sum = 0;
            for (int j = 0; j < k; j++)
                sum += s[i][j] * s[k][j]*d[j][j];
            if (i == k)
            {
                s[i][k] = sqrt(abs(rab[i][i] - sum));
                d[i][k] = sign(abs(rab[i][i] - sum));
            }
            else
                s[i][k]= (1.0 / s[k][k] * (rab[i][k] - sum));
        }
    }
    s = transp(s, n, n);
    cout << "Mатриця S:" << endl;
    print(n, n, s);
    cout << "Вектор Y: (";
    for (i = 0; i < n; i++)
    {
        double sum = 0;
        for (int k = 0; k <= i - 1; k++)
            sum += y[k] * s[k][i];
        y[i] = (b1[i] - sum) / s[i][i];
        cout << y[i];
        if (i != n - 1)
            cout << ", ";
        else
            cout << ")";
    }
    cout << endl;
    cout << "Вектор X: (";
    for (i = n - 1; i >= 0; i--)
    {
        double sum = 0;
        for (int k = i + 1; k <= n - 1; k++)
            sum += s[i][k] * x[k];
        x[i] = (y[i] - sum) / s[i][i];
    }
    for (i = 0; i < n; i++)
    {
        opred *= s[i][i] * s[i][i] * d[i][i];
        cout << x[i];
        if (i != n - 1)
            cout << ", ";
        else
            cout << ")";
    }
    cout << endl;
    cout << "Визначник матриці: "<<opred<<endl;
    cout << "Вектор нев'язки: (";
    for (i = 0; i < n; i++)
    {
        nev[i] = -b[i];
        for (k = 0; k < n; k++)
            nev[i] += a[i][k] * x[k];
        cout << nev[i];
        if (i != n - 1)
            cout << ", ";
        else
            cout << ")" << endl;
    }
    double n1 = Norm(a);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
            rab[i][j] = a[i][j];
    }
    inversion(rab);
    double n2 = Norm(rab);
    cout << "Матриця A^-1:" << endl;
    print(n, n, rab);
    cout << "Норма матриці A:" << endl;
    cout << n1 << endl;
    cout << "Норма матриці A^-1:" << endl;
    cout << n2 << endl;
    cout << "Число обумовленості матриці A:" << endl;
    cout << n1 * n2 << endl;
}

int main()
{
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
    setlocale (LC_CTYPE, "ukr");
    double **Matrix, *b, *y, *x;
    input(Matrix, b);

    cout << "1) Метод квадратного кореня: " << endl;

    SquareRoot(Matrix, b);

    cout << "2) Метод Зейделя: " << endl;
    if(checkMatrix(Matrix, b))
    {
        x = new double[n]; /// Початкові оцінки
        for (int i = 0; i < n; i++)
            x[i] = 1.0;
        y = new double[n];
        Zeydel(Matrix, b, x);
        cout << "Результат: ";
        Print(x);
    }
    else
    {
        cout << "Достатня умова збіжності методу Зейделя не виконується" << endl;
    }
    system("pause");
    return 0;
}
