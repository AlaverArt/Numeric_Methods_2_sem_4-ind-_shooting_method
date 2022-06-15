#include<iostream>
#include<fstream>
#include<vector>
#include<functional>
#include<iomanip>
using namespace std;

//ex1 linear n = 2, f(t)=0
//constexpr int n = 2;
//double n0 = 0.1;//y(0.5)=n0=0.1//расчёт в двух точках
//double n1 = 0.5;//y(0.5)=n1=0.5
//double y_diff_0 = 0.3;
//double y_end = 0.3;
//constexpr double t_0 = 0.5;//start
//constexpr double t_end = 3.5;//end

//ex2 linear system n = 2, f(t)!=0
//constexpr int n = 2;
//double n0 = 0.1;//y(0.5)=n0=0.1//расчёт в двух точках
//double n1 = 0.5;//y(0.5)=n1=0.5
//double y_diff_0 = 0.73418997228521;
//double y_end = -9.25;
//constexpr double t_0 = -8.0;//start
//constexpr double t_end = -3.0;//end

//ex3 linear system n = 2, f(t)=0
constexpr int n = 2;
double n0 = 0.1;//y(0.5)=n0=0.1//расчёт в двух точках
double n1 = 0.5;//y(0.5)=n1=0.5
double y_diff_0 = -8.3;
double y_end = 13.4851671;
constexpr double t_0 = 1.0;//start
constexpr double t_end = 7.0;//end

constexpr double tau = 0.05;//step
constexpr double b0 = 0.5;
constexpr double b1 = 0.5;

constexpr double epsSimpleIter = 0.0000000001;
constexpr double EPS = 0.00001;

vector<double> f(double t, vector<double> u) {
	vector<double> res(n);

	//ex1 695 FilipovDU
	//res[0] = u[1];//y'=z
	//res[1] = ((2 * t + 4) * u[1] - 2 * u[0])/(t*(t+4));//z'= ((2x+4)z-2y) / (x(x+4))

	//ex2 704 FilipovDU
	//res[0] = u[1];//y'=z
	//res[1] = (6 * t - 4 * t * u[1] - 2 * u[0]) / (t * t - 1);//z' = (2x - 4xz - 2y) / (x^2 - 1)

	//ex3 702 FilipovDU
	res[0] = u[1];//y'=z
	res[1] = ( t + 1/t - (t+2)*u[1] + u[0] ) / ( t*(t+1) );//z' = ( x + 1/x -(x+2)z + y ) / ( x(x+1) )

	return res;
}

vector<double> operator*(double a, vector<double> vect) {
	vector<double> res(vect.size());

	for (int i = 0; i < vect.size(); i++)
		res[i] = vect[i] * a;

	return res;
}

vector<double> operator+(vector<double> vect1, vector<double> vect2) {
	vector<double> res;
	if (vect1.size() != vect2.size()) return res;
	res.resize(vect1.size());

	for (int i = 0; i < vect1.size(); i++)
		res[i] = vect1[i] + vect2[i];

	return res;
}

vector<double> operator-(vector<double> vect1, vector<double> vect2) {
	vector<double> res;
	if (vect1.size() != vect2.size()) return res;
	res.resize(vect1.size());

	for (int i = 0; i < vect1.size(); i++)
		res[i] = vect1[i] - vect2[i];

	return res;
}

double getNorm(vector<double> vect) {
	//Norm l_infinity
	double norm = abs(vect[0]);

	for (int i = 1; i < vect.size(); i++)
		if (abs(vect[i]) > norm)
			norm = abs(vect[i]);

	return norm;
}

vector<vector<double>> adamsMethod_simpleIters(std::function<vector<double>(double, vector<double>)> f, vector<double> u0, double t0, double end_t) {
	int numb_un = (int)(abs(end_t - t0) / tau) + 1;
	vector<vector<double>> u(numb_un, vector<double>(n));

	for (int i = 0; i < n; i++)
	{
		u[0][i] = u0[i];
	}

	vector<double> F(u.size());
	vector<double> u_pred(n);
	vector<double> fn_sub_1(n);
	vector<double> fn(n);
	fn = f(t0, u[0]);
	double tn = t0;
	for (int k = 1; k < numb_un; k++) {
		tn += tau;
		//computing yn
		fn_sub_1 = f(tn - tau, u[k - 1]);
		u[k] = tau * b1 * fn_sub_1 + u[k - 1];
		int nkk = 0;
		do {
			u_pred = u[k];
			fn_sub_1 = fn;
			fn = f(tn, u_pred);
			u[k] = tau * (b0 * fn + b1 * fn_sub_1) + u[k - 1];
			//F = u[k] - tau * (b0 * fn + b1 * fn_sub_1) - u[k - 1];
		} while (getNorm(u[k] - u_pred) > epsSimpleIter /* || getNorm(F) > epsSimpleIter*/);
	}

	return u;
}

vector<vector<double>> shootingMethod(std::function<vector<double>(double, vector<double>)> func, double t0,
	double end_t, double y_diff_0, double n0, double n1) {
	int numb_un = (int)(abs(end_t - t0) / tau) + 1;
	vector<vector<double>> u(numb_un, vector<double>(n));
	vector<vector<double>> y_find_history(2, vector<double>(numb_un));
	vector<double> y(numb_un);
	double y_0 = n0;
	double y_0_prev;
	double y_0_prev_prev;
	double bigF;
	double bigF_prev;
	double bigF_prev_prev;

	double shoots_count = 0;

	y_0_prev = n0;
	vector<double> u0 = { y_0_prev, y_diff_0 };//start u(t)=(u(t), u'(t))
	u = adamsMethod_simpleIters(func, u0, t0, end_t);
	bigF_prev = u[numb_un - 1][0] - y_end;
	shoots_count++;

	cout << "(" << shoots_count << ") " << y_0_prev << endl;
	//y_find_history.resize(y_find_history.size() + 1);
	for (int i = 0; i < u.size(); i++)
		y_find_history[y_find_history.size() - 2][i] = u[i][0];

	y_0 = n1;
	u0 = { y_0, y_diff_0 };//start u(t)=(u(t), u'(t))
	u = adamsMethod_simpleIters(func, u0, t0, end_t);
	bigF = u[numb_un - 1][0] - y_end;
	shoots_count++;

	cout << "(" << shoots_count << ") " << y_0 << endl;
	//y_find_history.resize(y_find_history.size() + 1);
	for (int i = 0; i < u.size(); i++)
		y_find_history[y_find_history.size() - 1][i] = u[i][0];

	//метод простой итерации для функции bigF(y_0)
	do {
		bigF_prev_prev = bigF_prev;
		bigF_prev = bigF;

		y_0_prev_prev = y_0_prev;
		y_0_prev = y_0;

		y_0 = y_0_prev - (y_0_prev - y_0_prev_prev) * bigF_prev / (bigF_prev - bigF_prev_prev);

		u0 = { y_0, y_diff_0 };//start u(t)=(u(t), u'(t))
		u = adamsMethod_simpleIters(func, u0, t0, end_t);
		
		bigF = u[numb_un - 1][0] - y_end;
		shoots_count++;

		cout << "(" << shoots_count << ") " << setprecision(10) << y_0 << endl;
		y_find_history.resize(y_find_history.size() + 1, vector<double>(numb_un));
		for (int i = 0; i < u.size(); i++)
			y_find_history[y_find_history.size() - 1][i] = u[i][0];

	} while (abs(y_0 - y_0_prev) > EPS && abs(bigF) > EPS);//Добавил условие bigF > EPS

	for (int i = 0; i < u.size(); i++) {
		y[i] = u[i][0];
	}

	return y_find_history;
}

int main() {
	ofstream fout("output.txt");
	ofstream fout_history("history.txt");
	vector<vector<double>> y_find_history = shootingMethod(f, t_0, t_end, y_diff_0, n0, n1);
	for (int i = 0; i < y_find_history[y_find_history.size() - 1].size(); i++) {

		fout << y_find_history[y_find_history.size() - 1][i] << endl;
		//cout << y_find_history[y_find_history.size() - 1][i] << endl;
	}

	for (int j = 0; j < y_find_history.size(); j++) {
		fout_history << "y" << j + 1;
		if (j != y_find_history.size() - 1)
			fout_history << "\t";
	}
	fout_history << endl;

	for (int i = 0; i < y_find_history[0].size(); i++) {
		for (int j = 0; j < y_find_history.size(); j++) {
			if (j == 0)
				fout_history << y_find_history[j][i];
			else
				fout_history << "\t" << y_find_history[j][i];
		}
		fout_history << endl;
	}

	return 0;
}
