#include<iostream>
#include<vector>
#include<functional>
using namespace std;

//ex1 linear system n = 2, f(t)=0
constexpr int n = 2;
double n0 = 0.1;//y(0.5)=n0=0.1//расчёт в двух точках
double n1 = 0.5;//y(0.5)=n1=0.5
double y_diff_0 = 0.3;
double y_end = 0.3;
constexpr double t_0 = 0.5;//start
constexpr double t_end = 3.5;//end

constexpr double tau = 0.05;//step
constexpr double b0 = 0.5;
constexpr double b1 = 0.5;

constexpr double epsSimpleIter = 0.0000000001;

vector<double> f(double t, vector<double> u) {
	vector<double> res(n);

	//ex1 695 FilipovDU
	res[0] = u[1];//y'=z
	res[1] = ((2 * t + 4) * u[1] - 2 * u[0])/(t*(t+4));//z'= ((2x+4)z-2y) / (x(x+4))

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
	double norm = vect[0];

	for (int i = 1; i < vect.size(); i++)
		if (vect[i] > norm)
			norm = vect[i];

	return norm;
}

vector<vector<double>> adamsMethod_simpleIters(std::function<vector<double>(double, vector<double>)> f, vector<double> u0, double t0, double end_t) {
	int numb_un = (int)(abs(end_t - t0) / tau) + 1;
	vector<vector<double>> u(numb_un, vector<double>(n));

	for (int i = 0; i < n; i++)
	{
		u[0][i] = u0[i];
	}

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
		} while (getNorm(u[k] - u_pred) > epsSimpleIter);
	}

	return u;
}

vector<double> shootingMethod(std::function<vector<double>(double, vector<double>)> func, double t0, double end_t, double y_diff_0, double n0, double n1) {
	int numb_un = (int)(abs(end_t - t0) / tau) + 1;
	vector<vector<double>> u(numb_un, vector<double>(n));
	vector<double> y(numb_un);


	double y_0 = n0;
	vector<double> u0 = { y_0, y_diff_0 };//start u(t)=(u(t), u'(t))
	u = adamsMethod_simpleIters(func, u0, t0, end_t);
	double bigF_n0 = u[numb_un-1][0];
	cout << bigF_n0 << endl;

	y_0 = n1;
	u0 = { y_0, y_diff_0 };//start u(t)=(u(t), u'(t))
	u = adamsMethod_simpleIters(func, u0, t0, end_t);
	double bigF_n1 = u[numb_un - 1][0];
	cout << bigF_n1 << endl;

	return y;
}

int main() {
	
	vector<double> y = shootingMethod(f, t_0, t_end, y_diff_0, n0, n1);

	return 0;
}
