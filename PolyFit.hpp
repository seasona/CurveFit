//@note:	利用最小二乘法进行多项式曲线拟合，给出几个关键点坐标即可求出近似函数
//@Author:	Season
//@date:	2019/2/26

#ifndef POLYFIT_H_
#define POLYFIT_H_
#include <vector>

namespace curvefit {

//坐标点表示
template<typename T>
struct Point {
	T x;
	T y;
};

class PolyFit {
public:
	PolyFit() = default;
	~PolyFit() = default;

	template<typename T>
	void Fit(const std::vector<Point<T>> &key_points, int poly_n);

	//@berif 利用最小二乘法进行多项式拟合
	//@param key_points 关键点坐标数组
	//@param length		关键点个数
	//@param poly_n		拟合次数，次数越高越精确
	template<typename T>
	void Fit(const Point<T> *key_points, size_t length, int poly_n);

	std::vector<double> GetFactor();
private:
	std::vector<double> factor;	//拟合后的方程系数

	template<typename T>
	void GaussSolve(int n, std::vector<T> &a, std::vector<T> &x, std::vector<T> &b);

	//@berif 高斯消元法求线性方程组解,即AX=B这种形式
	//@param n 矩阵的阶数
	template<typename T>
	void GaussSolve(int n, T *a, T *x, T *b);

};

template<typename T>
void PolyFit::GaussSolve(int n, std::vector<T> &a, std::vector<T> &x, std::vector<T> &b) {
	GaussSolve(n, &a[0], &x[0], &b[0]);
}

//TODO 这个函数需要处理特殊情况，如0
//利用列主消元法求解
template<typename T>
void PolyFit::GaussSolve(int n, T *a, T *x, T *b) {
	int i, j, k;
	double max;
	int index;
	for (k = 0; k < n - 1; k++) {
		max = std::fabs(a[k*n + k]);
		index = k;
		for (i = k; i <= n; i++) {
			if (std::fabs(a[i*n + i]) > max) {
				max = fabs(a[i*n + i]);
				index = i;
			}
		}
		if (index != k) {
			for (i = 0; i < n; i++) {
				max = a[k*n + i];
				a[k*n + i] = a[index*n + i];
				a[index*n + i] = max;
			}
		}
		max = b[k];
		b[k] = b[index];
		b[index] = max;
		for (i = k + 1; i < n; i++) {
			for (j = k + 1; j < n; j++) {
				a[i*n + j] -= a[i*n + k] * a[k*n + j] / a[k*n + k];
			}
			b[i] -= a[i*n + k] * b[k] / a[k*n + k];
		}
	}
	for (i = n - 1; i >= 0; x[i] /= a[i*n + i], i--) {
		for (j = i + 1, x[i] = b[i]; j < n; j++) {
			x[i] -= a[i*n + j] * x[j];
		}
	}
}

template<typename T>
void PolyFit::Fit(const std::vector<Point<T>> &key_points, int poly_n) {
	Fit(&key_points[0], key_points.size(), poly_n);
}

template<typename T>
void PolyFit::Fit(const Point<T> *key_points, size_t length, int poly_n) {
	factor.resize(poly_n + 1, 0);
	int i, j;
	std::vector<double> tempx(length, 1.0);
	std::vector<double> tempy;
	std::vector<double> sumxx(poly_n * 2 + 1);
	std::vector<double> ata((poly_n + 1)*(poly_n + 1));
	std::vector<double> sumxy(poly_n + 1);

	for (i = 0; i < length; i++) {
		tempy.push_back(key_points[i].y);
	}

	for (i = 0; i < 2 * poly_n + 1; i++) {
		for (sumxx[i] = 0, j = 0; j < length; j++)
		{
			sumxx[i] += tempx[j];
			tempx[j] *= key_points[j].x;
		}
	}
	for (i = 0; i < poly_n + 1; i++) {
		for (sumxy[i] = 0, j = 0; j < length; j++)
		{
			sumxy[i] += tempy[j];
			tempy[j] *= key_points[j].x;
		}
	}
	for (i = 0; i < poly_n + 1; i++)
		for (j = 0; j < poly_n + 1; j++)
			ata[i*(poly_n + 1) + j] = sumxx[i + j];
	GaussSolve(poly_n + 1, ata, factor, sumxy);
}

std::vector<double> PolyFit::GetFactor() {
	return factor;
}

}

#endif
