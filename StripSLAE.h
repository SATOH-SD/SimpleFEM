#pragma once

#include <iostream>
#include <fstream>
#include <vector>

//Ленточная СЛАУ
//(матрица хранится и обрабатывается только в области ширины ленты)
class StripSLAE {

public:

	size_t N = 0; //размер СЛАУ

	size_t W = 1; //ширина ленты

	double* matrix = nullptr;  //массив элементов матрицы

	double* rp = nullptr;      //вектор правой части

	StripSLAE() = default;

	StripSLAE(size_t size, size_t width)
		: N(size), W(width) {
		rp = new double[N];
		matrix = new double[N * (2 * W + 1)];
	};

	StripSLAE(const StripSLAE& slae)
		: N(slae.N), W(slae.W) {
		rp = new double[N];
		matrix = new double[N * (2 * W + 1)];
		memcpy(matrix, slae.matrix, N * (2 * W + 1) * sizeof(double));
		memcpy(rp, slae.rp, N * sizeof(double));
	}

	~StripSLAE() {
		delete[] matrix;
		delete[] rp;
	}

	void resize(size_t size, size_t width) {
		N = size;
		W = width;
		delete[] matrix;
		delete[] rp;
		rp = new double[N];
		matrix = new double[N * (2 * W + 1)];
	}

	//Обращение к элементу матрицы
	//(сохраняется логика индексации полной матрицы)
	inline double& operator()(size_t i, size_t j) {
#ifdef _DEBUG
		if (i > N || j > i + W || long(j) < long(i) - long(W)) {
			std::cout << "!!! " << i << " " << j << "\n";
			return matrix[0];
		}
#endif
		return matrix[i * (2 * W + 1) + j + W - i];
	}

	//Константное обращение к элементу матрицы
	inline const double& operator()(size_t i, size_t j) const {
		return matrix[i * (2 * W + 1) + j + W - i];
	}

	double* data() {
		return matrix;
	}

	double* bData() {
		return rp;
	}

	//Вернуть размер СЛАУ
	size_t size() const {
		return N;
	}

	//Ширина ленты (в одну сторону)
	size_t width() const {
		return W;
	}

	//Обращение к правой части
	double& b(size_t i) {
		return rp[i];
	}

	//Константное обращение к правой части
	const double& b(size_t i) const {
		return rp[i];
	}

	//Очистить систему
	void clear() {
		/*for (size_t i = 0; i < N * (2 * W + 1); ++i)
			matrix[i] = 0.;
		for (size_t i = 0; i < N; ++i)
			rp[i] = 0.;*/

		memset(matrix, 0, N * (2 * W + 1) * sizeof(double));
		memset(rp, 0, N * sizeof(double));
	}

	//Вычислить невязку
	double calcDiscrepancy(const std::vector<double>& solution) const {
		std::vector<double> r(N);
		std::copy(rp, rp + N, r.begin());
		for (size_t i = 0; i < N; ++i) {
			double sum = 0.;
			size_t lineBegin = W > i ? W - i : 0;
			size_t lineEnd = 2 * W < W + N - i - 1 ? 2 * W : W + N - i - 1;
			for (size_t j = lineBegin; j <= lineEnd; ++j)
				sum += matrix[i * (2 * W + 1) + j] * solution[i - W + j];
			r[i] = sum - b(i);
		}
		double norm = 0.;
		for (size_t i = 0; i < N; ++i) {
			norm += r[i] * r[i];
		}
		return sqrt(norm);
	}

	//Вывод на экран для отладки
	void print() const {
		for (size_t i = 0; i < N; ++i)
			std::cout << rp[i] << "\n";
		for (size_t i = 0; i < N; ++i) {
			for (size_t j = 0; j <= 2 * W; ++j)
				std::cout << matrix[i * (2 * W + 1) + j] << " ";
			std::cout << "\n";
		}
	}

	bool saveAsMatrix(const std::string& fileName) const;

	bool saveAsStrip(const std::string& fileName) const {
		std::ofstream file(fileName, std::ios_base::out);
		for (size_t i = 0; i < N; ++i) {
			for (size_t j = 0; j <= 2 * W; ++j)
				file << matrix[i * (2 * W + 1) + j] << " ";
			file << "\n";
		}
		file.close();
		return true;
	}

	//Метод Гаусса, при решении меняет саму СЛАУ
	void Gauss(double* solution) {
		//прямой ход
		for (int k = 0; k < N; ++k) {
			double leadEl = operator()(k, k);
			size_t lineEnd = std::min((size_t)N, k + W + 1);
			//деление на ведущий элемент
#pragma omp parallel for
			for (int j = k + 1; j < lineEnd; ++j)
				operator()(k, j) /= leadEl;
			b(k) /= leadEl;
#pragma omp parallel for
			for (int i = k + 1; i < lineEnd; ++i) {
				for (int j = k + 1; j < lineEnd; ++j)
					operator()(i, j) -= operator()(k, j) * operator()(i, k);
				b(i) -= b(k) * operator()(i, k);
			}
		}
		//обратный ход
		for (int i = (int)N - 2; i >= 0; --i) {
			size_t lineEnd = std::min(N, i + W + 1);
			for (int j = (int)lineEnd - 1; j > i; --j)
				b(i) -= b(j) * operator()(i, j);
		}
		//std::vector<double> solution(N);
		
		//std::copy(rp, rp + N, solution);
		memcpy(solution, rp, N * sizeof(double));
	}


	double matrixHash() const {
		double sum = 0.;
		for (int i = 0; i < N * (2 * W + 1); ++i)
			sum += matrix[i];
		return sum;
	}

	double rpHash() const {
		double sum = 0.;
		for (int i = 0; i < N; ++i)
			sum += rp[i];
		return sum;
	}

};