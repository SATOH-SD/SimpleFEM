#pragma once

#include <iostream>

template<int s1, int s2, typename fp = double>
class StaticMatrix {

public:

	fp data[s1 * s2] = {};

	//Вернуть размер матрицы по первому индексу
	int size1() const {
		return s1;
	}

	//Вернуть размер матрицы по второму индексу
	int size2() const {
		return s2;
	}

	int dataSize() const {
		return s1 * s2;
	}

	inline fp& operator()(int i, int j) {
		return data[i * s2 + j];
	}

	//Константное обращение к элементу матрицы
	inline const fp& operator()(int i, int j) const {
		return data[i * s2 + j];
	}

	//Вывод на экран для отладки
	void print() const;
};


//Вывод на экран для отладки
template<int s1, int s2, typename fp>
void StaticMatrix<s1, s2, fp>::print() const {
	printf("\n");
	for (int i = 0; i < s1; ++i) {
		for (int j = 0; j < s2; ++j) {
			//std::cout.width(width);
			printf("%12.5e ", data[i * s2 + j]);
		}
		printf("\n");
	}
	//std::cout << std::endl;
}