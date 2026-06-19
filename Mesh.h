#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <omp.h>

#include "vec2.h"


// Сетка треугольных конечных элементов
class Mesh {

private:

	double refSize = 1.;  // Характерный размер

	void realloc(unsigned nodeCount, unsigned elemCount) {
		Mesh::nodeCount = nodeCount;
		Mesh::elemCount = elemCount;
		delete[] node;
		delete[] elem3;
		node = new vec2[nodeCount];
		elem3 = new unsigned[3 * elemCount];
	}

	static void deleteCommas(std::string& str) {
		for (int i = 0; i < str.size(); ++i)
			if (str[i] == ',') str[i] = ' ';
	}

	void calcRefSize() {
		double max = 0.;
		for (unsigned i = 0; i < nodeCount; ++i) {
			if (fabs(node[i].x) > max)
				max = fabs(node[i].x);
			if (fabs(node[i].y) > max)
				max = fabs(node[i].y);
		}
		refSize = max;
	}

	char* elemColor = nullptr;

	bool colored = false;

public:

	unsigned nodeCount = 0;    // Количество узлов

	unsigned elemCount = 0;    // Количество элементов

	vec2* node = nullptr;      // Массиы узлов

	unsigned* elem3 = nullptr; // Элементы

	unsigned borderCount = 0;

	unsigned* borderOffset = nullptr;

	unsigned* borders = nullptr;


	// Конструктор по умолчанию
	Mesh() = default;

	// Конструктор по количеству узлов и элементов
	Mesh(unsigned nodeCount, unsigned elemCount)
		: nodeCount(nodeCount), elemCount(elemCount) {
		node = new vec2[nodeCount];
		elem3 = new unsigned[3 * elemCount];
	}

	// Конструктор копирования
	Mesh(const Mesh& mesh)
		: nodeCount(mesh.nodeCount), elemCount(mesh.elemCount), borderCount(mesh.borderCount) {
		node = new vec2[nodeCount];
		elem3 = new unsigned[3 * elemCount];
		borderOffset = new unsigned[borderCount + 1];
		memcpy(node, mesh.node, nodeCount * sizeof(vec2));
		memcpy(elem3, mesh.elem3, 3 * elemCount * sizeof(int));

	}

	// Конструктор перемещения
	Mesh(Mesh&& mesh) noexcept
		: nodeCount(mesh.nodeCount), elemCount(mesh.elemCount), borderCount(mesh.borderCount) {
		mesh.nodeCount = 0;
		mesh.elemCount = 0;
		mesh.borderCount = 0;
		node = mesh.node;
		elem3 = mesh.elem3;
		borderOffset = mesh.borderOffset;
		borders = mesh.borders;
		mesh.node = nullptr;
		mesh.elem3 = nullptr;
		mesh.borderOffset = nullptr;
		mesh.borders = nullptr;
	}

	// Деструктор
	~Mesh() {
		delete[] node; node = nullptr;
		delete[] elem3; elem3 = nullptr;
		delete[] borderOffset; borderOffset = nullptr;
		delete[] borders; borders = nullptr;
	}

	// Обращение к координатам узла элемента
	vec2& elemNode(unsigned element, unsigned localNode) {
		return node[elem3[3 * element + localNode]];
	}

	// Константное обращение к координатам узла элемента
	const vec2& elemNode(unsigned element, unsigned localNode) const {
		return node[elem3[3 * element + localNode]];
	}

	// Номер узла конечного элемента
	unsigned elemNodeId(unsigned element, unsigned localNode) const {
		return elem3[3 * element + localNode];
	}

	unsigned borderLength(unsigned borderId) const {
		return borderOffset[borderId + 1] - borderOffset[borderId];
	}


	// Найти наибольшую разность номеров узлов в конечном элементе
	int findMaxIndexDiff() const {
		int max = 0;
		for (int e = 0; e < elemCount; ++e)
			for (int i = 3 * e; i < 3 * (e + 1); ++i)
				for (int j = i + 1; j < 3 * (e + 1); ++j) {
					int diff = abs((long)elem3[i] - (long)elem3[j]);
					if (diff > max)
						max = diff;
				}
		return max;
	}

	// Получить характерный размер модели
	double getRefSize() const {
		return refSize;
	}

	// Сгенерировать сетку в прямоугольной области
	void genRectangle(double x1, double x2, double y1, double y2, unsigned nx, unsigned ny);

	// Перенумеровать узлы вдоль направления
	void renumByDirection(vec2 direction = {1., 0.01});

	// Загрузить сетку из файла
	// Поддерживаются только файлы .inp для Abaqus
	bool loadFromFile(const std::string& fileName);

	// Сохранить сетку в файл .vtk для визуализации в Paraview
	void saveAsVtk(const std::string& fileName) const;

	// Сохранить сетку в файл .vtu для визуализации в Paraview
	void saveAsVtu(const std::string& fileName) const;

	void makeMultiColorSlow();
	void makeMultiColor();

};