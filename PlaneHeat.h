#pragma once

#include <omp.h>
#include <functional>

#include "Mesh.h"
#include "Material.h"
#include "StaticMatrix.h"
#include "StripSLAE.h"


class PlaneHeat {

	Material m;

	const Mesh* mesh = nullptr;

	StaticMatrix<2, 3>* B = nullptr;

	double* S = nullptr;               // Массив площадей конечных элементов

	StripSLAE slae;                    // Система уравнений

	double* T = nullptr;

	bool meshLoaded = false;           // Загружена ли в класс сетка
	bool assembled = false;            // Собрана ли матрица жёсткости
	bool solved = false;               // Решена ли система уравнений

	// Вычислить матрицы градиентов
	void calcBs();

	// Вставить матрицу жёсткости элемента в глобальную
	void putKe(const StaticMatrix<3, 3>& Ke, unsigned e);

	// Готовы ли данные к заданию условий
	bool condReady() const;


public:

	PlaneHeat() = default;

	PlaneHeat(const Material& material, const Mesh* mesh)
		: m(material), mesh(mesh) {
		meshLoaded = true;
	}

	~PlaneHeat() {
		delete[] B;
		delete[] S;
		delete[] T;
	}

	// Задать сетку
	void setMesh(const Mesh* mesh);

	// Задать материал
	void setMaterial(const Material& material);

	// Собрать глобальную матрицу жёсткости
	// (если указатель аргумент, то она будет сохранена в файл)
	void assembly(const std::string & savingFile = "");

	// Решить задачу
	void solve();

	// Сохранить результаты в файл .vtk для визуализации в ParaView
	void saveAsVtk(const std::string& fileName, bool polarCoords = false) const;

	// Сохранить результаты в файл .vtu для визуализации в ParaView
	void saveAsVtu(const std::string& fileName, bool polarCoords = false) const;


	bool setBorderT(unsigned borderId, double T);
	
	bool setPointT(vec2 point, double T);

	bool setNodeT(unsigned nodeId, double T);

	bool setVertLineT(double x, double T);

	bool setHoriLineT(double y, double T);

	
	bool setBorderFlow(unsigned borderId, double T);
	
	bool setPointFlow(vec2 point, double flow);

	bool setNodeFlow(unsigned nodeId, double flow);

	bool setVertLineFlow(double x, double flow);

	bool setHoriLineFlow(double );

	bool setEdgeFlow();


	friend class PlaneElast;

};