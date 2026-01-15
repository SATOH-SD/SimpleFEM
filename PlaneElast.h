#pragma once

#include <omp.h>
#include <functional>

#include "Mesh.h"
#include "Material.h"
#include "StaticMatrix.h"
#include "StripSLAE.h"


// Двухмерный тензор
struct Tensor2D {

	double xx, yy, xy;

	Tensor2D operator+=(const Tensor2D& tensor) {
		xx += tensor.xx;
		yy += tensor.yy;
		xy += tensor.xy;
		return *this;
	}

	Tensor2D operator*=(double k) {
		xx *= k;
		yy *= k;
		xy *= k;
		return *this;
	}

	// Поворот тензора в полярную систему координат в заданной точке
	Tensor2D polar(vec2 point) const {
		double r = point.norm();
		double cos_ = point.x / r, cos12 = point.y / r;
		Tensor2D polarTensor = {
			cos_* cos_* xx + cos12 * cos12 * yy + 2. * cos_ * cos12 * xy,   // rr
			cos12 * cos12 * xx + cos_ * cos_ * yy - 2. * cos_ * cos12 * xy, // phiphi
			cos_ * cos12 * (-xx + yy) + xy * (cos_ * cos_ - cos12 * cos12)  // rphi
		};
		return polarTensor;
	}

	// Интенсивность тензора
	double intensity() const {
		double d = xx - yy;
		return sqrt(0.5 * (xx * xx + yy * yy + d * d + 6. * xy * xy));
	}

};

// Плоское напряжённое или деформиованное состояние
enum class planeCond { stress, strain };


// Решение плоской задачи упругости
class PlaneElast {

private:

	Material m;                    // Материал

	const Mesh* mesh = nullptr;    // Указатель на сетку

	planeCond pc = planeCond::stress;  // Плоское состояние

	StaticMatrix<3, 6>* B = nullptr;   // Массив матриц градиентов

	double* S = nullptr;               // Массив площадей конечных элементов

	StaticMatrix<3, 3> C;              // Матрица задачи

	StripSLAE slae;                    // Система уравнений

	double* uv = nullptr;              // Вектор перемещения (решение системы)

	double* intE = nullptr;            // Интенсивности деформаций (в элементах)
	double* intS = nullptr;            // Интенсивности напряжений (в элементах)

	Tensor2D* strain = nullptr;        // Деформации (в элементах)
	Tensor2D* stress = nullptr;        // Напряжения (в элементах)

	double* nodalIntS = nullptr;       // Интенсивности напряжений (в узлах)
	Tensor2D* nodalStress = nullptr;   // Напряжения (в узлах)
	
	bool meshLoaded = false;           // Загружена ли в класс сетка
	bool assembled = false;            // Собрана ли матрица жёсткости
	bool solved = false;               // Решена ли система уравнений

	// Вычислить матрицы градиентов
	void calcBs();

	// Вставить матрицу жёсткости элемента в глобальную
	void putKe(const StaticMatrix<6, 6>& Ke, unsigned e);

	// Готовы ли данные к заданию условий нагружения
	bool condReady() const;

	// Вычислить деформации и напряжения
	void calcStress();

public:

	// Сохранение напряжений в узлах
	bool nodalStressSaving = false;

	// Конструктор по умолчанию
	PlaneElast() = default;

	// Конструктор по материалу и указателю на сетку
	PlaneElast(const Material& material, const Mesh* mesh)
		: m(material), mesh(mesh) {
		meshLoaded = true;
	}

	// Деструктор
	~PlaneElast() {
		delete[] B;
		delete[] S;
		delete[] uv;
		delete[] strain;
		delete[] stress;
		delete[] intE;
		delete[] intS;
		delete[] nodalStress;
		delete[] nodalIntS;
	}

	// Задать сетку
	void setMesh(const Mesh* mesh);

	// Задать материал
	void setMaterial(const Material& material);

	// Собрать глобальную матрицу жёсткости
	// (если указатель аргумент, то она будет сохранена в файл)
	void assembly(const std::string& savingFile = "");

	// Решить задачу
	void solve();

	// Сохранить результаты в файл .vtk для визуализации в ParaView
	void saveAsVtk(const std::string& fileName, bool polarCoords = false) const;


	// Задать перемещения в точке
	bool setPointUV(vec2 point, vec2 uv = { 0., 0. });

	// Задать горизонтальное перемещение в точке
	bool setPointU(vec2 point, double u = 0.);

	// Задать вертикальное перемещение в точке
	bool setPointV(vec2 point, double v = 0.);

	// Задать перемещения в узле
	bool setNodeUV(unsigned nodeId, vec2 uv = { 0., 0. });

	// Задать горизонтальное перемещение в узле
	bool setNodeU(unsigned nodeId, double u = 0.);

	// Задать вертикальное перемещение в узле
	bool setNodeV(unsigned nodeId, double v = 0.);

	// Задать перемещения на вертикальной линии
	bool setVertLineUV(double x, vec2 uv = { 0., 0. });

	// Задать горизонтальное перемещение на вертикальной линии
	bool setVertLineU(double x, double u = 0.);

	// Задать вертикальное перемещение на вертикальной линии
	bool setVertLineV(double x, double v = 0.);
	
	// Задать перемещения на горизонтальной линии
	bool setHoriLineUV(double y, vec2 uv = { 0., 0. });

	// Задать горизонтальное перемещение на горизонтальной линии
	bool setHoriLineU(double y, double u = 0.);

	// Задать вертикальное перемещение на горизонтальной линии
	bool setHoriLineV(double y, double v = 0.);

	// Задать силу в точке
	bool setPointForce(vec2 point, vec2 force);

	// Задать силу в узле
	bool setNodeForce(unsigned nodeId, vec2 force);

	// Задать силу на вертикальной линии
	bool setVertLineForce(double x, vec2 force);

	// Задать силу на горизонтальной линии
	bool setHoriLineForce(double y, vec2 force);

	// Задать силу на ребре элемента
	bool setEdgeForce(unsigned nodeId1, unsigned nodeId2, vec2 force);


	// Задать объёмные силы при вращении вокруг точки
	void setRotationPoint(vec2 point, double omega);

	// Задать объёмные силы при вращении вокруг оси
	void setRotationAxis(vec2 point, vec2 direction, double omega);
};