#pragma once

// Свойства материала
class Material {

public:

	double E = 1.;   // Модуль упругости
	double mu = 0.;  // Коэффициент Пуассона

	double lambda = 1.;  // Коэффициент теплопроводности
	double alpha = 0.;   // Коэффициент температурного расширения
};