#include "PlaneHeat.h"

// Готовы ли данные к заданию условий нагружения
bool PlaneHeat::condReady() const {
	if (!meshLoaded) {
		std::cerr << "ERROR! Mesh is not loaded to solver! Condition skipped\n";
		return false;
	}
	if (!assembled) {
		std::cerr << "ERROR! Matrix is not assambled! Condition skipped\n";
		return false;
	}
	return true;
}

bool PlaneHeat::setBorderT(unsigned borderId, double T) {
	if (!condReady())
		return false;
	if (borderId >= mesh->borderCount) {
		std::cerr << "WARNING! There no border with ID " << borderId << "\n";
		return false;
	}
	unsigned begin = mesh->borderOffset[borderId];
	unsigned end = mesh->borderOffset[borderId + 1] - (mesh->borders[begin] == mesh->borderOffset[borderId + 1] - 1 ? 1 : 0);
	for (unsigned j = begin; j < end; ++j) {
		unsigned i = mesh->borders[j];
		slae(i, i) *= 1e20;
		slae.b(i) = slae(i, i) * T;
	}
	return true;
}

bool PlaneHeat::setPointT(vec2 point, double T) {
	if (!condReady())
		return false;
	double relZero = 1e-14 * mesh->getRefSize();
	for (int i = 0; i < mesh->nodeCount; ++i)
		if ((mesh->node[i] - point).norm() <= relZero) {
			slae(i, i) *= 1e20;
			slae.b(i) = slae(i, i) * T;
			std::cout << "Temperature " << T \
				<< " in point (" << point.x << ", " << point.y << ") applied\n";
			return true;
		}
	std::cerr << "WARNING! There is no node with coordinates (" \
		<< point.x << ", " << point.y << ") for temperature condition " \
		<< T << "\n";
	return false;
}

//bool PlaneHeat::setNodeT(unsigned nodeId, double T);

bool PlaneHeat::setVertLineT(double x, double T) {
	if (!condReady())
		return false;
	double relZero = 1e-14 * mesh->getRefSize();
	bool found = false;
	for (int i = 0; i < mesh->nodeCount; ++i)
		if (fabs(mesh->node[i].x - x) <= relZero) {
			slae(i, i) *= 1e20;
			slae.b(i) = slae(i, i) * T;
			found = true;
		}
	if (found) {
		std::cout << "Temperature " << T \
			<< " on line x = " << x << " applied\n";
		return true;
	}
	else {
		std::cerr << "WARNING! There are no nodes with x = " \
			<< x << " for temperature condition " \
			<< T << "\n";
		return false;
	}
}

bool PlaneHeat::setHoriLineT(double y, double T) {
	if (!condReady())
		return false;
	double relZero = 1e-14 * mesh->getRefSize();
	bool found = false;
	for (int i = 0; i < mesh->nodeCount; ++i)
		if (fabs(mesh->node[i].y - y) <= relZero) {
			slae(i, i) *= 1e20;
			slae.b(i) = slae(i, i) * T;
			found = true;
		}
	if (found) {
		std::cout << "Temperature " << T \
			<< " on line y = " << y << " applied\n";
		return true;
	}
	else {
		std::cerr << "WARNING! There are no nodes with y = " \
			<< y << " for temperature condition " \
			<< T << "\n";
		return false;
	}
}


//bool PlaneHeat::setPointFlow(vec2 point, double flow);

//bool PlaneHeat::setNodeFlow(unsigned nodeId, double flow);

//bool PlaneHeat::setVertLineFlow(double x, double flow);

//bool PlaneHeat::setHoriLineFlow(double);

//bool PlaneHeat::setEdgeFlow();