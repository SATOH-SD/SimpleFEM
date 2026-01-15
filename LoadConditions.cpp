#include "PlaneElast.h"


// Готовы ли данные к заданию условий нагружения
bool PlaneElast::condReady() const {
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


// Задать перемещения в точке
bool PlaneElast::setPointUV(vec2 point, vec2 uv) {
	if (!condReady())
		return false;
	double relZero = 1e-14 * mesh->getRefSize();
	for (int i = 0; i < mesh->nodeCount; ++i)
		if ((mesh->node[i] - point).norm() <= relZero) {
			unsigned i2 = 2 * i;
			slae(i2, i2) *= 1e20;
			slae.b(i2) = slae(i2, i2) * uv.x;
			slae(i2 + 1, i2 + 1) *= 1e20;
			slae.b(i2 + 1) = slae(i2 + 1, i2 + 1) * uv.y;
			std::cout << "Displacement (" << uv.x << ", " << uv.y \
				<< ") in point (" << point.x << ", " << point.y << ") applied\n";
			return true;
		}
	std::cerr << "WARNING! There is no node with coordinates (" \
		<< point.x << ", " << point.y << ") for displacement condition (" \
		<< uv.x << ", " << uv.y << ")\n";
	return false;
}


// Задать горизонтальное перемещение в точке
bool PlaneElast::setPointU(vec2 point, double u) {
	if (!condReady())
		return false;
	double relZero = 1e-14 * mesh->getRefSize();
	for (int i = 0; i < mesh->nodeCount; ++i)
		if ((mesh->node[i] - point).norm() <= relZero) {
			unsigned i2 = 2 * i;
			slae(i2, i2) *= 1e20;
			slae.b(i2) = slae(i2, i2) * u;
			std::cout << "Horizontal displacement " << u << " in point (" \
				<< point.x << ", " << point.y << ") applied\n";
			return true;
		}
	std::cerr << "WARNING! There is no node with coordinates (" \
		<< point.x << ", " << point.y << ") for horizontal displacement condition " \
		<< u << "\n";
	return false;
}


// Задать вертикальное перемещение в точке
bool PlaneElast::setPointV(vec2 point, double v) {
	if (!condReady())
		return false;
	double relZero = 1e-14 * mesh->getRefSize();
	for (int i = 0; i < mesh->nodeCount; ++i)
		if ((mesh->node[i] - point).norm() <= relZero) {
			unsigned i2 = 2 * i;
			slae(i2 + 1, i2 + 1) *= 1e20;
			slae.b(i2 + 1) = slae(i2 + 1, i2 + 1) * v;
			std::cout << "Vertical displacement " << v << " in point (" \
				<< point.x << ", " << point.y << ") applied\n";
			return true;
		}
	std::cerr << "WARNING! There is no node with coordinates (" \
		<< point.x << ", " << point.y << ") for vertical displacement condition " \
		<< v << "\n";
	return false;
}

// Задать перемещения в узле
//bool PlaneElast::setNodeUV(unsigned nodeId, vec2 uv = { 0., 0. }) {}

// Задать горизонтальное перемещение в узле
//bool PlaneElast::setNodeU(unsigned nodeId, double u = 0.) {}

// Задать вертикальное перемещение в узле
//bool PlaneElast::setNodeV(unsigned nodeId, double v = 0.) {}


// Задать перемещения на вертикальной линии
bool PlaneElast::setVertLineUV(double x, vec2 uv) {
	if (!condReady())
		return false;
	double relZero = 1e-14 * mesh->getRefSize();
	bool found = false;
	for (int i = 0; i < mesh->nodeCount; ++i)
		if (fabs(mesh->node[i].x - x) <= relZero) {
			unsigned i2 = 2 * i;
			slae(i2, i2) *= 1e20;
			slae.b(i2) = slae(i2, i2) * uv.x;
			slae(i2 + 1, i2 + 1) *= 1e20;
			slae.b(i2 + 1) = slae(i2 + 1, i2 + 1) * uv.y;
			found = true;
		}
	if (found) {
		std::cout << "Displacement (" << uv.x << ", " << uv.y \
			<< ") on line x = " << x << " applied\n";
		return true;
	}
	else {
		std::cerr << "WARNING! There are no nodes with x = " \
			<< x << " for displacement condition (" \
			<< uv.x << ", " << uv.y << ")\n";
		return false;
	}
}


// Задать горизонтальное перемещение на вертикальной линии
bool PlaneElast::setVertLineU(double x, double u) {
	if (!condReady())
		return false;
	double relZero = 1e-14 * mesh->getRefSize();
	bool found = false;
	for (int i = 0; i < mesh->nodeCount; ++i)
		if (fabs(mesh->node[i].x - x) <= relZero) {
			unsigned i2 = 2 * i;
			slae(i2, i2) *= 1e20;
			slae.b(i2) = slae(i2, i2) * u;
			found = true;
		}
	if (found) {
		std::cout << "Horizontal displacement " << u \
			<< " on line x = " << x << " applied\n";
		return true;
	}
	else {
		std::cerr << "WARNING! There are no nodes with x = " \
			<< x << " for horizontal displacement condition " \
			<< u << "\n";
		return false;
	}
}


// Задать вертикальное перемещение на вертикальной линии
bool PlaneElast::setVertLineV(double x, double v) {
	if (!condReady())
		return false;
	double relZero = 1e-14 * mesh->getRefSize();
	bool found = false;
	for (int i = 0; i < mesh->nodeCount; ++i)
		if (fabs(mesh->node[i].x - x) <= relZero) {
			unsigned i2 = 2 * i;
			slae(i2 + 1, i2 + 1) *= 1e20;
			slae.b(i2 + 1) = slae(i2 + 1, i2 + 1) * v;
			found = true;
		}
	if (found) {
		std::cout << "Vertical displacement " << v \
			<< " on line x = " << x << " applied\n";
		return true;
	}
	else {
		std::cerr << "WARNING! There are no nodes with x = " \
			<< x << " for vertical displacement condition " \
			<< v << "\n";
		return false;
	}
}


// Задать перемещения на горизонтальной линии
bool PlaneElast::setHoriLineUV(double y, vec2 uv) {
	if (!condReady())
		return false;
	double relZero = 1e-14 * mesh->getRefSize();
	bool found = false;
	for (int i = 0; i < mesh->nodeCount; ++i)
		if (fabs(mesh->node[i].y - y) <= relZero) {
			unsigned i2 = 2 * i;
			slae(i2, i2) *= 1e20;
			slae.b(i2) = slae(i2, i2) * uv.x;
			slae(i2 + 1, i2 + 1) *= 1e20;
			slae.b(i2 + 1) = slae(i2 + 1, i2 + 1) * uv.y;
			found = true;
		}
	if (found) {
		std::cout << "Displacement (" << uv.x << ", " << uv.y \
			<< ") on line y = " << y << " applied\n";
		return true;
	}
	else {
		std::cerr << "WARNING! There are no nodes with y = " \
			<< y << " for displacement condition (" \
			<< uv.x << ", " << uv.y << ")\n";
		return false;
	}
}


// Задать горизонтальное перемещение на горизонтальной линии
bool PlaneElast::setHoriLineU(double y, double u) {
	if (!condReady())
		return false;
	double relZero = 1e-14 * mesh->getRefSize();
	bool found = false;
	for (int i = 0; i < mesh->nodeCount; ++i)
		if (fabs(mesh->node[i].y - y) <= relZero) {
			unsigned i2 = 2 * i;
			slae(i2, i2) *= 1e20;
			slae.b(i2) = slae(i2, i2) * u;
			found = true;
		}
	if (found) {
		std::cout << "Horizontal displacement " << u \
			<< " on line y = " << y << " applied\n";
		return true;
	}
	else {
		std::cerr << "WARNING! There are no nodes with y = " \
			<< y << " for horizontal displacement condition " \
			<< u << "\n";
		return false;
	}
}


// Задать вертикальное перемещение на горизонтальной линии
bool PlaneElast::setHoriLineV(double y, double v) {
	if (!condReady())
		return false;
	double relZero = 1e-14 * mesh->getRefSize();
	bool found = false;
	for (int i = 0; i < mesh->nodeCount; ++i)
		if (fabs(mesh->node[i].y - y) <= relZero) {
			unsigned i2 = 2 * i;
			slae(i2 + 1, i2 + 1) *= 1e20;
			slae.b(i2 + 1) = slae(i2 + 1, i2 + 1) * v;
			found = true;
		}
	if (found) {
		std::cout << "Vertical displacement " << v \
			<< " on line y = " << y << " applied\n";
		return true;
	}
	else {
		std::cerr << "WARNING! There are no nodes with y = " \
			<< y << " for vertical displacement condition " \
			<< v << "\n";
		return false;
	}
}


// Задать силу в точке
bool PlaneElast::setPointForce(vec2 point, vec2 force) {
	if (!condReady())
		return false;
	double relZero = 1e-14 * mesh->getRefSize();
	for (int i = 0; i < mesh->nodeCount; ++i)
		if ((mesh->node[i] - point).norm() <= relZero) {
			unsigned i2 = 2 * i;
			slae.b(i2) += force.x;
			slae.b(i2 + 1) += force.y;
			std::cout << "Force (" << force.x << ", " << force.y \
				<< ") in point (" << point.x << ", " << point.y << ") applied\n";
			return true;
		}
	std::cerr << "WARNING! There is no node with coordinates (" \
		<< point.x << ", " << point.y << ") for force condition (" \
		<< force.x << ", " << force.y << ")\n";
	return false;
}

// Задать силу в узле
//bool PlaneElast::setNodeForce(unsigned nodeId, vec2 force) {}


// Задать силу на вертикальной линии
bool PlaneElast::setVertLineForce(double x, vec2 force) {
	if (!condReady())
		return false;
	double relZero = 1e-14 * mesh->getRefSize();
	unsigned edgeCount = 0;
	for (unsigned e = 0; e < mesh->elemCount; ++e) {
		unsigned edgeNodesCount = 0;
		unsigned edgeNodes[3];
		for (unsigned i = 0; i < 3; ++i)
			if (fabs(mesh->elemNode(e, i).x - x) <= relZero)
				edgeNodes[edgeNodesCount++] = mesh->elemNodeId(e, i);
		if (edgeNodesCount == 2) {
			vec2 forceValue = force * (mesh->node[edgeNodes[0]] - mesh->node[edgeNodes[1]]).norm();
			slae.b(2 * edgeNodes[0]) += forceValue.x;
			slae.b(2 * edgeNodes[0] + 1) += forceValue.y;
			slae.b(2 * edgeNodes[1]) += forceValue.x;
			slae.b(2 * edgeNodes[1] + 1) += forceValue.y;
			++edgeCount;
		}
	}
	if (edgeCount == 0) {
		std::cerr << "WARNING! There are no edges on line with coordinate x = " \
			<< x << " for force condition (" << force.x << ", " << force.y << ")\n";
		return false;
	}
	else {
		std::cout << "Force (" << force.x << ", " << force.y \
			<< ") on line x = " << x << " applied\n";
		return true;
	}
}


// Задать силу на горизонтальной линии
bool PlaneElast::setHoriLineForce(double y, vec2 force) {
	if (!condReady())
		return false;
	double relZero = 1e-14 * mesh->getRefSize();
	unsigned edgeCount = 0;
	for (unsigned e = 0; e < mesh->elemCount; ++e) {
		unsigned edgeNodesCount = 0;
		unsigned edgeNodes[3];
		for (unsigned i = 0; i < 3; ++i)
			if (fabs(mesh->elemNode(e, i).y - y) <= relZero)
				edgeNodes[edgeNodesCount++] = mesh->elemNodeId(e, i);
		if (edgeNodesCount == 2) {
			vec2 forceValue = force * (mesh->node[edgeNodes[0]] - mesh->node[edgeNodes[1]]).norm();
			slae.b(2 * edgeNodes[0]) += forceValue.x;
			slae.b(2 * edgeNodes[0] + 1) += forceValue.y;
			slae.b(2 * edgeNodes[1]) += forceValue.x;
			slae.b(2 * edgeNodes[1] + 1) += forceValue.y;
			++edgeCount;
		}
	}
	if (edgeCount == 0) {
		std::cerr << "WARNING! There are no edges on line with coordinate y = " \
			<< y << " for force condition (" << force.x << ", " << force.y << ")\n";
		return false;
	}
	else {
		std::cout << "Force (" << force.x << ", " << force.y \
			<< ") on line y = " << y << " applied\n";
		return true;
	}
}


// Задать силу на ребре элемента
//bool PlaneElast::setEdgeForce(unsigned nodeId1, unsigned nodeId2, vec2 force) {}


// Задать объёмные силы при вращении вокруг точки
//void PlaneElast::setRotationPoint(vec2 point, double omega) {}

// Задать объёмные силы при вращении вокруг оси
//void PlaneElast::setRotationAxis(vec2 point, vec2 direction, double omega) {}