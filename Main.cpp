#include "PlaneElast.h"

void meshTest() {
	Mesh mesh;

	//mesh.genRectangle(-3., 3., -2., 2., 30, 20);
	//mesh.saveAsVtk("data/mesh1.vtk");

	mesh.loadFromFile("data/RectHole1.inp");
	mesh.saveAsVtk("data/RectHole.vtk");
}

//Задача Кирша (только отверстие кривое)
void test1() {
	Mesh mesh;
	mesh.loadFromFile("data/RectHole1.inp");
	mesh.renumByDirection();
	//mesh.makeMultiColor();
	//mesh.saveAsVtk("data/ColoredMesh.vtk");
	//mesh.saveAsVtk("data/RectHole.vtk");
	//mesh.genRectangle(-0.3, 0.3, -0.2, 0.2, 3, 2);

	Material m;
	m.E = 220e3;
	m.mu = 0.3;

	PlaneElast solver(m, &mesh);
	solver.assembly();
	//solver.assembly("data/matrix.txt");

	//solver.setVertLineU(-0.3, 0.);
	//solver.setPointV({ -0.3, 0. }, 0.);
	//solver.setVertLineU(0.3, 0.1);
	//solver.setPointV({ 0.3, 0. }, 0.);
	solver.setHoriLineV(0.);
	solver.setVertLineU(0.);
	solver.setVertLineForce(0.3, { 100., 0. });
	solver.setVertLineForce(-0.3, { -100., 0. });

	solver.nodalStressSaving = true;
	solver.solve();
	solver.saveAsVtk("data/RectHole.vtk", true);
}

//Балка
void test2() {
	double c = 0.1, h = 0.2, b = 0.01, l = 3.;
	Mesh mesh;
	mesh.genRectangle(0., l, -c, c, 600, 40);
	//mesh.makeMultiColor();

	double E = 2e11, mu = 0.3;
	Material m;
	m.E = E;
	m.mu = mu;

	PlaneElast solver(m, &mesh);
	solver.assembly();

	solver.setPointUV({ 0., 0. });
	solver.setPointUV({ l, 0. });
	//solver.setVertLineUV(0.);
	//solver.setVertLineUV(l);
	solver.setPointForce({ 0.5 * l, c }, { 0., -1e6 });

	solver.nodalStressSaving = true;
	solver.solve();
	solver.saveAsVtk("data/Beam1.vtk");
}

//Изогнутая балка
void test3() {
	Mesh mesh;
	mesh.loadFromFile("data/CurvedBeam.inp");
	mesh.renumByDirection({ 1., -1. });

	double E = 2e11, mu = 0.3;
	Material m;
	m.E = E;
	m.mu = mu;

	PlaneElast solver(m, &mesh);
	solver.assembly();

	solver.setVertLineU(0.);
	solver.setPointV({ 0., 2. });
	solver.setHoriLineForce(0., { -1e6, 0. });

	solver.nodalStressSaving = true;
	solver.solve();
	solver.saveAsVtk("data/CurvedBeam.vtk", true);
}


int main() {

	//meshTest();

	test1();

}