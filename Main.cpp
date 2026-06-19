#include "PlaneElast.h"
#include "PlaneHeat.h"

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
	//mesh.saveAsVtu("data/BeamMesh.vtu");
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
	solver.saveAsVtu("data/CurvedBeam.vtu", true);
}


void pointForceTest() {
	Mesh mesh;
	mesh.genRectangle(-6, 6, 0, 3, 120, 30);

	double E = 2e11, mu = 0.3;
	Material m;
	m.E = E;
	m.mu = mu;

	PlaneElast solver(m, &mesh);
	solver.assembly();

	solver.setHoriLineUV(0.);
	solver.setVertLineU(-6.);
	solver.setVertLineU(6.);
	solver.setPointForce({ 0., 3. }, { 0., -1e8 });

	solver.nodalStressSaving = true;
	solver.solve();
	solver.saveAsVtu("data/PointContact.vtu");
}


void colorTest() {
	Mesh mesh;
	//mesh.loadFromFile("data/RectHole1.inp");
	//mesh.genRectangle(-6, 6, 0, 3, 120 * 30, 30 * 30);
	mesh.genRectangle(-6, 6, 0, 3, 120, 30);

	//mesh.renumByDirection();
	mesh.makeMultiColor();
	mesh.saveAsVtk("data/ColoredMesh.vtu");
}


void heatTest() {
	Mesh mesh;
	mesh.genRectangle(0, 1, 0, 1, 100, 100);

	double E = 2e11, mu = 0.3;
	Material m;
	m.E = E;
	m.mu = mu;
	m.lambda = 1.;
	m.alpha = 15e-6;

	PlaneHeat heat(m, &mesh);
	heat.assembly();
	heat.setHoriLineT(0., 300);
	heat.setPointT({ 0.5, 1. }, 400);
	//heat.setHoriLineT(1., 400);

	heat.solve();
	//heat.saveAsVtk("data/heat1.vtk");

	PlaneElast elast(m, &mesh);
	elast.assembly();
	elast.setHoriLineV(0.);
	elast.setVertLineU(0.);
	elast.setHoriLineV(1.);
	elast.setVertLineU(1.);
	elast.setTemperature(heat, 300);

	elast.nodalStressSaving = true;
	elast.solve();
	elast.saveAsVtu("data/TempStress.vtu");
}


void bladeTest1() {
	Mesh mesh;
	mesh.loadFromFile("data/Blade1.inp");
	mesh.renumByDirection();
	//mesh.saveAsVtu("data/Blade1.vtu");

	double E = 2e11, mu = 0.3;
	Material m;
	m.E = E;
	m.mu = mu;
	m.lambda = 1.;
	m.alpha = 15e-6;

	PlaneHeat heat(m, &mesh);
	heat.assembly();
	heat.setBorderT(0, 500);
	heat.setBorderT(1, 300);

	heat.solve();
	//heat.saveAsVtu("data/Blade1.vtu");

	PlaneElast elast(m, &mesh);
	elast.assembly();
	/*elast.setPointUV({ 0., 0. });
	elast.setPointV({ 60., 0. });*/
	elast.setNodeUV(0);
	elast.setNodeV(2000);
	elast.setTemperature(heat, 300);

	elast.nodalStressSaving = true;
	elast.solve();
	elast.saveAsVtu("data/Blade1.vtu");
}


void bladeTest2() {
	Mesh mesh;
	mesh.loadFromFile("data/Blade3.inp");
	mesh.renumByDirection({ 1., -1. });
	//mesh.saveAsVtu("data/Blade2.vtu");

	Material m;
	m.E = 10730;
	m.mu = 0.3;
	m.lambda = 9;
	m.alpha = 15e-6;

	PlaneHeat heat(m, &mesh);
	heat.assembly();
	heat.setBorderT(0, 904);
	heat.setBorderT(1, 972);
	heat.setBorderT(2, 870.4);
	heat.setBorderT(3, 877.8);
	heat.setBorderT(4, 877.3);
	heat.setBorderT(5, 872.8);
	heat.setBorderT(6, 848.3);
	heat.setBorderT(7, 817.6);
	heat.setBorderT(8, 812.6);
	heat.setBorderT(9, 807.6);
	heat.setBorderT(10, 800.6);

	/*heat.setBorderT(2, -170);
	heat.setBorderT(3, -170);
	heat.setBorderT(4, -170);
	heat.setBorderT(5, -170);
	heat.setBorderT(6, -170);
	heat.setBorderT(7, -170);
	heat.setBorderT(8, -170);
	heat.setBorderT(9, -170);*/

	heat.solve();
	//heat.saveAsVtu("data/Blade2.vtu");

	PlaneElast elast(m, &mesh);
	elast.assembly();
	//elast.setPointUV({ 0., 0. });
	//elast.setPointV({ 60., 0. });
	elast.setNodeUV(0);
	elast.setNodeV(30000);
	elast.setTemperature(heat, 0);

	elast.nodalStressSaving = true;
	elast.solve();
	elast.saveAsVtu("data/Blade3.vtu");
}


int main() {

	//meshTest();

	//test1();

	//pointForceTest();

	//colorTest();

	//heatTest();

	bladeTest2();

}