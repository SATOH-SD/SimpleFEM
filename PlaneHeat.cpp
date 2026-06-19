#include "PlaneHeat.h"

// Вычислить матрицу градиентов для элемента с координатами узлов в массиве x
static void calcB(StaticMatrix<2, 3>& B, double& S, const vec2* x) {
	double x1 = x[0].x, y1 = x[0].y,
		x2 = x[1].x, y2 = x[1].y,
		x3 = x[2].x, y3 = x[2].y;
	//std::cout << "\n" << x1 << " " << y1 << "  " << x2 << " " << y2 << "  " << x3 << " " << y3;
	double
		a1 = x2 * y3 - x3 * y2,
		a2 = x3 * y1 - x1 * y3,
		a3 = x1 * y2 - x2 * y1;
	double
		b1 = y2 - y3, c1 = x3 - x2,
		b2 = y3 - y1, c2 = x1 - x3,
		b3 = y1 - y2, c3 = x2 - x1;
	S = 0.5 * (a1 + a2 + a3);
	double S2_ = 2. / S; // 2 обратные площади

	B(0, 0) = b1 * S2_;
	B(0, 1) = b2 * S2_;
	B(0, 2) = b3 * S2_;
	B(1, 0) = c1 * S2_;
	B(1, 1) = c2 * S2_;
	B(1, 2) = c3 * S2_;
}


// Вычислить матрицы градиентов
void PlaneHeat::calcBs() {
	delete[] B;
	B = new StaticMatrix<2, 3>[mesh->elemCount];
	delete[] S;
	S = new double[mesh->elemCount];
	for (unsigned e = 0; e < mesh->elemCount; ++e) {
		vec2 x[3] = { mesh->elemNode(e, 0), mesh->elemNode(e, 1), mesh->elemNode(e, 2) };
		calcB(B[e], S[e], x);
		//B[e].print();
	}
}


// Вставить матрицу жёсткости элемента в глобальную
void PlaneHeat::putKe(const StaticMatrix<3, 3>& Ke, unsigned e) {
	for (unsigned i = 0; i < 3; ++i)
		for (unsigned j = 0; j < 3; ++j) {
			unsigned gi = mesh->elem3[3 * e + i];  //глобальне индексы
			unsigned gj = mesh->elem3[3 * e + j];
			slae(gi, gj) += Ke(i, j);
		}
}


// Собрать глобальную матрицу жёсткости
// (если указатель аргумент, то она будет сохранена в файл)
void PlaneHeat::assembly(const std::string& savingFile) {
	double t = -omp_get_wtime();

	slae.resize(mesh->nodeCount, mesh->findMaxIndexDiff() + 1);
	std::cout << "\nK size: " << slae.size();
	std::cout << "\nStrip width: " << slae.width() << "\n";

	calcBs();
	slae.clear();
	for (unsigned e = 0; e < mesh->elemCount; ++e) {
		StaticMatrix<2, 3>& B = PlaneHeat::B[e];
		StaticMatrix<3, 2> BC;
		for (unsigned k = 0; k < 2; ++k)
			for (unsigned i = 0; i < 3; ++i)
				for (unsigned j = 0; j < 2; ++j)
					BC(i, j) += B(k, i) * (k == j ? m.lambda : 0.);
		StaticMatrix<3, 3> Ke;
		for (unsigned i = 0; i < 3; ++i)
			for (unsigned k = 0; k < 2; ++k)
				for (unsigned j = 0; j < 3; ++j)
					Ke(i, j) += BC(i, k) * B(k, j);
		for (unsigned i = 0; i < 9; ++i)
			Ke.data[i] *= S[e];
		//Ke.print();
		//std::cout << "\n" << e << "\n";
		putKe(Ke, e);
	}

	t += omp_get_wtime();
	std::cout << "Assembled [" << t << " sec]\n";
	assembled = true;

	if (savingFile != "") {
		std::cout << "Saving...";
		slae.saveAsStrip(savingFile);
		std::cout << "\rMatrix saved in " << savingFile << "\n";
	}
	std::cout << std::endl;
}


// Решить задачу
void PlaneHeat::solve() {
	std::cout << "\nSolving...\n";
	double t = -omp_get_wtime();

	T = new double[mesh->nodeCount];
	slae.Gauss(T);

	t += omp_get_wtime();
	std::cout << "Solved [" << t << " sec]\n";
	solved = true;
}


// Сохранить результаты в файл .vtk для визуализации в ParaView
void PlaneHeat::saveAsVtk(const std::string& fileName, bool polarCoords) const {
	std::cout << "\nSaving...\n";
	double t = -omp_get_wtime();
	std::ofstream file(fileName, std::ios_base::out);
	file << "# vtk DataFile Version 2.0\n";
	file << "Solution\n";
	file << "ASCII\n";
	file << "DATASET POLYDATA\n";
	file << "POINTS " << mesh->nodeCount << " float\n";
	for (int i = 0; i < mesh->nodeCount; ++i)
		file << mesh->node[i].x << " " << mesh->node[i].y << " " << 0 << "\n";
	file << "POLYGONS " << mesh->elemCount << " " << 4 * mesh->elemCount;
	for (unsigned i = 0; i < mesh->elemCount; ++i) {
		file << "\n3 ";
		for (unsigned j = 3 * i; j < 3 * (i + 1); ++j)
			file << mesh->elem3[j] << " ";
	}

	file << "\nPOINT_DATA " << mesh->nodeCount << \
		"\nSCALARS NodeID int 1\nLOOKUP_TABLE my_table";
	for (unsigned i = 0; i < mesh->nodeCount; ++i)
		file << "\n" << i;
	file << "\nFIELD FieldData2 " << 1 << "\nTemperature 1 " << mesh->nodeCount << " float\n";
	for (unsigned i = 0; i < mesh->nodeCount; ++i)
		file << T[i] << " ";

	file.close();
	t += omp_get_wtime();
	std::cout << "Solution saved in file " << fileName << " [" << t << " sec]\n";
}


// Сохранить результаты в файл .vtu для визуализации в Paraview
void PlaneHeat::saveAsVtu(const std::string& fileName, bool polarCoords) const {
	std::cout << "\nSaving...\n";
	double t = -omp_get_wtime();
	std::ofstream file(fileName, std::ios_base::out);
	file << "<?xml version=\"1.0\"?>" \
		<< "\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\" >" \
		<< "\n\t<UnstructuredGrid>" \
		<< "\n\t\t<Piece NumberOfPoints=\"" << mesh->nodeCount << "\" NumberOfCells=\"" << mesh->elemCount << "\">";

	file << "\n\t\t\t<PointData>";  // Данные в узлах
	// Номера узлов
	file << "\n\t\t\t\t<DataArray type=\"UInt32\" Name=\"NodeID\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
	for (unsigned i = 0; i < mesh->nodeCount; ++i)
		file << i << " ";
	file << "\n\t\t\t\t</DataArray>";
	// Перемещения
	file << "\n\t\t\t\t<DataArray type=\"Float64\" Name=\"Temperature\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
	for (unsigned i = 0; i < mesh->nodeCount; ++i)
		file << T[i] << " ";
	file << "\n\t\t\t\t</DataArray>";
	file << "\n\t\t\t</PointData>";

	file << "\n\t\t\t<CellData>";  // Данные в элементах
	// Номера элементов
	file << "\n\t\t\t\t<DataArray type=\"UInt32\" Name=\"ElementID\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
	for (unsigned i = 0; i < mesh->elemCount; ++i)
		file << i << " ";
	file << "\n\t\t\t\t</DataArray>";
	file << "\n\t\t\t</CellData>";

	file << "\n\t\t\t<Points>";  // Узлы
	file << "\n\t\t\t\t<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" >\n";
	for (unsigned i = 0; i < mesh->nodeCount; ++i)
		file << mesh->node[i].x << " " << mesh->node[i].y << " 0 ";
	file << "\n\t\t\t\t</DataArray>";
	file << "\n\t\t\t</Points>";

	file << "\n\t\t\t<Cells>";  // Элементы
	file << "\n\t\t\t\t<DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\" >\n";
	for (int i = 0; i < mesh->elemCount * 3; ++i)
		file << mesh->elem3[i] << " ";
	file << "\n\t\t\t\t</DataArray>";
	file << "\n\t\t\t\t<DataArray type=\"UInt32\" Name=\"offsets\" format=\"ascii\" >\n";
	for (int i = 3; i <= mesh->elemCount * 3; i += 3)
		file << i << " ";
	file << "\n\t\t\t\t</DataArray>";
	file << "\n\t\t\t\t<DataArray type=\"UInt32\" Name=\"types\" format=\"ascii\" >\n";
	for (int i = 0; i < mesh->elemCount; ++i)
		file << "5 ";
	file << "\n\t\t\t\t</DataArray>";
	file << "\n\t\t\t</Cells>";

	file << "\n\t\t</Piece>";
	file << "\n\t</UnstructuredGrid>";
	file << "\n</VTKFile>";
	file.close();
	t += omp_get_wtime();
	std::cout << "Solution saved in file " << fileName << " [" << t << " sec]\n";
}