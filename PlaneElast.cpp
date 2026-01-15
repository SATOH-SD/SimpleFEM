#include "PlaneElast.h"


// Вычислить матрицу градиентов для элемента с координатами узлов в массиве x
static void calcB(StaticMatrix<3, 6>& B, double& S, const vec2* x) {
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

	B(0, 0) = B(2, 1) = b1 * S2_;
	B(1, 1) = B(2, 0) = c1 * S2_;
	B(0, 2) = B(2, 3) = b2 * S2_;
	B(1, 3) = B(2, 2) = c2 * S2_;
	B(0, 4) = B(2, 5) = b3 * S2_;
	B(1, 5) = B(2, 4) = c3 * S2_;
}


// Вычислить матрицы градиентов
void PlaneElast::calcBs() {
	delete[] B;
	B = new StaticMatrix<3, 6>[mesh->elemCount];
	delete[] S;
	S = new double[mesh->elemCount];
	for (unsigned e = 0; e < mesh->elemCount; ++e) {
		vec2 x[3] = { mesh->elemNode(e, 0), mesh->elemNode(e, 1), mesh->elemNode(e, 2) };
		calcB(B[e], S[e], x);
		//B[e].print();
	}
}


// Вставить матрицу жёсткости элемента в глобальную
void PlaneElast::putKe(const StaticMatrix<6, 6>& Ke, unsigned e) {
	for (unsigned i = 0; i < 3; ++i)
		for (unsigned j = 0; j < 3; ++j) {
			unsigned gi2 = 2 * mesh->elem3[3 * e + i];  //глобальне индексы
			unsigned gj2 = 2 * mesh->elem3[3 * e + j];
			
			unsigned i2 = 2 * i, j2 = 2 * j;
			slae(gi2, gj2) += Ke(i2, j2);
			slae(gi2 + 1, gj2) += Ke(i2 + 1, j2);
			slae(gi2, gj2 + 1) += Ke(i2, j2 + 1);
			slae(gi2 + 1, gj2 + 1) += Ke(i2 + 1, j2 + 1);
		}
}

// Собрать глобальную матрицу жёсткости
// (если указатель аргумент, то она будет сохранена в файл)
void PlaneElast::assembly(const std::string& savingFile) {
	double t = -omp_get_wtime();

	slae.resize(2 * mesh->nodeCount, 2 * mesh->findMaxIndexDiff() + 1);
	std::cout << "\nK size: " << slae.size();
	std::cout << "\nStrip width: " << slae.width() << "\n";

	double coef;
	switch (pc) {
	case planeCond::stress:
		coef = m.E / (1. - m.mu * m.mu);
		C(0, 0) = C(1, 1) = coef;
		C(1, 0) = C(0, 1) = coef * m.mu;
		C(2, 2) = coef * 0.5 * (1. - m.mu);
		break;
	case planeCond::strain:
		coef = m.E / ((1. + m.mu) * (1. - 0.5 * m.mu));
		C(0, 0) = C(1, 1) = coef * (1. - m.mu);
		C(1, 0) = C(0, 1) = coef * m.mu;
		C(2, 2) = m.E * 0.5 / (1. + m.mu);
		break;
	}
	//C.print();

	calcBs();
	slae.clear();
	for (unsigned e = 0; e < mesh->elemCount; ++e) {
		StaticMatrix<3, 6>& B = PlaneElast::B[e];
		StaticMatrix<6, 3> BC;
		for (unsigned k = 0; k < 3; ++k)
			for (unsigned i = 0; i < 6; ++i)
				for (unsigned j = 0; j < 3; ++j)
					BC(i, j) += B(k, i) * C(k, j);
		StaticMatrix<6, 6> Ke;
		for (unsigned i = 0; i < 6; ++i)
			for (unsigned k = 0; k < 3; ++k)
				for (unsigned j = 0; j < 6; ++j)
					Ke(i, j) += BC(i, k) * B(k, j);
		for (unsigned i = 0; i < 36; ++i)
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


// Вычислить деформации и напряжения
void PlaneElast::calcStress() {
	delete[] strain;
	delete[] stress;
	delete[] intE;
	delete[] intS;
	strain = new Tensor2D[mesh->elemCount];
	stress = new Tensor2D[mesh->elemCount];
	intE = new double[mesh->elemCount];
	intS = new double[mesh->elemCount];
	
	auto sqr = [](double x) { return x * x; };

	for (unsigned e = 0; e < mesh->elemCount; ++e) {
		unsigned begin = 3 * e;
		/*double uve[6] = {
			uv[2 * mesh->elem3[begin]],
			uv[2 * mesh->elem3[begin] + 1],
			uv[2 * mesh->elem3[begin + 1]],
			uv[2 * mesh->elem3[begin + 1] + 1],
			uv[2 * mesh->elem3[begin + 2]],
			uv[2 * mesh->elem3[begin + 2] + 1]
		};*/
		double
			u1 = uv[2 * mesh->elem3[begin]],
			v1 = uv[2 * mesh->elem3[begin] + 1],
			u2 = uv[2 * mesh->elem3[begin + 1]],
			v2 = uv[2 * mesh->elem3[begin + 1] + 1],
			u3 = uv[2 * mesh->elem3[begin + 2]],
			v3 = uv[2 * mesh->elem3[begin + 2] + 1];
		strain[e] = {};
		StaticMatrix<3, 6>& Be = B[e];

		strain[e].xx = u1 * Be(0, 0) + u2 * Be(0, 2) + u3 * Be(0, 4);
		strain[e].yy = v1 * Be(1, 1) + v2 * Be(1, 3) + v3 * Be(1, 5);
		strain[e].xy = u1 * Be(2, 0) + v1 * Be(2, 1) + u2 * Be(2, 2) + \
			v2 * Be(2, 3) + u3 * Be(2, 4) + v3 * Be(2, 5);

		stress[e].xx = C(0, 0) * strain[e].xx + C(0, 1) * strain[e].yy;
		stress[e].yy = C(0, 1) * strain[e].xx + C(1, 1) * strain[e].yy;
		stress[e].xy = C(2, 2) * strain[e].xy;

		intS[e] = stress[e].intensity();
		intE[e] = intS[e] / m.E;
	}
	if (nodalStressSaving) {
		unsigned* links = new unsigned[mesh->nodeCount];
		memset(links, 0, mesh->nodeCount * sizeof(unsigned));
		nodalStress = new Tensor2D[mesh->nodeCount];
		memset(nodalStress, 0, mesh->nodeCount * sizeof(Tensor2D));
		nodalIntS = new double[mesh->nodeCount];
		memset(nodalIntS, 0, mesh->nodeCount * sizeof(double));
		for (unsigned e = 0; e < mesh->elemCount; ++e)
			for (unsigned i = 0; i < 3; ++i) {
				unsigned gi = mesh->elemNodeId(e, i);
				++links[gi];
				nodalStress[gi] += stress[e];
				nodalIntS[gi] += intS[e];
			}
		for (unsigned i = 0; i < mesh->nodeCount; ++i) {
			double coef = 1. / links[i];
			nodalStress[i] *= coef;
			nodalIntS[i] *= coef;
		}
		delete[] links;
	}
}


// Решить задачу
void PlaneElast::solve() {
	std::cout << "\nSolving...\n";
	double t = -omp_get_wtime();

	uv = new double[2 * mesh->nodeCount];
	slae.Gauss(uv);
	calcStress();

	t += omp_get_wtime();
	std::cout << "Solved [" << t << " sec]\n";
	solved = true;
}


// Сохранить результаты в файл .vtk для визуализации в ParaView
void PlaneElast::saveAsVtk(const std::string& fileName, bool polarCoords) const {
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
	file << "\nCELL_DATA " << mesh->elemCount << "\n";
	file << "SCALARS ElementID int 1\n";
	file << "LOOKUP_TABLE default\n";
	for (unsigned i = 0; i < mesh->elemCount; ++i) {
		file << i << "\n";
	}
	file << "FIELD FieldData " << 8 + polarCoords * 6;
	file << "\nEps_xx 1 " << mesh->elemCount << " float\n";
	for (unsigned i = 0; i < mesh->elemCount; ++i)
		file << strain[i].xx << " ";
	file << "\nEps_yy 1 " << mesh->elemCount << " float\n";
	for (unsigned i = 0; i < mesh->elemCount; ++i)
		file << strain[i].yy << " ";
	file << "\nGamma_xy 1 " << mesh->elemCount << " float\n";
	for (unsigned i = 0; i < mesh->elemCount; ++i)
		file << strain[i].xy << " ";
	file << "\nEps_intensity 1 " << mesh->elemCount << " float\n";
	for (unsigned i = 0; i < mesh->elemCount; ++i)
		file << intE[i] << " ";
	file << "\nSigma_xx 1 " << mesh->elemCount << " float\n";
	for (unsigned i = 0; i < mesh->elemCount; ++i)
		file << stress[i].xx << " ";
	file << "\nSigma_yy 1 " << mesh->elemCount << " float\n";
	for (unsigned i = 0; i < mesh->elemCount; ++i)
		file << stress[i].yy << " ";
	file << "\nTau_xy 1 " << mesh->elemCount << " float\n";
	for (unsigned i = 0; i < mesh->elemCount; ++i)
		file << stress[i].xy << " ";
	file << "\nSigma_intensity 1 " << mesh->elemCount << " float\n";
	for (unsigned i = 0; i < mesh->elemCount; ++i)
		file << intS[i] << " ";
	if (polarCoords) {
		Tensor2D* polarTensor = new Tensor2D[mesh->elemCount];
		for (unsigned i = 0; i < mesh->elemCount; ++i) {
			strain[i].xy *= 0.5;
			vec2 r = (mesh->elemNode(i, 0) + mesh->elemNode(i, 0) + mesh->elemNode(i, 0)) / 3.;
			polarTensor[i] = strain[i].polar(r);
			polarTensor[i].xy *= 2.;
			strain[i].xy *= 2.;
		}
		file << "\nEps_rr 1 " << mesh->elemCount << " float\n";
		for (unsigned i = 0; i < mesh->elemCount; ++i)
			file << polarTensor[i].xx << " ";
		file << "\nEps_phiphi 1 " << mesh->elemCount << " float\n";
		for (unsigned i = 0; i < mesh->elemCount; ++i)
			file << polarTensor[i].yy << " ";
		file << "\nGamma_rphi 1 " << mesh->elemCount << " float\n";
		for (unsigned i = 0; i < mesh->elemCount; ++i)
			file << polarTensor[i].xy << " ";
		for (unsigned i = 0; i < mesh->elemCount; ++i) {
			vec2 r = (mesh->elemNode(i, 0) + mesh->elemNode(i, 0) + mesh->elemNode(i, 0)) / 3.;
			polarTensor[i] = stress[i].polar(r);
		}
		file << "\nSigma_rr 1 " << mesh->elemCount << " float\n";
		for (unsigned i = 0; i < mesh->elemCount; ++i)
			file << polarTensor[i].xx << " ";
		file << "\nSigma_phiphi 1 " << mesh->elemCount << " float\n";
		for (unsigned i = 0; i < mesh->elemCount; ++i)
			file << polarTensor[i].yy << " ";
		file << "\nTau_rphi 1 " << mesh->elemCount << " float\n";
		for (unsigned i = 0; i < mesh->elemCount; ++i)
			file << polarTensor[i].xy << " ";
		delete[] polarTensor;
	}

	file << "\nPOINT_DATA " << mesh->nodeCount << \
		"\nSCALARS NodeID int 1\nLOOKUP_TABLE my_table";
	for (unsigned i = 0; i < mesh->nodeCount; ++i)
		file << "\n" << i;
	file << "\nFIELD FieldData2 " << 1 + nodalStressSaving * 4 + polarCoords * 3 << "\nDisplacement 3 " << mesh->nodeCount << " float\n";
	for (unsigned i = 0; i < mesh->nodeCount; ++i)
		file << uv[2 * i] << " " << uv[2 * i + 1] << " 0 ";
	if (nodalStressSaving) {
		file << "\nNodal_sigma_xx 1 " << mesh->nodeCount << " float\n";
		for (unsigned i = 0; i < mesh->nodeCount; ++i)
			file << nodalStress[i].xx << " ";
		file << "\nNodal_sigma_yy 1 " << mesh->nodeCount << " float\n";
		for (unsigned i = 0; i < mesh->nodeCount; ++i)
			file << nodalStress[i].yy << " ";
		file << "\nNodal_tau_xy 1 " << mesh->nodeCount << " float\n";
		for (unsigned i = 0; i < mesh->nodeCount; ++i)
			file << nodalStress[i].xy << " ";
		file << "\nNodal_sigma_intensity 1 " << mesh->nodeCount << " float\n";
		for (unsigned i = 0; i < mesh->nodeCount; ++i)
			file << nodalIntS[i] << " ";
		if (polarCoords) {
			Tensor2D* polarTensor = new Tensor2D[mesh->nodeCount];
			for (unsigned i = 0; i < mesh->nodeCount; ++i)
				polarTensor[i] = nodalStress[i].polar(mesh->node[i]);
			file << "\nNodal_sigma_rr 1 " << mesh->nodeCount << " float\n";
			for (unsigned i = 0; i < mesh->nodeCount; ++i)
				file << polarTensor[i].xx << " ";
			file << "\nNodal_sigma_phiphi 1 " << mesh->nodeCount << " float\n";
			for (unsigned i = 0; i < mesh->nodeCount; ++i)
				file << polarTensor[i].yy << " ";
			file << "\nNodal_tau_rphi 1 " << mesh->nodeCount << " float\n";
			for (unsigned i = 0; i < mesh->nodeCount; ++i)
				file << polarTensor[i].xy << " ";
			delete[] polarTensor;
		}
	}

	file.close();
	t += omp_get_wtime();
	std::cout << "Solution saved in file " << fileName << " [" << t << " sec]\n";
}