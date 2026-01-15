#include "Mesh.h"


// Сгенерировать сетку в прямоугольной области
void Mesh::genRectangle(double x1, double x2, double y1, double y2, unsigned nx, unsigned ny) {
	realloc((nx + 1) * (ny + 1), 2 * nx * ny);
	double hx = (x2 - x1) / nx, hy = (y2 - y1) / ny;
	for (unsigned i = 0; i <= nx; ++i)
		for (unsigned j = 0; j <= ny; ++j) {
			unsigned k = i * (ny + 1) + j;
			node[k] = { x1 + hx * i, y1 + hy * j };
		}
	for (unsigned i = 0; i < nx; ++i)
		for (unsigned j = 0; j < ny; ++j) {
			unsigned begin = 6 * (i * ny + j);
			elem3[begin] = i * (ny + 1) + j;
			elem3[begin + 1] = (i + 1) * (ny + 1) + j;
			elem3[begin + 2] = i * (ny + 1) + j + 1;

			elem3[begin + 3] = (i + 1) * (ny + 1) + j;
			elem3[begin + 4] = (i + 1) * (ny + 1) + j + 1;
			elem3[begin + 5] = i * (ny + 1) + j + 1;
		}
	colored = false;
	std::cout << "Mesh generated:\n";
	std::cout << nodeCount << " nodes, " << elemCount << " elements\n";
}


// Сортировка слиянием для перенумерации вдоль направления
static void renumSort(double* rBegin, double* rEnd, int* nBegin, double* bufferR, int* bufferN) {
	int dist = rEnd - rBegin;
	//std::cout << dist << "\n";
	if (dist != 1) {
		double* rMiddle = rBegin + dist / 2;
		renumSort(rBegin, rMiddle, nBegin, bufferR, bufferN);
		int half = rMiddle - rBegin;
		renumSort(rMiddle, rEnd, nBegin + half, bufferR + half, bufferN + half);
		//merge
		double* ir1 = rBegin, * ir2 = rMiddle, * ibr = bufferR;
		int* in1 = nBegin, * in2 = nBegin + half, * ibn = bufferN;
		while (ir1 != rMiddle && ir2 != rEnd) {
			if (*ir1 < *ir2) {
				*ibr++ = *ir1++;
				*ibn++ = *in1++;
			}
			else {
				*ibr++ = *ir2++;
				*ibn++ = *in2++;
			}
		}
		//copy rest to buffers
		memcpy(ibr, ir1, (rMiddle - ir1) * sizeof(double));
		memcpy(ibr, ir2, (rEnd - ir2) * sizeof(double));
		memcpy(ibn, in1, (rMiddle - ir1) * sizeof(int));
		memcpy(ibn, in2, (rEnd - ir2) * sizeof(int));
		//copy from buffers
		memcpy(rBegin, bufferR, (rEnd - rBegin) * sizeof(double));
		memcpy(nBegin, bufferN, (rEnd - rBegin) * sizeof(int));
	}
}


// Перенумеровать узлы вдоль направления
void Mesh::renumByDirection(vec2 direction) {
	std::cout << "\nMesh renumeration by direction (" << direction.x << ", " << direction.y << ")\n";
	double t = -omp_get_wtime();

	/*int keyNode = 0;
	vec2 r = node[keyNode];
	for (int i = 1; i < nodeCount; ++i) {
		if (r.x <= node[i].x && r.y <= node[i].y) {
			keyNode = i;
			r = node[i];
		}
	}*/
	int* newNodes = new int[nodeCount];
	int* oldNodes = new int[nodeCount];
	int* bufferN = new int[nodeCount];
	double* ranges = new double[nodeCount];
	double* bufferR = new double[nodeCount];

	for (int i = 0; i < nodeCount; ++i) {
		oldNodes[i] = i;
		//ranges[i] = (node[i] - node[keyNode]).norm();
		ranges[i] = node[i] * direction;
	}

	renumSort(ranges, ranges + nodeCount, oldNodes, bufferR, bufferN);

	for (int i = 0; i < nodeCount; ++i)
		newNodes[oldNodes[i]] = i;

	for (int i = 0; i < 3 * elemCount; ++i)
		elem3[i] = newNodes[elem3[i]];
	vec2* nodes = new vec2[nodeCount];
	for (size_t i = 0; i < nodeCount; ++i)
		nodes[newNodes[i]] = node[i];
	memcpy(node, nodes, nodeCount * sizeof(vec2));

	delete[] newNodes;
	delete[] oldNodes;
	delete[] bufferN;
	delete[] ranges;
	delete[] bufferR;
	delete[] nodes;
	t += omp_get_wtime();
	std::cout << "Renumerated [" << t << " sec]\n";
}


// Загрузить сетку из файла
// Поддерживаются только файлы .inp для Abaqus
bool Mesh::loadFromFile(const std::string& fileName) {
	if (fileName.substr(fileName.size() - 4) != ".inp") {
		std::cerr << "Unknown mesh file format\n";
		return false;
	}
	std::ifstream file(fileName);
	if (!file.is_open()) {
		std::cerr << "There is no such file\n";
		return false;
	}

	elemCount = nodeCount = 0;
	std::string line;
	while (std::getline(file, line)) {
		if (line == "*NODE")
			while (std::getline(file, line)) {
				if (line.size() != 0 && line[0] != '*')
					++nodeCount;
				else break;
			}
		else if (line.substr(0, 8) == "*ELEMENT") {
			std::string type;
			for (int i = 0; i < line.length(); ++i) {
				if (line.substr(i, 5) == "TYPE=") {
					i += 5;
					//std::clog << line.substr(i, 3) << "\n";
					if (line.substr(i, 3) == "S3R") {
						while (std::getline(file, line)) {
							if (line.size() != 0 && line[0] != '*')
								++elemCount;
							else break;
						}
					}
				}
			}
		}
	}
	file.clear();
	file.seekg(0);
	realloc(nodeCount, elemCount);
	int curElem3 = 0, input = 0;
	while (std::getline(file, line)) {
		if (line == "*NODE")
			while (std::getline(file, line)) {
				if (line.size() != 0 && line[0] != '*') {
					deleteCommas(line);
					//std::clog << line << "\n";
					std::istringstream ss(line);
					ss >> input;
					--input;
					if (input < nodeCount)
						ss >> node[input].x >> node[input].y;
					//std::clog << input << " " << node[input].x << " " << node[input].y << "\n";
				}
				else break;
			}
		else if (line.substr(0, 8) == "*ELEMENT") {
			std::string type;
			for (int i = 0; i < line.length(); ++i) {
				if (line.substr(i, 5) == "TYPE=") {
					i += 5;
					if (line.substr(i, 3) == "S3R") {
						while (std::getline(file, line)) {
							if (line.size() != 0 && line[0] != '*') {
								deleteCommas(line);
								std::istringstream ss(line);
								ss >> input;
								int begin = 3 * curElem3;
								for (int i = 0; i < 3; ++i) {
									ss >> elem3[begin + i];
									--elem3[begin + i];
								}
								++curElem3;
							}
							else break;
						}
					}
				}
			}
		}
	}
	file.close();
	colored = false;
	std::cout << "Mesh loaded:\n";
	std::cout << nodeCount << " nodes, " << elemCount << " elements\n";
}


// Сохранить сетку в файл .vtk для визуализации в Paraview
void Mesh::saveAsVtk(const std::string& fileName) const {
	std::ofstream file(fileName, std::ios_base::out);
	file << "# vtk DataFile Version 2.0\n";
	file << "Mesh\n";
	file << "ASCII\n";
	file << "DATASET POLYDATA\n";
	file << "POINTS " << nodeCount << " float\n";
	for (int i = 0; i < nodeCount; ++i)
		file << node[i].x << " " << node[i].y << " " << 0 << "\n";
	file << "POLYGONS " << elemCount << " " << 4 * elemCount;
	for (size_t i = 0; i < elemCount; ++i) {
		file << "\n3 ";
		for (size_t j = 3 * i; j < 3 * (i + 1); ++j)
			file << elem3[j] << " ";
	}
	file << "\nCELL_DATA " << elemCount << "\n";
	file << "SCALARS ElementID int 1\n";
	file << "LOOKUP_TABLE default\n";
	for (size_t i = 0; i < elemCount; ++i) {
		file << i << "\n";
	}
	if (colored) {
		file << "FIELD FieldData " << 1;
		file << "\nColor 1 " << elemCount << " int\n";
		for (size_t i = 0; i < elemCount; ++i)
			file << (int)elemColor[i] << " ";
	}
	file << "\nPOINT_DATA " << nodeCount << \
		"\nSCALARS NodeID int 1\nLOOKUP_TABLE my_table";
	for (size_t i = 0; i < nodeCount; ++i)
		file << "\n" << i;
}


void Mesh::makeMultiColor() {
	std::cout << "Mesh coloring...\n";
	double t = -omp_get_wtime();

	delete[] elemColor;
	elemColor = new char[elemCount];
	memset(elemColor, 0, elemCount * sizeof(char));

	char colorCount = 1;
	elemColor[0] = 1;

	for (unsigned i = 1; i < elemCount; ++i) {
		bool colors[10] = {};
		for (unsigned j = 0; j < elemCount; ++j) {
			for (unsigned ni = 0; ni < 3; ++ni)
				for (unsigned nj = 0; nj < 3; ++nj)
					if (elemNodeId(i, ni) == elemNodeId(j, nj)) {
						colors[elemColor[j]] = true;
						goto doneFindIntersection;
					}
		doneFindIntersection: {}
		}
		char newColor = 1;
		for (; newColor <= colorCount; ++newColor)
			if (!colors[newColor])
				break;
		if (newColor > colorCount)
			++colorCount;
		elemColor[i] = newColor;
	}

	t += omp_get_wtime();
	std::cout << "Colored [" << t << " sec]\n\n";
	colored = true;
}