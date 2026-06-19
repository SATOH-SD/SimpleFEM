#include "Mesh.h"

#include <iomanip>
#include <list>

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

	for (unsigned i = 0; i < nodeCount; ++i)
		newNodes[oldNodes[i]] = i;

	for (unsigned i = 0; i < 3 * elemCount; ++i)
		elem3[i] = newNodes[elem3[i]];
	for (unsigned i = 0; i < borderOffset[borderCount]; ++i)
		borders[i] = newNodes[borders[i]];
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

	delete[] borderOffset;
	borderOffset = new unsigned[1];
	borderOffset[0] = 0;

	unsigned maxIndex = 0, input = 0;
	elemCount = nodeCount = 0;
	std::string line;
	while (std::getline(file, line)) {
		if (line == "*NODE")
			while (std::getline(file, line)) {
				if (line.size() != 0 && line[0] != '*') {
					++nodeCount;
					std::istringstream ss(line);
					ss >> input;
					if (input > maxIndex)
						maxIndex = input;
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
							if (line.size() != 0 && line[0] != '*')
								++elemCount;
							else break;
						}
					}
				}
			}
		}
		else if (line == "*MPC")
			for (;;) {
				++borderCount;
				unsigned* tempBorderOffset = new unsigned[borderCount + 1];
				memcpy(tempBorderOffset, borderOffset, borderCount * sizeof(int));
				std::swap(tempBorderOffset, borderOffset);
				delete[] tempBorderOffset;
				borderOffset[borderCount] = borderCount ? borderOffset[borderCount - 1] : 0;
				while (std::getline(file, line)) {
					if (line.size() != 0 && line[0] != '*')
						++borderOffset[borderCount];
					else break;
				}
				if (line != "*MPC") break;
			}
	}

	//std::cout << "maxIndex " << maxIndex << "\n";

	file.clear();
	file.seekg(0);
	realloc(nodeCount, elemCount);
	borders = new unsigned[borderOffset[borderCount]];
	unsigned curNode = 0, curElem3 = 0, curBorder = 0;
	unsigned* indexes = new unsigned[maxIndex + 1];
	std::list<unsigned> superNodes;

	while (std::getline(file, line)) {
		if (line == "*NODE")
			while (std::getline(file, line)) {
				if (line.size() != 0 && line[0] != '*') {
					deleteCommas(line);
					std::istringstream ss(line);
					ss >> input;
					indexes[input] = curNode;
					ss >> node[curNode].x >> node[curNode].y;
					++curNode;
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
								}
								++curElem3;
							}
							else break;
						}
					}
				}
			}
		}
		else if (line == "*MPC")
			for (;;) {
				unsigned curPoint = 0;
				unsigned curOffset = borderOffset[curBorder];
				while (std::getline(file, line)) {
					if (line.size() != 0 && line[0] != '*') {
						deleteCommas(line);
						line[3] = line[2] = line[1] = line[0] = ' ';
						std::istringstream ss(line);
						ss >> borders[curOffset + curPoint];
						if (curPoint == 0) {
							int super;
							ss >> super;
							superNodes.push_back(super);
						}
						borders[curOffset + curPoint];
						++curPoint;
					}
					else break;
				}
				++curBorder;
				if (line != "*MPC") break;
			}
	}

	for (unsigned i = 0; i < 3 * elemCount; ++i)
		elem3[i] = indexes[elem3[i]];
	for (unsigned i = 0; i < borderOffset[borderCount]; ++i)
		borders[i] = indexes[borders[i]];
	for (unsigned& super : superNodes)
		super = indexes[super];
	for (unsigned super : superNodes) {
		for (unsigned i = super + 1; i < nodeCount; ++i)
			node[i - 1] = node[i];
		--nodeCount;
		for (unsigned i = 0; i < 3 * elemCount; ++i)
			if (elem3[i] > super)
				--elem3[i];
		for (unsigned i = 0; i < borderOffset[borderCount]; ++i)
			if (borders[i] > super)
				--borders[i];
		for (unsigned& rest : superNodes)
			if (rest > super)
				--rest;
	}

	file.close();
	colored = false;
	std::cout << "Mesh loaded:\n";
	std::cout << nodeCount << " nodes, " << elemCount << " elements\n";
}


// Сохранить сетку в файл .vtk для визуализации в Paraview
void Mesh::saveAsVtk(const std::string& fileName) const {
	std::cout << "\nSaving...\n";
	double t = -omp_get_wtime();
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
	file.close();
	t += omp_get_wtime();
	std::cout << "Mesh saved in file " << fileName << " [" << t << " sec]\n";
}


// Сохранить сетку в файл .vtu для визуализации в Paraview
void Mesh::saveAsVtu(const std::string& fileName) const {
	std::cout << "\nSaving...\n";
	double t = -omp_get_wtime();
	std::ofstream file(fileName, std::ios_base::out);
	file << "<?xml version=\"1.0\"?>" \
		<< "\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\" >" \
		<< "\n\t<UnstructuredGrid>" \
		<< "\n\t\t<Piece NumberOfPoints=\"" << nodeCount << "\" NumberOfCells=\"" << elemCount << "\">";
	
	file << "\n\t\t\t<PointData>";
	file << "\n\t\t\t\t<DataArray type=\"UInt32\" Name=\"NodeID\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
	for (unsigned i = 0; i < nodeCount; ++i)
		file << i << " ";
	file << "\n\t\t\t\t</DataArray>";
	file << "\n\t\t\t</PointData>";
	
	file << "\n\t\t\t<CellData>";
	file << "\n\t\t\t\t<DataArray type=\"UInt32\" Name=\"ElementID\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
	for (unsigned i = 0; i < elemCount; ++i)
		file << i << " ";
	file << "\n\t\t\t\t</DataArray>";
	file << "\n\t\t\t</CellData>";

	file << "\n\t\t\t<Points>";
	file << "\n\t\t\t\t<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" >\n";
	for (unsigned i = 0; i < nodeCount; ++i)
		file << node[i].x << " " << node[i].y << " 0 ";
	file << "\n\t\t\t\t</DataArray>";
	file << "\n\t\t\t</Points>";

	file << "\n\t\t\t<Cells>";
	file << "\n\t\t\t\t<DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\" >\n";
	for (int i = 0; i < elemCount * 3; ++i)
		file << elem3[i] << " ";
	file << "\n\t\t\t\t</DataArray>";
	file << "\n\t\t\t\t<DataArray type=\"UInt32\" Name=\"offsets\" format=\"ascii\" >\n";
	for (int i = 3; i <= elemCount * 3; i += 3)
		file << i << " ";
	file << "\n\t\t\t\t</DataArray>";
	file << "\n\t\t\t\t<DataArray type=\"UInt32\" Name=\"types\" format=\"ascii\" >\n";
	for (int i = 0; i < elemCount; ++i)
		file << "5 ";
	file << "\n\t\t\t\t</DataArray>";
	file << "\n\t\t\t</Cells>";

	file << "\n\t\t</Piece>";
	file << "\n\t</UnstructuredGrid>";
	file << "\n</VTKFile>";
	file.close();
	t += omp_get_wtime();
	std::cout << "Mesh saved in file " << fileName << " [" << t << " sec]\n";
}


// Сохранить сетку в файл .vtu для визуализации в Paraview
//void Mesh::saveAsVtu_(const std::string& fileName) const {
//	std::ofstream file(fileName, std::ios_base::out);
//	size_t offset = 0;
//	file << "<?xml version=\"1.0\"?>" \
//		<< "\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\" >" \
//		<< "\n\t<UnstructuredGrid>" \
//		<< "\n\t\t<Piece NumberOfPoints=\"" << nodeCount << "\" NumberOfCells=\"" << elemCount << "\">";
//
//	file << "\n\t\t\t<PointData>";
//	file << "\n\t\t\t\t<DataArray type=\"UInt32\" Name=\"NodeID\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" \
//		<< offset << "\" />\n";
//	offset += nodeCount;
//	file << "\n\t\t\t</PointData>";
//
//	file << "\n\t\t\t<CellData>";
//	file << "\n\t\t\t\t<DataArray type=\"UInt32\" Name=\"ElementID\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" \
//		<< offset << "\" />\n";
//	offset += elemCount;
//	file << "\n\t\t\t</CellData>";
//
//	file << "\n\t\t\t<Points>";
//	file << "\n\t\t\t\t<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" \
//		<< offset << "\" />\n";
//	offset += nodeCount * 3;
//	file << "\n\t\t\t</Points>";
//
//	file << "\n\t\t\t<Cells>";
//	file << "\n\t\t\t\t<DataArray type=\"UInt32\" Name=\"connectivity\" format=\"appended\" offset=\"" \
//		<< offset << "\" />\n";
//	offset += elemCount * 3;
//	file << "\n\t\t\t\t<DataArray type=\"UInt32\" Name=\"offsets\" format=\"appended\" offset=\"" \
//		<< offset << "\" />\n";
//	offset += elemCount;
//	file << "\n\t\t\t\t<DataArray type=\"UInt32\" Name=\"types\" format=\"appended\" offset=\"" \
//		<< offset << "\" />\n";
//	offset += elemCount;
//	file << "\n\t\t\t</Cells>";
//
//	file << "\n\t\t</Piece>";
//	file << "\n\t</UnstructuredGrid>";
//
//	file << "\n\t<AppendedData encoding=\"raw\">\n_";
//	double zero = 0;
//	for (unsigned i = 0; i < nodeCount; ++i)
//		file << reinterpret_cast<const char*>(&i);
//	for (unsigned i = 0; i < elemCount; ++i)
//		file << reinterpret_cast<const char*>(&i);
//	for (unsigned i = 0; i < nodeCount; ++i)
//		file << reinterpret_cast<const char*>(&node[i].x) << reinterpret_cast<const char*>(&node[i].y) << reinterpret_cast<const char*>(&zero);
//	for (int i = 0; i < elemCount * 3; ++i)
//		file << reinterpret_cast<const char*>(&elem3[i]);
//	for (int i = 3; i <= elemCount * 3; i += 3)
//		file << reinterpret_cast<const char*>(&i);
//	unsigned type = 5;
//	for (int i = 0; i < elemCount; ++i)
//		file << reinterpret_cast<const char*>(&type);
//
//	file << "\n\t</AppendedData>";
//	file << "\n</VTKFile>";
//	file.close();
//}


void Mesh::makeMultiColorSlow() {
	std::cout << "Mesh coloring... (slow)\n";
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


void Mesh::makeMultiColor() {
	std::cout << "Mesh coloring...\n";
	double t = -omp_get_wtime();

	delete[] elemColor;
	elemColor = new char[elemCount];
	memset(elemColor, 0, elemCount * sizeof(char));

	char colorCount = 1;
	elemColor[0] = 1;

	char* nodeInCount = new char[nodeCount];
	memset(nodeInCount, 0, nodeCount * sizeof(char));
	
	for (unsigned i = 0; i < elemCount; ++i)
		for (unsigned ni = 0; ni < 3; ++ni)
			++nodeInCount[elemNodeId(i, ni)];

	unsigned dataWidth = 0;
	for (unsigned i = 0; i < nodeCount; ++i) {
		unsigned count = nodeInCount[i];
		if (count > dataWidth)
			dataWidth = count;
	}
	std::cout << "dataWidth = " << dataWidth << "\n";

	unsigned* nodeInElems = new unsigned[nodeCount * dataWidth];
	memset(nodeInCount, 0, nodeCount * sizeof(char));

	for (unsigned i = 0; i < elemCount; ++i) {
		bool colors[10] = {};
		for (unsigned ni = 0; ni < 3; ++ni) {
			unsigned nodeId = elemNodeId(i, ni);
			for (unsigned j = 0; j < nodeInCount[nodeId]; ++j)
				colors[elemColor[nodeInElems[nodeId * dataWidth + j]]] = true;
			nodeInElems[nodeId * dataWidth + nodeInCount[nodeId]] = i;
			++nodeInCount[nodeId];
		}
		char newColor = 1;
		for (; newColor <= colorCount; ++newColor)
			if (!colors[newColor])
				break;
		if (newColor > colorCount)
			++colorCount;
		elemColor[i] = newColor;
	}

	std::cout << "colorCount = " << (int)colorCount << "\n";

	/*for (unsigned i = 0; i < elemCount; ++i)
		std::cout << i << " " << (int)elemColor[i] << "\n";*/

	delete[] nodeInCount;
	delete[] nodeInElems;

	unsigned* colorCounts = new unsigned[colorCount + 1];
	memset(colorCounts, 0, (colorCount + 1) * sizeof(unsigned));
	for (unsigned i = 0; i < elemCount; ++i)
		++colorCounts[elemColor[i]];
	for (unsigned i = 1; i <= colorCount; ++i)
		std::cout << i << ": " << colorCounts[i] << "\n";
	delete[] colorCounts;

	t += omp_get_wtime();
	std::cout << "Colored [" << t << " sec]\n\n";
	colored = true;
}