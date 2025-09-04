#include <iostream>
#include <string>
#include <string.h>
#include <omp.h>
#include <math.h>
#include <fstream>


//Красивый вывод IA, JA
void print_csr(int* IA, int* JA, int N) {
	std::cout << "\n\nCSR FORMAT MATRIX \nN: " << N << "\nJA Matrix: \n{";
	
	for (int i = 0; i < N; ++i) {
		for (int col = IA[i]; col < IA[i + 1]; ++col) {
			
			if (col == IA[i + 1] - 1)
				std::cout << JA[col];
			else
				std::cout << JA[col] << ", ";

		}
		if ( i != N - 1) std::cout << " | ";
	}
	std::cout << "}\nIA Matrix:\n{";
	for (int i = 0; i < N + 1; ++i) {
		if (i == N)
			std::cout << IA[i];
		else
			std::cout << IA[i] << ", ";
	}
	std::cout << "}\n";
}
//Постный вывод IA, JA
//(вдруг потом в файл надо будет вывести и дальше как-то с данными работать)
void print_clear_csr(int* IA, int* JA, int N) {
	std::cout << "\n\nCSR FORMAT MATRIX \nN: " << N << "\nJA Matrix:\n";

	for (int i = 0; i < N; ++i)
		for (int col = IA[i]; col < IA[i + 1]; ++col) 
			std::cout << JA[col] << " ";

	std::cout << "\nIA Matrix:\n";
	for (int i = 0; i < N + 1; ++i)
		std::cout << IA[i] << " ";
	std::cout << "\n";
}
//Красивый вывод матрицы, хранящейся в CSR формате с правой частью для решения СЛАУ
void print_full_csr(int* IA, int* JA, int N, double* A, double* b) {
	std::cout << "\n\nCSR FORMAT MATRIX \nN: " << N << "\n\nA Matrix: {";
	
	for (int i = 0; i < N; ++i) {
		for (int col = IA[i]; col < IA[i + 1]; ++col) {

			if (col == IA[i + 1] - 1)
				std::cout << A[col];
			else
				std::cout << A[col] << ", ";

		}
		if (i != N - 1) std::cout << " | ";
	}
	std::cout << "}\n\nJA Matrix: {";

	for (int i = 0; i < N; ++i) {
		for (int col = IA[i]; col < IA[i + 1]; ++col) {

			if (col == IA[i + 1] - 1)
				std::cout << JA[col];
			else
				std::cout << JA[col] << ", ";

		}
		if (i != N - 1) std::cout << " | ";
	}
	std::cout << "}\n\nIA Matrix:\n{";
	for (int i = 0; i < N + 1; ++i) {
		if (i == N)
			std::cout << IA[i];
		else
			std::cout << IA[i] << ", ";
	}
	std::cout << "}\n\nb:\n{";
	for (int i = 0; i < N; ++i) {
		std::cout << b[i];
		if (i != N - 1) std::cout << ", ";
	}
	std::cout << "}\n\n";
}
//Постный вывод матрицы, хранящейся в CSR формате с правой частью для решения СЛАУ
//(вдруг потом в файл надо будет вывести и дальше как-то с данными работать)
void print_clear_full_csr(int* IA, int* JA, int N, double* A, double* b) {
	std::cout << "\n\nCSR FORMAT MATRIX \nN: " << N << "\n\nA Matrix: ";
	for (int i = 0; i < N; ++i)
		for (int col = IA[i]; col < IA[i + 1]; ++col)
			std::cout << A[col] << " ";
	std::cout << "\n\nJA Matrix: ";
	for (int i = 0; i < N; ++i)
		for (int col = IA[i]; col < IA[i + 1]; ++col)
			std::cout << JA[col] << " ";
	std::cout << "\n\nIA Matrix:\n";
	for (int i = 0; i < N + 1; ++i)
		std::cout << IA[i] << " ";
	std::cout << "\n\nb:\n";
	for (int i = 0; i < N; ++i)
		std::cout << b[i] << " ";
	std::cout << "\n\n";
}

//Вывод информации о расходе памяти
int memoryUsage() {
    std::ifstream statusFile("/proc/self/status");
    std::string line;
    while (std::getline(statusFile, line)) {
        if (line.find("VmRSS") != std::string::npos) {
			return stoi(line.substr(6));
        }
    }
	return -1;
}

//Вывод информации о существовании памяти
int memoryExist() {
    std::ifstream statusFile("/proc/meminfo");
    std::string line;
    while (std::getline(statusFile, line)) {
        if (line.find("MemTotal") != std::string::npos) {
			return stoi(line.substr(9));
        }
    }
	return -1;
}
//Генерация JA, IA массивов по заданным параметрам сетки
void generate(int Nx, int Ny, int k1, int k2, int& N, int*& IA, int*& JA) {

	int SIZE = (Nx + 1) * (Ny + 1) + //самовхождения (число ячеек_x+1)*(число ячеек_y+1) +...
		2 * (Nx * (Ny + 1) + //горизонтальные рёбра (число ячеек_x)*(число ячеек_y+1) +...
			(Nx + 1) * Ny + //вертикальные рёбра (число ячеек_x+1)*(число ячеек_y) +...
			(Nx * Ny) / (k1 + k2) * k2 + //количество целых (полных по k2) косых рёбер
			((Nx * Ny) % (k1 + k2) > k1 ? (Nx * Ny) % (k1 + k2) - k1 : 0)); //количество последних косых рёбер 

	JA = new int[SIZE];//номера столбцов ненулевых элементов
	
	if(memoryUsage() > (memoryExist() / 2)){
		delete[] JA;
		std::cout << "Too much memory allocated during malloc JA\n";
		exit(1);
	}
	N = (Ny + 1) * (Nx + 1);
	IA = new int[N + 1]; //информация о позиции начала списка cтолбцов данной строки
	if(memoryUsage() > (memoryExist() / 2)){
		delete[] JA;
		delete[] IA;
		std::cout << "Too much memory allocated during malloc IA\n";
		exit(1);
	}
	IA[N] = SIZE;

	int I_node, I;
	int indecies = 0, kk = k1 + k2;
	for (int i = 0; i <= Ny; ++i) {
		for (int j = 0; j <= Nx; ++j) {
			I_node = i * (Nx+1) + j; //Индекс узла в сетке
			I = i * Nx + j; //Абсолютный индекс ячейки в сетке
			IA[I_node] = indecies;
			if (i && j && (I - Nx - 1) % kk >= k1) 
				JA[indecies++] = I_node - Nx - 2;// ↖ влево-вверх
			if (i)			//Верхний сосед есть у всех узлов, кроме первого слоя сетки по i
				JA[indecies++] = I_node - Nx - 1; // ↑ вверх
			if (j)			//Левый сосед есть у всех узлов, кроме первого слоя сетки по j
				JA[indecies++] = I_node - 1; // ← влево
			JA[indecies++] = I_node; // самовхождение
			if (Nx - j)		//Правый сосед есть у всех узлов, кроме последнего слоя сетки по j
				JA[indecies++] = I_node + 1; // → вправо
			if (Ny - i)		//Нижний сосед есть у всех узлов, кроме последнего слоя сетки по i
				JA[indecies++] = I_node + Nx + 1; // ↓ вниз
			if (Nx - j && Ny - i && (I % kk >= k1))
				JA[indecies++] = I_node + Nx + 2; // ↘ вправо-вниз
		}
	}
}

//Заполнение матрицы A и правой части b по заданным формулам в CSR формате
void fill(int* IA, int* JA, int N, double*& A, double*& b) {
	A = new double[IA[N]];

	if(memoryUsage() > (memoryExist() / 2)){
		delete[] JA;
		delete[] IA;
		delete[] A;
		std::cout << "Too much memory allocated during malloc A\n";
		exit(1);
	}
	
	b = new double[N];
	if(memoryUsage() > (memoryExist() / 2)){
		delete[] JA;
		delete[] IA;
		delete[] A;
		delete[] b;
		std::cout << "Too much memory allocated during malloc b\n";
		exit(1);
	}
	double sum;
	int indecies = 0; // Абсолютный индекс для прохождения по массиву A
	int index_i; // Переменная, хранящая диагональный индекс
	for (int i = 0; i < N; ++i) {
		sum = 0;
		for (int col = IA[i]; col < IA[i + 1]; ++col) {
			if (JA[col] - i) { //i!=j
				A[indecies] = cos(i * JA[col] + i + JA[col]);
				sum += fabs(A[indecies++]); //∑|aij|, i!=j
			}
			else//i=j
				index_i = indecies++;
		}
		A[index_i] = 1.234 * sum;
		b[i] = sin(i);//можно раскрутить
	}
}


int main(int argc, char** argv) {
	
	int Nx, Ny, k1, k2;
	char* p_end{};
	if(argc == 1){
		std::cout << "_____________________________________________\n"
				  << "|        ULTIMATE CSR MAKER VIA MESH        |\n"
				  << "¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
				  << "\nHow to launch me correctly?\n"
		 		  << "\nType command like:\n" 
				  << "./program Nx Ny k1 k2\n"
				  << "where Nx, Ny, k1, k2 are integers\n\n"
				  << "Nx is size of mesh in horizontal\n"
				  << "Ny is size of mesh in vertical\n"
				  << "k1 is number of full cells in a row\n"
				  << "k2 is number of split cells in a row\n"
				  << "\nIf you wonna out CSR data, please add \"out\" like:\n"
				  << "./program Nx Ny k1 k2 out\n"
				  << "\nPlease, launch me again with correct data :)\n";
		return 0;
	}
	if(argc < 5){
		std::cout << "Incorrect input data. Please, use me correctly\n"; 
		return 0;
	}
	try{
		Nx = std::stoi(argv[1]);
		Ny = std::stoi(argv[2]);
		k1 = std::stoi(argv[3]);
		k2 = std::stoi(argv[4]);
	} catch(std::invalid_argument& ex)
	{
		std::cout << "std::invalid_argument::what(): " << ex.what() << '\n';
		return 1;
	}
	catch (std::out_of_range const& ex)
	{
		std::cout << "std::out_of_range::what(): " << ex.what() << '\n';
		return 1;
	}
	bool is_out = false;
	if (argc >= 6)
		if (!strcmp(argv[5], "out"))
			is_out=true;
	
	
	if(Nx < 1 || Ny < 1 || k1==0 && k2==0 || k1 < 0 || k2 < 0){
		std::cout << "Incorrect input data. Please, use me correctly\n";
		return 0;
	}
	
	int* JA = NULL, *IA = NULL;
	double* A = NULL, *b = NULL;
	int N;
	double time;
	std::cout << "Memory usage before generating: " << memoryUsage() << " kB\n";
	time = omp_get_wtime();
	generate(Nx, Ny, k1, k2, N, IA, JA);
	time = omp_get_wtime() - time;
	//print_csr(IA, JA, N);
	std::cout << "Generation step time: " << time << std::endl;
	std::cout << "Memory usage after generating: " << memoryUsage() << " kB\n";
	time = omp_get_wtime();
	fill(IA, JA, N, A, b);
	time = omp_get_wtime() - time;
	
	std::cout << "Fill matrix step time: " << time << std::endl;
	std::cout << "Memory usage after filling matrix A: " << memoryUsage() << " kB\n";
	if(is_out)
		print_full_csr(IA, JA, N, A, b);

	delete[] JA;
	delete[] IA;
	delete[] A;
	delete[] b;
	std::cout << "Memory usage after freeing memory: " << memoryUsage() << " kB\n";
	return 0;
}