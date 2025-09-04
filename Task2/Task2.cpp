#include <iostream>
#include <string>
#include <string.h>
#include <omp.h>
#include <math.h>
#include <fstream>


//Красивый вывод IA, JA
void print_csr(const int* IA, const int* JA, const int N) {
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
void print_clear_csr(const int* IA, const int* JA, const int N) {
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
void print_full_csr(const int* IA, const int* JA, const int N, const double* A, const double* b , const double* x) {
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
	std::cout << "}\n\nx:{";
	for (int i = 0; i < N; ++i) {
		std::cout << x[i];
		if (i != N - 1) std::cout << ", ";
	}
	std::cout << "}\n";

}
//Постный вывод матрицы, хранящейся в CSR формате с правой частью для решения СЛАУ
//(вдруг потом в файл надо будет вывести и дальше как-то с данными работать)
void print_clear_full_csr(int* IA, int* JA, const int N, const double* A, const double* b) {
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

//Вывод информации о расходе памяти
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
void generate(const int Nx, const int Ny, const int k1, const int k2, int& N, int*& IA, int*& JA) {

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
void fill(const int* IA, const int* JA, const int N, double*& A, double*& b) {
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
	#pragma omp parallel for schedule(dynamic, 100)
	for (int i = 0; i < N; ++i) {
		double sum = 0;
		int index_i; // Переменная, хранящая диагональный индекс
		for (int col = IA[i]; col < IA[i + 1]; ++col) {
			if (JA[col] - i) { //i!=j
				A[col] = cos(i * JA[col] + i + JA[col]);
				sum += fabs(A[col]); //∑|aij|, i!=j
			}
			else//i=j
				index_i = col;
		}
		A[index_i] = 1.234 * sum;
		b[i] = sin(i);//можно раскрутить
	}
}


double dot(const double* a, const double* b, const int N){
	double sum = .0;
	#pragma omp for reduction(+:sum) default(shared) schedule(static)
	for (int i = 0; i < N; ++i){
		sum+=a[i]*b[i];
	}
	return sum;
}


void axpy(const double* a, const double beta, const double* b, double*& res, const int N){

	#pragma omp for schedule(static) default(none) shared(N, a, res, beta, b)
	for (int i = 0; i < N; ++i){
		res[i]=a[i]+beta*b[i];
	}

}

//Умножение матрицы на вектор с применением OpenMP
void SpmV(const int* IA, const int* JA, const double*A, const double* v, double*& res, const int N){
	

		#pragma omp for schedule(dynamic, 100)
		for (int i = 0; i < N; ++i) {
			double sum = 0.0;
			for (int col = IA[i]; col < IA[i + 1]; ++col) {
				sum+=A[col]*v[JA[col]];
			}
			res[i] = sum;
		}
}

void SpmV_M(const double*M, const double* v, double*& res, const int N){
	

		#pragma omp for schedule(static)
		for (int i = 0; i < N; ++i) {
			res[i] = M[i]*v[i];
		}
}



//Решение СЛАУ методом сопряжённых гардиентов с применением OpenMP
void solve(const int* IA, const int* JA, const double* A, const double* b, const int N, const double eps, const int maxit, const bool norms, double*& x, int& step, double& norm){
	double* r = new double[N];
	double* p = new double[N];
	double* q = new double[N];
	double* z = new double[N];
	double* M = new double[N];

	#pragma omp parallel for 
    for (int i = 0; i < N; ++i) {
        x[i] = 0;
    }
	#pragma omp parallel for 
    for (int i = 0; i < N; ++i) {
        r[i] = b[i];
    }
	
	#pragma omp parallel for schedule(dynamic, 100)
	for (int i = 0; i < N; ++i) {
		for (int col = IA[i]; col < IA[i + 1]; ++col) {
			if (!(JA[col] - i)) { //i=j
				M[i] = 1.0 / A[col];
				break;
			}
		}
	}

	step = 0;
	double po, beta, po_temp, alpha;
	#pragma omp parallel
	do{
		++step;
		SpmV_M(M,r,z,N);
		po=dot(r, z, N);
		if (step == 1) {
			#pragma omp for 
    			for (int i = 0; i < N; ++i) {
        			p[i] = z[i];
    			}
		} else {
			beta = po / po_temp;
			axpy(z, beta, p, p, N);
			}
		SpmV(IA, JA, A, p, q, N);
		alpha = po / dot(p, q, N);
		axpy(x, alpha, p, x, N);
		axpy(r, -alpha, q, r, N);
		po_temp = po;
		if (norms){
			norm = sqrt(dot(r, r, N));
			std::cout << "solve[" << step << "]: l2 norm=" << norm << std::endl;
		}
	}while (po > eps*eps && step<maxit);

	norm = sqrt(dot(r,r,N));
	delete[] z;
	delete[] p;
	delete[] r;
	delete[] q;
	delete[] M;
}

void report(double norm, double timeG, double timeF, double timeS, double dott, double axpyt, double SpmVt, int N, int nnz){

	std::cout << "Number of threads: " << omp_get_max_threads() << std::endl;
	std::cout << "L2 norm of solving Ax-b: " << norm << std::endl;
	std::cout << "Generation step time: " << timeG << std::endl;
	std::cout << "Fill matrix step time: " << timeF << std::endl;
	std::cout << "Solve matrix step time: " << timeS << std::endl;
	printf("dot[%d] time = %g GFLOPS = %g GB/s = %g \n",
		N, dott, 2.0*N/(dott*1E9), (2.0*N*sizeof(double))/(dott*1E9));

	printf("axpy[%d] time = %g GFLOPS = %g GB/s = %g \n",
		N, axpyt, 2.0*N/(axpyt*1E9), (3.0*N*sizeof(double))/(axpyt*1E9));
	
	printf("SpmV[%d] time = %g GFLOPS = %g \n",
		N, SpmVt, (2.0*nnz)/(SpmVt*1E9));
	printf("Solve GFLOPS = %g \n",
		(14.0*N+2.0*nnz)/(timeS*1E9));//2N -l2, 2*N - spmv-m, 2nnz - spmv-a, 3*2N - axpy, 2*2N - dot
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
				  <<  argv[0]<<" Nx Ny k1 k2\n"
				  << "where Nx, Ny, k1, k2 are integers\n\n"
				  << "Nx is size of mesh in horizontal\n"
				  << "Ny is size of mesh in vertical\n"
				  << "k1 is number of full cells in a row\n"
				  << "k2 is number of split cells in a row\n"
				  << "\nIf you wonna out CSR data, please add \"out\" like:\n"
				  <<  argv[0]<<" Nx Ny k1 k2 out\n"
				  << "\nIf you wonna out l2 norms by step while solving, please add \"l2\" like:\n"
				  <<  argv[0]<<" Nx Ny k1 k2 l2\n"
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
	bool is_out = false, norms = false;
	if (argc == 6){
		if (!strcmp(argv[5], "out"))
			is_out=true;
		else if (!strcmp(argv[5], "l2"))
			norms=true;
	}
	if (argc > 6){
		if (!strcmp(argv[5], "out") || !strcmp(argv[6], "out"))
			is_out=true;
		if (!strcmp(argv[5], "l2") || !strcmp(argv[6], "l2"))
			norms=true;
	}
	if(Nx < 1 || Ny < 1 || k1==0 && k2==0 || k1 < 0 || k2 < 0){
		std::cout << "Incorrect input data. Please, use me correctly\n";
		return 0;
	}
	int step;
	int* JA = NULL, *IA = NULL;
	double* A = NULL, *b = NULL, *x = NULL;
	double norm;
	int N;
	double time, timeG, timeF, timeS;
	
	time = omp_get_wtime();
	generate(Nx, Ny, k1, k2, N, IA, JA);
	timeG = omp_get_wtime() - time;	
	
	time = omp_get_wtime();
	fill(IA, JA, N, A, b);
	timeF = omp_get_wtime() - time;
	x = new double[N];
	time = omp_get_wtime();
	for (int i = 10; i--;)
		solve(IA, JA, A, b, N, 1e-7, 100000, norms, x, step, norm);
	timeS = (omp_get_wtime() - time) / 10 / step;
	
	
	double* a, * c;

	a = new double[N];
	c = new double[N];
	#pragma omp parallel for
	for(int i = 0; i<N; ++i){
		a[i] = i%456;
		c[i] = (i*i)%134;
	}

	time = omp_get_wtime();
#pragma omp parallel
	for (int i = 100; i--;)
		dot(a, c, N);
	double dott = (omp_get_wtime() - time) / 100;
	
	time = omp_get_wtime();
#pragma omp parallel
	for (int i = 100; i--;)
		axpy(a, 1.4, c, a, N);
	double axpyt = (omp_get_wtime() - time) / 100;
	
	time = omp_get_wtime();
#pragma omp parallel
	for (int i = 100; i--;)
		SpmV(IA, JA, A, c, c, N);
	double SpmVt = (omp_get_wtime() - time) / 100;
	
	
	delete[] a;
	delete[] c;

	report(norm, timeG, timeF, timeS, dott, axpyt, SpmVt, N, IA[N]);

	if(is_out)
		print_full_csr(IA, JA, N, A, b, x);

	delete[] JA;
	delete[] IA;
	delete[] A;
	delete[] b;
	delete[] x;
	
	return 0;
}
