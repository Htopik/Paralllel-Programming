#include <iostream>
#include <string>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <mpi.h>

int ranksize = 1;
int rank = 0;
int mpi_initialized = 0; // MPI initialization flag
MPI_Comm MCW = MPI_COMM_WORLD; // default communicator


#define crash(...) exit(Crash( __VA_ARGS__)) // via exit macro, so a static analyzer knows that it is an exit point
static int Crash(const char* fmt, ...) { // termination of the program in case of error
	va_list ap;
	if (mpi_initialized) fprintf(stderr, "\nEpic fail: MyID = %d\n", rank);
	else fprintf(stderr, "\nEpic fail: \n");

	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);

	fprintf(stderr, "\n");
	fflush(stderr);

	if (mpi_initialized) MPI_Abort(MCW, -1);
	return -1;
}
// wrapper for MPI barrier
static inline void barrier() {
	if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) crash("Base lib barrier: MPI_Barrier failed! \n");
}

void print_csr(const int* IA, const int* JA, const int N) {
	std::cout << "\n\nCSR FORMAT MATRIX \nN: " << N << "\nJA Matrix: \n{";

	for (int i = 0; i < N; ++i) {
		for (int col = IA[i]; col < IA[i + 1]; ++col) {

			if (col == IA[i + 1] - 1)
				std::cout << JA[col];
			else
				std::cout << JA[col] << ", ";

		}
		if (i != N - 1) std::cout << " | ";
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

void print_full_csr(const int* IA, const int* JA, const int N, const int sizeN, const double* A, const double* b, const double* x) {
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
//Получение глобального индекса по локальному внутреннему
int Global(const int index, const int Nx_loc, const int Nx, const int i_0, const int j_0) {
	return Nx * (i_0 + index / Nx_loc) + j_0 + index % Nx_loc;
}
//Массив взаимосвязи всех внутренних-интерфейсных-гало элементов с их глобальными индексами
void Local2Global(const int Nx, const int Ny, const int k1, const int k2, const int px, const int py,
	const int Nx_loc, const int Ny_loc, int& N, int& No, int*& L2G, int*& Part) {

	int up = py > 1 && rank / px ? Nx_loc : 0;
	int left = px > 1 && (rank % px) ? Ny_loc : 0;
	int right = px > 1 && ((rank + 1) % px) ? Ny_loc : 0;
	int down = py > 1 && rank / px != py - 1 ? Nx_loc : 0;
	bool is_l_up = up && left ?
		((Nx-1) * ((rank / px) * (Ny / py) + (Ny % py ? std::min(rank / px, Ny % py) : 0) - 1) + //номер строки ячейки
			(Nx % px ? std::min(rank % px, Nx % px) : 0) + (rank % px) * (Nx / px) - 1) //номер столбца ячейки
		% (k1 + k2) >= k1 : 0;

	bool is_r_down = down && right ?
		((Nx-1) * (Ny_loc - 1 + (rank / px) * (Ny / py) + (Ny % py ? std::min(rank / px, Ny % py) : 0)) +//номер строки ячейки
			Nx_loc - 1 + ((Nx % px) ? std::min(rank % px, Nx % px) : 0) + (rank % px) * (Nx / px))//номер столбца ячейки
		% (k1 + k2) >= k1 : 0;

	int HALO = up + down + left + right + is_r_down + is_l_up;
	No = Nx_loc * Ny_loc;
	N = No + HALO;

	L2G = new int[N];
	Part = new int[N];
	//Индекс строки глобальной позиции локального первого узла
	int i_0 = (rank / px) * (Ny / py) + ((Ny % py) ? std::min(rank / px, Ny % py) : 0);
	//Индекс столбца глобальной позиции локального первого узла
	int j_0 = ((Nx % px) ? std::min((rank % px), (Nx % px)) : 0) + (rank % px) * (Nx / px);
	for (int i = 0; i < Ny_loc + (down > 0); ++i) { //inner+interface+down
		L2G[i * Nx_loc] = Global(i * Nx_loc, Nx_loc, Nx, i_0, j_0);

		Part[i * Nx_loc] = i < Ny_loc ? rank : rank + px;

		for (int j = 1; j < Nx_loc ; ++j) {
			L2G[i * Nx_loc + j] = L2G[i * Nx_loc + j - 1] + 1;
			Part[i * Nx_loc + j] = i < Ny_loc ? rank : rank + px;
		}
	}


	int index = No + down;

	for (int i = 0; i < left; ++i) { //цикл не сработает, если left===0
		L2G[index] = L2G[i * Nx_loc] - 1;
		Part[index++] = rank - 1;
	}
	for (int i = 0; i < up; ++i) { //цикл не сработает, если up===0
		L2G[index] = L2G[i] - Nx;
		Part[index++] = rank - px;
	}
	for (int i = 0; i < right; ++i) { //цикл не сработает, если right===0
		L2G[index] = L2G[(i + 1) * Nx_loc - 1] + 1;
		Part[index++] = rank + 1;
	}
	if (is_l_up) {
		L2G[index] = L2G[0] - Nx - 1;
		Part[index++] = rank - px - 1;
	}
	if (is_r_down) {
		L2G[index] = L2G[No + down - 1] + 1;
		Part[index++] = rank + px + 1;
	}
}
//Обратный для Local2Global массив, чтобы узнать локальные позиции по глобальным индексам
void Global2Local(const int* L2G, const int N, int& Np, const int Nx, const int Ny, int*& G2L) {
	Np = (Nx + 1) * (Ny + 1);
	G2L = new int[Np];
	for (int i = 0; i < Np; ++i) {
		G2L[i] = -1;
	}
	for (int i = 0; i < N; ++i) {
		G2L[L2G[i]] = i;
	}
}

//Генерация локальных JA, IA массивов по заданным параметрам сетки и рангу
void generate(const int Nx, const int Ny, const int k1, const int k2, const int px, const int py,
	const int Nx_loc, const int Ny_loc, int& N, int& Np, int& No, int*& L2G, int*& Part, int*& G2L, int*& IA, int*& JA) {

	Local2Global(Nx, Ny, k1, k2, px, py, Nx_loc, Ny_loc, N, No, L2G, Part);

	Global2Local(L2G, N, Np, Nx, Ny, G2L);

	int kk = k1 + k2;

	int SIZE = Nx_loc * Ny_loc + //самовхождения (число ячеек_x+1)*(число ячеек_y+1) +...
		2 * (Nx_loc * (Ny_loc - 1) + //вертикальные рёбра (число ячеек_x+1)*(число ячеек_y) +...
			(Nx_loc - 1) * Ny_loc) + //горизонтальные рёбра (число ячеек_x)*(число ячеек_y+1) +...
		+N - No; //гало вертикальные+горизонтальные+граничные-диагональные рёбра.
	//осталось учесть диагональные + гало-диагональные рёбра
		
	int I; // глобальный индекс ячейки
	for (int i = 0; i < No - 1; ++i) {

		I = L2G[i] / Nx * (Nx-1) + (L2G[i] % Nx);
		if (px > 1) {
			if ((rank % px) && !(i % Nx_loc) && i < (Ny_loc-1) * Nx_loc && (I - 1) % kk >= k1) //левый
				SIZE++;
			if (((rank + 1) % px) && (i % Nx_loc == Nx_loc - 1) && I % kk >= k1) //правый
				SIZE++;
		}
		if (py > 1) {
			if (rank / px && i < Nx_loc - 1 && (I - Nx + 1) % kk >= k1) //верхний
				SIZE++;
			if (rank / px != py - 1 && i >= (Ny_loc-1) * Nx_loc && I % kk >= k1) //нижний
				SIZE++;
		}
		if ((i % Nx_loc) < Nx_loc - 1 && i < (Ny_loc-1) * Nx_loc && I % kk >= k1) //внутренние диагонали
			SIZE += 2;
	}	

	JA = new int[SIZE];//номера столбцов ненулевых элементов

	IA = new int[No + 1]; //информация о позиции начала списка cтолбцов данной строки
	IA[No] = SIZE;
	
	int I_node;
	int indecies = 0;
	int i_g, j_g;

	for (int i = 0; i < Ny_loc; ++i) {
		for (int j = 0; j < Nx_loc; ++j) {
			I_node = L2G[i * Nx_loc + j]; //Индекс узла в сетке
			I = I_node / Nx * (Nx-1) + I_node % Nx; //Абсолютный индекс ячейки в сетке
			i_g = I_node / Nx;
			j_g = I_node % Nx;
			IA[G2L[I_node]] = indecies;
			if (i_g && j_g && (I - Nx) % kk >= k1)
				JA[indecies++] = G2L[I_node - Nx - 1];// ↖ влево-вверх
			if (i_g)			//Верхний сосед есть у всех узлов, кроме первого слоя сетки по i
				JA[indecies++] = G2L[I_node - Nx]; // ↑ вверх
			if (j_g)			//Левый сосед есть у всех узлов, кроме первого слоя сетки по j
				JA[indecies++] = G2L[I_node - 1]; // ← влево
			JA[indecies++] = G2L[I_node]; // самовхождение
			if (Nx - j_g - 1)		//Правый сосед есть у всех узлов, кроме последнего слоя сетки по j
				JA[indecies++] = G2L[I_node + 1]; // → вправо
			if (Ny - i_g - 1)		//Нижний сосед есть у всех узлов, кроме последнего слоя сетки по i
				JA[indecies++] = G2L[I_node + Nx]; // ↓ вниз
			if (Nx - j_g - 1 && Ny - i_g - 1 && (I % kk >= k1))
				JA[indecies++] = G2L[I_node + Nx + 1]; // ↘ вправо-вниз
		}
	}
}

//Заполнение матрицы A и правой части b по заданным формулам в CSR формате
void fill(const int* IA, const int* JA, const int* L2G, const int No, const int N, double*& A, double*& b) {
	for (int i = 0; i < No; ++i) {
		double sum = 0;
		int index_i; // Переменная, хранящая диагональный индекс
		for (int col = IA[i]; col < IA[i + 1]; ++col) {
			if (L2G[JA[col]] - L2G[i]) { //i!=j
				A[col] = cos(L2G[i] * L2G[JA[col]] + L2G[i] + L2G[JA[col]]);
				sum += fabs(A[col]); //∑|aij|, i!=j
			}
			else//i=j
				index_i = col;
		}
		A[index_i] = 1.234 * sum;
	}
	for (int i = 0; i < N; ++i) {
		b[i] = sin(L2G[i]);
	}
}

//Схема обмена граничных элементов для обновления гало
void Comm(const int N, const int No, const int* IA, const int* JA, const int* Part, const int* L2G, const int* G2L,
	const int px, const int py, const int k1, const int k2, const int Nx, const int Nx_loc,
	std::vector<int>& Neighbours, int*& RecvOffset, int*& SendOffset, int*& Recv, int*& Send) {

	for (int i = 0; i < N; ++i) {
		if (Part[i] != rank && Part[i] != Part[i - 1])
			Neighbours.push_back(Part[i]);
	}

	int* P2N = new int[ranksize];
	for (int i = 0; i < ranksize; ++i) {
		P2N[i] = -1;
	}
	for (int i = 0; i < Neighbours.size(); ++i) {
		P2N[Neighbours[i]] = i;
	}

	std::vector<std::vector<int>> RecvFromProcess(Neighbours.size());
	std::vector<std::vector<int>> SendToProcess(Neighbours.size());
	for (int i = 0; i < No; ++i) {
		for (int col = IA[i]; col < IA[i + 1]; ++col) {
			if (JA[col] >= No) {
				RecvFromProcess[P2N[Part[JA[col]]]].push_back(JA[col]);
				SendToProcess[P2N[Part[JA[col]]]].push_back(i);
			}
		}
	}

	for (auto& vec : RecvFromProcess) {
		auto last = std::unique(vec.begin(), vec.end());
		vec.erase(last, vec.end());
	}
	for (auto& vec : SendToProcess) {
		auto last = std::unique(vec.begin(), vec.end());
		vec.erase(last, vec.end());
	}

	RecvOffset = new int[Neighbours.size() + 1];
	SendOffset = new int[Neighbours.size() + 1];
	SendOffset[0] = 0;	RecvOffset[0] = 0;
	int total_size1 = 0;
	for (const auto& row : RecvFromProcess) {
		total_size1 += row.size();
	}
	int total_size2 = 0;
	for (const auto& row : SendToProcess) {
		total_size2 += row.size();
	}
	Recv = new int[total_size1];
	Send = new int[total_size2];
	int index = 0;
	for (int i = 0; i < Neighbours.size(); ++i) {
		RecvOffset[i + 1] = RecvOffset[i] + RecvFromProcess[i].size();
		SendOffset[i + 1] = SendOffset[i] + SendToProcess[i].size();
		for (int j = 0; j < RecvFromProcess[i].size(); ++j) {
			Recv[index] = RecvFromProcess[i][j];
			Send[index++] = SendToProcess[i][j];
		}
	}	
	delete[] P2N;
}

void UpdateInit(const std::vector<int>& Neighbours, const int* RecvOffset, const int* SendOffset,
	const int* Recv, const int* Send, MPI_Request*& reqs, MPI_Status*& stats, std::vector<double>& SENDBUF, std::vector<double>& RECVBUF){
	int B = Neighbours.size();

	if (B == 0) return; 
	int sendCount = SendOffset[B]; // размер общего списка на отправку по всем соседям
	int recvCount = RecvOffset[B]; // размер общего списка на прием по всем соседям

	// MPI данные - сделаем статиками, поскольку это высокочастотная функция,
	// чтобы каждый раз не реаллокать (так делать вовсе не обязательно).

	
	// ресайзим, надо
	if (sendCount > (int)SENDBUF.size()) SENDBUF.resize(sendCount);
	if (recvCount > (int)RECVBUF.size()) RECVBUF.resize(recvCount);

	reqs = new MPI_Request[2*B];// реквесты для неблокирующих обменов
	stats = new MPI_Status[2*B];// статусы для неблокирующих обменов

}

void UpdateClose(MPI_Request*& reqs, MPI_Status*& stats){
	if (reqs) delete[] reqs;
	if (stats) delete[] stats;
}

//Обновление гало по схеме Comm
void Update(double*& V, const std::vector<int>& Neighbours, const int* RecvOffset, const int* SendOffset,
	const int* Recv, const int* Send, MPI_Request*& reqs, MPI_Status*& stats, std::vector<double>& SENDBUF, std::vector<double>& RECVBUF) {

	int B = Neighbours.size();

	if (B == 0) return; // нет соседей - нет проблем
	int sendCount = SendOffset[B]; // размер общего списка на отправку по всем соседям
	int recvCount = RecvOffset[B]; // размер общего списка на прием по всем соседям


	int nreq = 0; // сквозной счетчик реквестов сообщений

	for (int p = 0; p < B; p++) {
		int SZ = (RecvOffset[p + 1] - RecvOffset[p]); // размер сообщения
		if (SZ <= 0) continue; // если нечего слать - пропускаем соседа
		int NB_ID = Neighbours[p]; // узнаем номер процесса данного соседа
		int mpires = MPI_Irecv(&RECVBUF[RecvOffset[p]], SZ, MPI_DOUBLE, NB_ID, 0, MCW, &reqs[nreq]);
		if (mpires != MPI_SUCCESS)
			crash("MPI_Irecv failed");
		nreq++;
	}

	// пакуем данные с интерфейса по единому списку сразу по всем соседям
	for (int i = 0; i < sendCount; ++i) SENDBUF[i] = V[Send[i]/*номер ячейки на отправку*/];

	// инициируем отправку сообщений
	for (int p = 0; p < B; p++) {
		int SZ = (SendOffset[p + 1] - SendOffset[p]); // размер сообщения
		if (SZ <= 0) continue; // если нечего принимать - пропускаем соседа
		int NB_ID = Neighbours[p]; // узнаем номер процесса данного соседа
		int mpires = MPI_Isend(&SENDBUF[SendOffset[p]], SZ, MPI_DOUBLE, NB_ID, 0, MCW, &reqs[nreq]);
		if (mpires != MPI_SUCCESS)
			crash("MPI_Isend failed");
		nreq++;
	}

	if (nreq > 0) { // ждем завершения всех обменов
		int mpires = MPI_Waitall(nreq, reqs, stats);
		if (mpires != MPI_SUCCESS)
			crash("MPI_Waitall failed");
	}
	// разбираем данные с гало ячеек по единому списку сразу по всем соседям

	for (int i = 0; i < recvCount; ++i)
		V[Recv[i]] = RECVBUF[i]; //Recv[i]/*номер ячейки на прием*/
}

void axpy(const double* a, const double beta, const double* b, double*& res, const int N) {
	for (int i = 0; i < N; ++i) {
		res[i] = a[i] + beta * b[i];
	}
}

void SpmV_M(const double* M, const double* v, double*& res, const int N) {
	for (int i = 0; i < N; ++i) {
		res[i] = M[i] * v[i];
	}
}

//Умножение матрицы на вектор
void SpmV(const int* IA, const int* JA, const double* A, const double* v, double*& res, const int N) {
	for (int i = 0; i < N; ++i) {
		double sum = 0.0;
		for (int col = IA[i]; col < IA[i + 1]; ++col) {
			sum += A[col] * v[JA[col]];
		}
		res[i] = sum;
	}
}

double dot(const double* a, const double* b, const int N) {
	double lsum = .0;
	double gsum;
	for (int i = 0; i < N; ++i)
		lsum += a[i] * b[i];
	int mpires = MPI_Allreduce(&lsum, &gsum, 1, MPI_DOUBLE, MPI_SUM, MCW);
	if (mpires != MPI_SUCCESS)
		crash("MPI_Allreduce failed");
	return gsum;
}


//Решение СЛАУ методом сопряжённых гардиентов с применением OpenMP
void solve(const int* IA, const int* JA, const double* A, double* b, const int N, const int No, const double eps,
	const int maxit, const bool norms, double*& x, int& step, double& norm,
	const std::vector<int>& Neighbours, const int* RecvOffset, const int* SendOffset, const int* Recv, const int* Send, 
MPI_Request*& reqs, MPI_Status*& stats, std::vector<double>& SENDBUF, std::vector<double>& RECVBUF)
{
	double* r = new double[N];
	double* p = new double[N];
	double* q = new double[N];
	double* z = new double[N];
	double* M = new double[No];

	for (int i = 0; i < N; ++i) {
		x[i] = 0;
	}

	for (int i = 0; i < N; ++i) {
		r[i] = b[i];
	}

	for (int i = 0; i < No; ++i) {
		for (int col = IA[i]; col < IA[i + 1]; ++col) {
			if (!(JA[col] - i)) { //i=j
				M[i] = 1.0 / A[col];
				break;
			}
		}
	}

	step = 0;
	double po, beta, po_temp, alpha;
	do {
		++step;
		SpmV_M(M, r, z, No);
		po = dot(r, z, No);
		if (step == 1) {

			for (int i = 0; i < N; ++i) {
				p[i] = z[i];
			}
		}
		else {
			beta = po / po_temp;
			axpy(z, beta, p, p, No);
		}
		Update(p, Neighbours, RecvOffset, SendOffset, Recv, Send, reqs, stats, SENDBUF, RECVBUF);
		SpmV(IA, JA, A, p, q, No);
		alpha = po / dot(p, q, No);
		axpy(x, alpha, p, x, No);
		axpy(r, -alpha, q, r, No);
		po_temp = po;
		if (norms) {
			norm = sqrt(dot(r, r, No));
			std::cout << "solve[" << step << "]: l2 norm=" << norm << std::endl;
		}
	} while (po > eps * eps && step < maxit);

	norm = sqrt(dot(r, r, No));
	delete[] z;
	delete[] p;
	delete[] r;
	delete[] q;
	delete[] M;
}

void report(const double norm, const double timeG, const double timeF, const double timeS, const double dott, const double axpyt, const double SpmVt, const int N) {
	std::cout << std::endl;
	std::cout << "Number of threads: " << ranksize << std::endl;
	std::cout << "L2 norm of solving Ax-b: " << norm << std::endl;
	std::cout << "Generation step time: " << timeG << std::endl;
	std::cout << "Fill matrix step time: " << timeF << std::endl;
	std::cout << "Solve matrix step time: " << timeS << std::endl;
	printf("dot[%d] time = %g  \n", N, dott);
	printf("axpy[%d] time = %g \n", N, axpyt);
	printf("SpmV[%d] time = %g \n", N, SpmVt);

}

int main(int argc, char** argv) {
	int mpi_res;

	mpi_res = MPI_Init(&argc, &argv);
	if (mpi_res != MPI_SUCCESS)
		crash("MPI_Init failed (code %d)\n", mpi_res);
	mpi_res = MPI_Comm_rank(MCW, &rank);
	if (mpi_res != MPI_SUCCESS)
		crash("MPI_Comm_rank failed (code %d)\n", mpi_res);
	mpi_res = MPI_Comm_size(MCW, &ranksize);
	if (mpi_res != MPI_SUCCESS)
		crash("MPI_Comm_size failed (code %d)\n", mpi_res);
	mpi_initialized = 1;

	int Nx, Ny, k1, k2, px, py;
	bool is_out = 0, norms = 0, is_check = 0;
	if (rank == 0) {
		if (argc == 1) {
			std::cout << "___________________________________________________\n"
				<< "|           ULTIMATE CSR MAKER VIA MESH           |\n"
				<< "___________________________________________________\n"
				<< "\nHow to launch me correctly?\n"
				<< "\nType command like:\n"
				<< argv[0] << " Nx Ny k1 k2 px py\n"
				<< "where Nx, Ny, k1, k2 are integers\n\n"
				<< "Nx is size of mesh in horizontal\n"
				<< "Ny is size of mesh in vertical\n"
				<< "k1 is number of full cells in a row\n"
				<< "k2 is number of split cells in a row\n"
				<< "px is number of procs in horizontal\n"
				<< "py is number of procs in vertical\n"
				<< "\nIf you wonna out CSR data, please add \"out\" like:\n"
				<< argv[0] << " Nx Ny k1 k2 px py out\n"
				<< "\nIf you wonna out l2 norms by step while solving, please add \"l2\" like:\n"
				<< argv[0] << " Nx Ny k1 k2 px py l2\n"
				<< "\nIf you wonna out checks, please add \"check\" like:\n"
				<< argv[0] << " Nx Ny k1 k2 px py check\n";
			crash("Please, launch me again with correct data :)\n");
		}
	}
	barrier();
	if (argc < 7)
		crash("Incorrect input data. Please, use me correctly\n");
	
	try {
		Nx = std::stoi(argv[1]);
		Ny = std::stoi(argv[2]);
		k1 = std::stoi(argv[3]);
		k2 = std::stoi(argv[4]);
		px = std::stoi(argv[5]);
		py = std::stoi(argv[6]);
	}
	catch (std::invalid_argument& ex) { crash("std::invalid_argument::what(): %s\n", ex.what()); }
	catch (std::out_of_range const& ex) { crash("std::out_of_range::what(): %s\n", ex.what()); }


	if (argc >= 8) {
		if (!strcmp(argv[7], "out"))
			is_out = true;
		else
			if (argc >= 9 && !strcmp(argv[8], "out"))
				is_out = true;
			else
				if (argc >= 10 && !strcmp(argv[9], "out"))
					is_out = true;
		if (!strcmp(argv[7], "l2"))
		norms = true;
		else
			if (argc >= 9 && !strcmp(argv[8], "l2"))
				norms = true;
			else
				if (argc >= 10 && !strcmp(argv[9], "l2"))
					norms = true;
		if (!strcmp(argv[7], "check"))
		is_check = true;
		else
			if (argc >= 9 && !strcmp(argv[8], "check"))
				is_check = true;
			else
				if (argc >= 10 && !strcmp(argv[9], "check"))
					is_check = true;
	}
	

	if (rank == 0) {
		if (Nx < 1 || Ny < 1 || k1 == 0 && k2 == 0 || k1 < 0 || k2 < 0 || px < 0 || py < 0 || Nx < px || Ny < py)
			crash("Incorrect value input data. Please, use me correctly\n");
	}
	int step, soseds;
	int* L2G = NULL, * Part = NULL, * G2L = NULL;
	std::vector<int> Neighbours;
	int* RecvOffset = NULL, * SendOffset = NULL, * Recv = NULL, * Send = NULL;
	int* JA = NULL, * IA = NULL;
 	MPI_Request* reqs = NULL; 
	MPI_Status* stats = NULL;
	static std::vector<double> SENDBUF, RECVBUF; // буферы на отправку и прием по всем соседям
	double norm;
	int Np, No, N;
	double time, timeG, timeF, timeS;
	int Nx_loc, Ny_loc;
	Nx++; Ny++;
	Nx_loc = (Nx % px) ? Nx / px + ((rank % px) < (Nx % px)) : Nx / px;
	Ny_loc = (Ny % py) ? Ny / py + ((rank / px) < (Ny % py)) : Ny / py;

	barrier();
	time = MPI_Wtime();
	generate(Nx, Ny, k1, k2, px, py, Nx_loc, Ny_loc, N, Np, No, L2G, Part, G2L, IA, JA);
	barrier();
	timeG = MPI_Wtime() - time;
	double* A = new double[IA[No]];
	double* b = new double[N];
	barrier();
	time = MPI_Wtime();
	fill(IA, JA, L2G, No, N, A, b);

	barrier();
	timeF = MPI_Wtime() - time;

	Comm(N, No, IA, JA, Part, L2G, G2L, px, py, k1, k2, Nx, Nx_loc, Neighbours, RecvOffset, SendOffset, Recv, Send);

	double* x  = new double[N];

	barrier();
	time = MPI_Wtime();
	UpdateInit(Neighbours, RecvOffset, SendOffset, Recv, Send, reqs, stats, SENDBUF, RECVBUF);
	for (int i = 0; i < 10; ++i)
		solve(IA, JA, A, b, N, No, 1e-7, 10000, norms, x, step, norm, Neighbours, RecvOffset, SendOffset, Recv, Send, reqs, stats, SENDBUF, RECVBUF);
	UpdateClose(reqs, stats);
	barrier();
	timeS = (MPI_Wtime() - time) / step / 10;

	if (is_out && rank == 0) {
		print_full_csr(IA, JA, No, N, A, b, x);
	}
	
	double* a, * c;
	a = new double[N];
	c = new double[N];
	for (int i = 0; i < N; ++i) {
		a[i] = i % 456;
		c[i] = (i * i) % 134;
	}

	barrier();
	time = MPI_Wtime();
	for (int i = 50; i--;)
		dot(a, a, No);
	barrier();
	double dott = (MPI_Wtime() - time) / 50;
	barrier();
	time = MPI_Wtime();
	for (int i = 50; i--;)
		axpy(a, 1.4, c, a, No);
	barrier();
	double axpyt = (MPI_Wtime() - time) / 50;
	barrier();
	time = MPI_Wtime();
	for (int i = 50; i--;)
		SpmV(IA, JA, A, c, c, No);
	barrier();
	double SpmVt = (MPI_Wtime() - time) / 50;

	
	if (rank == 0) {
		report(norm, timeG, timeF, timeS, dott, axpyt, SpmVt, No);
	}
	if (is_check) {
		double check1 = dot(x, x, No);//l2 norm
		double lcheck2 = 0, lcheck3 = x[0], check2, check3;
		for (int k = 0; k < No; ++k)
			lcheck2 += x[k];
		MPI_Reduce(&lcheck2, &check2, 1, MPI_DOUBLE, MPI_SUM, 0, MCW);
		for (int k = 1; k < No; ++k)
			lcheck3 = std::min(x[k], lcheck3);
		MPI_Reduce(&lcheck3, &check3, 1, MPI_DOUBLE, MPI_MIN, 0, MCW);
		if (rank == 0) {
			std::cout << std::endl
				<< "l2(x) MPI: " << check1 << std::endl
				<< "sum(x) MPI: " << check2 << std::endl
				<< "min(x) MPI: " << check3 << std::endl;
		}

	}
	
	
	if (IA) delete[] IA;
	if (JA) delete[] JA;
	if (A) delete[] A;
	if (b) delete[] b;
	if (x) delete[] x;
	if (RecvOffset) delete[] RecvOffset;
	if (SendOffset) delete[] SendOffset;
	if (Recv) delete[] Recv;
	if (Send) delete[] Send;
	if (L2G) delete[] L2G;
	if (Part) delete[] Part;
	if (G2L) delete[] G2L;
	delete[] a;
	delete[] c;



	MPI_Finalize();
	return 0;
}
