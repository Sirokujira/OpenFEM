/*
mpiwrap.c

MPI の薄い包み (任意依存、既定 OFF)。

並列化の粒度は**周波数掃引の点** (solve() が点をプロセスに round-robin で
配り、rank 0 が結果を集めて書く)。各点の演算は直列と完全に同一なので、
**mpirun -np N の出力 (ofe_sweep.csv / rlc.csv / ofe.out) は直列と
バイト単位で一致する** — これが検証の恒等式になる (rlc_check.sh)。

mpi.h をこのファイルの外に出さないこと。WITH_MPI=OFF のビルドは
mpi.h を include すらせず、MSVC を含む 3 OS のビルドが一切変わらない。

**無効ビルドで mpirun 配下に置かれたら入力エラーにする。** 黙って直列に
落ちると、同じ作業ディレクトリで N 個の直列プロセスが同じ出力ファイルに
書いて壊し合う (書いたつもりでファイルが壊れているのが一番困る)。
mpirun が立てる環境変数 (OpenMPI / MPICH / Intel MPI / MS-MPI) で検出する。
*/

#include "fem.h"
#include "fem_prototype.h"

#ifdef OFE_MPI

#include <mpi.h>

int mpi_init(int *argc, char ***argv)
{
	MPI_Init(argc, argv);

	return 0;
}


int mpi_rank(void)
{
	int r = 0;

	MPI_Comm_rank(MPI_COMM_WORLD, &r);

	return r;
}


int mpi_size(void)
{
	int s = 1;

	MPI_Comm_size(MPI_COMM_WORLD, &s);

	return s;
}


void mpi_finalize(void)
{
	MPI_Finalize();
}


void mpi_send_dbl(const double *buf, int n, int dst, int tag)
{
	MPI_Send(buf, n, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD);
}


void mpi_recv_dbl(double *buf, int n, int src, int tag)
{
	MPI_Recv(buf, n, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

#else		// !OFE_MPI : スタブ (直列)

int mpi_init(int *argc, char ***argv)
{
	(void)argc;
	(void)argv;

	// mpirun 配下かどうかを主要な MPI 実装の環境変数で検出する
	static const char *marks[] = {
		"OMPI_COMM_WORLD_SIZE",		// OpenMPI
		"PMI_SIZE",					// MPICH / Intel MPI / MS-MPI
		"PMIX_RANK",				// OpenMPI 5 (PMIx)
		"SLURM_NTASKS",				// slurm srun
	};
	for (size_t i = 0; i < sizeof(marks) / sizeof(marks[0]); i++) {
		const char *v = getenv(marks[i]);
		if ((v != NULL) && (atoi(v) > 1)) {
			printf("*** this binary was built without MPI (WITH_MPI=OFF); "
				"rebuild with -DWITH_MPI=ON, or run it without mpirun "
				"(N serial processes would corrupt each other's output files)\n");
			return 1;
		}
	}

	return 0;
}


int mpi_rank(void)
{
	return 0;
}


int mpi_size(void)
{
	return 1;
}


void mpi_finalize(void)
{
}


void mpi_send_dbl(const double *buf, int n, int dst, int tag)
{
	(void)buf; (void)n; (void)dst; (void)tag;
}


void mpi_recv_dbl(double *buf, int n, int src, int tag)
{
	(void)buf; (void)n; (void)src; (void)tag;
}

#endif		// OFE_MPI
