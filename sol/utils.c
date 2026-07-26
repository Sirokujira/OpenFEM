/*
utils.c

CPU 時間・実行条件の表示・コマンドライン引数。
*/

#include "fem.h"
#include "fem_prototype.h"

double cputime(void)
{
#ifdef _WIN32
	return (double)clock() / CLOCKS_PER_SEC;
#else
	struct timespec ts;
	clock_gettime(CLOCK_REALTIME, &ts);
	return (ts.tv_sec + (ts.tv_nsec * 1e-9));
#endif
}


// 標準出力とログの両方に出す
void monitor1(FILE *fp, const char *msg)
{
	printf("%s\n", msg);
	fflush(stdout);
	if (fp != NULL) {
		fprintf(fp, "%s\n", msg);
		fflush(fp);
	}
}


// 解析条件
void monitor2(FILE *fp, int nthread)
{
	time_t now;
	time(&now);

	const int64_t nnode = num_node();
	const int64_t ncell = (int64_t)Nx * Ny * Nz;

	// CRS 行列 1 個分のおおよそのメモリ (27 点ステンシル)
	const double mem = ((nnode * 27.0 * (sizeof(double) + sizeof(int32_t)))
	                 + (nnode * 6.0 * sizeof(double))) / (1024 * 1024);

	fprintf(fp, "%s", ctime(&now));
	fprintf(fp, "Title = %s\n", Title);
	fprintf(fp, "Threads = %d\n", nthread);
	fprintf(fp, "Cells = %d x %d x %d = %lld\n", Nx, Ny, Nz, (long long)ncell);
	fprintf(fp, "Nodes = %d x %d x %d = %lld\n", Nx + 1, Ny + 1, Nz + 1, (long long)nnode);
	fprintf(fp, "No. of Materials  = %d\n", NMaterial);
	fprintf(fp, "No. of Geometries = %d\n", NGeometry);
	fprintf(fp, "No. of Conductors = %d (ports = %d)\n", NConductor, NPort);
	fprintf(fp, "Analysis =%s%s%s\n",
		((Analysis & ANALYSIS_C) ? " C" : ""),
		((Analysis & ANALYSIS_L) ? " L" : ""),
		((Analysis & ANALYSIS_R) ? " R" : ""));
	if (Tline) {
		fprintf(fp, "Transmission line axis = %c (length = %.6e [m])\n", Tline, TlineLength);
	}
	fprintf(fp, "Excitation voltage = %.6e [V]\n", Volt);
	fprintf(fp, "Solver = %d %d %.3e\n", Solver.maxiter, Solver.nout, Solver.converg);
	fprintf(fp, "Memory size = %.1f [MB]\n", mem);
	fflush(fp);
}


static void usage(void)
{
	printf("Usage: ofe [-n <thread>] [--help] <file.ofe>\n");
	printf("  -n <thread> : number of OpenMP threads (default 1)\n");
}


void comline(int argc, char *argv[], int *nthread, int *prompt, char *fn_in)
{
	*nthread = 1;
	*prompt = 0;
	strcpy(fn_in, "");

	int n = 1;
	while (n < argc) {
		if      (!strcmp(argv[n], "-n")) {
			if (n + 1 < argc) {
				*nthread = atoi(argv[++n]);
				if (*nthread < 1) *nthread = 1;
			}
		}
		else if (!strcmp(argv[n], "--prompt")) {
			*prompt = 1;
		}
		else if (!strcmp(argv[n], "--help") || !strcmp(argv[n], "-h")) {
			usage();
			exit(0);
		}
		else if (!strcmp(argv[n], "--version") || !strcmp(argv[n], "-v")) {
			printf("%s Ver.%d.%d.%d\n", PROGRAM, VERSION_MAJOR, VERSION_MINOR, VERSION_BUILD);
			exit(0);
		}
		else {
			strcpy(fn_in, argv[n]);
		}
		n++;
	}

	if (strlen(fn_in) == 0) {
		usage();
		exit(1);
	}
}
