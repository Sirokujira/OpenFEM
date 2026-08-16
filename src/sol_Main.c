/*
OpenFEM (CPU + OpenMP)

準静的 FEM ソルバー ofe のエントリ。
	入力読込 (.ofe) -> セットアップ -> solve() -> 回路パラメータ出力
*/

#define MAIN
#include "fem.h"
#undef MAIN

#include "fem_prototype.h"

static void error_check(int ierr, int prompt)
{
	if (ierr) {
		if (prompt) {
			printf("Press ENTER key.\n");
			getchar();
		}
		exit(1);
	}
}


int main(int argc, char *argv[])
{
	const char errfmt[] = "*** file %s open error.\n";
	char str[BUFSIZ];
	int ierr = 0;
	FILE *fp_in = NULL, *fp_out = NULL, *fp_log = NULL;

	// MPI (WITH_MPI=ON のときだけ実体がある。スタブは「無効ビルドで
	// mpirun 配下に置かれた」ことを検出してエラーにする)
	if (mpi_init(&argc, &argv)) {
		exit(1);
	}
	const int mrank = mpi_rank();

	int nthread = 1;
	int prompt = 0;
	char fn_in[BUFSIZ];
	comline(argc, argv, &nthread, &prompt, fn_in);

#ifdef _OPENMP
	omp_set_num_threads(nthread);
#endif

	const double t0 = cputime();

	// 入力
	if ((fp_in = fopen(fn_in, "r")) == NULL) {
		printf(errfmt, fn_in);
		error_check(1, prompt);
	}
	ierr = input_data(fp_in);
	fclose(fp_in);
	error_check(ierr, prompt);

	// ログ。rank 0 以外は tmpfile に書く (全 rank が同じ ofe.log に
	// 書くと壊し合う。遠隔の点の失敗は solve() が rank 0 のログに集める)
	fp_log = ((mrank == 0) ? fopen(FN_log, "w") : tmpfile());
	if (fp_log == NULL) {
		printf(errfmt, FN_log);
		error_check(1, prompt);
	}

	if (mrank == 0) {
		sprintf(str, "<<< %s (CPU+OpenMP) Ver.%d.%d.%d >>>",
			PROGRAM, VERSION_MAJOR, VERSION_MINOR, VERSION_BUILD);
		monitor1(fp_log, str);
	}

	// 入力解釈の警告は ofe.log を開く前に出るので、ここでログにも残す
	for (int w = 0; w < NInputWarn; w++) {
		fprintf(fp_log, "*** warning : %s\n", InputWarn[w]);
	}

	// セットアップ
	ierr = setup();
	if (ierr) {
		fprintf(fp_log, "*** setup error\n");
		fclose(fp_log);
		error_check(ierr, prompt);
	}
	if (mrank == 0) monitor2(fp_log, nthread);

	const double t1 = cputime();

	// 求解
	ierr = solve(fp_log);
	if (ierr) {
		fprintf(fp_log, "*** solver error\n");
		fclose(fp_log);
		error_check(ierr, prompt);
	}

	const double t2 = cputime();

	// 出力 (rank 0 だけが書く)
	if (mrank == 0) {
		outputRLC(fp_log);

		if ((fp_out = fopen(FN_out, "wb")) == NULL) {
			printf(errfmt, FN_out);
			fclose(fp_log);
			error_check(1, prompt);
		}
		writeout(fp_out);
		fclose(fp_out);

		const double t3 = cputime();

		fprintf(fp_log, "\ncpu time [sec] : setup = %.3f, solve = %.3f, output = %.3f, total = %.3f\n",
			t1 - t0, t2 - t1, t3 - t2, t3 - t0);
		fprintf(fp_log, "output files : %s, %s\n", FN_log, FN_out);

		monitor1(fp_log, "=== normal end ===");
	}

	fclose(fp_log);
	memfree();
	mpi_finalize();

	return 0;
}
