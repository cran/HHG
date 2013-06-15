//============================================================================
// Name        : HHG.cpp
// Author      : Shachar Kaufman
// Version     :
// Copyright   : TAU
// Description : An implementation of the Heller-Heller-Gorfine test
//============================================================================

// BUILD NOTE:
// - The standalone debugging variant should be built with -DDEBUG_INTERFACE -lpthread
// - On Windows, the R variant of this needs to be built with (don't forget to define R_INTERFACE):
//   R CMD SHLIB HHG.cpp SequentialTest.cpp StatsComputer.cpp -pthread
// - The standalone GWAS variant should be built with -DGWAS_INTERFACE -lpthread

#ifdef WIN32
#include <Windows.h>
#else
#include <unistd.h>
#include <sys/time.h>
#endif

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <pthread.h>

#include "HHG.h"

#include "SequentialTest.h"

using namespace std;

static int get_available_nr_threads(void);
#ifndef R_INTERFACE
static void read_data_file(ifstream& ifs, int skp_first_lines, int skp_first_fields, int n, int m, double* buf);
static void compute_01_distances(int n, int p, double* v, double *dv);
static void compute_hamming_distances(int n, int p, double* v, double *dv);
static void compute_lp_distances(int n, int p, int norm_p, double* v, double *dv);
static void compute_rm_distances(int n, int p, double* v, double *dv);
static string portable_double2str(double);
#endif

#ifdef DEBUG_INTERFACE

static void test_hhg(string test_type_str, int nr_reps, int nr_perm, int nr_threads);
static long long get_time_ms(void);

// Entry point
// ============================================================================

int main(int argc, char** argv) {
	if (argc != 4) {
		cerr << "Usage: " << argv[0] << " test_type nr_reps nr_perm" << endl;
		exit(1);
	}

	string test_type_str = argv[1];
	int nr_reps = atoi(argv[2]);
	int nr_perm = atoi(argv[3]);

	int nr_threads = get_available_nr_threads();

	test_hhg(test_type_str, nr_reps, nr_perm, nr_threads);

#ifndef NO_THREADS
	pthread_exit(NULL); // recommended by pthreads documentation
#endif

	return 0;
}

static void test_hhg(string test_type_str, int nr_reps, int nr_perm, int nr_threads) {
	TestType tt;
	string x_filename, y_filename;

#ifdef WIN32
	x_filename = "C:\\Library\\Work with Ruth Heller\\HHG\\";
	y_filename = "C:\\Library\\Work with Ruth Heller\\HHG\\";
#else
	x_filename = "../";
	y_filename = "../";
#endif

	if (test_type_str == "g") {
		tt = GENERAL_TEST;
		x_filename += "tmp.g.data.Dx";
		y_filename += "tmp.g.data.Dy";
	} else if (test_type_str == "k") {
		tt = K_SAMPLE_TEST;
		x_filename += "tmp.k.data.Dx";
		y_filename += "tmp.k.data.y";
	} else if (test_type_str == "2") {
		tt = TWO_SAMPLE_TEST;
		x_filename += "tmp.2.data.Dx";
		y_filename += "tmp.2.data.y";
	} else if (test_type_str == "1") {
		tt = UDF_PPR_33_OBS;
		x_filename += "tmp.1.data.x";
		y_filename += "tmp.1.data.y";
	}

	double alpha = 0.002;
	double alpha0 = 0.05;
	double beta0 = 0.000025;
	double eps = 0.01;
	double w_sum = 0;
	double w_max = 2;

	ifstream x_ifs(x_filename.c_str());
	ifstream y_ifs(y_filename.c_str());
	if (!x_ifs || !y_ifs) {
		cerr << "where is the data?" << endl;
		exit(1);
	}

	string buf;
	int n = 0;
	while (getline(y_ifs, buf)) ++n;
	y_ifs.clear();
	y_ifs.seekg(0);

	// TODO determine from input files
	int y_ncol = 1;

	double* y = new double[n * y_ncol];
	double* dx = new double[n * n];
	double* dy = new double[n * n];

	cout << "Reading data from debug files" << endl;

	switch (tt) {
		case TWO_SAMPLE_TEST:
		case K_SAMPLE_TEST:
			read_data_file(x_ifs, 0, 0, n, n, dx);
			read_data_file(y_ifs, 0, 0, n, 1, y);
		break;

		case GENERAL_TEST:
			read_data_file(x_ifs, 0, 0, n, n, dx);
			read_data_file(y_ifs, 0, 0, n, n, dy);
		break;

		case UDF_SPR_OBS:
		case UDF_SPR_ALL:
			read_data_file(x_ifs, 0, 0, n, 1, dx);
			read_data_file(y_ifs, 0, 0, n, 1, y);
		break;

		default:
			cerr << "What test is this?" << endl;
			exit(1);
		break;
	}

	cout << "Computing HHG" << endl;

	double sum_chi, sum_like, max_chi, max_like, ht, edist;
	double p_sum_chi, p_sum_like, p_max_chi, p_max_like, p_ht, p_edist;

	long long ts_start = get_time_ms();
	for (int i = 0; i < nr_reps; ++i) {
		cout << "  Iteration " << i << "..." << endl;

		SequentialTest seq(tt, n, y_ncol, dx, dy, y, w_sum, w_max, true, alpha, alpha0, beta0, eps, nr_perm, nr_threads, 0, false, false);

		seq.get_observed_stats(sum_chi, sum_like, max_chi, max_like, ht, edist);
		seq.get_pvalues(p_sum_chi, p_sum_like, p_max_chi, p_max_like, p_ht, p_edist);
	}

	cout << "Milliseconds elapsed during computation: " << get_time_ms() - ts_start << endl;

	cout << "Computed statistics (and p-values):" << endl;
	cout << "SC: " << sum_chi << " (" << p_sum_chi << ")" << endl;
	cout << "SL: " << sum_like << " (" <<  p_sum_like << ")"  << endl;
	cout << "MC: " << max_chi << " (" <<  p_max_chi << ")"  << endl;
	cout << "ML: " << max_like << " (" <<  p_max_like << ")"  << endl;
	cout << "HT: " << ht << " (" <<  p_ht << ")"  << endl;
	cout << "ED: " << edist << " (" <<  p_edist << ")"  << endl;

	delete[] y;
	delete[] dx;
	delete[] dy;
}

static long long get_time_ms(void) {
#ifdef WIN32
	FILETIME ft;
	LARGE_INTEGER li;

	/* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC) and copy it
	 * to a LARGE_INTEGER structure. */
	GetSystemTimeAsFileTime(&ft);
	li.LowPart = ft.dwLowDateTime;
	li.HighPart = ft.dwHighDateTime;

	unsigned long long ret = li.QuadPart;
	ret -= 116444736000000000LL; /* Convert from file time to UNIX epoch time. */
	ret /= 10000; /* From 100 nano seconds (10^-7) to 1 millisecond (10^-3) intervals */

	return ret;
#else // => Linux
	struct timeval tv;

	gettimeofday(&tv, NULL);

	unsigned long long ret = tv.tv_usec;
	/* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
	ret /= 1000;

	/* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
	ret += (tv.tv_sec * 1000);

	return ret;
#endif // WIN32/Linux
}

#endif // DEBUG_INTERFACE

#ifdef R_INTERFACE

extern "C" {

SEXP HHG(SEXP R_test_type, SEXP R_dx, SEXP R_dy, SEXP R_y,
		 SEXP R_w_sum, SEXP R_w_max, SEXP R_extra_params,
		 SEXP R_is_sequential, SEXP R_alpha, SEXP R_alpha0, SEXP R_beta0, SEXP R_eps, SEXP R_nr_perm,
		 SEXP R_nr_threads, SEXP R_tables_wanted, SEXP R_perm_stats_wanted)
{
	try {
		int test_type = *INTEGER(R_test_type);
		double* dx = REAL(R_dx);
		double* dy = REAL(R_dy);
		double* y = REAL(R_y);
		double w_sum = *REAL(R_w_sum);
		double w_max = *REAL(R_w_max);
		double* extra_params = REAL(R_extra_params);
		bool is_sequential = (*INTEGER(R_is_sequential) != 0);
		double alpha = *REAL(R_alpha);
		double alpha0 = *REAL(R_alpha0);
		double beta0 = *REAL(R_beta0);
		double eps = *REAL(R_eps);
		int nr_perm = *INTEGER(R_nr_perm);
		int nr_threads = *INTEGER(R_nr_threads);
		bool tables_wanted = (*INTEGER(R_tables_wanted) != 0);
		bool perm_stats_wanted = (*INTEGER(R_perm_stats_wanted) != 0);

		if (nr_threads == 0) {
			nr_threads = get_available_nr_threads();
		}

		SEXP Rdim = getAttrib(R_y, R_DimSymbol);
		int n = INTEGER(Rdim)[0];
		int y_ncol = INTEGER(Rdim)[1];

		double sum_chi, sum_like, max_chi, max_like, ht, edist;
		double p_sum_chi, p_sum_like, p_max_chi, p_max_like, p_ht, p_edist;
		SequentialTest seq((TestType)test_type, n, y_ncol, dx, dy, y, w_sum, w_max, extra_params, is_sequential, alpha, alpha0, beta0, eps, nr_perm, nr_threads, 0, tables_wanted, perm_stats_wanted);

		SEXP R_res;
		PROTECT(R_res = allocMatrix(REALSXP, 12 + 4 * n * n * tables_wanted + 6 * nr_perm * perm_stats_wanted, 1));
		double* res = REAL(R_res);

		seq.get_observed_stats(sum_chi, sum_like, max_chi, max_like, ht, edist);
		seq.get_pvalues(p_sum_chi, p_sum_like, p_max_chi, p_max_like, p_ht, p_edist);

		res[0] = p_sum_chi;
		res[1] = p_sum_like;
		res[2] = p_max_chi;
		res[3] = p_max_like;
		res[4] = p_ht;
		res[5] = p_edist;

		res[ 6] = sum_chi;
		res[ 7] = sum_like;
		res[ 8] = max_chi;
		res[ 9] = max_like;
		res[10] = ht;
		res[11] = edist;

		if (tables_wanted) {
			seq.get_observed_tables(res + 12);
		}

		if (perm_stats_wanted) {
			seq.get_perm_stats(res + 12 + 4 * n * n * tables_wanted);
		}

		UNPROTECT(1);

		// FIXME hmm, I can't call the recommended pthread_exit(). It doesn't seem to be critical, unless joining encountered a problem.

		return (R_res);
	} catch (exception& e) {
		Rprintf(e.what());
		SEXP R_res;
		PROTECT(R_res = allocMatrix(REALSXP, 12, 1));
		UNPROTECT(1);
		return (R_res);
	}
}

}

#endif // R_INTERFACE

#ifdef GWAS_INTERFACE

// Entry point
// ============================================================================

/*
 *
 * This program is invoked on a cluster node, with command line arguments that
 * identify a genotype file, a phenotype file, and a genotype-ranges file that
 * determines windows or genes to be tested as groups [and other parameters
 * that I'm sure I'll add later].
 *
 * Each SNP group is tested for association with the (possibly quantitative,
 * multivariate) phenotype using the Heller-Heller-Gorfine test with a
 * Wald-sequential permutation test so that p-values at a resolution that
 * facilitates an implicit multiplicity correction.
 *
 * Details about inputs:
 * - The genotype file is in transposed plink tped format (output with plink
 *   --recode12 --transpose)
 * - The phenotype file is an ASCII file where lines are individuals and in
 *   each line there is a tab delimited list of phenotypes. There is also a
 *   header line with tab delimited phenotype names.
 *
 */

//#define PHENO_DIST_EXTERNAL
//#define BINARY_PHENO // make sure PHENO_DIST_EXTERNAL is undefined

static void read_tped_12_snps(ifstream& ifs, int cohort_size, int nr_snps_in_range, double* x);
static char* get_strtime(void) { time_t now = time(0); char* strtime = asctime(localtime(&now)); strtime[strlen(strtime) - 1] = '\0'; return strtime; }

int main(int argc, char** argv) {
	if (argc != 16) {
		cerr << "Usage: " << endl;
		cerr << argv[0] << " phenotypes_filename genotypes_filename snp_ranges_filename output_filename_prefix w_sum w_max is_sequential alpha alpha0 beta0 eps nr_perm nr_threads base_seed global_nr_snp_ranges" << endl;
		cerr << endl;
		cerr << "Windows debugging example:" << endl;
		cerr << argv[0] << " \"C:\\Library\\Work with Ruth Heller\\Data\\KORA\\phenotypes_dist.txt\" \"C:\\Library\\Work with Ruth Heller\\Data\\KORA\\genotypes.chr22.tped\" \"C:\\Library\\Work with Ruth Heller\\Data\\KORA\\chr22-gene-ranges\" \"C:\\Library\\Work with Ruth Heller\\Data\\KORA\\hhg.out\" 0 2 1 0.05 0.05 0.01 0.01 100 0 0 0" << endl;
		cerr << endl;
		cerr << "Linux debugging example:" << endl;
		cerr << argv[0] << " dy.txt x.tped ranges.txt out 0 2 1 0.00015 0.05 0.0000015 0.01 1000000 0 0 0" << endl;
		cerr << endl;
		exit(1);
	}

	istringstream buffer;

	string executable = argv[0];

	string phenotypes_filename    = argv[1];
	string genotypes_filename     = argv[2];
	string snp_ranges_filename    = argv[3];
	string output_filename_prefix = argv[4];

	double w_sum; 				buffer.str(argv[ 5]); buffer >> w_sum; buffer.clear();
	double w_max; 				buffer.str(argv[ 6]); buffer >> w_max; buffer.clear();
	bool is_sequential; 		buffer.str(argv[ 7]); buffer >> is_sequential; buffer.clear();
	double alpha; 				buffer.str(argv[ 8]); buffer >> alpha; buffer.clear();
	double alpha0; 				buffer.str(argv[ 9]); buffer >> alpha0; buffer.clear();
	double beta0; 				buffer.str(argv[10]); buffer >> beta0; buffer.clear();
	double eps; 				buffer.str(argv[11]); buffer >> eps; buffer.clear();
	int nr_perm; 				buffer.str(argv[12]); buffer >> nr_perm; buffer.clear();
	int nr_threads; 			buffer.str(argv[13]); buffer >> nr_threads; buffer.clear();
	int base_seed; 				buffer.str(argv[14]); buffer >> base_seed; buffer.clear();
	int global_nr_snp_ranges; 	buffer.str(argv[15]); buffer >> global_nr_snp_ranges; buffer.clear();

	if (nr_threads <= 0) {
		nr_threads = get_available_nr_threads();
	}

	cout << "Gene-GWAS running with the following options:" << endl << endl;
	cout << "phenotypes_filename ...... " << phenotypes_filename << endl;
	cout << "genotypes_filename ....... " << genotypes_filename << endl;
	cout << "snp_ranges_filename ...... " << snp_ranges_filename << endl;
	cout << "output_filename_prefix ... " << output_filename_prefix << endl;
	cout << "w_sum .................... " << w_sum << endl;
	cout << "w_max .................... " << w_max << endl;
	cout << "is_sequential ............ " << is_sequential << endl;
	cout << "alpha .................... " << alpha << endl;
	cout << "alpha0 ................... " << alpha0 << endl;
	cout << "beta0 .................... " << beta0 << endl;
	cout << "eps ...................... " << eps << endl;
	cout << "nr_perm .................. " << nr_perm << endl;
	cout << "nr_threads ............... " << nr_threads << endl;
	cout << "base_seed ................ " << base_seed << endl;
	cout << "global_nr_snp_ranges ..... " << global_nr_snp_ranges << endl << endl;

	string line;

	// Open input files
	ifstream ifs_phenotypes(phenotypes_filename.c_str());
	if (!ifs_phenotypes.is_open()) {
		cerr << "Could not open phenotypes file " << phenotypes_filename << " for reading." << endl;
		exit(1);
	}

	ifstream ifs_snp_ranges(snp_ranges_filename.c_str());
	if (!ifs_snp_ranges.is_open()) {
		cerr << "Could not open groups file " << snp_ranges_filename << " for reading." << endl;
		exit(1);
	}

	ifstream ifs_genotypes(genotypes_filename.c_str());
	if (!ifs_genotypes.is_open()) {
		cerr << "Could not open genotypes file " << genotypes_filename << " for reading." << endl;
		exit(1);
	}

	string output_filename_observed = output_filename_prefix + "-observed.txt";
	ofstream ofs_observed(output_filename_observed.c_str());
	if (!ofs_observed.is_open()) {
		cerr << "Could not open output file " << output_filename_observed << " for writing." << endl;
		exit(1);
	}
	ofs_observed << "sum_chi sum_like max_chi max_like start_idx stop_idx start_pos end_pos range_name" << endl; // header

	string output_filename_pvals = output_filename_prefix + "-pvals.txt";
	ofstream ofs_pvals(output_filename_pvals.c_str());
	if (!ofs_pvals.is_open()) {
		cerr << "Could not open output file " << output_filename_pvals << " for writing." << endl;
		exit(1);
	}
	ofs_pvals << "range nr_perms p_sum_chi p_sum_like p_max_chi p_max_like" << endl; // header

	// Probe input dimensions
	int cohort_size = -1;
	int nr_phenotypes = -1;

#ifdef PHENO_DIST_EXTERNAL
	nr_phenotypes = 1; // will have no effect anyway
#else
	getline(ifs_phenotypes, line);
	++cohort_size;
	buffer.clear();
	buffer.str(line);
	string foo;
	nr_phenotypes = 0;
	while (buffer.good()) {
		buffer >> foo;
		++nr_phenotypes;
	}
#endif

	while (ifs_phenotypes.good()) {
		getline(ifs_phenotypes, line);
		++cohort_size;
#ifdef DEBUG_LIMIT_COHORT_SIZE
		if (cohort_size == 10) {
			break;
		}
#endif
	}
	ifs_phenotypes.clear();
	ifs_phenotypes.seekg(0);

	int nr_snp_ranges = -1;
	while (ifs_snp_ranges.good()) {
		getline(ifs_snp_ranges, line);
		++nr_snp_ranges;
#ifdef DEBUG_LIMIT_NR_RANGES
		if (nr_snp_ranges == 10) {
			break;
		}
#endif
	}
	ifs_snp_ranges.clear();
	ifs_snp_ranges.seekg(0);

#ifndef PHENO_DIST_EXTERNAL
	cout << "Dataset contains " << nr_snp_ranges << " SNP ranges and " << nr_phenotypes << " phenotypes on " << cohort_size << " individuals." << endl;
#else
	cout << "Dataset contains " << nr_snp_ranges << " SNP ranges and an unknown number of phenotypes on " << cohort_size << " individuals." << endl;
#endif

	if (nr_snp_ranges <= 0 || cohort_size <= 0 || nr_phenotypes <= 0) {
		cout << "Nothing to do, finishing." << endl;
		exit(0);
	}

	double* y = new double[cohort_size * nr_phenotypes];
	double* dy = new double[cohort_size * cohort_size];

#ifdef DATAIN_DEBUG_PRINTS
	cout << "Reading phenotypes:" << endl;
#endif

#ifdef PHENO_DIST_EXTERNAL
	// Read precomputed phenotype distances
	read_data_file(ifs_phenotypes, 0, 0, cohort_size, cohort_size, dy);
#else
	// Read the phenotypes
	read_data_file(ifs_phenotypes, 1, 0, cohort_size, nr_phenotypes, y);

	// Compute dy according to a robust Mahalanobis distance
#ifdef BINARY_PHENO
	compute_01_distances(cohort_size, nr_phenotypes, y, dy);
#else
	compute_lp_distances(cohort_size, nr_phenotypes, 2, y, dy);
#endif
#endif

	ifs_phenotypes.close();

	// For each hypothesis (i.e. SNP range)
	// NOTE: SNP ranges are not assumed to be ordered and in any case may overlap

	for (int g = 0; g < nr_snp_ranges; ++g) {
		// Read SNP range
		int snp1, snp2, nr_snps_in_range, start_pos, end_pos;
		string range_name;
		getline(ifs_snp_ranges, line);
		istringstream iss(line);
		iss >> snp1 >> snp2 >> start_pos >> end_pos >> range_name;
		nr_snps_in_range = snp2 - snp1 + 1;

		cout << get_strtime() << ": Working on range " << (g + 1) << " of " << nr_snp_ranges << ": " << range_name << ", SNPs " << snp1 << " to " << snp2 << endl;

		double* x = new double[cohort_size * nr_snps_in_range];
		double* dx = new double[cohort_size * cohort_size];

		// Read SNPs
		ifs_genotypes.clear();
		ifs_genotypes.seekg(0);
		for (int snpi = 0; snpi < snp1; ++snpi) {
			ifs_genotypes.ignore(1024 + cohort_size * 4, '\n');
		}
		read_tped_12_snps(ifs_genotypes, cohort_size, nr_snps_in_range, x);

		// Compute dx according to the L2 norm
		compute_lp_distances(cohort_size, nr_snps_in_range, 2, x, dx);

		// Perform the sequential test on dx and dy
#ifdef BINARY_PHENO
		SequentialTest seq(TWO_SAMPLE_TEST, cohort_size, nr_phenotypes, dx, dy, y, w_sum, w_max, is_sequential, alpha, alpha0, beta0, eps, nr_perm, nr_threads, base_seed, false, false);
#else
		SequentialTest seq(GENERAL_TEST, cohort_size, nr_phenotypes, dx, dy, y, w_sum, w_max, is_sequential, alpha, alpha0, beta0, eps, nr_perm, nr_threads, base_seed, false, false);
#endif

		double sum_chi, sum_like, max_chi, max_like, ht, edist;
		seq.get_observed_stats(sum_chi, sum_like, max_chi, max_like, ht, edist);
		ofs_observed << portable_double2str(sum_chi) << " " << portable_double2str(sum_like) << " "
				<< portable_double2str(max_chi) << " " << portable_double2str(max_like) << " "
				<< snp1 << " " << snp2 << " " << start_pos << " " << end_pos << " " << range_name << endl;

		// Write p-values to output file
		double p_sum_chi, p_sum_like, p_max_chi, p_max_like, p_ht, p_edist;
		seq.get_pvalues(p_sum_chi, p_sum_like, p_max_chi, p_max_like, p_ht, p_edist);
		ofs_pvals << range_name << " " << nr_perm << " " << portable_double2str(p_sum_chi) << " " << portable_double2str(p_sum_like) << " " << portable_double2str(p_max_chi) << " " << portable_double2str(p_max_like) << endl;

		delete[] x;
		delete[] dx;

		if ((global_nr_snp_ranges > 0) && (p_sum_like * global_nr_snp_ranges / (g + 1) > alpha0)) {
			cout << "FDR-like p-value threshold exceeded at nominal " << alpha0 << ", giving up on further ranges." << endl;
			break;
		}
	}

	ifs_genotypes.close();
	ifs_snp_ranges.close();
	ofs_pvals.close();

	delete[] y;
	delete[] dy;

	cout << get_strtime() << ": " << argv[0] << " finished. " << endl;

	return 0;
}

static void read_tped_12_snps(ifstream& ifs, int cohort_size, int nr_snps_in_range, double* x) {
	string line, foo;
	int a1, a2;

	for (int i = 0; i < nr_snps_in_range; ++i) {
		getline(ifs, line);
		istringstream iss(line);
		iss >> foo >> foo >> foo >> foo;
		for (int j = 0; j < cohort_size; ++j) {
			iss >> a1 >> a2;
			x[i * cohort_size + j] = a1 + a2 - 2;
		}
	}
}

#endif // GWAS_INTERFACE

// Misc. utility functions
// ============================================================================

// FIXME this could probably be improved, in particular made more general and accurate

static int get_available_nr_threads(void) {
#ifdef NO_THREADS
	return (1);
#else
#ifdef _WIN32
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);
	return (sysinfo.dwNumberOfProcessors);
#else
	return (sysconf(_SC_NPROCESSORS_ONLN));
#endif
#endif
}

#ifndef R_INTERFACE
static void read_data_file(ifstream& ifs, int skp_first_lines, int skp_first_fields, int n, int m, double* buf) {
	string line;
	double fld;

#ifdef DATAIN_DEBUG_PRINTS
	cout << "Reading data:" << endl;
#endif

	for (int i = 0; i < skp_first_lines; ++i) {
		getline(ifs, line);
	}

	for (int i = 0; i < n; ++i) {
		getline(ifs, line);
		istringstream iss(line);
		for (int j = 0; j < skp_first_fields; ++j) {
			iss >> fld;
		}
		for (int j = 0; j < m; ++j) {
			iss >> fld;
#ifdef DATAIN_DEBUG_PRINTS
			cout << fld << " ";
#endif
			buf[j * n + i] = fld;
		}
#ifdef DATAIN_DEBUG_PRINTS
		cout << endl;
#endif
	}

	ifs.close();
}

static void compute_01_distances(int n, int p, double* v, double *dv) {
	for (int k = 0; k < n; ++k) {
		for (int l = 0; l < n; ++l) {
			bool d = false;
			for (int j = 0; j < p; ++j) {
				d |= (v[j * n + k] != v[j * n + l]);
			}

			dv[l * n + k] = d;
		}
	}
}

static void compute_hamming_distances(int n, int p, double* v, double *dv) {
	for (int k = 0; k < n; ++k) {
		for (int l = 0; l < n; ++l) {
			int d = 0;
			for (int j = 0; j < p; ++j) {
				d += (v[j * n + k] != v[j * n + l]);
			}

			dv[l * n + k] = d;
		}
	}
}

static void compute_lp_distances(int n, int p, int norm_p, double* v, double *dv) {
	double rcp_norm_p = 1.0 / norm_p;
	double d;

	for (int k = 0; k < n; ++k) {
		for (int l = 0; l < n; ++l) {
			d = 0;
			for (int j = 0; j < p; ++j) {
				d += pow(abs(v[j * n + k] - v[j * n + l]), norm_p);
			}

			dv[l * n + k] = pow(d, rcp_norm_p);
		}
	}
}

// Robust Mahalanobis distance as suggested by Rosenbaum
static void compute_rm_distances(int n, int p, double* v, double *dv) {
	compute_lp_distances(n, p, 2, v, dv);

	// Well, since this is kind of painful to do in C, and not a computational bottleneck
	// I think I'll just do it in R/Matlab and load the distance matrix instead.

	// TODO:
	// for every column of v, compute ranks of the elements => r
	// compute covariance matrix of r => cv
	// compute variance of [1, 2, ..., n] => vuntied
	// rat <- sqrt(vuntied/diag(cv))
    // cv <- diag(rat) %*% cv %*% diag(rat)
    // out <- matrix(NA, m, n - m)
    // Xc <- X[z == 0, ]
    // Xt <- X[z == 1, ]
    // rownames(out) <- rownames(X)[z == 1]
    // colnames(out) <- rownames(X)[z == 0]
    // icov <- ginv(cv)
    // for (i in 1:m) {
	// 	out[i, ] <- mahalanobis(Xc, Xt[i, ], icov, inverted = T)
	// }
}

static string portable_double2str(double x) {
	if (isnan(x)) {
		return "nan";
	} else if (isinf(x)) {
		return "inf";
	} else {
		ostringstream oss;
		oss << x;
		return oss.str();
	}

	return "";
}
#endif
