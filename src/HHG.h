#ifndef HHG_H_
#define HHG_H_

#include <vector>

// Define exactly one of these (I sometimes define them in project files, makefiles, or build commands)
//#define DEBUG_INTERFACE
#define R_INTERFACE
//#define GWAS_INTERFACE

//#define NO_THREADS
//#define DEBUG_THREADS
//#define DEBUG_CHECKS
//#define DEBUG_PRINTS // NOTE: don't use this with R_INTERFACE, and you probably should define NO_THREADS (I don't print to the R console, and I don't synchronize printing)
//#define ST_DEBUG_PRINTS
//#define DATAIN_DEBUG_PRINTS
//#define DEBUG_LIMIT_COHORT_SIZE
//#define DEBUG_LIMIT_NR_RANGES

#ifdef R_INTERFACE
#undef ERROR
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#endif

typedef enum {
	TWO_SAMPLE_TEST 		= 0,
	K_SAMPLE_TEST 			= 1,
	NO_TIES_TEST 			= 2,
	GENERAL_TEST 			= 3,
	UDF_SPR_OBS      		= 4,
	UDF_SPR_ALL      		= 5,
	UDF_PPR_22_OBS   		= 6,
	UDF_PPR_22_ALL   		= 7,
	UDF_PPR_33_OBS   		= 8,
	UDF_PPR_33_ALL   		= 9,
	UDF_TPR_OBS      		= 10,
	UDF_TPR_ALL      		= 11,
	UDF_SPPR_OBS     		= 12,
	UDF_SPPR_ALL     		= 13,
	UDF_DDP_OBS     		= 14,
	UDF_DDP_ALL     		= 15,
	CI_UVZ_NN  				= 16,
	CI_UVZ_GAUSSIAN			= 17,
	CI_MVZ_NN  				= 18,
	CI_MVZ_GAUSSIAN			= 19,
	CI_UDF_ADP_MVZ_NN		= 20,
	CI_MVZ_NN_GRID_BW		= 21,
	K_SAMPLE_DDP_M2			= 22,
	K_SAMPLE_DDP_M3			= 23,
	K_SAMPLE_DDP			= 24,
	GOF_DDP_M2		        = 25,
	GOF_DDP_M3		        = 26,
	GOF_DDP			        = 27,
	K_SAMPLE_EXISTING		= 100,
	GOF_EXISTING		    = 101
} TestType;

#define IS_UDF_TEST(tt) ((tt) == UDF_SPR_OBS       || (tt) == UDF_SPR_ALL     || (tt) == UDF_PPR_22_OBS  || \
		                 (tt) == UDF_PPR_22_ALL    || (tt) == UDF_PPR_33_OBS  || (tt) == UDF_PPR_33_ALL  || \
		                 (tt) == UDF_TPR_OBS       || (tt) == UDF_TPR_ALL     || (tt) == UDF_SPPR_OBS    || \
		                 (tt) == UDF_SPPR_ALL      || (tt) == UDF_DDP_OBS     || (tt) == UDF_DDP_ALL     || \
		                 (tt) == CI_UDF_ADP_MVZ_NN)

#define IS_CI_MVZ_TEST(tt) ((tt) == CI_MVZ_NN         || (tt) == CI_MVZ_GAUSSIAN || \
		                    (tt) == CI_UDF_ADP_MVZ_NN || (tt) == CI_MVZ_NN_GRID_BW)

#define IS_K_SAMPLE_DDP(tt) ((tt) == K_SAMPLE_DDP_M2 || (tt) == K_SAMPLE_DDP_M3 || (tt) == K_SAMPLE_DDP)

#define IS_GOF_DDP(tt) ((tt) == GOF_DDP_M2 || (tt) == GOF_DDP_M3 || (tt) == GOF_DDP)

//#define UDF_ALLOW_DEGENERATE_PARTITIONS
#define UDF_NORMALIZE

typedef std::vector< std::vector<double> > matrix;
typedef std::pair<double, int> dbl_int_pair;
typedef std::vector<dbl_int_pair> dbl_int_pair_vector;
typedef std::vector< std::vector<dbl_int_pair> > dbl_int_pair_matrix;

struct ExtraParams {
	double w_sum;
	double w_max;
	int K; // number of unique y values in K-sample test, also used for HHGCI NN kernel width, and for order of ADP/DDP partitions
	int M; // order of ADP/DDP partition for the K-sample ADP/DDP test
	int* y_counts; // counts observed for each unique y value, sorted by y value
	bool correct_mi_bias;
	double sig;
	int nnh; // NN kernel width used for our statistic in the CI test
	int nnh_lsb; // NN kernel width used for locally smoothed bootstrap (for computing p-values in the CI test)
	int nnh_grid_cnt;
	int* nnh_grid;

	ExtraParams() {
		w_sum = 0;
		w_max = 0;
		K = 0;
		M = 0;
		y_counts = NULL;
		correct_mi_bias = false;
		sig = 0;
		nnh = 0;
		nnh_lsb = 0;
		nnh_grid_cnt = 0;
		nnh_grid = NULL;
	}
};

#endif /* HHG_H_ */
