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

// Please note: some test variants are still work in progress; only tests with exposed R interface are working and usable.
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
	CI_NN     				= 16
} TestType;

#define IS_UDF_TEST(tt) ((tt) == UDF_SPR_OBS    || (tt) == UDF_SPR_ALL    || (tt) == UDF_PPR_22_OBS || \
		                 (tt) == UDF_PPR_22_ALL || (tt) == UDF_PPR_33_OBS || (tt) == UDF_PPR_33_ALL || \
		                 (tt) == UDF_TPR_OBS    || (tt) == UDF_TPR_ALL    || (tt) == UDF_SPPR_OBS   || \
		                 (tt) == UDF_SPPR_ALL   || (tt) == UDF_DDP_OBS    || (tt) == UDF_DDP_ALL)

//#define UDF_ALLOW_DEGENERATE_PARTITIONS
#define UDF_NORMALIZE

typedef std::vector< std::vector<double> > matrix;
typedef std::pair<double, int> dbl_int_pair;
typedef std::vector<dbl_int_pair> dbl_int_pair_vector;
typedef std::vector< std::vector<dbl_int_pair> > dbl_int_pair_matrix;

struct ExtraParams {
	double w_sum;
	double w_max;
	int K; // number of unique y values, also used for NN kernel width, and for order of ADP/DDP partitions
	int* y_counts; // counts observed for each unique y value, sorted by y value
	bool correct_mi_bias;
};

#endif /* HHG_H_ */
