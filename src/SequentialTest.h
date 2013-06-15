/*
 * SequentialTest.h
 *
 * A sequential testing procedure.
 *
 * Currently using Wald's sequential procedure, on the a specialized HHG test,
 * but maybe I'll derive and override in the future with other procedures.
 *
 */

#ifndef SEQUENTIALTEST_H_
#define SEQUENTIALTEST_H_

#include "StatsComputer.h"

struct Compute_permutations_thread_arg; // see below

class SequentialTest {
public:
	SequentialTest(TestType tt, int xy_nrow, int y_ncol, double* dx, double* dy, double* y,
			double w_sum, double w_max, double* extra_params,
			bool is_sequential, double alpha, double alpha0, double beta0, double eps, int nr_perm,
			int nr_threads, int base_seed, bool tables_wanted, bool perm_stats_wanted);
	virtual ~SequentialTest();

	bool is_null_rejected(void);
	void get_pvalues(double& p_sum_chi, double& p_sum_like, double& p_max_chi, double& p_max_like, double& p_ht, double& p_edist);
	void get_observed_stats(double& sum_chi, double& sum_like, double& max_chi, double& max_like, double& ht, double& edist);
	void get_observed_tables(double* tbls); // returning double for passing back to R in a single structure (could save some memory if smarter about this)
	void get_perm_stats(double* ps);

protected:
	void run();
	bool update_sequential(int statistic_idx, bool is_null_more_extreme);
	bool update_sequential_all(double perm_sum_chi, double perm_sum_like, double perm_max_chi, double perm_max_like, double perm_ht, double perm_edist);

	TestType tt;

	int nr_perm;
	int nr_threads;
	int base_seed;

	double alpha;
	double alpha0;
	double beta0;
	double eps;

	bool tables_wanted;
	bool perm_stats_wanted;

	double lA, lB;
	double exp1, exp2;

	int nr_statistics;
	double* llr;
	int* pvalc;
	bool* stopped_high;
	bool* stopped_low;
	int* perm_counter;

    double obs_sum_chi, obs_sum_like, obs_max_chi, obs_max_like, obs_ht, obs_edist;
    int* obs_tbls;
    double *perm_sum_chi_v, *perm_sum_like_v, *perm_max_chi_v, *perm_max_like_v, *perm_ht_v, *perm_edist_v;
    int orig_nr_perm;

protected:
	double *dx, *dy, *y;
	int xy_nrow, y_ncol;
	bool is_sequential;

	static inline bool dbl_int_pair_comparator(const dbl_int_pair& l, const dbl_int_pair& r) { return l.first < r.first; }
	static inline bool int_comparator(const int& l, const int& r) { return l < r; }

	dbl_int_pair_matrix sorted_dx, sorted_dy;

	void count_unique_y(void);
	void compute_distances(void);
	void sort_x_distances_per_row(void);
	void sort_y_distances_per_row(void);

	StatsComputer** scs;
	volatile bool stop_all_flag;
	pthread_mutex_t mutex;

	ExtraParams extra_params;

public: // this shouldn't be public but is made so due to pthread/c++ limitations
	void compute_one_permutation(int t);
	void compute_permutations(Compute_permutations_thread_arg* carg);
};

// Helper object for use with threading
struct Compute_permutations_thread_arg {
	SequentialTest* seq;
	int t; // thread number
	bool done_flag;

	Compute_permutations_thread_arg(SequentialTest* seq, int t) : seq(seq), t(t), done_flag(false) {}
};

#endif /* SEQUENTIALTEST_H_ */
