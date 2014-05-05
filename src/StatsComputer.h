/*
 * StatsComputer.h
 *
 */

#ifndef STATSCOMPUTER_H_
#define STATSCOMPUTER_H_

#include "HHG.h"

class StatsComputer {
public:
	StatsComputer(TestType tt, int xy_nrow, int y_ncol, double* dx, double* dy, double* y,
			dbl_int_pair_matrix* sorted_dx, dbl_int_pair_matrix* sorted_dy, dbl_int_pair_matrix* sorted_dz,
			ExtraParams* extra_params);
	virtual ~StatsComputer();

	void compute(void);
	void compute_and_store_tbls(int* tbls);
	void permute_and_compute(void);
	void get_stats(double& sum_chi, double& sum_like, double& max_chi, double& max_like, double& ht, double& edist);
	void get_tables(double* tbl);
	void get_grid_stats(double *sum_chi_grid, double *sum_like_grid, double *max_chi_grid, double *max_like_grid);

protected:
	void resample_univariate(void);
	void resample_multivariate(void);
	void resample_uvz_ci(void);
	void resample_mvz_ci(void);
	void resample_dummy(void);

	void sort_xy_distances_per_row(void);

	void hhg_two_sample(void);
	void hhg_k_sample(void);
	void hhg_no_ties(void);
	void hhg_general(void);
	void hhg_udf_spr_obs(void);
	void hhg_udf_spr_all(void);
	void hhg_udf_ppr_22_obs(void);
	void hhg_udf_ppr_22_all(void);
	void hhg_udf_ppr_33_obs(void);
	void hhg_udf_ppr_33_all(void);
	void hhg_udf_tpr_obs(void);
	void hhg_udf_tpr_all(void);
	void hhg_udf_sppr_obs(void);
	void hhg_udf_sppr_all(void);
	void hhg_udf_ddp(void);
	void hhg_udf_adp(void);
	void hhg_ci_uvz_nn(void);
	void hhg_ci_uvz_gaussian(void);
	void hhg_ci_mvz_nn(void);
	void hhg_ci_mvz_gaussian(void);
	void hhg_ci_udf_adp_mvz_nn(void);
	void hhg_ci_mvz_nn_grid(void);
	void hhg_k_sample_existing(void);
	void hhg_k_sample_ddp_m2(void);
	void hhg_k_sample_ddp_m3(void);
	void hhg_k_sample_ddp(void);
	void hhg_gof_existing(void);
	void hhg_gof_ddp_m2(void);
	void hhg_gof_ddp_m3(void);
	void hhg_gof_ddp(void);

	void other_stats_two_sample(void);
	void other_stats_k_sample(void);
	void other_stats_general(void);
	void other_stats_univar_dist_free(void);
	void other_stats_ci(void);

	void accumulate_2x2_contingency_table(double a00, double a01, double a10, double a11, double nrmlz, double reps);
	void accumulate_local_stats(double chi, double like, double emin);
	void hhg_gen_inversions(int *permutation, int *source, int *inversion_count, int dim);
	void hhg_gen_merge(int *permutation, int *source, int *inversion_count, int dim);
	void compute_single_integral(int n, double* xx, double* yy);
	void compute_double_integral(int n, double* xx, double* yy);
	int count_sample_points_in_rect(int xl, int xh, int yl, int yh);
	double count_ddp_with_given_cell(int xl, int xh, int yl, int yh);
	double count_adp_with_given_cell(int xl, int xh, int yl, int yh);
	void precompute_adp(void);
	void precompute_adp_k_sample(void);
	double my_choose(int n, int k);
	int my_rand(int lo, int hi);

	void compute_ht(void);
	void compute_edist(void);

	void compute_spr_obs(int xi, int yi, int n, int pn, int nm1, double nm1d);
	void compute_spr_all(int xi, int yi, int n, int pn, double nd);
	void compute_ppr_22(int xr_lo, int xr_hi, int yr_lo, int yr_hi, int pn, int nm2, double nm2s);
	void compute_ppr_33(int xr_lo, int xr_hi, int yr_lo, int yr_hi, int n, int pn, double nm2);
	void compute_tpr(int xl, int xm, int xh, int yl, int ym, int yh, int n, int pn, double nm3);

	TestType tt;
	int xy_nrow;
	int y_ncol;
	double* dx;
	double* dy;
	double* y;
	double* z;
	dbl_int_pair_matrix *sorted_dx, *sorted_dy, *sorted_dz;
	int K; // k-sample test: number of unique y values; DDP/ADP: partition order
	int M; // in the k-sample ADP/DDP test, this is the partition order and K is the number of unique y values
	int* y_counts; // counts observed for each unique y value, sorted by y value
	double w_sum;
	double w_max;
	double min_w;
	double *x_ordered_by_y, *y_ordered_by_x;
	double sig;
	int nnh; // NN kernel width used for our statistic in the CI test
	int nnh_lsb; // NN kernel width used for locally smoothed bootstrap (for computing p-values in the CI test)
	int nnh_grid_cnt;
	int* nnh_grid;

	void (StatsComputer::*stats_func)(void);
	void (StatsComputer::*stats_func2)(void);
	void (StatsComputer::*perm_y_func)(void);

    double sum_chi, sum_like, max_chi, max_like, ht, edist;
    double *sum_chi_grid, *sum_like_grid, *max_chi_grid, *max_like_grid;

	int *y0_idx, *y1_idx; // indices of samples with y_i == 0 and 1 respectively
	int *idx_1_to_n;
	int *idx_perm, *idx_perm_inv;
	int *hhg_gen_inversion_count, *hhg_gen_source, *hhg_gen_xy_perm, *hhg_gen_xy_perm_temp, *hhg_gen_y_rev;
	int *hhg_gen_left_buffer, *hhg_gen_right_buffer, *hhg_gen_left_source_buffer, *hhg_gen_right_source_buffer;
	bool should_randomize;

	double kahan_c_chi, kahan_c_like;
	int* double_integral;
	int dintegral_zero_based_idxs;
	int dintegral_pn;
	double* adp, *adp_l, *adp_r;
	int ng_chi, ng_like;
	bool correct_bias;
	dbl_int_pair_vector nn_sorted_x, nn_sorted_y;
	double* null_dist;

	struct dbl_dbl_int {
		double x;
		double y;
		int i;
	};

	typedef std::vector< std::vector<dbl_dbl_int> > dbl_dbl_int_matrix;

	static inline bool dbl_int_pair_comparator(const dbl_int_pair& l, const dbl_int_pair& r) {
		return l.first < r.first;
	}

	static inline bool dbl_dbl_int_pair_comparator_xy(const dbl_dbl_int& l, const dbl_dbl_int& r) {
		return ((l.x < r.x) || ((l.x == r.x) && (l.y > r.y)));
	}

	dbl_dbl_int_matrix sorted_dx_gen;
	int* tbls;
	bool store_tables;
};

#endif /* STATSCOMPUTER_H_ */
