/*
 * StatsComputer.cpp
 *
 */

#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <ctime>
#include <cstring>

#ifdef WIN32
#include <stddef.h>
#endif

#include "StatsComputer.h"

using namespace std;

StatsComputer::StatsComputer(TestType tt, int xy_nrow, int y_ncol, double* dx, double* dy, double* y,
		dbl_int_pair_matrix* sorted_dx, dbl_int_pair_matrix* sorted_dy, dbl_int_pair_matrix* sorted_dz,
		ExtraParams* extra_params)
{
	this->tt = tt;
	this->xy_nrow = xy_nrow;
	this->y_ncol = y_ncol;

	// I will assume it is safe to use these constructor arguments for the
	// lifetime of the object without creating a copy
	this->dx = dx;
	this->dy = dy;
	this->z = y; // in the case of some CI tests, y is actually z or dz

	if (tt == TWO_SAMPLE_TEST || tt == K_SAMPLE_TEST || IS_UDF_TEST(tt) || IS_K_SAMPLE_DDP(tt) || tt == K_SAMPLE_EXISTING) {
		// A copy of y is created which will be permuted
		this->y = new double[xy_nrow * y_ncol];
		memcpy(this->y, y, sizeof(double) * xy_nrow * y_ncol);
	} else {
		this->y = NULL;
	}

	this->sorted_dx = sorted_dx;
	this->sorted_dy = sorted_dy;
	this->sorted_dz = sorted_dz;

	// FIXME now it starts making sense to inherit from extra_params
	w_sum = extra_params->w_sum;
	w_max = extra_params->w_max;
	min_w = min(w_sum, w_max);

	tbls = NULL;
	store_tables = false;
	nnh_grid_cnt = 0;

	if (tt == TWO_SAMPLE_TEST || tt == K_SAMPLE_TEST) {
		K = extra_params->K;
		y_counts = extra_params->y_counts;
	} else if (tt == UDF_DDP_OBS || tt == UDF_DDP_ALL) {
		x_ordered_by_y = dy;
		y_ordered_by_x = dy + xy_nrow;
		K = extra_params->K;
		correct_bias = extra_params->correct_mi_bias;
	} else if (tt == CI_UVZ_NN || tt == CI_MVZ_NN) {
		nnh = extra_params->nnh;
		nnh_lsb = extra_params->nnh_lsb;
	} else if (tt == CI_UVZ_GAUSSIAN || tt == CI_MVZ_GAUSSIAN) {
		sig = extra_params->sig;
		nnh_lsb = extra_params->nnh_lsb;
	} else if (tt == CI_UDF_ADP_MVZ_NN) {
		K = extra_params->K;
		correct_bias = extra_params->correct_mi_bias;
		nnh = extra_params->nnh;
		nnh_lsb = extra_params->nnh_lsb;
	} else if (tt == CI_MVZ_NN_GRID_BW) {
		nnh_grid_cnt = extra_params->nnh_grid_cnt;
		nnh_grid = extra_params->nnh_grid; // NOTE: not copying, assuming reference will persist
		nnh_lsb = extra_params->nnh_lsb;
		sum_chi_grid  = new double[nnh_grid_cnt];
		sum_like_grid = new double[nnh_grid_cnt];
		max_chi_grid  = new double[nnh_grid_cnt];
		max_like_grid = new double[nnh_grid_cnt];
	} else if (IS_K_SAMPLE_DDP(tt)) {
		K = extra_params->K; // number of sub-samples (as in "K-sample problem")
		M = extra_params->M; // order of partition
		y_counts = extra_params->y_counts;
	} else if (tt == K_SAMPLE_EXISTING) {
		K = extra_params->K; // number of sub-samples (as in "K-sample problem")
		M = 2; // order of partition (all classical tests use dichotomies)
		y_counts = extra_params->y_counts;
	} else if (IS_GOF_DDP(tt)) {
		K = 1;
		M = extra_params->M;
	} else if (tt == GOF_EXISTING) {
		K = 1;
		M = 2;
	}

	// Allocate temporary buffers
	if (tt == NO_TIES_TEST || tt == GENERAL_TEST || tt == CI_UVZ_GAUSSIAN || tt == CI_MVZ_GAUSSIAN) {
		idx_1_to_n = new int[xy_nrow];
		idx_perm = new int[xy_nrow];
		idx_perm_inv = new int[xy_nrow];

		for (int i = 0; i < xy_nrow; ++i) {
			idx_perm[i] = idx_perm_inv[i] = idx_1_to_n[i] = i;
		}
	} else if (tt == CI_UVZ_NN) {
		idx_perm = new int[xy_nrow];
		idx_perm_inv = new int[xy_nrow];

		for (int i = 0; i < xy_nrow; ++i) {
			idx_perm[i] = idx_perm_inv[i] = i;
		}
	} else if (tt == CI_MVZ_NN || tt == CI_UDF_ADP_MVZ_NN || tt == CI_MVZ_NN_GRID_BW) {
		idx_perm = new int[xy_nrow];
		idx_perm_inv = new int[xy_nrow];

		for (int i = 0; i < xy_nrow; ++i) {
			idx_perm[i] = idx_perm_inv[i] = i;
		}
	}

	if (tt == NO_TIES_TEST || tt == GENERAL_TEST) {
		hhg_gen_inversion_count = new int[xy_nrow];
		hhg_gen_source = new int[xy_nrow];
		hhg_gen_xy_perm = new int[xy_nrow];
		hhg_gen_xy_perm_temp = new int[xy_nrow];
		hhg_gen_y_rev = new int[xy_nrow];
		hhg_gen_left_buffer = new int[xy_nrow / 2];
		hhg_gen_right_buffer = new int[xy_nrow / 2 + (xy_nrow & 1)];
		hhg_gen_left_source_buffer = new int[xy_nrow / 2];
		hhg_gen_right_source_buffer = new int[xy_nrow / 2 + (xy_nrow & 1)];
	}

	if (IS_K_SAMPLE_DDP(tt) || tt == K_SAMPLE_EXISTING) {
		// OK, so in this case it is not a *double* integral; it is K single integrals
		// the last row will hold the overall x integral (ignoring y)
		dintegral_pn = xy_nrow + 1;
		double_integral = new int[(K + 1) * dintegral_pn];
	} else if (IS_UDF_TEST(tt)) {
		if (tt == CI_UDF_ADP_MVZ_NN) {
			dintegral_pn = nnh + 2;
			nn_sorted_x.resize(nnh);
			nn_sorted_y.resize(nnh);
		} else {
			dintegral_pn = xy_nrow + 2;
			dintegral_zero_based_idxs = !(tt == UDF_DDP_OBS || tt == UDF_DDP_ALL);
		}
		double_integral = new int[dintegral_pn * dintegral_pn];
	} else if (IS_GOF_DDP(tt) || tt == GOF_EXISTING) {
		null_dist = new double[xy_nrow + 1]; // FIXME can just do this padding in R and then we won't have to create a copy
		null_dist[0] = 0;
		memcpy(null_dist + 1, y, sizeof(double) * xy_nrow);
		null_dist[xy_nrow] = 1;
	}

	if (tt == UDF_DDP_ALL || tt == CI_UDF_ADP_MVZ_NN) {
		adp = new double[xy_nrow];
		adp_l = new double[xy_nrow];
		adp_r = new double[xy_nrow];
		precompute_adp();
	} else if (IS_K_SAMPLE_DDP(tt) || IS_GOF_DDP(tt)) {
		adp = new double[xy_nrow];
		adp_l = new double[xy_nrow];
		precompute_adp_k_sample();
	}

	y0_idx = y1_idx = NULL;
	should_randomize = false;

	switch (tt) {
		case TWO_SAMPLE_TEST:
			assert(y_ncol == 1 && K == 2);

			y0_idx = new int[y_counts[0]];
			y1_idx = new int[y_counts[1]];

			stats_func  = &StatsComputer::hhg_two_sample;
			stats_func2 = &StatsComputer::other_stats_two_sample;
			perm_y_func = &StatsComputer::resample_univariate;
		break;

		case K_SAMPLE_TEST:
			stats_func  = &StatsComputer::hhg_k_sample;
			stats_func2 = &StatsComputer::other_stats_k_sample;
			perm_y_func = &StatsComputer::resample_univariate;
		break;

		case NO_TIES_TEST:
			stats_func  = &StatsComputer::hhg_no_ties;
			stats_func2 = &StatsComputer::other_stats_general;
			perm_y_func = &StatsComputer::resample_multivariate; // actually y can be univariate here too and I could optimize for such a case
		break;

		case GENERAL_TEST:
			stats_func  = &StatsComputer::hhg_general;
			stats_func2 = &StatsComputer::other_stats_general;
			perm_y_func = &StatsComputer::resample_multivariate; // actually y can be univariate here too and I could optimize for such a case

			sorted_dx_gen.resize(xy_nrow);
			for (int k = 0; k < xy_nrow; ++k) {
				sorted_dx_gen[k].resize(xy_nrow);
			}

			sort_xy_distances_per_row();
		break;

		case UDF_SPR_OBS:
			stats_func  = &StatsComputer::hhg_udf_spr_obs;
			stats_func2 = &StatsComputer::other_stats_univar_dist_free;
			perm_y_func = &StatsComputer::resample_univariate;
		break;

		case UDF_SPR_ALL:
			stats_func  = &StatsComputer::hhg_udf_spr_all;
			stats_func2 = &StatsComputer::other_stats_univar_dist_free;
			perm_y_func = &StatsComputer::resample_univariate;
		break;

		case UDF_PPR_22_OBS:
			stats_func  = &StatsComputer::hhg_udf_ppr_22_obs;
			stats_func2 = &StatsComputer::other_stats_univar_dist_free;
			perm_y_func = &StatsComputer::resample_univariate;
		break;

		case UDF_PPR_22_ALL:
			stats_func  = &StatsComputer::hhg_udf_ppr_22_all;
			stats_func2 = &StatsComputer::other_stats_univar_dist_free;
			perm_y_func = &StatsComputer::resample_univariate;
		break;

		case UDF_PPR_33_OBS:
			stats_func  = &StatsComputer::hhg_udf_ppr_33_obs;
			stats_func2 = &StatsComputer::other_stats_univar_dist_free;
			perm_y_func = &StatsComputer::resample_univariate;
		break;

		case UDF_PPR_33_ALL:
			stats_func  = &StatsComputer::hhg_udf_ppr_33_all;
			stats_func2 = &StatsComputer::other_stats_univar_dist_free;
			perm_y_func = &StatsComputer::resample_univariate;
		break;

		case UDF_TPR_OBS:
			stats_func  = &StatsComputer::hhg_udf_tpr_obs;
			stats_func2 = &StatsComputer::other_stats_univar_dist_free;
			perm_y_func = &StatsComputer::resample_univariate;
		break;

		case UDF_TPR_ALL:
			stats_func  = &StatsComputer::hhg_udf_tpr_all;
			stats_func2 = &StatsComputer::other_stats_univar_dist_free;
			perm_y_func = &StatsComputer::resample_univariate;
		break;

		case UDF_SPPR_OBS:
			stats_func  = &StatsComputer::hhg_udf_sppr_obs;
			stats_func2 = &StatsComputer::other_stats_univar_dist_free;
			perm_y_func = &StatsComputer::resample_univariate;
		break;

		case UDF_SPPR_ALL:
			stats_func  = &StatsComputer::hhg_udf_sppr_all;
			stats_func2 = &StatsComputer::other_stats_univar_dist_free;
			perm_y_func = &StatsComputer::resample_univariate;
		break;

		case UDF_DDP_OBS:
			stats_func  = &StatsComputer::hhg_udf_ddp;
			stats_func2 = &StatsComputer::other_stats_univar_dist_free;
			perm_y_func = &StatsComputer::resample_univariate;
		break;

		case UDF_DDP_ALL:
			stats_func  = &StatsComputer::hhg_udf_adp;
			stats_func2 = &StatsComputer::other_stats_univar_dist_free;
			perm_y_func = &StatsComputer::resample_univariate;
		break;

		case CI_UVZ_NN:
			stats_func  = &StatsComputer::hhg_ci_uvz_nn;
			stats_func2 = &StatsComputer::other_stats_ci;
			perm_y_func = &StatsComputer::resample_uvz_ci;
		break;

		case CI_UVZ_GAUSSIAN:
			stats_func  = &StatsComputer::hhg_ci_uvz_gaussian;
			stats_func2 = &StatsComputer::other_stats_ci;
			perm_y_func = &StatsComputer::resample_dummy;
		break;

		case CI_MVZ_NN:
			stats_func  = &StatsComputer::hhg_ci_mvz_nn;
			stats_func2 = &StatsComputer::other_stats_ci;
			perm_y_func = &StatsComputer::resample_mvz_ci;
		break;

		case CI_MVZ_GAUSSIAN:
			stats_func  = &StatsComputer::hhg_ci_mvz_gaussian;
			stats_func2 = &StatsComputer::other_stats_ci;
			perm_y_func = &StatsComputer::resample_dummy;
		break;

		case CI_UDF_ADP_MVZ_NN:
			stats_func  = &StatsComputer::hhg_ci_udf_adp_mvz_nn;
			stats_func2 = &StatsComputer::other_stats_ci;
			perm_y_func = &StatsComputer::resample_mvz_ci;
		break;

		case CI_MVZ_NN_GRID_BW:
			stats_func  = &StatsComputer::hhg_ci_mvz_nn_grid;
			stats_func2 = &StatsComputer::other_stats_ci;
			perm_y_func = &StatsComputer::resample_mvz_ci;
		break;

		case K_SAMPLE_DDP_M2:
			stats_func  = &StatsComputer::hhg_k_sample_ddp_m2;
			stats_func2 = &StatsComputer::other_stats_univar_dist_free;
			perm_y_func = &StatsComputer::resample_univariate;
		break;

		case K_SAMPLE_DDP_M3:
			stats_func  = &StatsComputer::hhg_k_sample_ddp_m3;
			stats_func2 = &StatsComputer::other_stats_univar_dist_free;
			perm_y_func = &StatsComputer::resample_univariate;
		break;

		case K_SAMPLE_DDP:
			stats_func  = &StatsComputer::hhg_k_sample_ddp;
			stats_func2 = &StatsComputer::other_stats_univar_dist_free;
			perm_y_func = &StatsComputer::resample_univariate;
		break;

		case GOF_DDP_M2:
			stats_func  = &StatsComputer::hhg_gof_ddp_m2;
			stats_func2 = &StatsComputer::other_stats_univar_dist_free;
			perm_y_func = &StatsComputer::resample_univariate;
		break;

		case GOF_DDP_M3:
			stats_func  = &StatsComputer::hhg_gof_ddp_m3;
			stats_func2 = &StatsComputer::other_stats_univar_dist_free;
			perm_y_func = &StatsComputer::resample_univariate;
		break;

		case GOF_DDP:
			stats_func  = &StatsComputer::hhg_gof_ddp;
			stats_func2 = &StatsComputer::other_stats_univar_dist_free;
			perm_y_func = &StatsComputer::resample_univariate;
		break;

		case GOF_EXISTING:
			stats_func  = &StatsComputer::hhg_gof_existing;
			stats_func2 = &StatsComputer::other_stats_univar_dist_free;
			perm_y_func = &StatsComputer::resample_univariate;
		break;

		case K_SAMPLE_EXISTING:
			stats_func  = &StatsComputer::hhg_k_sample_existing;
			stats_func2 = &StatsComputer::other_stats_univar_dist_free;
			perm_y_func = &StatsComputer::resample_univariate;
		break;

		default:
#ifdef DEBUG_CHECKS
			cerr << "Unexpected test type specified" << endl;
			exit(1);
#endif
		break;
	}

    sum_chi = sum_like = max_chi = max_like = ht = edist = 0;
	kahan_c_chi = kahan_c_like = 0;
}

StatsComputer::~StatsComputer() {
	if (y != NULL) {
		delete[] y;
	}

	if (y0_idx != NULL) {
		delete[] y0_idx;
		delete[] y1_idx;
	}

	if (tt == NO_TIES_TEST || tt == GENERAL_TEST || tt == CI_UVZ_GAUSSIAN || tt == CI_MVZ_GAUSSIAN) {
		delete[] idx_1_to_n;
		delete[] idx_perm;
		delete[] idx_perm_inv;
	} else if (tt == CI_UVZ_NN || tt == CI_MVZ_NN || tt == CI_UDF_ADP_MVZ_NN || tt == CI_MVZ_NN_GRID_BW) {
		delete[] idx_perm;
		delete[] idx_perm_inv;
	}

	if (tt == NO_TIES_TEST || tt == GENERAL_TEST) {
		delete[] hhg_gen_inversion_count;
		delete[] hhg_gen_source;
		delete[] hhg_gen_xy_perm;
		delete[] hhg_gen_xy_perm_temp;
		delete[] hhg_gen_y_rev;
		delete[] hhg_gen_left_buffer;
		delete[] hhg_gen_right_buffer;
		delete[] hhg_gen_left_source_buffer;
		delete[] hhg_gen_right_source_buffer;
	}

	if (IS_UDF_TEST(tt) || IS_K_SAMPLE_DDP(tt) || tt == K_SAMPLE_EXISTING) {
		delete[] double_integral;
	}
  
	if (tt == UDF_DDP_ALL || tt == CI_UDF_ADP_MVZ_NN) {
		delete[] adp;
		delete[] adp_l;
		delete[] adp_r;
	} else if (IS_K_SAMPLE_DDP(tt) || IS_GOF_DDP(tt)) {
		delete[] adp;
		delete[] adp_l;
	}

	if (IS_GOF_DDP(tt) || tt == GOF_EXISTING) {
		delete[] null_dist;
	}

	if (tt == CI_MVZ_NN_GRID_BW) {
		delete[] sum_chi_grid;
		delete[] sum_like_grid;
		delete[] max_chi_grid;
		delete[] max_like_grid;
	}
}

void StatsComputer::compute(void) {
	(this->*stats_func)();
	(this->*stats_func2)();
}

void StatsComputer::compute_and_store_tbls(int* tbls) {
	store_tables = true;
	this->tbls = tbls;
	compute();
	store_tables = false;
	this->tbls = NULL;
}

void StatsComputer::permute_and_compute(void) {
	// permute the y's (I'm doing it on top of the previous permutation, I don't see a problem with that)
	(this->*perm_y_func)();

	if (tt == GENERAL_TEST) {
		// re-sort first according to x ascending, then according to y descending
		sort_xy_distances_per_row();
	}

	// compute statistic
	compute();
}

void StatsComputer::get_stats(double& sum_chi, double& sum_like, double& max_chi, double& max_like, double& ht, double& edist) {
	sum_chi 	= this->sum_chi;
	sum_like 	= this->sum_like;
	max_chi 	= this->max_chi;
	max_like 	= this->max_like;
	ht 			= this->ht;
	edist 		= this->edist;
}

void StatsComputer::get_grid_stats(double *sum_chi_grid, double *sum_like_grid, double *max_chi_grid, double *max_like_grid) {
	// NOTE: this has to be safe to call even if not a grid test
	for (int i = 0; i < nnh_grid_cnt; ++i) {
		sum_chi_grid [i] = this->sum_chi_grid [i];
		sum_like_grid[i] = this->sum_like_grid[i];
		max_chi_grid [i] = this->max_chi_grid [i];
		max_like_grid[i] = this->max_like_grid[i];
	}
}

// Generate a new random permutation of the y vector (Fisher-Yates)
// FIXME this can probably be optimized and may take a large percent of the CPU time
// Also, rand() is not really thread safe. On Windows it is kind of ok to use but
// on Linux this is a particularly bad idea since the state is process-wide.

void StatsComputer::resample_univariate(void) {
	for (int i = xy_nrow - 1; i > 0; --i) {
		int j = rand() % (i + 1); // should use my_rand() implemented below

		double temp = y[j];
		y[j] = y[i];
		y[i] = temp;
	}
}

void StatsComputer::resample_multivariate(void) {
	// In this case we don't really care about y values, but rather their
	// per (permuted) row sorted (permuted) indices.

	// (FIXME I may have a problem with the other tests implemented here,
	// that we compare against: they do care about the values themselves,
	// but for now this has no effect)

	for (int i = 0; i < xy_nrow; ++i) {
		int j = rand() % (i + 1); // should use my_rand() implemented below

		idx_perm[i] = idx_perm[j];
		idx_perm[j] = i;
	}

	for (int i = 0; i < xy_nrow; ++i) {
		idx_perm_inv[idx_perm[i]] = i;
	}
}

void StatsComputer::resample_uvz_ci(void) {
	// A smoothed bootstrap. The smoothing is done over the values of z, so for each
	// sample i we sample uniformly (i.e. currently I am assuming a NN kernel for the
	// smoothing, and even assuming its width is the same as the HHG statistic's
	// kernel) from the z-neighborhood of i. We need to generate resampling indices
	// for x (stored in idx_perm_inv), and independent indices for y (stored in idx_perm).
	// The current implementation also assumes that the observations have been presorted
	// according to their z values (which is only meaningful for univariate z).

	// I guess this could be optimized...

	int bwh = (nnh_lsb >> 1);
	for (int i = 0; i < xy_nrow; ++i) {
		int i1 = max(0, i - bwh);
		int i2 = min(xy_nrow - 1, i + bwh);
		idx_perm    [i] = my_rand(i1, i2);
		idx_perm_inv[i] = my_rand(i1, i2);
	}
}

void StatsComputer::resample_mvz_ci(void) {
	// Same as the univariate version, but has to look in the sorted dz in order to
	// figure out the indices of neighbors

	for (int i = 0; i < xy_nrow; ++i) {
		int nn_x = my_rand(0, nnh_lsb - 1);
		int nn_y = my_rand(0, nnh_lsb - 1);
		idx_perm_inv[i] = (*sorted_dz)[i][nn_x].second;
		idx_perm    [i] = (*sorted_dz)[i][nn_y].second;
	}
}

void StatsComputer::resample_dummy(void) {
	should_randomize = true;
}

void StatsComputer::sort_xy_distances_per_row(void) {
	// FIXME this REALLY deserves parallelizing. Trivial to do per row.

	for (int k = 0; k < xy_nrow; ++k) {
		for (int l = 0; l < xy_nrow; ++l) {
			sorted_dx_gen[k][l].x = dx[l * xy_nrow + k];
			sorted_dx_gen[k][l].y = dy[idx_perm[l] * xy_nrow + idx_perm[k]];
			sorted_dx_gen[k][l].i = l;
		}

		sort(sorted_dx_gen[k].begin(), sorted_dx_gen[k].end(), dbl_dbl_int_pair_comparator_xy);
	}
}

void StatsComputer::hhg_two_sample(void) {
	int n = xy_nrow;
    int a00, a01, a10, a11;
    double nrmlz = 1.0 / (n - 2);
    int i, j, k;
    int y_i, total_same_y, total_same_y_count_so_far, curr_dx_same_y_count;

	sum_chi  = 0;
	sum_like = 0;
	max_chi  = 0;
	max_like = 0;

    for (i = 0; i < n; i++) {
        k = 0;
        y_i = (int)y[i];
        total_same_y = y_counts[y_i]; // NOTE: this assumes that y itself is an index in [0, 1]
        total_same_y_count_so_far = 0;
        curr_dx_same_y_count = 0;

    	for (j = 0; j < n - 1; j++) {
        	k += ((*sorted_dx)[i][k].second == i); // exclude d(i,i)
    		curr_dx_same_y_count += (y[(*sorted_dx)[i][k].second] == y_i);

        	if ((k == n - 1) || ((*sorted_dx)[i][k + 1 + ((*sorted_dx)[i][k+1].second == i)].first > (*sorted_dx)[i][k].first)) {
        		// found all duplicates at current dx

        		if (curr_dx_same_y_count > 0) {
        			// In the notation of the HHG test paper (describing the 2x2 table):
        			//
					// The value "A_1.": j+1 is the number of dx smaller or equal to
        			// the current dx, but we don't want to include the current point
        			// in the count so: j
        			//
					// The value "A_.1": in the two sample setup is
        			// necessarily total_same_y - 2
        			//
					// The value "A_11": total_same_y_count_so_far + curr_dx_same_y_count - 1
        			// where the minus one is for the current point again.
        			//
        			// All remaining table values can be computed from these.
        			//
					// When adding to the sum statistics, add the resulting chisq times-
					// curr_dx_same_y_count, since each such sample has the same table
					// (and all the others have degenerate tables that contribute nothing)

					a00 = total_same_y_count_so_far + curr_dx_same_y_count - 1;
					a01 = j - a00;
					a10 = total_same_y - 2 - a00;
					a11 = n - 2 - j - a10;

					// Note that it is expected this would only be necessary for the computing the
					// observed statistic (the permutation is identity).
					// FIXME it might be a better idea to create a separate copy of this function
					// and add this only in the copy.
					if (store_tables) {
						int row = i * n + (*sorted_dx)[i][k].second;
						tbls[        row] = a00;
						tbls[  n*n + row] = a01;
						tbls[2*n*n + row] = a10;
						tbls[3*n*n + row] = a11;
					}

#ifdef DEBUG_CHECKS
					if (!((a00 >= 0) && (a01 >= 0) && (a10 >= 0) && (a11 >= 0) && (a00 + a01 + a10 + a11 == n - 2))) {
						cout << "THIS IS NOT A VALID CONTINGENCY TABLE !!!" << endl;
						exit(1);
					}
#endif

					accumulate_2x2_contingency_table(a00, a01, a10, a11, nrmlz, curr_dx_same_y_count);
        		}

        		// update/reset the counters in preparation for the next unique dx
        		total_same_y_count_so_far += curr_dx_same_y_count;
        		curr_dx_same_y_count = 0;
        	}

        	++k;
		}
    }
}

// NOTE: Actually now this is exactly the same as the 2-sample implementation (just the y_counts here would be of length K)
void StatsComputer::hhg_k_sample(void) {
	int n = xy_nrow;
    int a00, a01, a10, a11;
    double nrmlz = 1.0 / (n - 2);
    int i, j, k;
    int y_i, total_same_y, total_same_y_count_so_far, curr_dx_same_y_count;

	sum_chi  = 0;
	sum_like = 0;
	max_chi  = 0;
	max_like = 0;

    for (i = 0; i < n; i++) {
#ifdef DEBUG_PRINTS
    	int pi = i; // FIXME

		cout << "Working on center point " << i << " (y-permuted to " << pi << ")" << endl;
		cout << "This point has y = " << y[i] << ", which is the same in " << y_counts[(int)y[i]] << " points." << endl;
		cout << "Distances dx, dy for this row:" << endl;
		for (j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += (k == i);
			cout << j << " (" << k << "): " << dx[k*n+i] << ", " << dy[k*n+i] << endl;
		}
		cout << "Marginally sorted distances dx (and src idx) for this row:" << endl;
		for (j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += ((*sorted_dx)[i][k].second == i);
			cout << j << ": " << (*sorted_dx)[i][k].first << " (" << (*sorted_dx)[i][k].second << ")" << endl;
		}
		cout << "Marginally sorted distances dy (and src idx) for this row:" << endl;
		for (j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += ((*sorted_dy)[pi][k].second == i);
			cout << j << ": " << (*sorted_dy)[pi][k].first << " (" << (*sorted_dy)[pi][k].second << ")" << endl;
		}
#endif

        k = 0;
        y_i = (int)y[i];
        total_same_y = y_counts[y_i]; // NOTE: this assumes that y itself is an index in [0, K-1]. Otherwise need to implement y_counts as an STL map or something
        total_same_y_count_so_far = 0;
        curr_dx_same_y_count = 0;

    	for (j = 0; j < n - 1; j++) {
        	k += ((*sorted_dx)[i][k].second == i); // exclude d(i,i)
    		curr_dx_same_y_count += (y[(*sorted_dx)[i][k].second] == y_i);

        	if ((k == n - 1) || ((*sorted_dx)[i][k + 1 + ((*sorted_dx)[i][k+1].second == i)].first > (*sorted_dx)[i][k].first)) {
        		// found all duplicates at current dx

        		if (curr_dx_same_y_count > 0) {
        			// In the notation of the HHG test paper (describing the 2x2 table):
        			//
					// The value "A_1.": j+1 is the number of dx smaller or equal to
        			// the current dx, but we don't want to include the current point
        			// in the count so: j
        			//
					// The value "A_.1": in the two sample setup is
        			// necessarily total_same_y - 2
        			//
					// The value "A_11": total_same_y_count_so_far + curr_dx_same_y_count - 1
        			// where the minus one is for the current point again.
        			//
        			// All remaining table values can be computed from these.
        			//
					// When adding to the sum statistics, add the resulting chisq times-
					// curr_dx_same_y_count, since each such sample has the same table
					// (and all the others have degenerate tables that contribute nothing)

					a00 = total_same_y_count_so_far + curr_dx_same_y_count - 1;
					a01 = j - a00;
					a10 = total_same_y - 2 - a00;
					a11 = n - 2 - j - a10;

					// Note that it is expected this would only be necessary for the computing the
					// observed statistic (the permutation is identity).
					// FIXME it might be a better idea to create a separate copy of this function
					// and add this only in the copy.
					if (store_tables) {
						int row = i * n + (*sorted_dx)[i][k].second;
						tbls[        row] = a00;
						tbls[  n*n + row] = a01;
						tbls[2*n*n + row] = a10;
						tbls[3*n*n + row] = a11;
					}

#ifdef DEBUG_PRINTS
					cout << "with point at position " << j << " of the sorted dx: ";
					cout << "a00 = " << a00 << ", a01 = " << a01 << ", a10 = " << a10 << ", a11 = " << a11 << endl;
#endif

#ifdef DEBUG_CHECKS
					if (!((a00 >= 0) && (a01 >= 0) && (a10 >= 0) && (a11 >= 0) && (a00 + a01 + a10 + a11 == n - 2))) {
						cout << "THIS IS NOT A VALID CONTINGENCY TABLE !!!" << endl;
						exit(1);
					}
#endif

					accumulate_2x2_contingency_table(a00, a01, a10, a11, nrmlz, curr_dx_same_y_count);
        		}

        		// update/reset the counters in preparation for the next unique dx
        		total_same_y_count_so_far += curr_dx_same_y_count;
        		curr_dx_same_y_count = 0;
        	}

        	++k;
		}

#ifdef DEBUG_PRINTS
		cout << "current stats: sum_like = " << sum_like << ", sum_chi = " << sum_chi << endl;
#endif
    }
}

void StatsComputer::hhg_no_ties(void) {
	int n = xy_nrow, src;
    int a00, a01, a10, a11;
    double nrmlz = 1.0 / (n - 2);

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	for (int i = 0; i < n; ++i) {
		int pi = idx_perm[i];

#ifdef DEBUG_PRINTS
		cout << "Working on center point " << i << " (y-permuted to " << pi << ")" << endl;
		cout << "Distances dx, dy for this row:" << endl;
		for (int j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += (k == i);
			cout << j << " (" << k << "): " << dx[k*n+i] << ", " << dy[k*n+i] << endl;
		}
		cout << "Marginally sorted distances dx (and src idx) for this row:" << endl;
		for (int j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += ((*sorted_dx)[i][k].second == i);
			cout << j << ": " << (*sorted_dx)[i][k].first << " (" << (*sorted_dx)[i][k].second << ")" << endl;
		}
		cout << "Marginally sorted distances dy (and src idx) for this row:" << endl;
		for (int j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += ((*sorted_dy)[pi][k].second == i);
			cout << j << ": " << (*sorted_dy)[pi][k].first << " (" << (*sorted_dy)[pi][k].second << ")" << endl;
		}
#endif

		// Use Yair's merge-sort-like implementation (assumes there are no ties)

		// NOTE: I am using the fact that the sorting of permuted y's is the same as
		// the sorting of original y's. I only have to go (a) to the permuted row instead
		// of the i'th row, and (b) pass the "sorted indices" through the permutation.

		for (int j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += ((*sorted_dy)[pi][k].second == pi);
			src = idx_perm_inv[(*sorted_dy)[pi][k].second]; // NOTE: k may be different than from the line above
			src -= (src > i);
			hhg_gen_y_rev[src] = j;
		}

		for (int j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += ((*sorted_dx)[i][k].second == i);
			src = (*sorted_dx)[i][k].second; // NOTE: k may be different than from the line above
			src -= (src > i);
			hhg_gen_xy_perm[j] = hhg_gen_y_rev[src];
			hhg_gen_source[j] = j;
			hhg_gen_inversion_count[j] = 0;
			hhg_gen_xy_perm_temp[j] = hhg_gen_xy_perm[j];
		}

		hhg_gen_inversions(hhg_gen_xy_perm_temp, hhg_gen_source, hhg_gen_inversion_count, n - 1);

		for (int j = 0, k = 0; j < n - 1; ++j, ++k) {
			a00 = j - hhg_gen_inversion_count[j];
			a01 = hhg_gen_inversion_count[j];
			a10 = hhg_gen_xy_perm[j] + hhg_gen_inversion_count[j] - j;
			a11 = n - hhg_gen_xy_perm[j] - hhg_gen_inversion_count[j] - 2;

			// Note that it is expected this would only be necessary for the computing the
			// observed statistic (the permutation is identity).
			// FIXME it might be a better idea to create a separate copy of this function
			// and add this only in the copy.
			if (store_tables) {
				k += ((*sorted_dx)[i][k].second == i);
				int row = i * n + (*sorted_dx)[i][k].second;
				tbls[        row] = a00;
				tbls[  n*n + row] = a01;
				tbls[2*n*n + row] = a10;
				tbls[3*n*n + row] = a11;
			}

#ifdef DEBUG_PRINTS
			cout << "with point at position " << j << " of the sorted dx: ";
			cout << "a00 = " << a00 << ", a01 = " << a01 << ", a10 = " << a10 << ", a11 = " << a11 << endl;
#endif

#ifdef DEBUG_CHECKS
			if (!((a00 >= 0) && (a01 >= 0) && (a10 >= 0) && (a11 >= 0) && (a00 + a01 + a10 + a11 == n - 2))) {
				cout << "THIS IS NOT A VALID CONTINGENCY TABLE !!!" << endl;
				//exit(1);
			}
#endif

			accumulate_2x2_contingency_table(a00, a01, a10, a11, nrmlz, 1);
		}
	}
}

void StatsComputer::hhg_general(void) {
	int n = xy_nrow, src;
    int a00, a01, a10, a11;
    double nrmlz = 1.0 / (n - 2);

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	for (int i = 0; i < n; ++i) {
		int pi = idx_perm[i];

#ifdef DEBUG_PRINTS
		cout << "Working on center point " << i << " (y-permuted to " << pi << ")" << endl;
		cout << "Distances dx, dy for this row:" << endl;
		for (int j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += (k == i);
			cout << j << " (" << k << "): " << dx[k*n+i] << ", " << dy[idx_perm[k]*n+pi] << endl;
		}
		cout << "Lexicographically [x asc, y desc] sorted distances dx (and dy, src idx) for this row:" << endl;
		for (int j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += (sorted_dx_gen[i][k].i == i);
			cout << j << ": " << sorted_dx_gen[i][k].x << " (dy = " << sorted_dx_gen[i][k].y << ", src = " << sorted_dx_gen[i][k].i << ")" << endl;
		}
		cout << "Marginally sorted distances dy (and src idx) for this row:" << endl;
		for (int j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += ((*sorted_dy)[pi][k].second == i);
			cout << j << ": " << (*sorted_dy)[pi][k].first << " (" << (*sorted_dy)[pi][k].second << ")" << endl;
		}
#endif

		// Use Yair's merge-sort-like implementation (assumes there are no ties)

		// NOTE: I am using the fact that the sorting of permuted y's is the same as
		// the sorting of original y's. I only have to go (a) to the permuted row instead
		// of the i'th row, and (b) pass the "sorted indices" through the permutation.

		double last_new_y = 0;
		int last_new_y_pos = -1;

		for (int j = n - 1, k = n - 1; j >= 1; --j, --k) {
			k -= ((*sorted_dy)[pi][k].second == pi);
			if (last_new_y_pos == -1 || (*sorted_dy)[pi][k].first != last_new_y) {
				last_new_y = (*sorted_dy)[pi][k].first;
				last_new_y_pos = j;
			}
			src = idx_perm_inv[(*sorted_dy)[pi][k].second]; // NOTE: k may be different than from the line above
			src -= (src > i);
			hhg_gen_y_rev[src] = last_new_y_pos;
		}

		for (int j = 0, k = 0; j < n - 1; ++j, ++k) {
			k += (sorted_dx_gen[i][k].i == i);
			src = sorted_dx_gen[i][k].i; // NOTE: k may be different than from the line above
			src -= (src > i);
			hhg_gen_xy_perm[j] = hhg_gen_y_rev[src];
			hhg_gen_source[j] = j;
			hhg_gen_inversion_count[j] = 0;
			hhg_gen_xy_perm_temp[j] = hhg_gen_xy_perm[j];
		}

#ifdef DEBUG_PRINTS
		cout << "Contents of y_rev: (should be: rank of y in marginal sorting, listed in original order)" << endl;
		for (int j = 0; j < n - 1; ++j) {
			cout << j << ": " << hhg_gen_y_rev[j] << endl;
		}
		cout << "Contents of xy_perm: (should be: rank of y of points in order of x)" << endl;
		for (int j = 0; j < n - 1; ++j) {
			cout << j << ": " << hhg_gen_xy_perm[j] << endl;
		}
#endif

		hhg_gen_inversions(hhg_gen_xy_perm_temp, hhg_gen_source, hhg_gen_inversion_count, n - 1);

#ifdef DEBUG_PRINTS
		cout << "Contents of inversion_count: (should be: number of inversions in y from left)" << endl;
		for (int j = 0; j < n - 1; ++j) {
			cout << j << ": " << hhg_gen_inversion_count[j] << endl;
		}
#endif

		// (This was an alternative, that I did not pursue)
		// When ties are possible, what we want is, for every point j in the sorted order of x's,
		// (going with j from last to first): how many points have smaller or equal x (that would
		// be the last positions we encountered while walking with j that had a new value) and
		// also have a larger y (this is the number of inversions in y from all those points.
		// this can be decomposed as the inversions in y to the left of j, plus the inversions
		// in y to the right-and-same-x-value. The former we already have, the latter I think we
		// can compute in nlog(n) with a very similar algorithm to "inversions()").

		// (This is what I did)
		// An alternative to the above is not to sort the x's but instead sort according to
		// "x ascending, then y descending" (and so this really has to be redone for every
		// permutation...). If you do this, then the remainder of the algorithm remains unchanged.

		double last_new_x = 0;
		int last_new_x_pos = -1;

		for (int j = n - 2, k = n - 1; j >= 0; --j, --k) {
			k -= ((*sorted_dx)[i][k].second == i);
			if (last_new_x_pos == -1 || (*sorted_dx)[i][k].first != last_new_x) {
				last_new_x = (*sorted_dx)[i][k].first;
				last_new_x_pos = j;
			}

			a00 = last_new_x_pos - hhg_gen_inversion_count[j];
			a01 = hhg_gen_inversion_count[j];
			a10 = hhg_gen_xy_perm[j] - 1 + hhg_gen_inversion_count[j] - last_new_x_pos;
			a11 = n - hhg_gen_xy_perm[j] - hhg_gen_inversion_count[j] - 1;

			// Note that it is expected this would only be necessary for the computing the
			// observed statistic (the permutation is identity).
			// FIXME it might be a better idea to create a separate copy of this function
			// and add this only in the copy.
			if (store_tables) {
				int row = i * n + (*sorted_dx)[i][k].second;
				tbls[        row] = a00;
				tbls[  n*n + row] = a01;
				tbls[2*n*n + row] = a10;
				tbls[3*n*n + row] = a11;
			}

#ifdef DEBUG_PRINTS
			cout << "with point at position " << j << " of the sorted dx: ";
			cout << "a00 = " << a00 << ", a01 = " << a01 << ", a10 = " << a10 << ", a11 = " << a11 << endl;
#endif

#ifdef DEBUG_CHECKS
			if (!((a00 >= 0) && (a01 >= 0) && (a10 >= 0) && (a11 >= 0) && (a00 + a01 + a10 + a11 == n - 2))) {
				cout << "THIS IS NOT A VALID CONTINGENCY TABLE !!!" << endl;
				//exit(1);
			}
#endif

			accumulate_2x2_contingency_table(a00, a01, a10, a11, nrmlz, 1);
		}

#ifdef DEBUG_PRINTS
		cout << "current stats: sum_like = " << sum_like << ", sum_chi = " << sum_chi << endl;
#endif
	}
}

// Conditional independence test using the nearest neighbor kernel (univariate z)
void StatsComputer::hhg_ci_uvz_nn(void) {
	int n = xy_nrow;
	int nnhh = nnh / 2;
	int nnht = nnhh * 2 - 1;
	double nrmlz = 1.0 / nnht;
	int count[2][2];

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	double dx_ij, dx_ik, dy_ij, dy_ik;

	for (int i = nnhh; i < n - nnhh; ++i) { // FIXME what about the edges? completely ignored??
		// FIXME: can the following be computed faster? if nnhh is expected to be a small number
		// it doesn't seem like a good idea to compute in a fancy way like hhg_general(). Maybe
		// better to exploit the symmetry in the distance comparison matrices.

		int pi_x = idx_perm_inv[i];
		int pi_y = idx_perm    [i];

		for (int j = i - nnhh; j <= i + nnhh; ++j) {
			// A smoothed bootstrap that resamples both x and y
			if (j != i) {
				count[0][0] = 0;
				count[0][1] = 0;
				count[1][0] = 0;
				count[1][1] = 0;

				int pj_x = idx_perm_inv[j];
				int pj_y = idx_perm    [j];

				dx_ij = dx[pj_x * n + pi_x];
				dy_ij = dy[pj_y * n + pi_y];

				for (int k = i - nnhh; k <= i + nnhh; ++k) {
					if ((k != j) && (k != i)) { // FIXME this condition can be avoided by unrolling to three loops
						int pk_x = idx_perm_inv[k];
						int pk_y = idx_perm    [k];

						dx_ik = dx[pk_x * n + pi_x];
						dy_ik = dy[pk_y * n + pi_y];

						++count[dx_ij > dx_ik][dy_ij > dy_ik];
					}
				}
			}

			accumulate_2x2_contingency_table(count[0][0], count[0][1], count[1][0], count[1][1], nrmlz, 1);
		}
	}
}

// Conditional independence test using a gaussian kernel (univariate z)
void StatsComputer::hhg_ci_uvz_gaussian(void) {
	int n = xy_nrow;
	double ebw = 3 * sig;
	double nrmlz = 1 / (2 * M_PI * sig * sig), enrmlz = -0.5 / (sig * sig), wijk;
	double count[2][2];

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	Rprintf("NOTE: THIS IS BROKEN\n");

	double dx_ij, dx_ik, dy_ij, dy_ik;
	int zli = 0, zhi = 0, ii, knn;

	for (int i = 0; i < n; ++i) {
		// Can do binary search, but if ebw is not too big that may be useless
		while ((zli < i) && (z[zli] < z[i] - ebw)) {
			++zli;
		}
		while ((zhi < n - 1) && (z[zhi] < z[i] + ebw)) {
			++zhi;
		}
		if (z[zhi] > z[i] + ebw) {
			--zhi;
		}
		knn = zhi - zli;

		for (ii = 0; ii < i - zli; ++ii) {
			idx_1_to_n[ii] = idx_perm[ii] = zli + ii;
		}
		for (; ii < knn; ++ii) {
			idx_1_to_n[ii] = idx_perm[ii] = zli + ii + 1;
		}
		if (should_randomize) {
			// FIXME smoothed bootstrap using a separate kernel
			for (ii = knn - 1; ii > 0; --ii) {
				int jj = rand() % (ii + 1); // this can be heavily biased since RAND_MAX is 32K and taking the modulo will give higher probability to smaller numbers

				int temp = idx_perm[jj];
				idx_perm[jj] = idx_perm[ii];
				idx_perm[ii] = temp;
			}
		}

		for (int j = 0; j < knn; ++j) {
			count[0][0] = 0;
			count[0][1] = 0;
			count[1][0] = 0;
			count[1][1] = 0;

			int jj = idx_1_to_n[j];
			dx_ij = dx[jj * n + i];
			dy_ij = dy[idx_perm[j] * n + i];

			for (int k = 0; k < knn; ++k) {
				if (k != j) {
					int kk = idx_1_to_n[k];
					dx_ik = dx[kk * n + i];
					dy_ik = dy[idx_perm[k] * n + i];

					wijk = nrmlz * exp(enrmlz * ((z[jj] - z[i]) * (z[jj] - z[i]) + (z[kk] - z[i]) * (z[kk] - z[i])));
					count[dx_ij > dx_ik][dy_ij > dy_ik] += wijk;
				}
			}

			accumulate_2x2_contingency_table(count[0][0], count[0][1], count[1][0], count[1][1], nrmlz, 1);
		}
	}

	should_randomize = false;
}

// Conditional independence test using the nearest neighbor kernel (multivariate z)
void StatsComputer::hhg_ci_mvz_nn(void) {
	int n = xy_nrow;
	double nrmlz = 1.0 / (nnh - 1);
	double dx_ij, dx_ik, dy_ij, dy_ik;
	double count[2][2];

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	// NOTE: here I am using nnh nearest neighbors (not necessarily nnh/2 on each "side")
	// because it is not clear to me what is an equivalent notion of a "side" high dimension

	for (int i = 0; i < n; ++i) {
		int pi_x = idx_perm_inv[i];
		int pi_y = idx_perm    [i];

		for (int j = 1; j <= nnh; ++j) {
			count[0][0] = 0;
			count[0][1] = 0;
			count[1][0] = 0;
			count[1][1] = 0;

			// We want the j'th z-neighbor of the i'th sample, in the bootstrapped data.
			int jj = (*sorted_dz)[i][j].second; // the z's are not bootstrapped
			int pj_x = idx_perm_inv[jj]; // and this gives us the original sample that is the j'th neighbor in the bootstrapped data
			int pj_y = idx_perm    [jj];

			dx_ij = dx[pj_x * n + pi_x]; // FIXME though this is more readable, working on the transposed matrix may have better locality or reference
			dy_ij = dy[pj_y * n + pi_y];

			for (int k = 1; k < j; ++k) {
				int kk = (*sorted_dz)[i][k].second;
				int pk_x = idx_perm_inv[kk];
				int pk_y = idx_perm    [kk];

				dx_ik = dx[pk_x * n + pi_x];
				dy_ik = dy[pk_y * n + pi_y];

				++count[dx_ij > dx_ik][dy_ij > dy_ik];
			}
			for (int k = j + 1; k <= nnh; ++k) {
				int kk = (*sorted_dz)[i][k].second;
				int pk_x = idx_perm_inv[kk];
				int pk_y = idx_perm    [kk];

				dx_ik = dx[pk_x * n + pi_x];
				dy_ik = dy[pk_y * n + pi_y];

				++count[dx_ij > dx_ik][dy_ij > dy_ik];
			}

			accumulate_2x2_contingency_table(count[0][0], count[0][1], count[1][0], count[1][1], nrmlz, 1);
		}
	}

	// FIXME it may be necessary to normalize the sum statistics only by the number of tables
	// passing the w_sum/w_max conditions
	sum_chi /= nnh;
	sum_like /= nnh;
}

// Conditional independence test using a gaussian kernel (multivariate z)
void StatsComputer::hhg_ci_mvz_gaussian(void) {
	int n = xy_nrow;
	double count[2][2];
	double dx_ij, dx_ik, dy_ij, dy_ik;
	int knn;

	// FIXME 1: I need to compensate for the truncation I am doing (probably best to
	// precompute the normalizing constant in R and send it here as another kernel
	// parameter)

	// FIXME 2: allow specification of a general Gaussian kernel? (with a general
	// covariance matrix). I don't see justification for this.

	double ebw = 3 * sig;
	double nrmlz = 1 / (2 * M_PI * sig * sig), enrmlz = -0.5 / (sig * sig), wijk;

	Rprintf("NOTE: THIS IS BROKEN\n");

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	for (int i = 0; i < n; ++i) {
		// Can do binary search, but if ebw is not too big that may be useless
		for (knn = 0; (knn < n - 1) && ((*sorted_dz)[i][knn + 1].first < ebw); ++knn) {
			idx_1_to_n[knn] = idx_perm[knn] = (*sorted_dz)[i][knn + 1].second;
		}

		if (should_randomize) {
			// FIXME smoothed bootstrap using a separate kernel
			for (int ii = knn - 1; ii > 0; --ii) {
				int jj = rand() % (ii + 1); // should use my_rand() implemented below

				int temp = idx_perm[jj];
				idx_perm[jj] = idx_perm[ii];
				idx_perm[ii] = temp;
			}
		}

		for (int j = 0; j < knn; ++j) {
			count[0][0] = 0;
			count[0][1] = 0;
			count[1][0] = 0;
			count[1][1] = 0;

			int jj = idx_1_to_n[j];
			dx_ij = dx[         jj * n + i];
			dy_ij = dy[idx_perm[j] * n + i];

			for (int k = 0; k < knn; ++k) {
				if (k != j) {
					int kk = idx_1_to_n[k];
					dx_ik = dx[         kk * n + i];
					dy_ik = dy[idx_perm[k] * n + i];

					wijk = nrmlz * exp(enrmlz * ((z[jj] - z[i]) * (z[jj] - z[i]) + (z[kk] - z[i]) * (z[kk] - z[i])));
					count[dx_ij > dx_ik][dy_ij > dy_ik] += wijk;
				}
			}

			accumulate_2x2_contingency_table(count[0][0], count[0][1], count[1][0], count[1][1], nrmlz, 1);
		}
	}

	should_randomize = false;
}

// Conditional independence test using the nearest neighbor kernel (multivariate z)
// In this version we go over a grid of bandwidth values and compute the regular
// test for each of them. The main set of statistics is then replaced with the
// maximum statistic across the grid.
void StatsComputer::hhg_ci_mvz_nn_grid(void) {
	double max_nnh_sum_chi  = 0;
	double max_nnh_max_chi  = 0;
	double max_nnh_sum_like = 0;
	double max_nnh_max_like = 0;

	for (int i = 0; i < nnh_grid_cnt; ++i) {
		nnh = nnh_grid[i];
		hhg_ci_mvz_nn();
		sum_chi_grid [i] = sum_chi;
		sum_like_grid[i] = sum_like;
		max_chi_grid [i] = max_chi;
		max_like_grid[i] = max_like;

		if (sum_chi > max_nnh_sum_chi) {
			max_nnh_sum_chi = sum_chi;
		}
		if (sum_like > max_nnh_sum_like) {
			max_nnh_sum_like = sum_like;
		}
		if (max_chi > max_nnh_max_chi) {
			max_nnh_max_chi = max_chi;
		}
		if (max_like > max_nnh_max_like) {
			max_nnh_max_like = max_like;
		}
	}

	sum_chi  = max_nnh_sum_chi;
	max_chi  = max_nnh_max_chi;
	sum_like = max_nnh_sum_like;
	max_like = max_nnh_max_like;
}

// General notes about the univariate distribution-free test (UDF) variants:
// -----------------------------------------------------------------------------
//
// In this test type, the relevant inputs are rx (stored in <dx>) and ry (stored
// in <y> so that it can be permuted the same way as in the 2/k-sample tests).
//
// NOTE: in this C implementation, I expect ranks to be in 0,1,2,...,n-1
//
// The SPR test uses single points in the paired rank-rank sample space as loci
// for examining the local dependence. This is equivalent to Hoeffding's test.
// The PPR test can be seen as direct extension of this idea, which uses pairs
// of points in order to examine local dependence. In this case we can,
// hopefully, model better the distribution of points relative to one another.
// The TPR test goes on to look at triplets of points. The SPPR variant is an
// attempt to combine the SPR and PPR scores.
//
// Given the sample size <n>, the inputs <xr> and <yr> (the ranks of some <x>
// and <y> on their original scales) are merely two permutations of the numbers
// 1:n. We thus only need to focus on this <n> over <n> grid. Our statistic
// concerns the number of (x,y) points in the sample that fall in the rectangle
// whose diagonal connects two points from the set and has sides parallel to the
// axes. There are n(n - 1) such rectangles, so merely going over them is O(n^2),
// defining a lower bound for complexity. An algorithm that achieves this bound
// first computes the cumulative sum over the grid, of the indicator of whether
// the point (i,j) is in the set or not. Subsequently we can go over the n(n - 1)
// pairs of points and compute the number of points inside the rectangle they
// define in O(1) using the "double integral" we have computed.
//
// The '33' (3x3) variant is almost the same as the '22' (2,2) variant. The
// difference is that here, for each pair of points from the set, the 3x3
// contingency table for the 3x3 partition of the (x,y) plane is computed,
// whose center cell is defined by the rectangle used in the (2,2) variant.
//
// The "all" variants are exactly the same as the "obs" variants, only they
// run over all rank-rank pairs, not just those observed.
//
// FIXME: I am replicating code like crazy instead of encapsulating it!
//
// -----------------------------------------------------------------------------

void StatsComputer::hhg_udf_spr_obs(void) {
	// First, compute the double integral (padded with a leading row and column of zeros)
	compute_double_integral(xy_nrow, dx, y);

	// Now compute the score using all points
	int n = xy_nrow;
	int nm1 = n - 1;
	double nm1d = n - 1;

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	ng_chi  = 0;
	ng_like = 0;

	// TODO: collect individual tables if wanted

	int yi, xi;

	for (int i = 0; i < n; ++i) {
		yi = y[i];
		xi = dx[i];

#ifndef UDF_ALLOW_DEGENERATE_PARTITIONS
		// FIXME points with extreme ranks can be filtered beforehand...
		if (xi == 0 || xi == nm1 || yi == 0 || yi == nm1) {
			continue;
		}
#endif

		compute_spr_obs(xi, yi, n, dintegral_pn, nm1, nm1d);
	}

#ifdef UDF_NORMALIZE
	ng_chi *= n;
	ng_like *= n;
	sum_chi /= ng_chi;
	sum_like /= ng_like;
#endif
}

void StatsComputer::hhg_udf_spr_all(void) {
	// First, compute the double integral (padded with a leading row and column of zeros)
	compute_double_integral(xy_nrow, dx, y);

	// Now compute the score using all points
	int n = xy_nrow;
	double nd = n;

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	ng_chi  = 0;
	ng_like = 0;

	// TODO: collect individual tables if wanted

#ifndef UDF_ALLOW_DEGENERATE_PARTITIONS
	for (int xi = 1; xi < n; ++xi) {
		for (int yi = 1; yi < n; ++yi) {
#else
	int np1 = n + 1;
	for (int xi = 0; xi < np1; ++xi) {
		for (int yi = 0; yi < np1; ++yi) {
#endif
			compute_spr_all(xi, yi, n, dintegral_pn, nd);
		}
	}

#ifdef UDF_NORMALIZE
	ng_chi *= n;
	ng_like *= n;
	sum_chi /= ng_chi;
	sum_like /= ng_like;
#endif
}

void StatsComputer::hhg_udf_ppr_22_obs(void) {
	// First, compute the double integral (padded with a leading row and column of zeros)
	compute_double_integral(xy_nrow, dx, y);

	// Now compute the score using all rectangles
	int n = xy_nrow;
#ifndef UDF_ALLOW_DEGENERATE_PARTITIONS
	int nm1 = n - 1;
#endif
	int nm2 = n - 2;
	double nm2s = nm2 * nm2;

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	ng_chi  = 0;
	ng_like = 0;

	// TODO: collect individual tables if wanted

	int yi, yj, xi, xj;
	int xr_lo, xr_hi, yr_lo, yr_hi;

	for (int i = 0; i < n; ++i) {
		for (int j = i + 1; j < n; ++j) {
			yi = y[i]; xi = dx[i];
			yj = y[j]; xj = dx[j];

		    if (xi < xj) {
		      xr_lo = xi;
		      xr_hi = xj;
		    } else {
		      xr_lo = xj;
		      xr_hi = xi;
		    }

		    if (yi < yj) {
		      yr_lo = yi;
		      yr_hi = yj;
		    } else {
		      yr_lo = yj;
		      yr_hi = yi;
		    }

#ifndef UDF_ALLOW_DEGENERATE_PARTITIONS
			if (xr_lo == 0 || xr_hi == nm1 || yr_lo == 0 || yr_hi == nm1 ||
				xr_hi - xr_lo == 1 || yr_hi - yr_lo == 1)
			{
				continue;
			}
#endif

			compute_ppr_22(xr_lo, xr_hi, yr_lo, yr_hi, dintegral_pn, nm2, nm2s);
		}
	}

#ifdef UDF_NORMALIZE
	ng_chi *= n;
	ng_like *= n;
	sum_chi /= ng_chi;
	sum_like /= ng_like;
#endif
}

void StatsComputer::hhg_udf_ppr_22_all(void) {
	// First, compute the double integral (padded with a leading row and column of zeros)
	compute_double_integral(xy_nrow, dx, y);

	// Now compute the score using all rectangles
	int n = xy_nrow;
	int nm1 = n - 1;
	int nm2 = n - 2;
#ifndef UDF_ALLOW_DEGENERATE_PARTITIONS
	int nm3 = n - 3;
#endif
	double nm2s = nm2 * nm2;

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	ng_chi  = 0;
	ng_like = 0;

	// TODO: collect individual tables if wanted

#ifndef UDF_ALLOW_DEGENERATE_PARTITIONS
	for (int xr_lo = 1; xr_lo < nm3; ++xr_lo) {
		for (int xr_hi = xr_lo + 2; xr_hi < nm1; ++xr_hi) {
			for (int yr_lo = 1; yr_lo < nm3; ++yr_lo) {
				for (int yr_hi = yr_lo + 2; yr_hi < nm1; ++yr_hi) {
#else
	for (int xr_lo = 0; xr_lo < nm1; ++xr_lo) {
		for (int xr_hi = xr_lo + 1; xr_hi < n; ++xr_hi) {
			for (int yr_lo = 0; yr_lo < nm1; ++yr_lo) {
				for (int yr_hi = yr_lo + 1; yr_hi < n; ++yr_hi) {
#endif
					// FIXME this must be replaced with a special version that works on half ranks, see spr_all
					compute_ppr_22(xr_lo, xr_hi, yr_lo, yr_hi, dintegral_pn, nm2, nm2s);
				}
			}
		}
	}

#ifdef UDF_NORMALIZE
	ng_chi *= n;
	ng_like *= n;
	sum_chi /= ng_chi;
	sum_like /= ng_like;
#endif
}

void StatsComputer::hhg_udf_ppr_33_obs(void) {
	// First, compute the double integral (padded with a leading row and column of zeros)
	compute_double_integral(xy_nrow, dx, y);

	// Now compute the score using all rectangles
	int n = xy_nrow;
#ifndef UDF_ALLOW_DEGENERATE_PARTITIONS
	int nm1 = n - 1;
#endif
	double nm2 = n - 2;

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	ng_chi  = 0;
	ng_like = 0;

	// TODO: collect individual tables if wanted

	int yi, yj, xi, xj;
	int xr_lo, xr_hi, yr_lo, yr_hi;

	for (int i = 0; i < n; ++i) {
		for (int j = i + 1; j < n; ++j) {
			yi = y[i]; xi = dx[i];
			yj = y[j]; xj = dx[j];

		    if (xi < xj) {
		      xr_lo = xi;
		      xr_hi = xj;
		    } else {
		      xr_lo = xj;
		      xr_hi = xi;
		    }

		    if (yi < yj) {
		      yr_lo = yi;
		      yr_hi = yj;
		    } else {
		      yr_lo = yj;
		      yr_hi = yi;
		    }

#ifndef UDF_ALLOW_DEGENERATE_PARTITIONS
			if (xr_lo == 0 || xr_hi == nm1 || yr_lo == 0 || yr_hi == nm1 ||
				xr_hi - xr_lo <= 1 || yr_hi - yr_lo <= 1)
			{
				continue;
			}
#endif

			compute_ppr_33(xr_lo, xr_hi, yr_lo, yr_hi, n, dintegral_pn, nm2);
		}
	}

#ifdef UDF_NORMALIZE
	ng_chi *= n;
	ng_like *= n;
	sum_chi /= ng_chi;
	sum_like /= ng_like;
#endif
}

void StatsComputer::hhg_udf_ppr_33_all(void) {
	// First, compute the double integral (padded with a leading row and column of zeros)
	compute_double_integral(xy_nrow, dx, y);

	// Now compute the score using all rectangles
	int n = xy_nrow;
	int nm1 = n - 1;
#ifndef UDF_ALLOW_DEGENERATE_PARTITIONS
	int nm3 = n - 3;
#endif
	double nm2 = n - 2;

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	ng_chi  = 0;
	ng_like = 0;

	// TODO: collect individual tables if wanted

#ifndef UDF_ALLOW_DEGENERATE_PARTITIONS
	for (int xr_lo = 1; xr_lo < nm3; ++xr_lo) {
		for (int xr_hi = xr_lo + 2; xr_hi < nm1; ++xr_hi) {
			for (int yr_lo = 1; yr_lo < nm3; ++yr_lo) {
				for (int yr_hi = yr_lo + 2; yr_hi < nm1; ++yr_hi) {
#else
	for (int xr_lo = 0; xr_lo < nm1; ++xr_lo) {
		for (int xr_hi = xr_lo + 1; xr_hi < n; ++xr_hi) {
			for (int yr_lo = 0; yr_lo < nm1; ++yr_lo) {
				for (int yr_hi = yr_lo + 1; yr_hi < n; ++yr_hi) {
#endif
					// FIXME this must be replaced with a special version that works on half ranks, see spr_all
					compute_ppr_33(xr_lo, xr_hi, yr_lo, yr_hi, n, dintegral_pn, nm2);
				}
			}
		}
	}

#ifdef UDF_NORMALIZE
	ng_chi *= n;
	ng_like *= n;
	sum_chi /= ng_chi;
	sum_like /= ng_like;
#endif
}

void StatsComputer::hhg_udf_tpr_obs(void) {
	// First, compute the double integral (padded with a leading row and column of zeros)
	compute_double_integral(xy_nrow, dx, y);

	// Now compute the score using all rectangles
	int n = xy_nrow;
#ifndef UDF_ALLOW_DEGENERATE_PARTITIONS
	int nm1 = n - 1;
#endif
	double nm3 = n - 3;

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	ng_chi  = 0;
	ng_like = 0;

	// TODO: collect individual tables if wanted

	int yi, xi, yj, xj, yk, xk, xl, xm, xh, yl, ym, yh;
	bool fij, fik, fjk;

	for (int i = 0; i < n; ++i) {
		for (int j = i + 1; j < n; ++j) {
			for (int k = j + 1; k < n; ++k) {
				yi = y [i]; xi = dx[i];
				yj = y [j]; xj = dx[j];
				yk = y [k]; xk = dx[k];

				fij = (xi < xj); fik = (xi < xk); fjk = (xj < xk);

				if (fij && fjk) {
					xl = xi; xm = xj; xh = xk;
				} else if (fik && !fjk) {
					xl = xi; xm = xk; xh = xj;
				} else if (!fij && fik) {
					xl = xj; xm = xi; xh = xk;
				} else if (fjk && !fik) {
					xl = xj; xm = xk; xh = xi;
				} else if (!fik && fij) {
					xl = xk; xm = xi; xh = xj;
				} else {
					xl = xk; xm = xj; xh = xi;
				}

				fij = (yi < yj); fik = (yi < yk); fjk = (yj < yk);

				if (fij && fjk) {
					yl = yi; ym = yj; yh = yk;
				} else if (fik && !fjk) {
					yl = yi; ym = yk; yh = yj;
				} else if (!fij && fik) {
					yl = yj; ym = yi; yh = yk;
				} else if (fjk && !fik) {
					yl = yj; ym = yk; yh = yi;
				} else if (!fik && fij) {
					yl = yk; ym = yi; yh = yj;
				} else {
					yl = yk; ym = yj; yh = yi;
				}

#ifndef UDF_ALLOW_DEGENERATE_PARTITIONS
				if (xl == 0 || xh == nm1 || yl == 0 || yh == nm1 ||
					xh - xm <= 1 || xm - xl <= 1 || yh - ym <= 1 || ym - yl <= 1) {
					// partition is degenerate (not sure if I should still allow it to contribute to the overall score
					continue;
				}
#endif

				compute_tpr(xl, xm, xh, yl, ym, yh, n, dintegral_pn, nm3);
			}
		}
	}

#ifdef UDF_NORMALIZE
	ng_chi *= n;
	ng_like *= n;
	sum_chi /= ng_chi;
	sum_like /= ng_like;
#endif
}

void StatsComputer::hhg_udf_tpr_all(void) {
	// First, compute the double integral (padded with a leading row and column of zeros)
	compute_double_integral(xy_nrow, dx, y);

	// Now compute the score using all rectangles
	int n = xy_nrow;
	int nm1 = n - 1;
	double nm3d = n - 3;

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	ng_chi  = 0;
	ng_like = 0;

	// TODO: collect individual tables if wanted

#ifndef UDF_ALLOW_DEGENERATE_PARTITIONS
	int nm3 = n - 3;
	int nm5 = n - 5;

	for (int xl = 1; xl < nm5; ++xl) {
		for (int xm = xl + 2; xm < nm3; ++xm) {
			for (int xh = xm + 2; xh < nm1; ++xh) {
				for (int yl = 1; yl < nm5; ++yl) {
					for (int ym = yl + 2; ym < nm3; ++ym) {
						for (int yh = ym + 2; yh < nm1; ++yh) {
#else
	int nm2 = n - 2;

	for (int xl = 0; xl < nm2; ++xl) {
		for (int xm = xl + 1; xm < nm1; ++xm) {
			for (int xh = xm + 1; xh < n; ++xh) {
				for (int yl = 0; yl < nm2; ++yl) {
					for (int ym = yl + 1; ym < nm1; ++ym) {
						for (int yh = ym + 1; yh < n; ++yh) {
#endif
							// FIXME this must be replaced with a special version that works on half ranks, see spr_all
							compute_tpr(xl, xm, xh, yl, ym, yh, n, dintegral_pn, nm3d);
						}
					}
				}
			}
		}
	}

#ifdef UDF_NORMALIZE
	ng_chi *= n;
	ng_like *= n;
	sum_chi /= ng_chi;
	sum_like /= ng_like;
#endif
}

void StatsComputer::hhg_udf_sppr_obs(void) {
	// TODO
}

void StatsComputer::hhg_udf_sppr_all(void) {
	// TODO
}

// This computes DDP ("DDP_OBS")
void StatsComputer::hhg_udf_ddp(void) {
	// First, compute the double integral (padded with a leading row and column of zeros)
	compute_double_integral(xy_nrow, dx, y);

	// Now compute the score using all rectangles
	int n = xy_nrow;

	sum_chi  = 0;
	sum_like = 0;

	// These can't be computed this way
	max_chi  = 0;
	max_like = 0;

	// FIXME: It's not important for now, since DDP is not used with permutations
	// in this C code (rather, data for null table simulations is permuted in R),
	// but if I even want to use the permutation feature here I would have to make
	// sure that not only y (and hence double_integral) is permuted, but also
	// x_ordered_by_y and y_ordered_by_x, needed in count_ddp_with_given_cell.

	int rect_o;
	double rect_e, rect_c, rect_l, cnt;
	double edenom = 1.0 / (n - K + 1); // I wonder if the denominator here makes sense for the degenerate partitions
	double kahan_c_chi = 0, kahan_c_like = 0, kahan_t;
	double nr_parts = 0;

	double nr_nonempty_cells = 0;

	for (int xl = 1; xl <= n; ++xl) {
		for (int xh = xl; xh <= n; ++xh) {
			for (int yl = 1; yl <= n; ++yl) {
				for (int yh = yl; yh <= n; ++yh) {
					cnt = count_ddp_with_given_cell(xl, xh, yl, yh);
					if (cnt > 0) {
						rect_o = count_sample_points_in_rect(xl, xh, yl, yh);
						rect_e = (xh - xl + 1) * (yh - yl + 1) * edenom;
						rect_c = ((rect_o - rect_e) * (rect_o - rect_e)) / rect_e * cnt - kahan_c_chi;
						rect_l = ((rect_o > 0) ? (rect_o * log(rect_o / rect_e)) : 0) * cnt - kahan_c_like;

						kahan_t = sum_chi + rect_c;
						kahan_c_chi = (kahan_t - sum_chi) - rect_c;
						sum_chi = kahan_t;

						kahan_t = sum_like + rect_l;
						kahan_c_like = (kahan_t - sum_like) - rect_l;
						sum_like = kahan_t;

						nr_parts += cnt;

						if (rect_o > 0) {
							nr_nonempty_cells += cnt;
						}
					}
				}
			}
		}
	}

#ifdef UDF_NORMALIZE
  // NOTE: since there are degenerate partitions, nr_parts is not a multiple of
  // K^2... this is irrelevant for testing but may be relevant for estimating 
  // mutual information (which DDP is not really intended for)

	nr_parts /= (K * K);

	if (correct_bias) {
		double mm_bias = ((2 * K - 1) * nr_parts - nr_nonempty_cells) / 2; // note this will also be normalized, then it will make sense
		sum_chi += mm_bias;
		sum_like += mm_bias;
	}

	double normalizer = nr_parts * n; // the *extra* n factor get us to the MI scale.

	sum_chi /= normalizer;
	sum_like /= normalizer;
#endif
}

// This computes ADP ("DDP_ALL")
void StatsComputer::hhg_udf_adp(void) {
	// First, compute the double integral (padded with a leading row and column of zeros)
	compute_double_integral(xy_nrow, dx, y);

	// Now compute the score using all rectangles
	int n = xy_nrow;

	sum_chi  = 0;
	sum_like = 0;

	// These can't be computed this way
	max_chi  = 0;
	max_like = 0;

	int rect_o, xh, yh;
	double cnt, rect_e, rect_c, rect_l;
	//double edenom = 1.0 / (n - K + 1); // I wonder if the denominator here makes sense for the degenerate partitions
	double edenom = 1.0 / n; // This makes much more sense (and is for partitions between ranks rather than at ranks)

	double nr_nonempty_cells = 0;

	// It is more stable numerically to accumulate the statistics in the order
	// of rectangle area, this is why I use the (xl, yl, w, h) enumeration.
	// While this helps stability, it may adversely affect memory locality...
	// Maybe with the Kahan summation, this is no longer necessary?
  
	double kahan_c_chi = 0, kahan_c_like = 0, kahan_t;

	for (int w = 1; w <= n; ++w) {
		for (int h = 1; h <= n; ++h) {
			for (int xl = 1; xl <= n - w + 1; ++xl) {
				for (int yl = 1; yl <= n - h + 1; ++yl) {
					xh = xl + w - 1;
					yh = yl + h - 1;
					cnt = count_adp_with_given_cell(xl, xh, yl, yh);
					if (cnt > 0) {
						// (in theory I could compute this for multiple K's at
						// once with relatively small overhead)
						rect_o = count_sample_points_in_rect(xl, xh, yl, yh);
						rect_e = w * h * edenom;
						rect_c = ((rect_o - rect_e) * (rect_o - rect_e)) / rect_e * cnt - kahan_c_chi;
						rect_l = ((rect_o > 0) ? (rect_o * log(rect_o / rect_e)) : 0) * cnt - kahan_c_like;

						kahan_t = sum_chi + rect_c;
						kahan_c_chi = (kahan_t - sum_chi) - rect_c;
						sum_chi = kahan_t;

						kahan_t = sum_like + rect_l;
						kahan_c_like = (kahan_t - sum_like) - rect_l;
						sum_like = kahan_t;
            
						if (rect_o > 0) {
							nr_nonempty_cells += cnt;
						}
					}
				}
			}
		}
	}

#ifdef UDF_NORMALIZE
	double nr_parts = choose(n - 1, K - 1); // this is actually only the sqrt of the nr parts
	nr_parts *= nr_parts;

	if (correct_bias) {
		double mm_bias = ((2 * K - 1) * nr_parts - nr_nonempty_cells) / 2; // note this will also be normalized, then it will make sense
		sum_chi += mm_bias;
		sum_like += mm_bias;
	}

	double normalizer = nr_parts * n; // gets us to the MI scale
	sum_chi /= normalizer;
	sum_like /= normalizer;
#endif
}

// A variant of the ADP statistic used for CI testing
// NOTE: this is currently NOT a distribution-free test!
void StatsComputer::hhg_ci_udf_adp_mvz_nn(void) {
	int n = xy_nrow;
	int rect_o, xh, yh;
	double cnt, rect_e, rect_c, rect_l;
	double edenom = 1.0 / nnh;

	sum_chi  = 0;
	sum_like = 0;
	max_chi  = 0; // NOTE:
	max_like = 0; // these will not be computed (the fast algorithm for ADP does not support max)

	for (int i = 0; i < n; ++i) {
		// 1. Start by computing the double integral of paired-sample indicators
		// NOTE: see further notes under compute_double_integral
		memset(double_integral, 0, sizeof(int) * dintegral_pn * dintegral_pn);

		// 1.1 Compute x- and y-ranks of the neighbors
		for (int j = 0; j < nnh; ++j) {
			// We want the j'th z-neighbor of the i'th sample, in the bootstrapped data.
			int jj = (*sorted_dz)[i][j + 1].second; // the z's are not bootstrapped
			int pj_x = idx_perm_inv[jj]; // and this gives us the original sample that is the j'th neighbor in the bootstrapped data
			int pj_y = idx_perm    [jj];
			nn_sorted_x[j].first = dx[pj_x];
			nn_sorted_y[j].first = dy[pj_y];
			nn_sorted_x[j].second = j;
			nn_sorted_y[j].second = j;
		}
		sort(nn_sorted_x.begin(), nn_sorted_x.end(), dbl_int_pair_comparator);
		sort(nn_sorted_y.begin(), nn_sorted_y.end(), dbl_int_pair_comparator);
		for (int j = 0; j < nnh; ++j) {
			nn_sorted_x[nn_sorted_x[j].second].first = j + 1;
			nn_sorted_y[nn_sorted_y[j].second].first = j + 1;
		}

		// 1.2. Populate the padded matrix with the indicator variables of whether there
		// is a point in the neighborhood
		for (int j = 0; j < nnh; ++j) {
			int rxj = nn_sorted_x[j].first;
			int ryj = nn_sorted_y[j].first;
			double_integral[ryj * dintegral_pn + rxj] = 1;
		}

		// 1.3. Then run linearly and compute the integral in one row-major pass
		int la = dintegral_pn;
		for (int j = 1; j < dintegral_pn; ++j) {
			int row_running_sum = 0;
			++la;
			for (int k = 1; k < dintegral_pn; ++k) {
				row_running_sum += double_integral[la];
				double_integral[la] = row_running_sum + double_integral[la - dintegral_pn];
				++la;
			}
		}

		// 2. Now use this to compute the ADP statistic for the neighborhood
		// NOTE: see further notes under hhg_udf_adp
		double nn_sum_chi = 0, nn_sum_like = 0;
		double nr_nonempty_cells = 0;
		double kahan_c_chi = 0, kahan_c_like = 0, kahan_t;

		for (int w = 1; w <= nnh; ++w) {
			for (int h = 1; h <= nnh; ++h) {
				for (int xl = 1; xl <= nnh - w + 1; ++xl) {
					for (int yl = 1; yl <= nnh - h + 1; ++yl) {
						xh = xl + w - 1;
						yh = yl + h - 1;
						cnt = count_adp_with_given_cell(xl, xh, yl, yh);
						if (cnt > 0) {
							rect_o = count_sample_points_in_rect(xl, xh, yl, yh);
							rect_e = w * h * edenom;
							rect_c = ((rect_o - rect_e) * (rect_o - rect_e)) / rect_e * cnt - kahan_c_chi;
							rect_l = ((rect_o > 0) ? (rect_o * log(rect_o / rect_e)) : 0) * cnt - kahan_c_like;

							kahan_t = nn_sum_chi + rect_c;
							kahan_c_chi = (kahan_t - nn_sum_chi) - rect_c;
							nn_sum_chi = kahan_t;

							kahan_t = nn_sum_like + rect_l;
							kahan_c_like = (kahan_t - nn_sum_like) - rect_l;
							nn_sum_like = kahan_t;

							if (rect_o > 0) {
								nr_nonempty_cells += cnt;
							}
						}
					}
				}
			}
		}

#ifdef UDF_NORMALIZE
		double nr_parts = choose(nnh - 1, K - 1); // this is actually only the sqrt of the nr parts
		nr_parts *= nr_parts;

		if (correct_bias) {
			double mm_bias = ((2 * K - 1) * nr_parts - nr_nonempty_cells) / 2; // note this will also be normalized, then it will make sense
			nn_sum_chi += mm_bias;
			nn_sum_like += mm_bias;
		}

		double normalizer = nr_parts * nnh; // gets us to the MI scale
		nn_sum_chi /= normalizer;
		nn_sum_like /= normalizer;
#endif

		sum_chi += nn_sum_chi;
		sum_like += nn_sum_like;
	}

#ifdef UDF_NORMALIZE
	sum_chi /= n;
	sum_like /= n;
#endif
}

// Existing K-sample distribution-free tests
// (used in simulations as a reference to compare the to performance of our tests)
void StatsComputer::hhg_k_sample_existing(void) {
	// NOTE: The tests implemented here also exist in versions that use at-ranks grid for
	// inducing dichotomies. This creates the question of how to estimate the
	// empirical CDF at the partitioning points, and some people split each
	// rank and count is as being left on each side of the partition. I avoid
	// this by performing partitions at half-ranks as in ADP and our K-sample
	// test. The effect of this is expected to be negligible.
	// Also, I do not use w_sum/w_max. Maybe that's not fair but the intention
	// is to compare to the "plain vanilla" classical tests.

	int n = xy_nrow;
	compute_single_integral(n, dx, y);

	// will be abused for storing the statistics:
	sum_chi  = 0; // Cramer-von Mises (Kiefer, 1959)
	max_chi  = 0; // Kolmogorov-Smirnov (Kiefer, 1959)
	sum_like = 0; // LR-based Cramer-von Mises (Zhang 2007)
	max_like = 0; // LR-based Kolmogorov-Smirnov (Zhang 2007)

	// FIXME: Could this be done faster?
	// FIXME: may need to address roundoff issues
	// TODO: collect individual tables if wanted

	// NOTE: actually, this is very confusing, since Zhang doesn't simply switch from Pearson chi-squared
	// to the G statistic, rather he also changes the weights when he defines his KS, AD, and CvM

	// NOTE: this uses a "between-ranks" grid, where we partition between xi and xi+1
	// M here is 2

	double Nki, Nk, Eki, dt, d, lr, nrcp = 1.0 / n;

	for (int i = 1; i < n; ++i) {
		d = 0;
		lr = 0;
		for (int k = 0; k < K; ++k) {
			Nk = y_counts[k];
			Eki = Nk * nrcp * i;
			Nki = double_integral[k * dintegral_pn + i];
			dt = Nki - Eki;
			d += dt * dt / Nk;
			lr += (Nki != 0 && Nki != Nk) ? (Nki * log(Nki / Eki) + (Nk - Nki) * log((Nk - Nki) / (Nk - Eki))) : 0;
		}

		sum_chi += d;
		max_chi = max(max_chi, d); // NOTE: sometimes the K-S statistic is defined with an L1 norm and sometimes with L2
		sum_like += lr;
		max_like = max(max_like, lr);
	}
}

// K-sample DDP (also ADP) Kx2 test
void StatsComputer::hhg_k_sample_ddp_m2(void) {
	compute_single_integral(xy_nrow, dx, y);

	int n = xy_nrow;
	double nrmlz = 1.0 / n;

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	ng_chi  = 0;
	ng_like = 0;

	// FIXME: Could this be done faster?
	// TODO: collect individual tables if wanted

	// NOTE: this uses a "between-ranks" grid, where we partition between xi and xi+1
	// M here is 2

	double* tbl_o = new double[K * M];
	double* tbl_e = new double[K * M];
	double chi, like, etmp, emin;
	int yiM;

	for (int xi = 1; xi < n; ++xi) {
		chi = like = 0;
		emin = n;

		for (int yi = 0; yi < K; ++yi) {
			yiM = yi * M;

			tbl_o[yiM + 0] = double_integral[yi * dintegral_pn + xi];
			tbl_o[yiM + 1] = y_counts[yi] - double_integral[yi * dintegral_pn + xi];
			tbl_e[yiM + 0] = y_counts[yi] * double_integral[K * dintegral_pn + xi] * nrmlz;
			tbl_e[yiM + 1] = y_counts[yi] * (n - double_integral[K * dintegral_pn + xi]) * nrmlz;

			chi += (  (tbl_o[yiM    ] - tbl_e[yiM    ]) * (tbl_o[yiM    ] - tbl_e[yiM    ]) / tbl_e[yiM    ]
			        + (tbl_o[yiM + 1] - tbl_e[yiM + 1]) * (tbl_o[yiM + 1] - tbl_e[yiM + 1]) / tbl_e[yiM + 1]);

			if (tbl_o[yiM] > 0) {
				like += tbl_o[yiM] * log(tbl_o[yiM] / tbl_e[yiM]);
			}
			if (tbl_o[yiM + 1] > 0) {
				like += tbl_o[yiM + 1] * log(tbl_o[yiM + 1] / tbl_e[yiM + 1]);
			}

			etmp = min(tbl_e[yiM], tbl_e[yiM + 1]);
			emin = min(emin, etmp);
		}

		accumulate_local_stats(chi, like, emin);
	}

	delete[] tbl_o;
	delete[] tbl_e;

	ng_chi *= n;
	ng_like *= n;
	sum_chi /= ng_chi;
	sum_like /= ng_like;
}

// K-sample DDP (also ADP) Kx3 test
void StatsComputer::hhg_k_sample_ddp_m3(void) {
	compute_single_integral(xy_nrow, dx, y);

	int n = xy_nrow;
	double nrmlz = 1.0 / n;

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	ng_chi  = 0;
	ng_like = 0;

	// FIXME: Could this be done in O(nlogn) time?
	// TODO: collect individual tables if wanted

	// NOTE: this uses a "between-ranks" grid, where we partition between xi and xi+1
	// M here is 3

	double* tbl_o = new double[K * M];
	double* tbl_e = new double[K * M];
	double chi, like, etmp, emin;
	int yiM;

	for (int xi = 1; xi < n - 1; ++xi) {
		for (int xj = xi + 1; xj < n; ++xj) {
			chi = like = 0;
			emin = n;

			for (int yi = 0; yi < K; ++yi) {
				yiM = yi * M;

				tbl_o[yiM + 0] =                 double_integral[yi * dintegral_pn + xi];
				tbl_o[yiM + 1] =                 double_integral[yi * dintegral_pn + xj] - double_integral[yi * dintegral_pn + xi];
				tbl_o[yiM + 2] =                 y_counts[yi]                            - double_integral[yi * dintegral_pn + xj];
				tbl_e[yiM + 0] = y_counts[yi] *  double_integral[K  * dintegral_pn + xi] * nrmlz;
				tbl_e[yiM + 1] = y_counts[yi] * (double_integral[K  * dintegral_pn + xj] - double_integral[K  * dintegral_pn + xi]) * nrmlz;
				tbl_e[yiM + 2] = y_counts[yi] * (n                                       - double_integral[K  * dintegral_pn + xj]) * nrmlz;

				chi += (  (tbl_o[yiM + 0] - tbl_e[yiM + 0]) * (tbl_o[yiM + 0] - tbl_e[yiM + 0]) / tbl_e[yiM + 0]
						+ (tbl_o[yiM + 1] - tbl_e[yiM + 1]) * (tbl_o[yiM + 1] - tbl_e[yiM + 1]) / tbl_e[yiM + 1]
						+ (tbl_o[yiM + 2] - tbl_e[yiM + 2]) * (tbl_o[yiM + 2] - tbl_e[yiM + 2]) / tbl_e[yiM + 2]);

				if (tbl_o[yiM + 0] > 0) {
					like += tbl_o[yiM + 0] * log(tbl_o[yiM + 0] / tbl_e[yiM + 0]);
				}
				if (tbl_o[yiM + 1] > 0) {
					like += tbl_o[yiM + 1] * log(tbl_o[yiM + 1] / tbl_e[yiM + 1]);
				}
				if (tbl_o[yiM + 2] > 0) {
					like += tbl_o[yiM + 2] * log(tbl_o[yiM + 2] / tbl_e[yiM + 2]);
				}

				etmp = min(tbl_e[yi * 2 + 0], tbl_e[yi * 2 + 1]);
				etmp = min(etmp, tbl_e[yi * 2 + 2]);
				emin = min(emin, etmp);
			}

			accumulate_local_stats(chi, like, emin);
		}
	}

	delete[] tbl_o;
	delete[] tbl_e;

	sum_chi /= (double(n) * ng_chi);
	sum_like /= (double(n) * ng_like);
}

// K-sample DDP (also ADP) KxM test
void StatsComputer::hhg_k_sample_ddp(void) {
	compute_single_integral(xy_nrow, dx, y);

	int n = xy_nrow;
	double nrmlz = 1.0 / n;

	sum_chi  = 0;
	sum_like = 0;

	// These can't be computed this way
	max_chi  = 0;
	max_like = 0;

	// TODO: collect individual tables if wanted

	// NOTE: this uses a "between-ranks" grid, where we partition between xi and xi+1
	// This algorithm works interval-wise rather than partition-wise.

	double normalized_cnt;
	double interval_o, interval_e, interval_c, interval_l;
	double kahan_c_chi = 0, kahan_c_like = 0, kahan_t;
	int xj, wmax;

	for (int xi = 0; xi < n; ++xi) {
		wmax = min(n - M - 1, n - xi);

		for (int w = 1; w <= wmax; ++w) {
			xj = xi + w;

			// TODO could be slightly optimized by going over edge intervals in a separate loop from mid intervals
			if (xi == 0 || xj == n) {
				// edge interval
				normalized_cnt = adp_l[w];
			} else {
				// mid interval
				normalized_cnt = adp[w];
			}

			for (int yi = 0; yi < K; ++yi) {
				interval_o =                 double_integral[yi * dintegral_pn + xj] - double_integral[yi * dintegral_pn + xi];
				interval_e = y_counts[yi] * (double_integral[K  * dintegral_pn + xj] - double_integral[K  * dintegral_pn + xi]) * nrmlz;

				interval_c = ((interval_o - interval_e) * (interval_o - interval_e)) / interval_e * normalized_cnt - kahan_c_chi;
				interval_l = ((interval_o > 0) ? (interval_o * log(interval_o / interval_e)) : 0) * normalized_cnt - kahan_c_like;

				kahan_t = sum_chi + interval_c;
				kahan_c_chi = (kahan_t - sum_chi) - interval_c;
				sum_chi = kahan_t;

				kahan_t = sum_like + interval_l;
				kahan_c_like = (kahan_t - sum_like) - interval_l;
				sum_like = kahan_t;
			}
		}
	}

	// the 1/n gets us to the MI scale
	sum_chi /= n;
	sum_like /= n;
}

// Existing GOF tests
// (used in simulations as a reference to compare the to performance of our tests)
void StatsComputer::hhg_gof_existing(void) {
	// See general notes under hhg_gof_ddp() and hhg_k_sample_existing().

	int n = xy_nrow;

	// will be abused for storing the statistics:
	sum_chi  = 0; // Cramer-von Mises (1928)
	max_chi  = 0; // Kolmogorov-Smirnov (1933)
	sum_like = 0; // LR-based Cramer-von Mises (Zhang 2006)
	max_like = 0; // LR-based Kolmogorov-Smirnov (Zhang 2006)

	// FIXME: Could this be done faster?
	// TODO: collect individual tables if wanted

	double Eki, dt, d, lr;

	for (int i = 1; i < n; ++i) {
		Eki = n * null_dist[i];
		dt = i - Eki;
		d = dt * dt / n;
		lr = (i != 0 && i != n) ? (i * log(i / Eki) + (n - i) * log((n - i) / (n - Eki))) : 0;

		sum_chi += d;
		max_chi = max(max_chi, d); // NOTE: sometimes the K-S statistic is defined with an L1 norm and sometimes with L2
		sum_like += lr;
		max_like = max(max_like, lr);
	}
}

// GOF DDP (also ADP) 1x2 test
void StatsComputer::hhg_gof_ddp_m2(void) {
	// See general notes under hhg_gof_ddp()

	int n = xy_nrow;

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	ng_chi  = 0;
	ng_like = 0;

	// FIXME: Could this be done faster?
	// TODO: collect individual tables if wanted

	// NOTE: this uses a "between-ranks" grid, where we partition between xi and xi+1
	// M here is 2

	double* tbl_o = new double[M];
	double* tbl_e = new double[M];
	double chi, like, emin;

	for (int xi = 1; xi < n; ++xi) {
		tbl_o[0] = xi;
		tbl_o[1] = n - xi;
		tbl_e[0] = n * null_dist[xi];
		tbl_e[1] = n * (1 - null_dist[xi]);

		chi = (  (tbl_o[0] - tbl_e[0]) * (tbl_o[0] - tbl_e[0]) / tbl_e[0]
			   + (tbl_o[1] - tbl_e[1]) * (tbl_o[1] - tbl_e[1]) / tbl_e[1]);

		like = 0;
		if (tbl_o[0] > 0) {
			like += tbl_o[0] * log(tbl_o[0] / tbl_e[0]);
		}
		if (tbl_o[1] > 0) {
			like += tbl_o[1] * log(tbl_o[1] / tbl_e[1]);
		}

		emin = min(tbl_e[0], tbl_e[1]);

		accumulate_local_stats(chi, like, emin);
	}

	delete[] tbl_o;
	delete[] tbl_e;

	ng_chi *= n;
	ng_like *= n;
	sum_chi /= ng_chi;
	sum_like /= ng_like;
}

// GOF DDP (also ADP) 1x3 test
void StatsComputer::hhg_gof_ddp_m3(void) {
	// See general notes under hhg_gof_ddp()

	int n = xy_nrow;

	sum_chi  = 0;
	max_chi  = 0;
	sum_like = 0;
	max_like = 0;

	ng_chi  = 0;
	ng_like = 0;

	// FIXME: Could this be done in O(nlogn) time?
	// TODO: collect individual tables if wanted

	// NOTE: this uses a "between-ranks" grid, where we partition between xi and xi+1
	// M here is 3

	double* tbl_o = new double[M];
	double* tbl_e = new double[M];
	double chi, like, etmp, emin;

	for (int xi = 1; xi < n - 1; ++xi) {
		for (int xj = xi + 1; xj < n; ++xj) {
			tbl_o[0] = xi;
			tbl_o[1] = xj - xi;
			tbl_o[2] = n - xj;
			tbl_e[0] = n * null_dist[xi];
			tbl_e[1] = n * (null_dist[xj] - null_dist[xi]);
			tbl_e[2] = n * (1 - null_dist[xj]);

			chi = (  (tbl_o[0] - tbl_e[0]) * (tbl_o[0] - tbl_e[0]) / tbl_e[0]
				   + (tbl_o[1] - tbl_e[1]) * (tbl_o[1] - tbl_e[1]) / tbl_e[1]
				   + (tbl_o[2] - tbl_e[2]) * (tbl_o[2] - tbl_e[2]) / tbl_e[2]);

			like = 0;
			if (tbl_o[0] > 0) {
				like += tbl_o[0] * log(tbl_o[0] / tbl_e[0]);
			}
			if (tbl_o[1] > 0) {
				like += tbl_o[1] * log(tbl_o[1] / tbl_e[1]);
			}
			if (tbl_o[2] > 0) {
				like += tbl_o[2] * log(tbl_o[2] / tbl_e[2]);
			}

			etmp = min(tbl_e[0], tbl_e[1]);
			emin = min(etmp, tbl_e[2]);

			accumulate_local_stats(chi, like, emin);
		}
	}

	delete[] tbl_o;
	delete[] tbl_e;

	sum_chi /= (double(n) * ng_chi);
	sum_like /= (double(n) * ng_like);
}

// GOF DDP (also ADP) 1xM test against a given null distribution
void StatsComputer::hhg_gof_ddp(void) {
	// It is assumed that y contains the null CDF, computed at the midpoints
	// between every two consecutive observations (so there are n-1 given
	// values between 0 and 1)

	int n = xy_nrow;

	sum_chi  = 0;
	sum_like = 0;

	// These can't be computed this way
	max_chi  = 0;
	max_like = 0;

	// TODO: collect individual tables if wanted

	// NOTE: this uses a "between-ranks" grid, where we partition between xi and xi+1
	// For GOF testing, we could opt to partition at any point in the interval
	// between every two consecutive observations, and this choice actually matters.
	// The choice of the midpoints is arbitrary, but will do for now.

	// This algorithm works interval-wise rather than partition-wise.

	double normalized_cnt;
	double interval_o, interval_e, interval_c, interval_l;
	double kahan_c_chi = 0, kahan_c_like = 0, kahan_t;
	int xj, wmax;

	for (int xi = 0; xi < n; ++xi) {
		wmax = min(n - M - 1, n - xi);

		for (int w = 1; w <= wmax; ++w) {
			xj = xi + w;

			// TODO could be slightly optimized by going over edge intervals in a separate loop from mid intervals
			if (xi == 0 || xj == n) {
				// edge interval
				normalized_cnt = adp_l[w];
			} else {
				// mid interval
				normalized_cnt = adp[w];
			}

			interval_o = xj - xi;
			interval_e = (null_dist[xj] - null_dist[xi]) * n;

			interval_c = ((interval_o - interval_e) * (interval_o - interval_e)) / interval_e * normalized_cnt - kahan_c_chi;
			interval_l = ((interval_o > 0) ? (interval_o * log(interval_o / interval_e)) : 0) * normalized_cnt - kahan_c_like;

			kahan_t = sum_chi + interval_c;
			kahan_c_chi = (kahan_t - sum_chi) - interval_c;
			sum_chi = kahan_t;

			kahan_t = sum_like + interval_l;
			kahan_c_like = (kahan_t - sum_like) - interval_l;
			sum_like = kahan_t;
		}
	}

	// the 1/n gets us to the MI scale
	sum_chi /= n;
	sum_like /= n;
}

// Some other statistics that we compare ourselves against
// ================================================================================================

void StatsComputer::other_stats_two_sample(void) {
	// This representation of the data simplifies the computations of edist and HT
	int j0 = 0, j1 = 0;
	for (int i = 0; i < xy_nrow; ++i) {
		if (y[i] == 0) {
			y0_idx[j0++] = i;
		} else {
			y1_idx[j1++] = i;
		}
	}

	// Compute edist (a version of distance correlation for the 2-sample problem)
	compute_edist();

	// Compute Hall & Tajvidi (2002) test
	compute_ht();
}

void StatsComputer::other_stats_k_sample(void) {
	// TODO
	edist = 0;
}

void StatsComputer::other_stats_general(void) {
	// TODO
	edist = 0;
}

void StatsComputer::other_stats_univar_dist_free(void) {
	// TODO, though I doubt this is meaningful
	edist = 0;
}

void StatsComputer::other_stats_ci(void) {
}

// Helper functions
// ================================================================================================

void StatsComputer::accumulate_2x2_contingency_table(double a00, double a01, double a10, double a11, double nrmlz, double reps) {
	double e00, e01, e10, e11, emin00_01, emin10_11, emin, current_chi, current_like;

	e00 = ((a00 + a01) * (a00 + a10)) * nrmlz;
	e01 = ((a00 + a01) * (a01 + a11)) * nrmlz;
	e10 = ((a10 + a11) * (a00 + a10)) * nrmlz;
	e11 = ((a10 + a11) * (a01 + a11)) * nrmlz;

	emin00_01 = min(e00, e01);
	emin10_11 = min(e10, e11);
	emin = min(emin00_01, emin10_11);

	if (emin > min_w) {
		current_chi = (a00 - e00) * (a00 - e00) / e00
					+ (a01 - e01) * (a01 - e01) / e01
					+ (a10 - e10) * (a10 - e10) / e10
					+ (a11 - e11) * (a11 - e11) / e11;
	} else {
		current_chi = 0;
	}

	if (emin > w_sum) {
		sum_chi += current_chi * reps;
	}

	if ((emin > w_max) && (current_chi > max_chi)) {
		max_chi = current_chi;
	}

	current_like = ((a00 > 0) ? (a00 * log(a00 / e00)) : 0)
				 + ((a01 > 0) ? (a01 * log(a01 / e01)) : 0)
				 + ((a10 > 0) ? (a10 * log(a10 / e10)) : 0)
				 + ((a11 > 0) ? (a11 * log(a11 / e11)) : 0);

	sum_like += current_like * reps;

	if (current_like > max_like) {
		max_like = current_like;
	}
}

void StatsComputer::accumulate_local_stats(double chi, double like, double emin) {
	double kahan_t, kahan_chi, kahan_like;

	if (emin > w_sum) {
		kahan_chi = chi - kahan_c_chi;
		kahan_t = sum_chi + kahan_chi;
		kahan_c_chi = (kahan_t - sum_chi) - kahan_chi;
		sum_chi = kahan_t;
		++ng_chi;
	}

	if ((emin > w_max) && (chi > max_chi)) {
		max_chi = chi;
	}

	kahan_like = like - kahan_c_like;
	kahan_t = sum_like + kahan_like;
	kahan_c_like = (kahan_t - sum_like) - kahan_like;
	sum_like = kahan_t;
	++ng_like;

	if (like > max_like) {
		max_like = like;
	}
}

void StatsComputer::hhg_gen_inversions(int *permutation, int *source, int *inversion_count, int dim) {
    if (dim <= 1) {
        return;
    } else {
    	hhg_gen_inversions(permutation, source, inversion_count, dim / 2);
    	hhg_gen_inversions(&permutation[dim / 2], &source[dim / 2], inversion_count, dim - dim / 2);
    	hhg_gen_merge(permutation, source, inversion_count, dim);
    }
}

void StatsComputer::hhg_gen_merge(int *permutation, int *source, int *inversion_count, int dim) {
    int left_index = 0, right_index = 0;
    int i, half_dim = dim / 2, nleft = half_dim, nright = dim - half_dim;

    int* left = hhg_gen_left_buffer;
    int* right = hhg_gen_right_buffer;
    int* left_source = hhg_gen_left_source_buffer;
    int* right_source = hhg_gen_right_source_buffer;

    for (i = 0; i < half_dim; i++) {
        left[i] = permutation[i];
        left_source[i] = source[i];
        right[i] = permutation[i + half_dim];
        right_source[i] = source[i + half_dim];
    }

    if (nleft < nright) {
        right[i] = permutation[i + half_dim];
        right_source[i] = source[i + half_dim];
    }

    for (i = 0; i < dim; i++) {
        if ((left_index < half_dim) && (right_index < dim - half_dim)) {
             if (left[left_index] <= right[right_index]) { // I added "=" in order to support ties
                permutation[i] = left[left_index];
                source[i] = left_source[left_index];
                left_index++;
            } else {
                permutation[i] = right[right_index];
                source[i] = right_source[right_index];
                inversion_count[source[i]] += (half_dim - left_index);
                right_index++;
            }
        } else {
            if (left_index < half_dim) {
                permutation[i] = left[left_index];
                source[i] = left_source[left_index];
                left_index++;
            }

            if (right_index < dim - half_dim) {
                permutation[i] = right[right_index];
                source[i] = right_source[right_index];
                right_index++;
            }
        }
    }
}

void StatsComputer::compute_single_integral(int n, double* xx, double* yy) {
	memset(double_integral, 0, sizeof(int) * (K + 1) * dintegral_pn);

	// Populate the padded matrix with the indicator variables of whether there
	// is a point in the set {(x_k, y_k)}_k=1^n in the grid slot (i, j) for i in 1,2,...,n and j in 1,2,...,K
	// NOTE: the last row holds the integral over all x, ignoring y
	for (int i = 0; i < n; ++i) {
		int yi = yy[i]; // assumed in 0..K-1, many ties, and won't be padded
		int xi = xx[i]; // assumed in 1..n, no ties, and will be padded
		double_integral[yi * dintegral_pn + xi] = 1;
		double_integral[K  * dintegral_pn + xi] = 1;
	}

	// Then run linearly and compute the integral row by row
	for (int k = 0; k <= K; ++k) {
		int row_running_sum = 0;
		for (int i = 1; i < dintegral_pn; ++i) {
			row_running_sum += double_integral[k * dintegral_pn + i];
			double_integral[k * dintegral_pn + i] = row_running_sum;
		}
	}
}

void StatsComputer::compute_double_integral(int n, double* xx, double* yy) {
	memset(double_integral, 0, sizeof(int) * dintegral_pn * dintegral_pn);

	// Populate the padded matrix with the indicator variables of whether there
	// is a point in the set {(x_k, y_k)}_k=1^n in the grid slot (i, j) for i,j in 1,2,...,n
	for (int i = 0; i < n; ++i) {
		int yi = yy[i] + dintegral_zero_based_idxs;
		int xi = xx[i] + dintegral_zero_based_idxs;
		double_integral[yi * dintegral_pn + xi] = 1;
	}

	// Then run linearly and compute the integral in one row-major pass
	int la = dintegral_pn;
	for (int i = 1; i < dintegral_pn; ++i) {
		int row_running_sum = 0;
		++la;
		for (int j = 1; j < dintegral_pn; ++j) {
			row_running_sum += double_integral[la];
			double_integral[la] = row_running_sum + double_integral[la - dintegral_pn];
			++la;
		}
	}
}

int StatsComputer::count_sample_points_in_rect(int xl, int xh, int yl, int yh) {
	return (  double_integral[ yh      * dintegral_pn + xh    ]
	        - double_integral[ yh      * dintegral_pn + xl - 1]
	        - double_integral[(yl - 1) * dintegral_pn + xh    ]
	        + double_integral[(yl - 1) * dintegral_pn + xl - 1]);
}

// This counts the number of data driven partitions of the data, that contain the given rectangle as a cell.
double StatsComputer::count_ddp_with_given_cell(int xl, int xh, int yl, int yh) {
	int n = xy_nrow, m;

	// This may not be necessary in practical usage scenarios
	if ((xl == 1 && xh == n) || (yl == 1 && yh == n)) {
		return (0);
	}

	if (xl == 1 && yl == 1) {
		// ll corner cell
		if ((y_ordered_by_x[xh] > yh) && (x_ordered_by_y[yh] > xh)) {
			m = count_sample_points_in_rect(xh + 2, n, yh + 2, n);

			if (x_ordered_by_y[yh] == xh + 1) {
				// The hh point is a sample point, it has to be chosen,
				// and we must choose the remaining K-2 points from m.
				return (choose(m, K - 2));
			} else {
				// We have to choose both the sample point with x = xh + 1, and the point
				// with y = yh + 1, which do not split the rectangle, and the remaining
				// K-3 points are chosen from m.
				return (choose(m, K - 3));
			}
		}
	} else if (xl == 1 && yh == n) {
		// lh corner cell
		if ((y_ordered_by_x[xh] < yl) && (x_ordered_by_y[yl - 2] > xh)) {
			m = count_sample_points_in_rect(xh + 2, n, 1, yl - 2);

			if (x_ordered_by_y[yl - 2] == xh + 1) {
				// The hl point is a sample point, it has to be chosen,
				// and we must choose the remaining K-2 points from m.
				return (choose(m, K - 2));
			} else {
				// We have to choose both the sample point with x = xh + 1, and the point
				// with y = yl - 1, which do not split the rectangle, and the remaining K-3
				// points are chosen from m.
				return (choose(m, K - 3));
			}
		}
	} else if (xh == n && yl == 1) {
		// hl corner cell
		if ((y_ordered_by_x[xl - 2] > yh) && (x_ordered_by_y[yh] < xl)) {
			m = count_sample_points_in_rect(1, xl - 2, yh + 2, n);

			if (x_ordered_by_y[yh] == xl - 1) {
				// The lh point is a sample point, it has to be chosen,
				// and we must choose the remaining K-2 points from m.
				return (choose(m, K - 2));
			} else {
				// We have to choose both the sample point with x = xl - 1, and the point
				// with y = yh + 1, which do not split the rectangle, and the remaining K-3
				// points are chosen from m.
				return (choose(m, K - 3));
			}
		}
	} else if (xh == n && yh == n) {
		// hh corner cell
		if ((y_ordered_by_x[xl - 2] < yl) && (x_ordered_by_y[yl - 2] < xl)) {
			m = count_sample_points_in_rect(1, xl - 2, 1, yl - 2);

			if (x_ordered_by_y[yl - 2] == xl - 1) {
				// The ll point is a sample point, it has to be chosen,
				// and we must choose the remaining K-2 points from m.
				return (choose(m, K - 2));
			} else {
				// We have to choose both the sample point with x = xl - 1, and the point
				// with y = yl - 1, which do not split the rectangle, and the remaining K-3
				// points are chosen from m.
				return (choose(m, K - 3));
			}
		}
	} else if (yl == 1) {
		// bottom edge cell
		if (  (x_ordered_by_y[yh    ] < xl || x_ordered_by_y[yh] > xh)
			&& y_ordered_by_x[xl - 2] > yh && y_ordered_by_x[xh] > yh)
		{
			m = count_sample_points_in_rect(1, xl - 2, yh + 2, n)
			  + count_sample_points_in_rect(xh + 2, n, yh + 2, n);

			if ((x_ordered_by_y[yh] == xl - 1) || (x_ordered_by_y[yh] == xh + 1)) {
				// One of the top corners is in the sample, so this point has to be
				// chosen, so does the point with the remaining x boundary, and we need to
				// further select the remaining K-3 points from m.
				return (choose(m, K - 3));
			} else {
				// We must choose the point with y = yh + 1, the point with x = xl - 1, and
				// the point with x = xh + 1, and the remaining K-4 are chosen from m.
				return (choose(m, K - 4));
			}
		}
	} else if (yh == n) {
		// top edge cell
		if (  (x_ordered_by_y[yl - 2] < xl || x_ordered_by_y[yl - 2] > xh)
			&& y_ordered_by_x[xl - 2] < yl && y_ordered_by_x[xh    ] < yl) {
			m = count_sample_points_in_rect(1, xl - 2, 1, yl - 2)
			  + count_sample_points_in_rect(xh + 2, n, 1, yl - 2);

			if ((x_ordered_by_y[yl - 2] == xl - 1) || (x_ordered_by_y[yl - 2] == xh + 1)) {
				// One of the bottom corners are in the sample, so this point has to be
				// chosen, so does the point with the remaining x boundary, and we need to
				// further select the remaining K-3 points from m.
				return (choose(m, K - 3));
			} else {
				// We must choose the point with y = yh + 1, the point with x = xl - 1, and
				// the point with x = xh + 1, and the remaining K-4 are chosen from m.
				return (choose(m, K - 4));
			}
		}
	} else if (xl == 1) {
		// left edge cell
		if (  (y_ordered_by_x[xh    ] < yl || y_ordered_by_x[xh] > yh)
			&& x_ordered_by_y[yl - 2] > xh && x_ordered_by_y[yh] > xh) {
			m = count_sample_points_in_rect(xh + 2, n, 1, yl - 2)
			  + count_sample_points_in_rect(xh + 2, n, yh + 2, n);

			if ((y_ordered_by_x[xh] == yl - 1) || (y_ordered_by_x[xh] == yh + 1)) {
				// One of the right corners are in the sample, so this point has to be
				// chosen, so does the point with the remaining y boundary, and we need to
				// further select the remaining K-3 points from m.
				return (choose(m, K - 3));
			} else {
				// We must choose the point with x = xh + 1, the point with y = yl - 1, and
				// the point with y = yh + 1, and the remaining K-4 are chosen from m.
				return (choose(m, K - 4));
			}
		}
	} else if (xh == n) {
		// right edge cell
		if (  (y_ordered_by_x[xl - 2] < yl || y_ordered_by_x[xl - 2] > yh)
			&& x_ordered_by_y[yl - 2] < xl && x_ordered_by_y[yh    ] < xl)
		{
			m = count_sample_points_in_rect(1, xl - 2, 1, yl - 2)
			  + count_sample_points_in_rect(1, xl - 2, yh + 2, n);

			if ((y_ordered_by_x[xl - 2] == yl - 1) || (y_ordered_by_x[xl - 2] == yh + 1)) {
				// One of the left corners are in the sample, so this point has to be
				// chosen, so does the point with the remaining y boundary, and we need to
				// further select the remaining K-3 points from m.
				return (choose(m, K - 3));
			} else {
				// We must choose the point with x = xl - 1, the point with y = yl - 1, and
				// the point with y = yh + 1, and the remaining K-4 are chosen from m.
				return (choose(m, K - 4));
			}
		}
	} else {
		// inner cell
		if (   (y_ordered_by_x[xl - 2] < yl || y_ordered_by_x[xl - 2] > yh)
			&& (y_ordered_by_x[xh    ] < yl || y_ordered_by_x[xh    ] > yh)
			&& (x_ordered_by_y[yl - 2] < xl || x_ordered_by_y[yl - 2] > xh)
			&& (x_ordered_by_y[yh    ] < xl || x_ordered_by_y[yh    ] > xh))
		{
			m = count_sample_points_in_rect(1, xl - 2, 1, yl - 2)
			  + count_sample_points_in_rect(1, xl - 2, yh + 2, n)
			  + count_sample_points_in_rect(xh + 2, n, 1, yl - 2)
			  + count_sample_points_in_rect(xh + 2, n, yh + 2, n);

			if (   ((x_ordered_by_y[yl - 2] == xl - 1) && (x_ordered_by_y[yh] == xh + 1))
				|| ((x_ordered_by_y[yl - 2] == xh + 1) && (x_ordered_by_y[yh] == xl - 1))) {
				// Two opposing corners are sample points. These must be chosen, and the
				// rest K-3 points are chosen from m.
				return (choose(m, K - 3));
			} else if (   (x_ordered_by_y[yl - 2] == xl - 1) || (x_ordered_by_y[yl - 2] == xh + 1)
					   || (x_ordered_by_y[yh    ] == xl - 1) || (x_ordered_by_y[yh    ] == xh + 1)) {
				// Exactly one corner is a sample point and has to be chosen. We also have
				// to choose the two points with x and y coordinates of the opposing
				// corner. The remaining K-4 points are chosen from m.
				return (choose(m, K - 4));
			} else {
				// No corner is a sample point, so we must choose the four points with
				// the bounding x's and y's (which we can). The remaining K-5 points are
				// chosen from m.
				return (choose(m, K - 5));
			}
		}
	}

	// If we reached here, then the specified rectangle can never be a cell in any
	// data driven partition of the given dataset
	return (0);
}

// This counts the number of unrestricted partitions of the data.
double StatsComputer::count_adp_with_given_cell(int xl, int xh, int yl, int yh) {
#if 0
	// FIXME this is quite silly, and I can probably speed this up considerably
	// by blocking in a way that takes cache hierarchy into account.
	return (adp[(xl - 1) * xy_nrow + xh - 1] * adp[(yl - 1) * xy_nrow + yh - 1]);
#else
	double cx, cy;

	if (xl == 1) {
		// left anchored interval
		cx = adp_l[xh - 1];
	} else if (xh == xy_nrow) {
		// right anchored interval
		cx = adp_r[xl - 1];
	} else {
		cx = adp[xh - xl];
	}

	if (yl == 1) {
		// left anchored interval
		cy = adp_l[yh - 1];
	} else if (yh == xy_nrow) {
		// right anchored interval
		cy = adp_r[yl - 1];
	} else {
		cy = adp[yh - yl];
	}

	return (cx * cy);
#endif
}

void StatsComputer::precompute_adp(void) {
	int n = xy_nrow;

#if 0 // partition at ranks
	for (int xl = 1; xl <= n; ++xl) {
		for (int xh = xl; xh <= n; ++xh) {
			int idx = (xl - 1) * xy_nrow + xh - 1;

			if (xl == 1) {
				// left anchored interval
				adp[idx] = my_choose(n - xh - 2 - (K - 2), K - 2);
			} else if (xh == n) {
				// right anchored interval
				adp[idx] = my_choose(xl - 3 - (K - 2), K - 2);
			} else if (xl == 2 || xh == n - 1 || K == 2) {
				adp[idx] = 0;
			} else if (K == 3) {
				adp[idx] = 1;
			} else {
				// nondegenerate regular interval
				adp[idx] = 0;
				for (int i = 0; i <= K - 3; ++i) {
					adp[idx] += my_choose(xl - 3 - i, i) * my_choose(n - xh - 2 - (K - 3 - i), K - 3 - i);
				}
			}
		}
	}
#else
	// left anchored interval (xl == 1)
	for (int xh = 1; xh <= n; ++xh) {
		adp_l[xh - 1] = my_choose(n - xh - 1, K - 2); // (the xh == n case is irrelevant)
	}

	// right anchored interval (xh == n)
	for (int xl = 1; xl <= n; ++xl) {
		adp_r[xl - 1] = my_choose(xl - 2, K - 2); // (the xl == 1 case is irrelevant)
	}

	// nondegenerate regular interval
	for (int xd = 0; xd < n; ++xd) {
		adp[xd] = my_choose(n - xd - 3, K - 3); // (the xd > n - 3 case is irrelevant)
	}
#endif
}

void StatsComputer::precompute_adp_k_sample(void) {
	int n = xy_nrow;

	// NOTE: in this test we would like to be able to use large samples, and the involved
	// binomial coefficients can get much larger than the DOUBLE_XMAX. We thus need to
	// take extra care to use logarithm binomial coefficients.
	double log_denom = lchoose(n - 1, M - 1);

	// edge-anchored interval (xi == 0)
	for (int w = 1; w < n; ++w) {
		adp_l[w] = exp(lchoose(n - w - 1, M - 2) - log_denom);
	}

	// mid interval
	for (int w = 1; w <= n - 3; ++w) {
		adp[w] = exp(lchoose(n - w - 3, M - 3) - log_denom);
	}
}

double StatsComputer::my_choose(int n, int k) {
    if (n < 0) {
        return (0);
    }
	return (choose(n, k));
}

void StatsComputer::compute_ht(void) {
	int n = xy_nrow;
    int x_c_i, total_other_y_count_so_far;
	double expected_total_other_y_count_so_far, nm1i = 1.0 / (n - 1);
	double ht0 = 0, ht1 = 0, tmp;
    int i, j;
	int n0 = y_counts[0];
	int n1 = y_counts[1];

	// NOTE: I'll compute the HT test in the version they refer to as "T", with
	// "gamma" being 2, and with "w" being all ones. They also have an "S"
	// version which is comparable maybe with the HHG "max" variants.
	// Also, for the moment, ties are not handled.

	// We already have the y1_idx and y0_idx from above, which define the
	// permuted partitioning of "Z" to "X" and "Y". We also already have sorted
	// dx's (these are the distances between the pooled sample "Z").

    for (i = 0; i < n0; i++) {
        x_c_i = y0_idx[i];
        total_other_y_count_so_far = 0; // this is "Mij" in the HT paper

    	for (j = 0; j < n1; j++) {
    		total_other_y_count_so_far += (y[(*sorted_dx)[x_c_i][j].second] == 1);
    		expected_total_other_y_count_so_far = n1 * j * nm1i;
    		tmp = (total_other_y_count_so_far - expected_total_other_y_count_so_far);
    		ht0 += tmp * tmp;
		}
    }

    for (i = 0; i < n1; i++) {
        x_c_i = y1_idx[i];
        total_other_y_count_so_far = 0; // this is "Nij" in the HT paper

    	for (j = 0; j < n0; j++) {
    		total_other_y_count_so_far += (y[(*sorted_dx)[x_c_i][j].second] == 0);
    		expected_total_other_y_count_so_far = n0 * j * nm1i;
    		tmp = (total_other_y_count_so_far - expected_total_other_y_count_so_far);
    		ht1 += tmp * tmp;
		}
    }

    ht = ht0 / n0 + ht1 / n1;
}

void StatsComputer::compute_edist(void) {
	int n0 = y_counts[0];
	int n1 = y_counts[1];

	double sum01 = 0;
	double sum00 = 0;
	double sum11 = 0;

	for (int i = 0; i < n0; ++i) {
		for (int j = 0; j < n1; ++j) {
			sum01 += dx[y1_idx[j] * xy_nrow + y0_idx[i]];
		}
	}

	for (int i = 0; i < n0; ++i) {
		for (int j = 0; j < n0; ++j) {
			sum00 += dx[y0_idx[j] * xy_nrow + y0_idx[i]];
		}
	}

	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n1; ++j) {
			sum11 += dx[y1_idx[j] * xy_nrow + y1_idx[i]];
		}
	}

	edist = (2.0 / (n0 * n1)) * sum01 - (1.0 / (n0 * n0)) * sum00 - (1.0 / (n1 * n1)) * sum11;
	edist *= ((double)(n0 * n1)) / (n0 + n1);
}

void StatsComputer::compute_spr_obs(int xi, int yi, int n, int pn, int nm1, double nm1d) {
	int a11, a12, a21, a22;
	double e11, e12, e21, e22, current_chi, current_like;
	double emin11_12, emin21_22, emin;

	a11 = double_integral[n  * pn + xi] - double_integral[(yi + 1) * pn + xi    ];
	a12 = double_integral[n  * pn + n ] - double_integral[ n       * pn + xi + 1] - double_integral[(yi + 1) * pn + n] + double_integral[(yi + 1) * pn + xi + 1];
	a21 = double_integral[yi * pn + xi];
	a22 = double_integral[yi * pn + n ] - double_integral[yi * pn + xi + 1];

	e11 =        xi  * (nm1 - yi) / nm1d;
	e12 = (nm1 - xi) * (nm1 - yi) / nm1d;
	e21 =        xi  *        yi  / nm1d;
	e22 = (nm1 - xi) *        yi  / nm1d;

	emin11_12 = min(e11, e12);
	emin21_22 = min(e21, e22);
	emin = min(emin11_12, emin21_22);

	if (emin > min_w) {
#ifdef UDF_ALLOW_DEGENERATE_PARTITIONS
		current_chi = ((e11 > 0) ? ((a11 - e11) * (a11 - e11) / e11) : 0)
					+ ((e21 > 0) ? ((a21 - e21) * (a21 - e21) / e21) : 0)
					+ ((e12 > 0) ? ((a12 - e12) * (a12 - e12) / e12) : 0)
					+ ((e22 > 0) ? ((a22 - e22) * (a22 - e22) / e22) : 0);
#else
		current_chi = (a11 - e11) * (a11 - e11) / e11
					+ (a21 - e21) * (a21 - e21) / e21
					+ (a12 - e12) * (a12 - e12) / e12
					+ (a22 - e22) * (a22 - e22) / e22;
#endif
	} else {
		current_chi = 0;
	}

	if (emin > w_sum) {
		sum_chi += current_chi;
		++ng_chi;
	}

	if ((emin > w_max) && (current_chi > max_chi)) {
		max_chi = current_chi;
	}

	current_like = ((a11 > 0) ? (a11 * log(a11 / e11)) : 0)
	             + ((a21 > 0) ? (a21 * log(a21 / e21)) : 0)
	             + ((a12 > 0) ? (a12 * log(a12 / e12)) : 0)
	             + ((a22 > 0) ? (a22 * log(a22 / e22)) : 0);

	sum_like += current_like;
	++ng_like;
	if (current_like > max_like) {
		max_like = current_like;
	}
}

void StatsComputer::compute_spr_all(int xi, int yi, int n, int pn, double nd) {
	int a11, a12, a21, a22;
	double e11, e12, e21, e22, current_chi, current_like;
	double emin11_12, emin21_22, emin;

	a11 = double_integral[n  * pn + xi] - double_integral[yi * pn + xi];
	a12 = double_integral[n  * pn + n ] - double_integral[ n * pn + xi] - double_integral[yi * pn + n] + double_integral[yi * pn + xi];
	a21 = double_integral[yi * pn + xi];
	a22 = double_integral[yi * pn + n ] - double_integral[yi * pn + xi];

	e11 =      xi  * (n - yi) / nd;
	e12 = (n - xi) * (n - yi) / nd;
	e21 =      xi  *      yi  / nd;
	e22 = (n - xi) *      yi  / nd;

	emin11_12 = min(e11, e12);
	emin21_22 = min(e21, e22);
	emin = min(emin11_12, emin21_22);

	if (emin > min_w) {
#ifdef UDF_ALLOW_DEGENERATE_PARTITIONS
		current_chi = ((e11 > 0) ? ((a11 - e11) * (a11 - e11) / e11) : 0)
					+ ((e21 > 0) ? ((a21 - e21) * (a21 - e21) / e21) : 0)
					+ ((e12 > 0) ? ((a12 - e12) * (a12 - e12) / e12) : 0)
					+ ((e22 > 0) ? ((a22 - e22) * (a22 - e22) / e22) : 0);
#else
		current_chi = (a11 - e11) * (a11 - e11) / e11
					+ (a21 - e21) * (a21 - e21) / e21
					+ (a12 - e12) * (a12 - e12) / e12
					+ (a22 - e22) * (a22 - e22) / e22;
#endif
	} else {
		current_chi = 0;
	}

	if (emin > w_sum) {
		sum_chi += current_chi;
		++ng_chi;
	}

	if ((emin > w_max) && (current_chi > max_chi)) {
		max_chi = current_chi;
	}

	current_like = ((a11 > 0) ? (a11 * log(a11 / e11)) : 0)
	             + ((a21 > 0) ? (a21 * log(a21 / e21)) : 0)
	             + ((a12 > 0) ? (a12 * log(a12 / e12)) : 0)
	             + ((a22 > 0) ? (a22 * log(a22 / e22)) : 0);

	sum_like += current_like;
	++ng_like;
	if (current_like > max_like) {
		max_like = current_like;
	}
}

void StatsComputer::compute_ppr_22(int xr_lo, int xr_hi, int yr_lo, int yr_hi, int pn, int nm2, double nm2s) {
	int pij_num, Aij, Bij;
	double pij, qij, current_chi, current_like;
	double emin;

    // in the notation of the grant document:
    pij_num = (yr_hi - yr_lo - 1) * (xr_hi - xr_lo - 1);

    Aij = double_integral[yr_hi * pn + xr_hi    ] + double_integral[(yr_lo + 1) * pn + xr_lo + 1]
        - double_integral[yr_hi * pn + xr_lo + 1] - double_integral[(yr_lo + 1) * pn + xr_hi    ];
    Bij = nm2 - Aij;

    pij = pij_num / nm2s;
    qij = 1 - pij;

	emin = min(pij, qij) * nm2s;

	if (emin > min_w) {
#ifndef UDF_ALLOW_DEGENERATE_PARTITIONS
		current_chi = ((Aij - nm2 * pij) * (Aij - nm2 * pij) / (nm2 * pij * (1 - pij)));
#else
		current_chi = (pij * (1 - pij) > 0) ? ((Aij - nm2 * pij) * (Aij - nm2 * pij) / (nm2 * pij * (1 - pij))) : 0;
#endif
	} else {
		current_chi = 0;
	}

	if (emin > w_sum) {
		sum_chi += current_chi;
		++ng_chi;
	}

	if ((emin > w_max) && (current_chi > max_chi)) {
		max_chi = current_chi;
	}

	current_like = ((Aij > 0) ? (Aij * log(Aij / (pij * nm2))) : 0)
		         + ((Bij > 0) ? (Bij * log(Bij / (qij * nm2))) : 0);

	sum_like += current_like;
	++ng_like;
	if (current_like > max_like) {
		max_like = current_like;
	}
}

void StatsComputer::compute_ppr_33(int xr_lo, int xr_hi, int yr_lo, int yr_hi, int n, int pn, double nm2) {
	int a11, a21, a31, a12, a22, a32, a13, a23, a33;
	double e11, e21, e31, e12, e22, e32, e13, e23, e33, current_chi, current_like;
	double emin1, emin2, emin3, emin4, emin;

	a11 = double_integral[(n    ) * pn + xr_lo] - double_integral[(yr_hi + 1) * pn + xr_lo    ];
	a21 = double_integral[(yr_hi) * pn + xr_lo] - double_integral[(yr_lo + 1) * pn + xr_lo    ];
	a31 = double_integral[(yr_lo) * pn + xr_lo];
	a12 = double_integral[(n    ) * pn + xr_hi] + double_integral[(yr_hi + 1) * pn + xr_lo + 1] - double_integral[(n    ) * pn + xr_lo + 1] - double_integral[(yr_hi + 1) * pn + xr_hi];
	a22 = double_integral[(yr_hi) * pn + xr_hi] + double_integral[(yr_lo + 1) * pn + xr_lo + 1] - double_integral[(yr_hi) * pn + xr_lo + 1] - double_integral[(yr_lo + 1) * pn + xr_hi];
	a32 = double_integral[(yr_lo) * pn + xr_hi] - double_integral[(yr_lo    ) * pn + xr_lo + 1];
	a13 = double_integral[(n    ) * pn + n    ] + double_integral[(yr_hi + 1) * pn + xr_hi + 1] - double_integral[(n    ) * pn + xr_hi + 1] - double_integral[(yr_hi + 1) * pn + n    ];
	a23 = double_integral[(yr_hi) * pn + n    ] + double_integral[(yr_lo + 1) * pn + xr_hi + 1] - double_integral[(yr_hi) * pn + xr_hi + 1] - double_integral[(yr_lo + 1) * pn + n    ];
	a33 = double_integral[(yr_lo) * pn + n    ] - double_integral[(yr_lo    ) * pn + xr_hi + 1];

#if 0
	int rs1, rs2, rs3, cs1, cs2, cs3;

	rs1 = a11 + a12 + a13;
	rs2 = a21 + a22 + a23;
	rs3 = a31 + a32 + a33;
	cs1 = a11 + a21 + a31;
	cs2 = a12 + a22 + a32;
	cs3 = a13 + a23 + a33;

	e11 = rs1 * cs1 / nm2;
	e21 = rs2 * cs1 / nm2;
	e31 = rs3 * cs1 / nm2;
	e12 = rs1 * cs2 / nm2;
	e22 = rs2 * cs2 / nm2;
	e32 = rs3 * cs2 / nm2;
	e13 = rs1 * cs3 / nm2;
	e23 = rs2 * cs3 / nm2;
	e33 = rs3 * cs3 / nm2;
#else
	e11 = (xr_lo            ) * (n     - 1 - yr_hi) / nm2;
	e21 = (xr_lo            ) * (yr_hi - 1 - yr_lo) / nm2;
	e31 = (xr_lo            ) * (yr_lo            ) / nm2;
	e12 = (xr_hi - 1 - xr_lo) * (n     - 1 - yr_hi) / nm2;
	e22 = (xr_hi - 1 - xr_lo) * (yr_hi - 1 - yr_lo) / nm2;
	e32 = (xr_hi - 1 - xr_lo) * (yr_lo            ) / nm2;
	e13 = (n     - 1 - xr_hi) * (n     - 1 - yr_hi) / nm2;
	e23 = (n     - 1 - xr_hi) * (yr_hi - 1 - yr_lo) / nm2;
	e33 = (n     - 1 - xr_hi) * (yr_lo            ) / nm2;
#endif

	emin1 = min(e11, e21);
	emin2 = min(e31, e12);
	emin3 = min(e22, e32);
	emin4 = min(e13, e23);
	emin1 = min(emin1, emin2);
	emin2 = min(emin3, emin4);
	emin1 = min(emin1, emin2);
	emin = min(emin1, e33);

	if (emin > min_w) {
#ifndef UDF_ALLOW_DEGENERATE_PARTITIONS
		current_chi = (a11 - e11) * (a11 - e11) / e11
					+ (a21 - e21) * (a21 - e21) / e21
					+ (a31 - e31) * (a31 - e31) / e31
					+ (a12 - e12) * (a12 - e12) / e12
					+ (a22 - e22) * (a22 - e22) / e22
					+ (a32 - e32) * (a32 - e32) / e32
					+ (a13 - e13) * (a13 - e13) / e13
					+ (a23 - e23) * (a23 - e23) / e23
					+ (a33 - e33) * (a33 - e33) / e33;
#else
		current_chi = ((e11 > 0) ? ((a11 - e11) * (a11 - e11) / e11) : 0)
					+ ((e21 > 0) ? ((a21 - e21) * (a21 - e21) / e21) : 0)
					+ ((e31 > 0) ? ((a31 - e31) * (a31 - e31) / e31) : 0)
					+ ((e12 > 0) ? ((a12 - e12) * (a12 - e12) / e12) : 0)
					+ ((e22 > 0) ? ((a22 - e22) * (a22 - e22) / e22) : 0)
					+ ((e32 > 0) ? ((a32 - e32) * (a32 - e32) / e32) : 0)
					+ ((e13 > 0) ? ((a13 - e13) * (a13 - e13) / e13) : 0)
					+ ((e23 > 0) ? ((a23 - e23) * (a23 - e23) / e23) : 0)
					+ ((e33 > 0) ? ((a33 - e33) * (a33 - e33) / e33) : 0);
#endif
	} else {
		current_chi = 0;
	}

	if (emin > w_sum) {
		sum_chi += current_chi;
		++ng_chi;
	}

	if ((emin > w_max) && (current_chi > max_chi)) {
		max_chi = current_chi;
	}

	current_like = ((a11 > 0) ? (a11 * log(a11 / e11)) : 0)
				 + ((a21 > 0) ? (a21 * log(a21 / e21)) : 0)
				 + ((a31 > 0) ? (a31 * log(a31 / e31)) : 0)
				 + ((a12 > 0) ? (a12 * log(a12 / e12)) : 0)
				 + ((a22 > 0) ? (a22 * log(a22 / e22)) : 0)
				 + ((a32 > 0) ? (a32 * log(a32 / e32)) : 0)
				 + ((a13 > 0) ? (a13 * log(a13 / e13)) : 0)
				 + ((a23 > 0) ? (a23 * log(a23 / e23)) : 0)
				 + ((a33 > 0) ? (a33 * log(a33 / e33)) : 0);

	sum_like += current_like;
	++ng_like;
	if (current_like > max_like) {
		max_like = current_like;
	}
}

void StatsComputer::compute_tpr(int xl, int xm, int xh, int yl, int ym, int yh, int n, int pn, double nm3) {
	int a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34, a41, a42, a43, a44;
	double e11, e12, e13, e14, e21, e22, e23, e24, e31, e32, e33, e34, e41, e42, e43, e44;
	double current_chi, current_like;
	double emin1, emin2, emin3, emin4, emin5, emin6, emin7, emin8, emin;

    a11 = double_integral[(n ) * pn + xl] - double_integral[(yh + 1) * pn + xl    ];
    a21 = double_integral[(yh) * pn + xl] - double_integral[(ym + 1) * pn + xl    ];
    a31 = double_integral[(ym) * pn + xl] - double_integral[(yl + 1) * pn + xl    ];
    a41 = double_integral[(yl) * pn + xl];

    a12 = double_integral[(n ) * pn + xm] - double_integral[(n     ) * pn + xl + 1] - double_integral[(yh + 1) * pn + xm] + double_integral[(yh + 1) * pn + xl + 1];
    a22 = double_integral[(yh) * pn + xm] - double_integral[(yh    ) * pn + xl + 1] - double_integral[(ym + 1) * pn + xm] + double_integral[(ym + 1) * pn + xl + 1];
    a32 = double_integral[(ym) * pn + xm] - double_integral[(ym    ) * pn + xl + 1] - double_integral[(yl + 1) * pn + xm] + double_integral[(yl + 1) * pn + xl + 1];
    a42 = double_integral[(yl) * pn + xm] - double_integral[(yl    ) * pn + xl + 1];

    a13 = double_integral[(n ) * pn + xh] - double_integral[(n     ) * pn + xm + 1] - double_integral[(yh + 1) * pn + xh] + double_integral[(yh + 1) * pn + xm + 1];
    a23 = double_integral[(yh) * pn + xh] - double_integral[(yh    ) * pn + xm + 1] - double_integral[(ym + 1) * pn + xh] + double_integral[(ym + 1) * pn + xm + 1];
    a33 = double_integral[(ym) * pn + xh] - double_integral[(ym    ) * pn + xm + 1] - double_integral[(yl + 1) * pn + xh] + double_integral[(yl + 1) * pn + xm + 1];
    a43 = double_integral[(yl) * pn + xh] - double_integral[(yl    ) * pn + xm + 1];

    a14 = double_integral[(n ) * pn + n ] - double_integral[(n     ) * pn + xh + 1] - double_integral[(yh + 1) * pn + n ] + double_integral[(yh + 1) * pn + xh + 1];
    a24 = double_integral[(yh) * pn + n ] - double_integral[(yh    ) * pn + xh + 1] - double_integral[(ym + 1) * pn + n ] + double_integral[(ym + 1) * pn + xh + 1];
    a34 = double_integral[(ym) * pn + n ] - double_integral[(ym    ) * pn + xh + 1] - double_integral[(yl + 1) * pn + n ] + double_integral[(yl + 1) * pn + xh + 1];
    a44 = double_integral[(yl) * pn + n ] - double_integral[(yl    ) * pn + xh + 1];

#if 0
	int rs1, rs2, rs3, rs4, cs1, cs2, cs3, cs4;

	rs1 = a11 + a12 + a13 + a14;
	rs2 = a21 + a22 + a23 + a24;
	rs3 = a31 + a32 + a33 + a34;
	rs4 = a41 + a42 + a43 + a44;
	cs1 = a11 + a21 + a31 + a41;
	cs2 = a12 + a22 + a32 + a42;
	cs3 = a13 + a23 + a33 + a43;
	cs4 = a14 + a24 + a34 + a44;

	e11 = rs1 * cs1 / nm3;
	e21 = rs2 * cs1 / nm3;
	e31 = rs3 * cs1 / nm3;
	e41 = rs4 * cs1 / nm3;
	e12 = rs1 * cs2 / nm3;
	e22 = rs2 * cs2 / nm3;
	e32 = rs3 * cs2 / nm3;
	e42 = rs4 * cs2 / nm3;
	e13 = rs1 * cs3 / nm3;
	e23 = rs2 * cs3 / nm3;
	e33 = rs3 * cs3 / nm3;
	e43 = rs4 * cs3 / nm3;
	e14 = rs1 * cs4 / nm3;
	e24 = rs2 * cs4 / nm3;
	e34 = rs3 * cs4 / nm3;
	e44 = rs4 * cs4 / nm3;
#else
	e11 = (xl         ) * (n  - 1 - yh) / nm3;
	e21 = (xl         ) * (yh - 1 - ym) / nm3;
	e31 = (xl         ) * (ym - 1 - yl) / nm3;
	e41 = (xl         ) * (yl         ) / nm3;
	e12 = (xm - 1 - xl) * (n  - 1 - yh) / nm3;
	e22 = (xm - 1 - xl) * (yh - 1 - ym) / nm3;
	e32 = (xm - 1 - xl) * (ym - 1 - yl) / nm3;
	e42 = (xm - 1 - xl) * (yl         ) / nm3;
	e13 = (xh - 1 - xm) * (n  - 1 - yh) / nm3;
	e23 = (xh - 1 - xm) * (yh - 1 - ym) / nm3;
	e33 = (xh - 1 - xm) * (ym - 1 - yl) / nm3;
	e43 = (xh - 1 - xm) * (yl         ) / nm3;
	e14 = (n  - 1 - xh) * (n  - 1 - yh) / nm3;
	e24 = (n  - 1 - xh) * (yh - 1 - ym) / nm3;
	e34 = (n  - 1 - xh) * (ym - 1 - yl) / nm3;
	e44 = (n  - 1 - xh) * (yl         ) / nm3;
#endif

	emin1 = min(e11, e21);
	emin2 = min(e31, e41);
	emin3 = min(e12, e22);
	emin4 = min(e32, e42);
	emin5 = min(e13, e23);
	emin6 = min(e33, e43);
	emin7 = min(e14, e24);
	emin8 = min(e34, e44);
	emin1 = min(emin1, emin2);
	emin2 = min(emin3, emin4);
	emin3 = min(emin5, emin6);
	emin4 = min(emin7, emin8);
	emin1 = min(emin1, emin2);
	emin2 = min(emin3, emin4);
	emin = min(emin1, emin2);

	if (emin > min_w) {
#ifndef UDF_ALLOW_DEGENERATE_PARTITIONS
		current_chi = (a11 - e11) * (a11 - e11) / e11
					+ (a21 - e21) * (a21 - e21) / e21
					+ (a31 - e31) * (a31 - e31) / e31
					+ (a41 - e41) * (a41 - e41) / e41
					+ (a12 - e12) * (a12 - e12) / e12
					+ (a22 - e22) * (a22 - e22) / e22
					+ (a32 - e32) * (a32 - e32) / e32
					+ (a42 - e42) * (a42 - e42) / e42
					+ (a13 - e13) * (a13 - e13) / e13
					+ (a23 - e23) * (a23 - e23) / e23
					+ (a33 - e33) * (a33 - e33) / e33
					+ (a43 - e43) * (a43 - e43) / e43
					+ (a14 - e14) * (a14 - e14) / e14
					+ (a24 - e24) * (a24 - e24) / e24
					+ (a34 - e34) * (a34 - e34) / e34
					+ (a44 - e44) * (a44 - e44) / e44;
#else
		current_chi = ((e11 > 0) ? ((a11 - e11) * (a11 - e11) / e11) : 0)
					+ ((e21 > 0) ? ((a21 - e21) * (a21 - e21) / e21) : 0)
					+ ((e31 > 0) ? ((a31 - e31) * (a31 - e31) / e31) : 0)
					+ ((e41 > 0) ? ((a41 - e41) * (a41 - e41) / e41) : 0)
					+ ((e12 > 0) ? ((a12 - e12) * (a12 - e12) / e12) : 0)
					+ ((e22 > 0) ? ((a22 - e22) * (a22 - e22) / e22) : 0)
					+ ((e32 > 0) ? ((a32 - e32) * (a32 - e32) / e32) : 0)
					+ ((e42 > 0) ? ((a42 - e42) * (a42 - e42) / e42) : 0)
					+ ((e13 > 0) ? ((a13 - e13) * (a13 - e13) / e13) : 0)
					+ ((e23 > 0) ? ((a23 - e23) * (a23 - e23) / e23) : 0)
					+ ((e33 > 0) ? ((a33 - e33) * (a33 - e33) / e33) : 0)
					+ ((e43 > 0) ? ((a43 - e43) * (a43 - e43) / e43) : 0)
					+ ((e14 > 0) ? ((a14 - e14) * (a14 - e14) / e14) : 0)
					+ ((e24 > 0) ? ((a24 - e24) * (a24 - e24) / e24) : 0)
					+ ((e34 > 0) ? ((a34 - e34) * (a34 - e34) / e34) : 0)
					+ ((e44 > 0) ? ((a44 - e44) * (a44 - e44) / e44) : 0);
#endif
	} else {
		current_chi = 0;
	}

	if (emin > w_sum) {
		sum_chi += current_chi;
		++ng_chi;
	}

	if ((emin > w_max) && (current_chi > max_chi)) {
		max_chi = current_chi;
	}

	current_like = ((a11 > 0) ? (a11 * log(a11 / e11)) : 0)
				 + ((a21 > 0) ? (a21 * log(a21 / e21)) : 0)
				 + ((a31 > 0) ? (a31 * log(a31 / e31)) : 0)
				 + ((a41 > 0) ? (a41 * log(a41 / e41)) : 0)
				 + ((a12 > 0) ? (a12 * log(a12 / e12)) : 0)
				 + ((a22 > 0) ? (a22 * log(a22 / e22)) : 0)
				 + ((a32 > 0) ? (a32 * log(a32 / e32)) : 0)
				 + ((a42 > 0) ? (a42 * log(a42 / e42)) : 0)
				 + ((a13 > 0) ? (a13 * log(a13 / e13)) : 0)
				 + ((a23 > 0) ? (a23 * log(a23 / e23)) : 0)
				 + ((a33 > 0) ? (a33 * log(a33 / e33)) : 0)
				 + ((a43 > 0) ? (a43 * log(a43 / e43)) : 0)
				 + ((a14 > 0) ? (a14 * log(a14 / e14)) : 0)
				 + ((a24 > 0) ? (a24 * log(a24 / e24)) : 0)
				 + ((a34 > 0) ? (a34 * log(a34 / e34)) : 0)
				 + ((a44 > 0) ? (a44 * log(a44 / e44)) : 0);

	sum_like += current_like;
	++ng_like;
	if (current_like > max_like) {
		max_like = current_like;
	}
}

inline int StatsComputer::my_rand(int lo, int hi) {
	// return a random number between lo and hi inclusive.

#if 0
	// perhaps more accurate, but slow and undeterministic.
	// rand() is not a high quality PRNG anyway

    int divisor = RAND_MAX / (hi + 1);
    int retval;

    do {
        retval = rand() / divisor;
    } while (retval > hi);

    return (retval + lo);
#else
    // slightly skewed but faster
    return (rand() % (hi - lo + 1) + lo);
#endif
}
