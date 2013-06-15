/*
 * SequentialTest.cpp
 *
 */

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <ctime>
#include <cstring>
#include <pthread.h>

#include "SequentialTest.h"

#ifdef GWAS_INTERFACE
static char* get_strtime(void) { time_t now = time(0); char* strtime = asctime(localtime(&now)); strtime[strlen(strtime) - 1] = '\0'; return strtime; }
#endif

using namespace std;

// Parallelization helper functions
// ================================================================================================

#ifndef NO_THREADS

extern "C" {

static void* compute_permutations_thread(void* arg) {
	Compute_permutations_thread_arg* carg = reinterpret_cast<Compute_permutations_thread_arg*>(arg);

	(carg->seq)->compute_permutations(carg); // runs permutations and updates the sequential test state through a mutex, eventually stopping early or after all permutations are completed

	carg->done_flag = true;
	return NULL; // implicitly calls pthread_exit()
}

} // extern "C"

#endif

// ================================================================================================

SequentialTest::SequentialTest(TestType tt, int xy_nrow, int y_ncol, double* dx, double* dy, double* y,
		double w_sum, double w_max, double* extra_params,
		bool is_sequential, double alpha, double alpha0, double beta0, double eps, int nr_perm,
		int nr_threads, int base_seed, bool tables_wanted, bool perm_stats_wanted)
{
	// NOTE: Everything is done in the constructor.
	this->tt = tt;

	this->xy_nrow = xy_nrow;
	this->y_ncol = y_ncol;
	this->dx = dx;
	this->dy = dy;
	this->y = y;

	this->extra_params.w_sum = w_sum;
	this->extra_params.w_max = w_max;

	if (tt == UDF_DDP_OBS || tt == UDF_DDP_ALL) {
		this->extra_params.K = extra_params[0];
		this->extra_params.correct_mi_bias = extra_params[1];
	} else if (tt == CI_NN) {
		this->extra_params.K = extra_params[0];
	}

	this->is_sequential = is_sequential;

	// Initialize Wald sequential test state per each statistic

	// FIXME This is not very elegant, with the identity of supported statistics hard-coded.
	// I imagine I'll generalized this at some point

	this->alpha = alpha;
	this->alpha0 = alpha0;
	this->beta0 = beta0;
	this->eps = eps;
	this->nr_perm = (nr_perm == 0) ? 1 : (nr_perm + (nr_threads - nr_perm % nr_threads) + 1);
	this->nr_threads = nr_threads;
	this->base_seed = base_seed;

#ifndef _WIN32
	srand(base_seed);
#endif

	double p0 = this->alpha * (1 + this->eps);
	double p1 = this->alpha * (1 - this->eps);

	exp1 = log((p1 / (1 - p1)) / (p0 / (1 - p0)));
	exp2 = log((1 - p1) / (1 - p0));

	lA = log((1 - this->beta0) / this->alpha0);
	lB = log(this->beta0 / (1 - this->alpha0));

	nr_statistics = 6;
	llr = new double[nr_statistics];
	pvalc = new int[nr_statistics];
	stopped_high = new bool[nr_statistics];
	stopped_low = new bool[nr_statistics];
	perm_counter = new int[nr_statistics];

	for (int i = 0; i < nr_statistics; ++i) {
		llr[i] = 0;
		pvalc[i] = 1;
		stopped_high[i] = false;
		stopped_low[i] = false;
		perm_counter[i] = 0;
	}

	this->extra_params.y_counts = NULL;

	this->tables_wanted = tables_wanted;
	if (tables_wanted) {
		obs_tbls = new int[4 * xy_nrow * xy_nrow];
		for (int i = 0; i < 4 * xy_nrow * xy_nrow; ++i) {
			obs_tbls[i] = -42; // will replace these by NA's in R
		}
	}

	this->perm_stats_wanted = (perm_stats_wanted && (this->nr_perm > 1));
	this->orig_nr_perm = nr_perm;
	if (perm_stats_wanted) {
		perm_sum_chi_v 	= new double[this->nr_perm - 1];
		perm_sum_like_v = new double[this->nr_perm - 1];
		perm_max_chi_v 	= new double[this->nr_perm - 1];
		perm_max_like_v = new double[this->nr_perm - 1];
		perm_ht_v 		= new double[this->nr_perm - 1];
		perm_edist_v 	= new double[this->nr_perm - 1];
	}

	pthread_mutex_init(&mutex, NULL);

	run();
}

SequentialTest::~SequentialTest() {
	delete[] llr;
	delete[] pvalc;
	delete[] stopped_high;
	delete[] stopped_low;
	delete[] perm_counter;
	pthread_mutex_destroy(&mutex);
}

void SequentialTest::run() {
	if (tt == TWO_SAMPLE_TEST || tt == K_SAMPLE_TEST) {
		count_unique_y();
	}

	if (!IS_UDF_TEST(tt)) {
		// Sort distances
		sort_x_distances_per_row();
		sort_y_distances_per_row();
	}

	// Instantiate one StatsComputer object per thread
	scs = new StatsComputer*[nr_threads];
	for (int t = 0; t < nr_threads; ++t) {
		scs[t] = new StatsComputer(tt, xy_nrow, y_ncol, dx, dy, y, &sorted_dx, &sorted_dy, extra_params);
	}

	// Compute observed statistics
	if (tables_wanted) {
		scs[0]->compute_and_store_tbls(obs_tbls);
	} else {
		scs[0]->compute();
	}
	scs[0]->get_stats(obs_sum_chi, obs_sum_like, obs_max_chi, obs_max_like, obs_ht, obs_edist);

#ifdef ST_DEBUG_PRINTS
	cout << "Observed stats: " << obs_sum_chi << " " << obs_sum_like << " " << obs_max_chi << " " << obs_max_like << " " << obs_ht << " " << obs_edist << endl;
#endif

	if (nr_perm > 1) {
		// Compute many null statistics and estimate the p-value
		stop_all_flag = false;

#ifdef NO_THREADS
		for (int k = 0; k < nr_perm - 1; ++k) {
			scs[0]->permute_and_compute();

			double perm_sum_chi, perm_sum_like, perm_max_chi, perm_max_like, perm_ht, perm_edist;
			scs[0]->get_stats(perm_sum_chi, perm_sum_like, perm_max_chi, perm_max_like, perm_ht, perm_edist);

#ifdef ST_DEBUG_PRINTS
			cout << "Perm " << k << " stats: " << perm_sum_chi << " " << perm_sum_like << " " << perm_max_chi << " " << perm_max_like << " " << perm_ht << " " << perm_edist << endl;
#endif

			if (perm_stats_wanted) {
				perm_sum_chi_v [perm_counter[0]] = perm_sum_chi;
				perm_sum_like_v[perm_counter[1]] = perm_sum_like;
				perm_max_chi_v [perm_counter[2]] = perm_max_chi;
				perm_max_like_v[perm_counter[3]] = perm_max_like;
				perm_ht_v      [perm_counter[4]] = perm_ht;
				perm_edist_v   [perm_counter[5]] = perm_edist;
			}

			// count more extreme values of the test statistics than that observed
			bool stp_all = update_sequential_all(perm_sum_chi, perm_sum_like, perm_max_chi, perm_max_like, perm_ht, perm_edist);

			if (stp_all) {
				break;
			}

#ifdef GWAS_INTERFACE
			int k_progress = (nr_perm - 1) / 100 + 1;
			if ((k % k_progress) == 0) {
				cout << get_strtime() << ": permutation test " << round((100.0 * k) / nr_perm) << "% complete" << endl;
			}
#endif
		}

#else // => use threads

		pthread_t* worker_threads = new pthread_t[nr_threads];

		pthread_attr_t attr;
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED); // due to an extremely annoying bug, I am not using joinable threads
		pthread_attr_setstacksize(&attr, 4 * 1024 * 1024); // if we have a maximal feasible xy_nrow = 10000, then this should leave plenty of room for recursive calls {1,2,...,log2(xy_nrow)}

		Compute_permutations_thread_arg** worker_args = new Compute_permutations_thread_arg*[nr_threads];
		int* pthread_status = new int[nr_threads];

		for (int t = 0; t < nr_threads; ++t) {
			worker_args[t] = new Compute_permutations_thread_arg(this, t);
			pthread_status[t] = pthread_create(&(worker_threads[t]), &attr, compute_permutations_thread, (void*)(worker_args[t]));
		}

		// Manually join the worker threads (might be unnecessarily slow for small datasets)
		bool all_done = false;
		while (!all_done) {
			all_done = true;
			for (int t = 0; t < nr_threads; ++t) {
				all_done &= worker_args[t]->done_flag;
			}

#ifdef _WIN32
			Sleep(100);
#else
			usleep(100);
#endif
		}

		for (int t = 0; t < nr_threads; ++t) {
			delete worker_args[t];
		}

		delete[] worker_args;
		pthread_attr_destroy(&attr);
		delete[] worker_threads;

#endif // NO_THREADS

#ifdef GWAS_INTERFACE
		if (stop_all_flag) {
			cout << get_strtime() << ": Sequential testing says stop on all statistics. Moving on." << endl;
		} else {
			cout << get_strtime() << ": Computed all permutations!" << endl;
			cout << "   p_sum_chi  = " << ((double)(pvalc[0])) / nr_perm << "\t(observed " << obs_sum_chi  << ")" << endl;
			cout << "   p_sum_like = " << ((double)(pvalc[1])) / nr_perm << "\t(observed " << obs_sum_like << ")" << endl;
			cout << "   p_max_chi  = " << ((double)(pvalc[2])) / nr_perm << "\t(observed " << obs_max_chi  << ")" << endl;
			cout << "   p_max_like = " << ((double)(pvalc[3])) / nr_perm << "\t(observed " << obs_max_like << ")" << endl;
		}
#endif
	}

	for (int t = 0; t < nr_threads; ++t) {
		delete scs[t];
	}
	delete[] scs;

	if (this->extra_params.y_counts != NULL) {
		delete[] this->extra_params.y_counts;
	}
}

// NOTE: This function runs in the context of a worker thread
void SequentialTest::compute_permutations(Compute_permutations_thread_arg* carg) {
	int t = carg->t;

	// NOTE: I made sure (nr_perm - 1) is a multiple of nr_threads by rounding it up
	int nr_perm_per_thread = (nr_perm - 1) / nr_threads;

	// rand() is not thread safe, and its behavior varies markedly between Linux (process state, usually) and Windows (thread state, usually)
	// The following is not a proper solution, but it'll have to do for now.
#ifdef _WIN32
	srand(base_seed + t); // FIXME in the future I might want to base this on a random draw, and allow an external seed to be provided and srand()'ed once at the main thread
#else
	// TODO use drand48_r() on Linux. For now it's not critical since I don't really use multiple
	// concurrent threads due to the CS cluster's limitations.
#endif

	for (int k = 0; k < nr_perm_per_thread; ++k) {
		scs[t]->permute_and_compute();

	    double perm_sum_chi, perm_sum_like, perm_max_chi, perm_max_like, perm_ht, perm_edist;
	    scs[t]->get_stats(perm_sum_chi, perm_sum_like, perm_max_chi, perm_max_like, perm_ht, perm_edist);

		pthread_mutex_lock(&mutex);

#ifdef ST_DEBUG_PRINTS
		if (t == 0) {
			cout << "Perm " << k + t << " stats: " << perm_sum_chi << " " << perm_sum_like << " " << perm_max_chi << " " << perm_max_like << " " << perm_ht << " " << perm_edist << endl;
		}
#endif

		if (perm_stats_wanted) {
			perm_sum_chi_v [perm_counter[0]] = perm_sum_chi;
			perm_sum_like_v[perm_counter[1]] = perm_sum_like;
			perm_max_chi_v [perm_counter[2]] = perm_max_chi;
			perm_max_like_v[perm_counter[3]] = perm_max_like;
			perm_ht_v      [perm_counter[4]] = perm_ht;
			perm_edist_v   [perm_counter[5]] = perm_edist;
		}

		// count more extreme values of the test statistics than that observed
		bool stp_all = update_sequential_all(perm_sum_chi, perm_sum_like, perm_max_chi, perm_max_like, perm_ht, perm_edist);

		if (stp_all) {
			pthread_mutex_unlock(&mutex);
			break;
		}

#ifdef GWAS_INTERFACE
		int k_progress = nr_perm_per_thread / 100 + 1;
		if ((t == 0) && (k % k_progress == 0)) {
			cout << get_strtime() << ": permutation test " << round((100.0 * k) / nr_perm_per_thread) << "% complete" << endl;
		}
#endif

		pthread_mutex_unlock(&mutex);
	}
}

bool SequentialTest::is_null_rejected(void) {
	return (((double)(pvalc[0])) / nr_perm < beta0); // FIXME allow user to choose statistic
}

void SequentialTest::get_observed_stats(double& sum_chi, double& sum_like, double& max_chi, double& max_like, double& ht, double& edist) {
	sum_chi 	= this->obs_sum_chi;
	sum_like 	= this->obs_sum_like;
	max_chi 	= this->obs_max_chi;
	max_like 	= this->obs_max_like;
	ht 			= this->obs_ht;
	edist 		= this->obs_edist;
}

void SequentialTest::get_pvalues(double& p_sum_chi, double& p_sum_like, double& p_max_chi, double& p_max_like, double& p_ht, double& p_edist) {
	if (nr_perm == 0) {
		p_sum_chi  = 1;
		p_sum_like = 1;
		p_max_chi  = 1;
		p_max_like = 1;
		p_ht       = 1;
		p_edist    = 1;
	} else {
		p_sum_chi  = ((double)(pvalc[0])) / nr_perm;
		p_sum_like = ((double)(pvalc[1])) / nr_perm;
		p_max_chi  = ((double)(pvalc[2])) / nr_perm;
		p_max_like = ((double)(pvalc[3])) / nr_perm;
		p_ht       = ((double)(pvalc[4])) / nr_perm;
		p_edist    = ((double)(pvalc[5])) / nr_perm;
	}
}

void SequentialTest::get_observed_tables(double* tbls) {
	if (tables_wanted) {
		for (int i = 0; i < 4 * xy_nrow * xy_nrow; ++i) {
			tbls[i] = obs_tbls[i];
		}
	}
}

void SequentialTest::get_perm_stats(double* ps) {
	if (perm_stats_wanted) {
		for (int i = 0; i < orig_nr_perm; ++i) {
			ps[0 * orig_nr_perm + i] = perm_sum_chi_v [i];
			ps[1 * orig_nr_perm + i] = perm_sum_like_v[i];
			ps[2 * orig_nr_perm + i] = perm_max_chi_v [i];
			ps[3 * orig_nr_perm + i] = perm_max_like_v[i];
			ps[4 * orig_nr_perm + i] = perm_ht_v      [i];
			ps[5 * orig_nr_perm + i] = perm_edist_v   [i];
		}
	}
}

void SequentialTest::count_unique_y(void) {
	int K = 0;
	for (int i = 0; i < xy_nrow; ++i) {
		K = max(K, (int)(y[i]));
	}
	K = max(2, K + 1);
	this->extra_params.K = K;

	this->extra_params.y_counts = new int[K];
	for (int k = 0; k < K; ++k) {
		this->extra_params.y_counts[k] = 0;
	}

	for (int i = 0; i < xy_nrow; ++i) {
		++this->extra_params.y_counts[(int)(y[i])];
	}
}

void SequentialTest::sort_x_distances_per_row(void) {
	sorted_dx.resize(xy_nrow);

	// FIXME this deserves parallelizing. Trivial to do per row.

	for (int k = 0; k < xy_nrow; ++k) {
		sorted_dx[k].resize(xy_nrow);

		for (int l = 0; l < xy_nrow; ++l) {
			sorted_dx[k][l].first = dx[l * xy_nrow + k];
			sorted_dx[k][l].second = l;
		}

		sort(sorted_dx[k].begin(), sorted_dx[k].end(), dbl_int_pair_comparator);
	}
}

void SequentialTest::sort_y_distances_per_row(void) {
	sorted_dy.resize(xy_nrow);

	// FIXME this deserves parallelizing. Trivial to do per row.

	for (int k = 0; k < xy_nrow; ++k) {
		sorted_dy[k].resize(xy_nrow);

		for (int l = 0; l < xy_nrow; ++l) {
			sorted_dy[k][l].first = dy[l * xy_nrow + k];
			sorted_dy[k][l].second = l;
		}

		sort(sorted_dy[k].begin(), sorted_dy[k].end(), dbl_int_pair_comparator);
	}
}

bool SequentialTest::update_sequential(int statistic_idx, bool is_null_more_extreme) {
	if (!is_sequential) {
		pvalc[statistic_idx] += is_null_more_extreme;
		++perm_counter[statistic_idx];
		return false;
	}

	if (stopped_high[statistic_idx]) {
		return true;
	}

	pvalc[statistic_idx] += is_null_more_extreme;
	llr[statistic_idx] += is_null_more_extreme * exp1 + exp2;
	++perm_counter[statistic_idx];

	if ((!stopped_low[statistic_idx]) && (llr[statistic_idx] <= lB)) {
		// mark as useless and stop permuting
		pvalc[statistic_idx] = nr_perm;
		stopped_high[statistic_idx] = true;
		return true;
	} else if (llr[statistic_idx] >= lA) {
		// we can probably stop permuting and conclude x is truly
		// associated with y. But we want the full p-value
		// in this case so that we can follow up with a data driven
		// multiplicity adjustment such as BH. So we have to just make
		// sure that we wont stop in the other direction (unlikely but
		// possible).
		stopped_low[statistic_idx] = true;
	}

	return false;
}

bool SequentialTest::update_sequential_all(double perm_sum_chi, double perm_sum_like, double perm_max_chi, double perm_max_like, double perm_ht, double perm_edist) {
	bool stp_all = true;

	stp_all &= update_sequential(0, perm_sum_chi  >= obs_sum_chi );
	stp_all &= update_sequential(1, perm_sum_like >= obs_sum_like);
	stp_all &= update_sequential(2, perm_max_chi  >= obs_max_chi );
	stp_all &= update_sequential(3, perm_max_like >= obs_max_like);
	stp_all &= update_sequential(4, perm_ht       >= obs_ht      );
	stp_all &= update_sequential(5, perm_edist    >= obs_edist   );

	stop_all_flag = stp_all;

	return (stp_all);
}
