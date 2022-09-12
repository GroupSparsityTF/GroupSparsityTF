/**
** @file Reidentification.cpp
** 
** @brief
** Functions in the re-identification (and evaluation) phase.
** 
** Copyright (c) 2017 National Institute of Advanced Industrial Science and Technology (AIST)
** 
** This software is released under the MIT License.
** http://opensource.org/licenses/mit-license.php
*/

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <omp.h>
#include "MemoryOperation.h"
#include "Prototype.h"
#include "Struct.h"

/**
** @brief 				Compute a log-likelihood using a personalized/population transition matrix.
** @param ttensor 		Transition probability tensor (User x From Region x To Region)
** @param ttensor_n 	Population transition matrix for each target user (User x From Region x To Region)
** @param prior			Personalized prior distribution (User x Region)
** @param prior_n		Population prior distribution (User x Region)
** @param user_num 		Number of users
** @param euser_num 	Number of target users
** @param user_rnd 		Randomized user IDs
** @param r_num			Number of regions
** @param v_num 		Number of vertical regions
** @param h_num 		Number of horizontal regions
** @param tetrace 		Traces
** @param ttrace_num	Number of training traces (per user)
** @param etrace_num	Number of testing traces
** @param loglikeli		Log-likelihood computed using a personalized transitoin matrix (Testing User x Testing Trace x Target User)
** @param loglikeli_n	Log-likelihood computed using a population transition matrix (Testing User x Testing Trace x Target User)
*/
void CalcLogLikeli(double ***ttensor, double ***ttensor_n, double **prior, double **prior_n, int user_num, int euser_num, int *user_rnd, int r_num, int v_num, int h_num, struct trace **tetrace, int *ttrace_num, int *etrace_num, double ***loglikeli, double ***loglikeli_n){
	int fr_no, tr_no;
	double **f;
	double *scale;

	int u, u2, i, j;
	double dn_avg = 0.0;

	double *avg_prior_t, *avg_prior_e;

	int *rnd_event;
	int eu, ei, event_num_now;

	// Malloc
	malloc1D(&avg_prior_t, r_num);
	malloc1D(&avg_prior_e, r_num);

	// u: testing user
	for (u = 0; u < user_num; u++){
		// i: anonymized trace of user u
		for (i = 0; i < etrace_num[user_rnd[u]]; i++){
			// eu: testing user no.
			eu = user_rnd[u];
			// ei: anonymized trace no.
			ei = ttrace_num[user_rnd[u]] + i;
			// Event num in the anonymized trace --> event_num_now
			event_num_now = tetrace[eu][ei].event_num;

			// Malloc
			malloc2D(&f, event_num_now, r_num);
			malloc1D(&scale, event_num_now);
			malloc1D(&rnd_event, event_num_now);

			// Initialization
			for (u2 = 0; u2 < euser_num; u2++) loglikeli[u][i][u2] = loglikeli_n[u][i][u2] = 0.0;

			// Compute a log-likelihood using a personalized/population transition matrix --> loglikeli[u][i][u2], loglikeli_n[u][i][u2]
			// j: event
			for (j = 0; j < event_num_now - 1; j++){
				fr_no = tetrace[eu][ei].r_no[j];
				tr_no = tetrace[eu][ei].r_no[j + 1];

				// Log-likelihood using a personalized/population transition matrix
				// u2: target user
				for (u2 = 0; u2 < euser_num; u2++){
					// Update a log-likelihood for target user u2 using a personalized transition matrix
					if (j == 0) loglikeli[u][i][u2] += log(prior[u2][fr_no]);
					loglikeli[u][i][u2] += log(ttensor[u2][fr_no][tr_no]);

					// Update a log-likelihood for a non-target user using a ppulation transition matrix
					if (j == 0) loglikeli_n[u][i][u2] += log(prior_n[u2][fr_no]);
					loglikeli_n[u][i][u2] += log(ttensor_n[u2][fr_no][tr_no]);
				}
			}

			// Free
			free2D(f, event_num_now);
			free1D(scale);
			free1D(rnd_event);
		}
	}

	// Free
	free1D(avg_prior_t);
	free1D(avg_prior_e);
}

/**
** @brief 						Compute FPR, TPR, and AUC, and output them to a file.
** @param outfile 				Output file
** @param e_loglikeli_ratio 	Log-likelihood ratios for target users
** @param e_llr_num			 	Number of Log-likelihood ratios for target users
** @param n_loglikeli_ratio 	Log-likelihood ratios for non-target users
** @param n_llr_num			 	Number of Log-likelihood ratios for non-target users
** @param max_sim			 	Maximum of the similarity value
** @param min_sim			 	Minimum of the similarity value
** @param thr_shift			 	Shift value of the threshold
*/
void CalcFPRTPR(char *outfile, double *e_loglikeli_ratio, int e_llr_num, double *n_loglikeli_ratio, int n_llr_num, double max_sim, double min_sim, double thr_shift){
	double thr;
	int fneg, fpos;		// Number of false negative (resp. positive)
	int gfreq, ifreq;	// Number of genuine (resp. impostor) scores that fall below the threshold
	int gp, ip;
	FILE *fp;
	double thr_1000[1001];
	double tpr_1000[1001];
	double auc;
	int i, j;

	// Open the output file
	FileOpen(&fp, outfile, "w");

	fileout(fp, "thr,count(non-target),prob(non-target),FP,FPR,count(target),prob(target),FN,FNR,TPR\n");

	// Sort log-likelihood ratios for target/non-target users
	qsort(e_loglikeli_ratio, e_llr_num, sizeof(double), (int(*) (const void *, const void *))compare_double);
	qsort(n_loglikeli_ratio, n_llr_num, sizeof(double), (int(*) (const void *, const void *))compare_double);

	// Initialization
	fpos = n_llr_num;
	fneg = 0;
	gp = ip = 0;

	// Compute FPR, FNR, and TPR while changing a threshold
	// thr: threshold
	for (thr = min_sim; thr <= max_sim; thr += thr_shift){
		gfreq = ifreq = 0;
		// Increase gfreq and fneg while e_loglikeli_ratio[gp] (genuine score) falls below the threshold
		while (gp < e_llr_num && e_loglikeli_ratio[gp] < thr){
			gfreq++;
			fneg++;
			gp++;
		}
		// Increase ifreq and decrease fpos while n_loglikeli_ratio[ip] (impostor score) falls below the threshold
		while (ip < n_llr_num && n_loglikeli_ratio[ip] < thr){
			ifreq++;
			fpos--;
			ip++;
		}

		// Output to the file
		fileout(fp, "%.2f,%d,%f,%d,%e,%d,%f,%d,%e,%e\n", thr, ifreq, (double)ifreq / (double)n_llr_num, fpos, (double)fpos / (double)n_llr_num, gfreq, (double)gfreq / (double)e_llr_num, fneg, (double)fneg / (double)e_llr_num, 1.0 - (double)fneg / (double)e_llr_num);
	}

	// Compute AUC (using the Trapezoidal rule)
	// Compute a threshold in the case when FPR is 100%, 99.9%, ..., 0.1%, or 0% --> thr_1000
	thr_1000[0] = -DBL_MAX;
	for (i = 1; i <= 1000; i++) thr_1000[i] = n_loglikeli_ratio[(n_llr_num / 1000 * i) - 1];

	// Compute TPR in the case when FPR is 100%, 99.9%, ..., 0.1%, or 0% --> tpr_1000
	for (i = 0; i <= 1000; i++){
		for (j = 0; j < e_llr_num; j++){
			if (e_loglikeli_ratio[j] > thr_1000[i]) break;
		}
		tpr_1000[i] = (double)(e_llr_num - j) / (double)e_llr_num;
	}

	// Compute AUC --> auc
	auc = 0.0;
	for (i = 0; i < 1000; i++){
		auc += (tpr_1000[i] + tpr_1000[i + 1]) / 2000.0;
	}

	// Output to the file
	for (i = 0; i <= 1000; i++) fileout(fp, "%d,%f,%f\n", i, thr_1000[i], tpr_1000[i]);

	fileout(fp, "\nAUC\n");
	fileout(fp, "%f\n", auc);
	fclose(fp);
}
