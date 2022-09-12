/**
** @file Train.cpp
** 
** @brief
** Functions in the training phase (except for TF and GSTF).
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
** @brief 				Decide region boundaries and assign region IDs to events in traces.
** @param tetrace 		Traces
** @param user_num 		Number of users
** @param trace_num		Number of traces (per user)
** @param v_num 		Number of regions in a vertical way
** @param h_num 		Number of regions in a horizontal way
** @param v_min 		Minimum latitude
** @param v_max 		Maximum latitude
** @param h_min			Minimum longitude
** @param h_max 		Maximum longitude
*/
void AssignRegionID(struct trace **tetrace, int user_num, int *trace_num, int v_num, int h_num, double v_min, double v_max, double h_min, double h_max){
	double *bou_v, *bou_h;
	int v_no, h_no, r_no;
	int u, i, j, k;

	// Malloc
	malloc1D(&bou_v, v_num);
	malloc1D(&bou_h, h_num);

	// Decide region boundaries --> bou_v, bou_h
	for (i = 0; i < v_num; i++){
		bou_v[i] = v_min + (v_max - v_min) * (double)i / (double)v_num;
	}
	for (i = 0; i < h_num; i++){
		bou_h[i] = h_min + (h_max - h_min) * (double)i / (double)h_num;
	}

	// Assign region IDs to events in traces
	for (u = 0; u < user_num; u++){
		for (i = 0; i < trace_num[u]; i++){
			for (j = 0; j < tetrace[u][i].event_num; j++){
				// Vertical region no., Horizontal region no., Region no. --> tetrace[u][i].v_no[j], tetrace[u][i].h_no[j], tetrace[u][i].r_no[j]
				for (k = 0; k < v_num - 1; k++){
					if (tetrace[u][i].v[j] >= bou_v[k] && tetrace[u][i].v[j] < bou_v[k + 1]) break;
				}
				v_no = k;
				for (k = 0; k < h_num - 1; k++){
					if (tetrace[u][i].h[j] >= bou_h[k] && tetrace[u][i].h[j] < bou_h[k + 1]) break;
				}
				h_no = k;

				tetrace[u][i].v_no[j] = v_no;
				tetrace[u][i].h_no[j] = h_no;
				r_no = v_no * h_num + h_no;
				tetrace[u][i].r_no[j] = r_no;
			}
		}
	}

	free1D(bou_v);
	free1D(bou_h);
}

/**
** @brief 				Compute a transition pattern structure from training traces.
** @param tetrace 		Traces
** @param user_num 		Number of users
** @param ttrace_num	Number of training traces (per user)
** @param tevent_num 	Number of training events (per user)
** @param ttrans 		Transition pattern
** @param ttrans_num 	Number of transition patterns (per user)
*/
void MakeTrans(struct trace **tetrace, int user_num, int *ttrace_num, int *tevent_num, struct trans **ttrans, int *ttrans_num){
	int fr_no, tr_no;
	int u, i, j, k;
	int tevent_num_tmp;

	for (u = 0; u < user_num; u++){
		ttrans_num[u] = 0;
		tevent_num_tmp = 0;
		for (i = 0; i < ttrace_num[u]; i++){
			// Update the number of training events (add the first event)
			tevent_num_tmp++;
			for (j = 0; j < tetrace[u][i].event_num - 1; j++){
				// Break if the number of training events is larger than or equal to tevent_num[u] (we use the first tevent_num[u] events for training)
				if (tevent_num_tmp >= tevent_num[u]) break;

				// Transition pattern --> fr_no, tr_no
				fr_no = tetrace[u][i].r_no[j];
				tr_no = tetrace[u][i].r_no[j + 1];

				if (fr_no == -1 || tr_no == -1) continue;

				// Find the same transition pattern as [fr_no --> tr_no] from ttrans[u]
				for (k = 0; k < ttrans_num[u]; k++){
					// If the same transition pattern is found, add one to the corresponding transition count.
					if (ttrans[u][k].fr_no == fr_no && ttrans[u][k].tr_no == tr_no){
						ttrans[u][k].num = ttrans[u][k].num + 1.0;
						break;
					}
				}
				// If [fr_no --> tr_no] is a new transition pattern, save it to ttrans[u]
				if (k == ttrans_num[u]){
					ttrans[u][k].fr_no = fr_no;
					ttrans[u][k].tr_no = tr_no;
					ttrans[u][k].num = 1.0;
					ttrans_num[u]++;
				}

				// Update the number of training events
				tevent_num_tmp++;
			}
		}
	}
}

/**
** @brief 				Compute a transition count tensor and transition count matrix for all users from a transition pattern structure.
** @param ttrans 		Transition pattern (User x Transition Pattern)
** @param ttrans_num 	Number of transition patterns (per user)
** @param user_rnd 		Randomized user IDs
** @param euser_num 	Number of target users
** @param r_num			Number of regions
** @param ttensor	 	Transition count tensor (User x From Region x To Region)
** @param tmatrix_all 	Transition count matrix for all users (From Region x To Region)
*/
void Trans2Tensor(struct trans **ttrans, int *ttrans_num, int *user_rnd, int euser_num, int r_num, double ***ttensor, double **tmatrix_all){
	int u, i, j;
	int fr_no, tr_no;
	int exst;
	double num;

	// Initialization (-1: unobserved element) --> ttensor
	for (u = 0; u < euser_num; u++){
		for (i = 0; i < r_num; i++){
			for (j = 0; j < r_num; j++){
				ttensor[u][i][j] = -1.0;
			}
		}
	}
	// Initialization (-1: unobserved element) --> tmatrix_all
	for (i = 0; i < r_num; i++){
		for (j = 0; j < r_num; j++){
			tmatrix_all[i][j] = -1.0;
		}
	}

	// Compute a transition count tensor and transition count matrix for all users from ttrans --> ttensor, tmatrix_all
	for (u = 0; u < euser_num; u++){
		for (i = 0; i < ttrans_num[user_rnd[u]]; i++){
			fr_no = ttrans[user_rnd[u]][i].fr_no;
			tr_no = ttrans[user_rnd[u]][i].tr_no;
			num = ttrans[user_rnd[u]][i].num;
			ttensor[u][fr_no][tr_no] = num;
			if (tmatrix_all[fr_no][tr_no] == -1.0) tmatrix_all[fr_no][tr_no] = num;
			else tmatrix_all[fr_no][tr_no] += num;
		}
	}

	// If there are transitions in a row, change -1 in the row into 0 --> ttensor
	for (u = 0; u < euser_num; u++){
		for (i = 0; i < r_num; i++){
			exst = 0;
			for (j = 0; j < r_num; j++){
				if (ttensor[u][i][j] > 0.0){
					exst = 1;
					break;
				}
			}
			if (exst){
				for (j = 0; j < r_num; j++){
					if (ttensor[u][i][j] == -1.0) ttensor[u][i][j] = 0.0;
				}
			}
		}
	}

	// If there are transitions in a row, change -1 in the row into 0 --> tmatrix_all
	for (i = 0; i < r_num; i++){
		exst = 0;
		for (j = 0; j < r_num; j++){
			if (tmatrix_all[i][j] > 0.0){
				exst = 1;
				break;
			}
		}
		if (exst){
			for (j = 0; j < r_num; j++){
				if (tmatrix_all[i][j] == -1.0) tmatrix_all[i][j] = 0.0;
			}
		}
	}
}

/**
** @brief 				Compute a transition probability tensor from model parameters in TF.
** @param v_u_f 		Model parameters (User x From Region)
** @param v_f_u 		Model parameters (From Region x User)
** @param v_f_t 		Model parameters (From Region x To Region)
** @param v_t_f 		Model parameters (To Region x From Region)
** @param v_u_t 		Model parameters (User x To Region)
** @param v_t_u 		Model parameters (To Region x User)
** @param k_u_f 		Number of feature dimensions in v_u_f and v_f_u
** @param k_f_t 		Number of feature dimensions in v_f_t and v_t_f
** @param k_u_t 		Number of feature dimensions in v_u_t and v_t_u
** @param euser_num 	Number of target users
** @param r_num			Number of regions
** @param ttensor 		Transition probability tensor (User x From Region x To Region)
*/
void TFParam2Tensor(double **v_u_f, double **v_f_u, double **v_f_t, double **v_t_f, double **v_u_t, double **v_t_u, int k_u_f, int k_f_t, int k_u_t, int euser_num, int r_num, double ***ttensor){
	double sum;
	int u, i, j, k;

	// Transition count --> ttensor
	for (u = 0; u < euser_num; u++){
		for (i = 0; i < r_num; i++){
			for (j = 0; j < r_num; j++){
				ttensor[u][i][j] = 0.0;
				for (k = 0; k < k_u_f; k++) ttensor[u][i][j] += v_u_f[u][k] * v_f_u[i][k];
				for (k = 0; k < k_f_t; k++) ttensor[u][i][j] += v_f_t[i][k] * v_t_f[j][k];
				for (k = 0; k < k_u_t; k++) ttensor[u][i][j] += v_u_t[u][k] * v_t_u[j][k];
				if (ttensor[u][i][j] < 0.0){
					printf("TFParam2Tensor: negative transition count (=%f) ([%d]-[%d]-[%d])D\n", ttensor[u][i][j], u, i, j);
					ttensor[u][i][j] = 0.0;
				}
			}
		}
	}

	// Normalize counts to probabilities --> ttensor
	for (u = 0; u < euser_num; u++){
		for (i = 0; i < r_num; i++){
			// Summation over To Region --> sum
			sum = 0.0;
			for (j = 0; j < r_num; j++) sum += ttensor[u][i][j];
			// Normalize counts to probabilities so that the summation over To Region is 1
			if (sum > 0){
				for (j = 0; j < r_num; j++) ttensor[u][i][j] = ttensor[u][i][j] / sum;
			}
			// If the summation over To Region is 0 (if there are no training data), assign 1 / r_num to each element (i.e. uniform distribution)
			else{
				for (j = 0; j < r_num; j++) ttensor[u][i][j] = 1.0 / r_num;
			}
		}
	}
}

/**
** @brief 				Train a transition probability tensor via the ML estimation method.
** @param ttensor 		Transition probability tensor (User x From Region x To Region)
** @param euser_num 	Number of target users
** @param r_num			Number of regions
*/
void ML(double ***ttensor, int euser_num, int r_num){
	double sum;
	int u, i, j;

	// Normalize counts to probabilities
	for (u = 0; u < euser_num; u++){
		for (i = 0; i < r_num; i++){
			// Summation over To Region --> sum
			sum = 0.0;
			for (j = 0; j < r_num; j++) sum += ttensor[u][i][j];
			// Normalize counts to probabilities so that the summation over To Region is 1
			if (sum > 0){
				for (j = 0; j < r_num; j++) ttensor[u][i][j] = ttensor[u][i][j] / sum;
			}
			// If the summation over To Region is less than 0, assign 1 / r_num to each element (i.e. uniform distribution)
			// (we assigned -1 to unobserved elements; see the Trans2Tensor function)
			else{
				for (j = 0; j < r_num; j++) ttensor[u][i][j] = 1.0 / r_num;
			}
		}
	}
}

/**
** @brief 				Compute a prior distribution (stationary distribution) from a transition probability tensor.
** @param ttensor 		Transition probability tensor (User x From Region x To Region)
** @param euser_num 	Number of target users
** @param r_num			Number of regions
** @param prior			Prior distribution (User x Region)
*/
void Tensor2Prior(double ***ttensor, int euser_num, int r_num, double **prior){
	double pre_prior[MAX_REGION_NUM];
	double chng_avg;
	int u, i, j, k;
	double r_inv = 1.0 / (double)r_num;

	double *** new_ttensor;
	malloc3D(&new_ttensor, euser_num, r_num, r_num);

	// Use new_ttensor for efficient caching
	for (u = 0; u < euser_num;u++)
		for (int i = 0; i < r_num;i++)
			for (int k = 0; k < r_num; k++)
				new_ttensor[u][i][k] = ttensor[u][k][i];

	// Initialization (uniform distribution) --> pre_prior
	for (i = 0; i < r_num; i++) pre_prior[i] = 1.0 / (double)r_num;

	for (u = 0; u < euser_num; u++){
		for (j = 0; j < 1000; j++){
			// Update a personalized prior distribution --> prior
			for (i = 0; i < r_num; i++)
				prior[u][i] = 0.0;
			for (i = 0; i < r_num; i++)
				for (k = 0; k < r_num; k++)
					prior[u][i] += pre_prior[k] * new_ttensor[u][i][k];

			// Average difference between prior and pre_prior --> chng_avg
			chng_avg = 0.0;
			for (i = 0; i < r_num; i++) 
				chng_avg += fabs(prior[u][i] - pre_prior[i]);
			chng_avg *= r_inv;
			
			if (chng_avg < PRIOR_CHNG_THR)
				break;

			// Personalized prior distribution --> pre_prior
			for (i = 0; i < r_num; i++) pre_prior[i] = prior[u][i];

		}
	}

	free3D(new_ttensor, euser_num, r_num);
}

/**
** @brief 				Compute a population transition matrix for each target user.
** @param ttensor 		Transition probability tensor (User x From Region x To Region)
** @param ttensor_n 	Population transition matrix for each target user (User x From Region x To Region)
** @param euser_num 	Number of target users
** @param r_num			Number of regions
*/
void CalcPopMat(double ***ttensor, double ***ttensor_n, int euser_num, int r_num){
	int u, u2, i, j;
	double avg_prob;

	double ***new_ttensor;
	malloc3D(&new_ttensor, r_num, r_num, euser_num);

	// Use new_ttensor for efficient caching
	for (u = 0; u < euser_num; u++){
		for (i = 0; i < r_num; i++){
			for (j = 0; j < r_num; j++){
				new_ttensor[i][j][u] = ttensor[u][i][j];
			}
		}
	}

	// Initalization
	for (u = 0; u < euser_num; u++){
		for (i = 0; i < r_num; i++){
			for (j = 0; j < r_num; j++){
				ttensor_n[u][i][j] = 0.0;
			}
		}
	}

	double euser_inv = 1.0 / (double)(euser_num - 1);
	// Compute a population transition matrix for each target user
	for (u = 0; u < euser_num; u++){
		for (i = 0; i < r_num; i++){
			for (j = 0; j < r_num; j++){
				// Average transition probabilities over target users except for u --> avg_prob
				avg_prob = 0.0;
				for (u2 = 0; u2 < u;u2++)
					avg_prob += new_ttensor[i][j][u2];
				for (u2 = u + 1; u2 < euser_num;u2++)
					avg_prob += new_ttensor[i][j][u2];
				avg_prob *= euser_inv;

				ttensor_n[u][i][j] = avg_prob;
			}
		}
	}

	free3D(new_ttensor, r_num, r_num);
}
