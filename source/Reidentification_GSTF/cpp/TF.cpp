/**
** @file TF.cpp
** 
** @brief
** Tensor Factorization (i.e. training method in [Murakami+, TIFS16]).
** 
** Copyright (c) 2017 National Institute of Advanced Industrial Science and Technology (AIST)
** 
** This software is released under the MIT License.
** http://opensource.org/licenses/mit-license.php
*/

#include <stdio.h>
#include <math.h>
#include "MemoryOperation.h"
#include "Prototype.h"
#include "parameter.h"

/**
** @brief 			Train model parameters (i.e. v_u_f, v_f_u, v_f_t, v_t_f, v_u_t, v_t_u) in TF (Tensor Factorization)
					from personalized transition matrices (i.e. ttensor).
** @param v_u_f 	Model parameters (User x From Region)
** @param v_f_u 	Model parameters (From Region x User)
** @param v_f_t 	Model parameters (From Region x To Region)
** @param v_t_f 	Model parameters (To Region x From Region)
** @param v_u_t 	Model parameters (User x To Region)
** @param v_t_u 	Model parameters (To Region x User)
** @param r_num 	Number of regions
** @param euser_num Number of target users
** @param k_u_f 	Number of feature dimensions in v_u_f and v_f_u
** @param k_f_t 	Number of feature dimensions in v_f_t and v_t_f
** @param k_u_t 	Number of feature dimensions in v_u_t and v_t_u
** @param ttensor 	Transition count tensor (User x From Region x To Region)
** @param alpha 	Regularization parameter
*/
void TF(double **v_u_f, double **v_f_u, double **v_f_t, double **v_t_f, double **v_u_t, double **v_t_u, int r_num, int euser_num, 
	int k_u_f, int k_f_t, int k_u_t, double ***ttensor, double alpha){

	double ***n_est;
	int tu, ti, tj, tk;
	double tn;
	double err = 0.0, err_pre = 0.0;
	double numerator, denominator;
	double v_tmp, v_pre;
	double chng, chng_max, chng_avg;
	int chng_num;
	int i, j, k;
	int converge = 0, converge_tmp;

	// Malloc
	malloc3D(&n_est, euser_num, r_num, r_num);

	// Initialization of model parameters (random initialization)
	for (i = 0; i < euser_num; i++){
		for (k = 0; k < k_u_f; k++){
			v_u_f[i][k] = getRandomValue();
			if (v_u_f[i][k] < EPSILON) v_u_f[i][k] = EPSILON;
		}
	}
	for (i = 0; i < r_num; i++){
		for (k = 0; k < k_u_f; k++){
			v_f_u[i][k] = getRandomValue();
			if (v_f_u[i][k] < EPSILON) v_f_u[i][k] = EPSILON;
		}
	}
	for (i = 0; i < r_num; i++){
		for (k = 0; k < k_f_t; k++){
			v_f_t[i][k] = getRandomValue();
			if (v_f_t[i][k] < EPSILON) v_f_t[i][k] = EPSILON;
		}
	}
	for (i = 0; i < r_num; i++){
		for (k = 0; k < k_f_t; k++){
			v_t_f[i][k] = getRandomValue();
			if (v_t_f[i][k] < EPSILON) v_t_f[i][k] = EPSILON;
		}
	}
	for (i = 0; i < euser_num; i++){
		for (k = 0; k < k_u_t; k++){
			v_u_t[i][k] = getRandomValue();
			if (v_u_t[i][k] < EPSILON) v_u_t[i][k] = EPSILON;
		}
	}
	for (i = 0; i < r_num; i++){
		for (k = 0; k < k_u_t; k++){
			v_t_u[i][k] = getRandomValue();
			if (v_t_u[i][k] < EPSILON) v_t_u[i][k] = EPSILON;
		}
	}

	// Estimate transition counts --> n_est
	for (tu = 0; tu < euser_num; tu++){
		for (ti = 0; ti < r_num; ti++){
			for (tj = 0; tj < r_num; tj++){
				n_est[tu][ti][tj] = 0.0;
				for (k = 0; k < k_u_f; k++) n_est[tu][ti][tj] += v_u_f[tu][k] * v_f_u[ti][k];
				for (k = 0; k < k_f_t; k++) n_est[tu][ti][tj] += v_f_t[ti][k] * v_t_f[tj][k];
				for (k = 0; k < k_u_t; k++) n_est[tu][ti][tj] += v_u_t[tu][k] * v_t_u[tj][k];
			}
		}
	}

	// Train model parameters until convergence
	j = 0;
	while (1){
		converge_tmp = 1;
		chng_max = chng_avg = 0.0;
		chng_num = 0;

		// Update v_u_f (User x From Region)
		// tk: Feature dimension
		for (tk = 0; tk < k_u_f; tk++){
			// tu: User
			for (tu = 0; tu < euser_num; tu++){
				numerator = denominator = 0.0;
				// ti: From Region
				for (ti = 0; ti < r_num; ti++){
					// continue if elements are unobserved
					tn = ttensor[tu][ti][0];
					if (tn == -1.0) continue;

					// tj: To Region
					for (tj = 0; tj < r_num; tj++){
						// Transition count --> tn
						tn = ttensor[tu][ti][tj];

						// Numerator --> numerator
						numerator += (tn - n_est[tu][ti][tj] + v_u_f[tu][tk] * v_f_u[ti][tk]) * v_f_u[ti][tk];
						// Denominator --> denominator
						denominator += v_f_u[ti][tk] * v_f_u[ti][tk];
					}
				}
				// Update v_u_f[tu][tk]
				if (denominator + alpha == 0.0) v_tmp = EPSILON;
				else{
					v_tmp = numerator / (denominator + alpha);
					if (v_tmp < EPSILON) v_tmp = EPSILON;	// NTF
				}

				chng = fabs(v_tmp - v_u_f[tu][tk]) / v_u_f[tu][tk];
				if (chng > ANLS_CHNG_THR) converge_tmp = 0;
				if (chng > chng_max) chng_max = chng;
				chng_avg += chng;
				chng_num++;
				v_pre = v_u_f[tu][tk];
				v_u_f[tu][tk] = v_tmp;

				// Estimate transition counts --> n_est
				for (ti = 0; ti < r_num; ti++){
					// continue if elements are unobserved
					tn = ttensor[tu][ti][0];
					if (tn == -1.0) continue;
					for (tj = 0; tj < r_num; tj++){
						n_est[tu][ti][tj] += v_u_f[tu][tk] * v_f_u[ti][tk] - v_pre * v_f_u[ti][tk];
					}
				}
			}
		}

		// Update v_f_u (From Region x User)
		// tk: Feature dimension
		for (tk = 0; tk < k_u_f; tk++){
			// ti: From Region
			for (ti = 0; ti < r_num; ti++){
				numerator = denominator = 0.0;
				// tu: User
				for (tu = 0; tu < euser_num; tu++){
					// continue if elements are unobserved
					tn = ttensor[tu][ti][0];
					if (tn == -1.0) continue;

					// tj: To Region
					for (tj = 0; tj < r_num; tj++){
						// Transition count --> tn
						tn = ttensor[tu][ti][tj];

						// Numerator --> numerator
						numerator += (tn - n_est[tu][ti][tj] + v_u_f[tu][tk] * v_f_u[ti][tk]) * v_u_f[tu][tk];
						// Denominator --> denominator
						denominator += v_u_f[tu][tk] * v_u_f[tu][tk];
					}
				}
				// Update v_f_u[ti][tk]
				if (denominator + alpha == 0.0) v_tmp = EPSILON;
				else{
					v_tmp = numerator / (denominator + alpha);
					if (v_tmp < EPSILON) v_tmp = EPSILON;	// NTF
				}

				chng = fabs(v_tmp - v_f_u[ti][tk]) / v_f_u[ti][tk];
				if (chng > ANLS_CHNG_THR) converge_tmp = 0;
				if (chng > chng_max) chng_max = chng;
				chng_avg += chng;
				chng_num++;
				v_pre = v_f_u[ti][tk];
				v_f_u[ti][tk] = v_tmp;

				// Estimate transition counts --> n_est
				for (tu = 0; tu < euser_num; tu++){
					// continue if elements are unobserved
					tn = ttensor[tu][ti][0];
					if (tn == -1.0) continue;
					for (tj = 0; tj < r_num; tj++){
						n_est[tu][ti][tj] += v_u_f[tu][tk] * (v_f_u[ti][tk] - v_pre);
					}
				}

			}
		}

		// Update v_f_t (From Region x To Region)
		// tk: Feature dimension
		for (tk = 0; tk < k_f_t; tk++){
			// ti: From Region
			for (ti = 0; ti < r_num; ti++){
				numerator = denominator = 0.0;
				// tu: User
				for (tu = 0; tu < euser_num; tu++){
					// continue if elements are unobserved
					tn = ttensor[tu][ti][0];
					if (tn == -1.0) continue;

					// tj: To Region
					for (tj = 0; tj < r_num; tj++){
						// Transition count --> tn
						tn = ttensor[tu][ti][tj];

						// Numerator --> numerator
						numerator += (tn - n_est[tu][ti][tj] + v_f_t[ti][tk] * v_t_f[tj][tk]) * v_t_f[tj][tk];
						// Denominator --> denominator
						denominator += v_t_f[tj][tk] * v_t_f[tj][tk];
					}
				}
				// Update v_f_t[ti][tk]
				if (denominator + alpha == 0.0) v_tmp = EPSILON;
				else{
					v_tmp = numerator / (denominator + alpha);
					if (v_tmp < EPSILON) v_tmp = EPSILON;	// NTF
				}

				chng = fabs(v_tmp - v_f_t[ti][tk]) / v_f_t[ti][tk];
				if (chng > ANLS_CHNG_THR) converge_tmp = 0;
				if (chng > chng_max) chng_max = chng;
				chng_avg += chng;
				chng_num++;
				v_pre = v_f_t[ti][tk];
				v_f_t[ti][tk] = v_tmp;

				// Estimate transition counts --> n_est
				for (tu = 0; tu < euser_num; tu++){
					// continue if elements are unobserved
					tn = ttensor[tu][ti][0];
					if (tn == -1.0) continue;
					for (tj = 0; tj < r_num; tj++){
						n_est[tu][ti][tj] += (v_f_t[ti][tk] - v_pre) * v_t_f[tj][tk];
					}
				}
			}
		}

		// Update v_t_f (To Region x From Region)
		// tk: Feature dimension
		for (tk = 0; tk < k_f_t; tk++){
			// tj: To Region
			for (tj = 0; tj < r_num; tj++){
				numerator = denominator = 0.0;
				// tu: User
				for (tu = 0; tu < euser_num; tu++){
					// ti: From Region
					for (ti = 0; ti < r_num; ti++){
						// continue if elements are unobserved
						tn = ttensor[tu][ti][0];
						if (tn == -1.0) continue;

						// Transition count --> tn
						tn = ttensor[tu][ti][tj];

						// Numerator --> numerator
						numerator += (tn - n_est[tu][ti][tj] + v_f_t[ti][tk] * v_t_f[tj][tk]) * v_f_t[ti][tk];
						// Denominator --> denominator
						denominator += v_f_t[ti][tk] * v_f_t[ti][tk];
					}
				}
				// Update v_t_f[tj][tk]
				if (denominator + alpha == 0.0) v_tmp = EPSILON;
				else{
					v_tmp = numerator / (denominator + alpha);
					if (v_tmp < EPSILON) v_tmp = EPSILON;	// NTF
				}

				chng = fabs(v_tmp - v_t_f[tj][tk]) / v_t_f[tj][tk];
				if (chng > ANLS_CHNG_THR) converge_tmp = 0;
				if (chng > chng_max) chng_max = chng;
				chng_avg += chng;
				chng_num++;
				v_pre = v_t_f[tj][tk];
				v_t_f[tj][tk] = v_tmp;

				// Estimate transition counts --> n_est
				for (tu = 0; tu < euser_num; tu++){
					for (ti = 0; ti < r_num; ti++){
						// continue if elements are unobserved
						tn = ttensor[tu][ti][0];
						if (tn == -1.0) continue;
						n_est[tu][ti][tj] += v_f_t[ti][tk] * (v_t_f[tj][tk] - v_pre);
					}
				}
			}
		}

		// Update v_u_t (User x To Region)
		// tk: Feature dimension
		for (tk = 0; tk < k_u_t; tk++){
			// tu: User
			for (tu = 0; tu < euser_num; tu++){
				numerator = denominator = 0.0;
				// ti: From Region
				for (ti = 0; ti < r_num; ti++){
					// continue if elements are unobserved
					tn = ttensor[tu][ti][0];
					if (tn == -1.0) continue;

					// tj: To Region
					for (tj = 0; tj < r_num; tj++){
						// Transition count --> tn
						tn = ttensor[tu][ti][tj];

						// Numerator --> numerator
						numerator += (tn - n_est[tu][ti][tj] + v_u_t[tu][tk] * v_t_u[tj][tk]) * v_t_u[tj][tk];
						// Denominator --> denominator
						denominator += v_t_u[tj][tk] * v_t_u[tj][tk];
					}
				}
				// Update v_u_t[tu][tk]
				if (denominator + alpha == 0.0) v_tmp = EPSILON;
				else{
					v_tmp = numerator / (denominator + alpha);
					if (v_tmp < EPSILON) v_tmp = EPSILON;	// NTF
				}

				chng = fabs(v_tmp - v_u_t[tu][tk]) / v_u_t[tu][tk];
				if (chng > ANLS_CHNG_THR) converge_tmp = 0;
				if (chng > chng_max) chng_max = chng;
				chng_avg += chng;
				chng_num++;
				v_pre = v_u_t[tu][tk];
				v_u_t[tu][tk] = v_tmp;

				// Estimate transition counts --> n_est
				for (ti = 0; ti < r_num; ti++){
					// continue if elements are unobserved
					tn = ttensor[tu][ti][0];
					if (tn == -1.0) continue;
					for (tj = 0; tj < r_num; tj++){
						n_est[tu][ti][tj] += (v_u_t[tu][tk] - v_pre) * v_t_u[tj][tk];
					}
				}
			}
		}

		// Update v_t_u (To Region x User)
		// tk: Feature dimension
		for (tk = 0; tk < k_u_t; tk++){
			// tj: To Region
			for (tj = 0; tj < r_num; tj++){
				numerator = denominator = 0.0;
				// tu: User
				for (tu = 0; tu < euser_num; tu++){
					// ti: From Region
					for (ti = 0; ti < r_num; ti++){
						// continue if elements are unobserved
						tn = ttensor[tu][ti][0];
						if (tn == -1.0) continue;

						// Transition count --> tn
						tn = ttensor[tu][ti][tj];

						// Numerator --> numerator
						numerator += (tn - n_est[tu][ti][tj] + v_u_t[tu][tk] * v_t_u[tj][tk]) * v_u_t[tu][tk];
						// Denominator --> denominator
						denominator += v_u_t[tu][tk] * v_u_t[tu][tk];
					}
				}
				// Update v_t_u[tj][tk]
				if (denominator + alpha == 0.0) v_tmp = EPSILON;
				else{
					v_tmp = numerator / (denominator + alpha);
					if (v_tmp < EPSILON) v_tmp = EPSILON;	// NTF
				}

				chng = fabs(v_tmp - v_t_u[tj][tk]) / v_t_u[tj][tk];
				if (chng > ANLS_CHNG_THR) converge_tmp = 0;
				if (chng > chng_max) chng_max = chng;
				chng_avg += chng;
				chng_num++;
				v_pre = v_t_u[tj][tk];
				v_t_u[tj][tk] = v_tmp;

				// Estimate transition counts --> n_est
				for (tu = 0; tu < euser_num; tu++){
					for (ti = 0; ti < r_num; ti++){
						// continue if elements are unobserved
						tn = ttensor[tu][ti][0];
						if (tn == -1.0) continue;
						n_est[tu][ti][tj] += v_u_t[tu][tk] * (v_t_u[tj][tk] - v_pre);
					}
				}
			}
		}

		chng_avg /= (double)chng_num;

		j++;

		if (chng_avg < ANLS_CHNG_THR) converge++;
		if (converge == 1 || j == MAX_ITERATION) break;
	}
	if (j == MAX_ITERATION) printf("TF: Reached the maximum number of iterations (=%d).\n", MAX_ITERATION);
	else printf("TF: Converged at the %d-th iterations\n", j);

	// Free
	free3D(n_est, euser_num, r_num);
}