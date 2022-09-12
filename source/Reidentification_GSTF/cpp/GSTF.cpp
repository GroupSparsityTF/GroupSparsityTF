/**
** @file GSTF.cpp
** 
** @brief
** Group Sparsity Tensor Factorization (i.e. training method in [Murakami+, TIFS17]).
** 
** Copyright (c) 2017 National Institute of Advanced Industrial Science and Technology (AIST)
** 
** This software is released under the MIT License.
** http://opensource.org/licenses/mit-license.php
*/

#include <stdio.h>
#include <math.h>
#include "MemoryOperation.h"
#include "parameter.h"
#include "Prototype.h"
#include "Struct.h"

/**
** @brief 			Train model parameters (i.e. v_u_f, v_f_u, v_f_t, v_t_f, v_u_t, v_t_u) in GSTF (Group Sparsity Tensor Factorization) 
					from personalized transition matrices (i.e. ttensor).
** @param v_u_f 	Model parameters (User x From Region)
** @param v_f_u 	Model parameters (From Region x User)
** @param v_f_t 	Model parameters (From Region x To Region)
** @param v_t_f 	Model parameters (To Region x From Region)
** @param v_u_t 	Model parameters (User x To Region)
** @param v_t_u 	Model parameters (To Region x User)
** @param r_num 	Number of regions
** @param euser_num Number of target users
** @param k_u_f 	Feature dimensionality of model parameters between USER and From Region (i.e. v_u_f, v_f_u)
** @param k_f_t 	Feature dimensionality of model parameters between From Region and To Region (i.e. v_f_t, v_t_f)
** @param k_u_t 	Feature dimensionality of model parameters between To Region and USER (i.e. v_u_t, v_t_u)
** @param ttensor 	Transition count tensor (User x From Region x To Region)
** @param alpha 	Regularization parameter alpha
** @param beta	 	Regularization parameter beta
** @param rg 		Region groups
** @param rg_num 	Number of region groups
*/
void GSTF(double **v_u_f, double **v_f_u, double **v_f_t, double **v_t_f, double **v_u_t, double **v_t_u, int r_num, int euser_num,
	int k_u_f, int k_f_t, int k_u_t, double ***ttensor, double alpha, double beta, struct rg_type *rg, int rg_num){

	double ***a_hat;
	int tu, ti, tj, tk;
	double a;
	double err = 0.0, err_pre = 0.0;
	double numerator, denominator;
	double v_tmp, v_pre;
	double chng_avg;
	int chng_num;
	int i, j, k;
	int converge = 0;
	int tg, tr;
	double *mu, *lambda, *s, *s_ast;
	double s_l2;

	// Malloc
	malloc3D(&a_hat, euser_num, r_num, r_num);
	malloc1D(&mu, r_num);
	malloc1D(&lambda, r_num);
	malloc1D(&s, r_num);
	malloc1D(&s_ast, r_num);

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

	// Estimate transition counts --> a_hat
	for (tu = 0; tu < euser_num; tu++){
		for (ti = 0; ti < r_num; ti++){
			for (tj = 0; tj < r_num; tj++){
				a_hat[tu][ti][tj] = 0.0;
				for (k = 0; k < k_u_f; k++) a_hat[tu][ti][tj] += v_u_f[tu][k] * v_f_u[ti][k];
				for (k = 0; k < k_f_t; k++) a_hat[tu][ti][tj] += v_f_t[ti][k] * v_t_f[tj][k];
				for (k = 0; k < k_u_t; k++) a_hat[tu][ti][tj] += v_u_t[tu][k] * v_t_u[tj][k];
			}
		}
	}

	// Train model parameters until convergence
	j = 0;
	while (1){
		chng_avg = 0.0;
		chng_num = 0;

		// Update v_u_f (User x From Region)
		// tk: Feature dimension
		for (tk = 0; tk < k_u_f; tk++){
			// tu: User
			for (tu = 0; tu < euser_num; tu++){
				numerator = denominator = 0.0;
				// ti: From Region
				for (ti = 0; ti < r_num; ti++){
					// Continue if elements are unobserved
					a = ttensor[tu][ti][0];
					if (a == -1.0) continue;

					// tj: To Region
					for (tj = 0; tj < r_num; tj++){
						// Transition count --> a
						a = ttensor[tu][ti][tj];

						// Numerator --> numerator
						numerator += (a - a_hat[tu][ti][tj] + v_u_f[tu][tk] * v_f_u[ti][tk]) * v_f_u[ti][tk];
						// Denominator --> denominator
						denominator += v_f_u[ti][tk] * v_f_u[ti][tk];
					}
				}

				// Update v_u_f[tu][tk] (= [numerator]+ / (denominator + alpha))
				if (denominator + alpha == 0.0) v_tmp = EPSILON;
				else{
					v_tmp = numerator / (denominator + alpha);
					if (v_tmp < EPSILON) v_tmp = EPSILON;	// NTF
				}
				chng_avg += fabs(v_tmp - v_u_f[tu][tk]) / v_u_f[tu][tk];
				chng_num++;
				v_pre = v_u_f[tu][tk];
				v_u_f[tu][tk] = v_tmp;

				// Estimate transition counts --> a_hat
				for (ti = 0; ti < r_num; ti++){
					// Continue if elements are unobserved
					a = ttensor[tu][ti][0];
					if (a == -1.0) continue;
					for (tj = 0; tj < r_num; tj++){
						a_hat[tu][ti][tj] += v_u_f[tu][tk] * v_f_u[ti][tk] - v_pre * v_f_u[ti][tk];
					}
				}
			}
		}

		// Update v_f_u (From Region x User)
		// tk: Feature dimension
		for (tk = 0; tk < k_u_f; tk++){
			// tg: group of From Region
			for (tg = 0; tg < rg_num; tg++){
				// Compute [lambda]+, mu, s, and the L2 norm of s --> lambda, mu, s, s_l2
				s_l2 = 0.0;
				for (tr = 0; tr < rg[tg].size; tr++){
					// ti: From Region
					ti = rg[tg].r_no[tr];

					numerator = denominator = 0.0;
					// tu: User
					for (tu = 0; tu < euser_num; tu++){
						// Continue if elements are unobserved
						a = ttensor[tu][ti][0];
						if (a == -1.0) continue;

						// tj: To Region
						for (tj = 0; tj < r_num; tj++){
							// Transition count --> a
							a = ttensor[tu][ti][tj];

							// Numerator --> numerator
							numerator += (a - a_hat[tu][ti][tj] + v_u_f[tu][tk] * v_f_u[ti][tk]) * v_u_f[tu][tk];
							// Denominator --> denominator
							denominator += v_u_f[tu][tk] * v_u_f[tu][tk];
						}
					}

					if (denominator == 0.0){
						mu[tr] = EPSILON;
						lambda[tr] = EPSILON;
					}
					else{
						mu[tr] = denominator;
						lambda[tr] = numerator / denominator;
						if (lambda[tr] < EPSILON) lambda[tr] = EPSILON;
					}
					s[tr] = 2.0 * mu[tr] * lambda[tr];
					s_l2 += s[tr] * s[tr];
				}
				s_l2 = sqrt(s_l2);

				for (tr = 0; tr < rg[tg].size; tr++){
					// ti: From Region
					ti = rg[tg].r_no[tr];

					// Update v_f_u[ti][tk] (group sparse solution)
					if (s_l2 <= beta){
						v_tmp = EPSILON;
					}
					// Update v_f_u[ti][tk] (= lambda[tr] - s_ast[tr] / (2.0 * mu[tr]))
					else{
						s_ast[tr] = (beta / s_l2) * s[tr];
						v_tmp = lambda[tr] - s_ast[tr] / (2.0 * mu[tr]);
						if (v_tmp < EPSILON) v_tmp = EPSILON;	// NTF
					}
					chng_avg += fabs(v_tmp - v_f_u[ti][tk]) / v_f_u[ti][tk];
					chng_num++;
					v_pre = v_f_u[ti][tk];
					v_f_u[ti][tk] = v_tmp;

					// Estimate transition counts --> a_hat
					for (tu = 0; tu < euser_num; tu++){
						// Continue if elements are unobserved
						a = ttensor[tu][ti][0];
						if (a == -1.0) continue;
						for (tj = 0; tj < r_num; tj++){
							a_hat[tu][ti][tj] += v_u_f[tu][tk] * (v_f_u[ti][tk] - v_pre);
						}
					}
				}
			}
		}

		// Update v_f_t (From Region x To Region)
		// tk: Feature dimension
		for (tk = 0; tk < k_f_t; tk++){
			// tg: group of From Region
			for (tg = 0; tg < rg_num; tg++){
				// Compute [lambda]+, mu, s, and the L2 norm of s --> lambda, mu, s, s_l2
				s_l2 = 0.0;
				for (tr = 0; tr < rg[tg].size; tr++){
					// ti: From Region
					ti = rg[tg].r_no[tr];

					numerator = denominator = 0.0;
					// tu: User
					for (tu = 0; tu < euser_num; tu++){
						// Continue if elements are unobserved
						a = ttensor[tu][ti][0];
						if (a == -1.0) continue;

						// tj: To Region
						for (tj = 0; tj < r_num; tj++){
							// Transition count --> a
							a = ttensor[tu][ti][tj];

							// Numerator --> numerator
							numerator += (a - a_hat[tu][ti][tj] + v_f_t[ti][tk] * v_t_f[tj][tk]) * v_t_f[tj][tk];
							// Denominator --> denominator
							denominator += v_t_f[tj][tk] * v_t_f[tj][tk];
						}
					}

					if (denominator == 0.0){
						mu[tr] = EPSILON;
						lambda[tr] = EPSILON;
					}
					else{
						mu[tr] = denominator;
						lambda[tr] = numerator / denominator;
						if (lambda[tr] < EPSILON) lambda[tr] = EPSILON;
					}
					s[tr] = 2.0 * mu[tr] * lambda[tr];
					s_l2 += s[tr] * s[tr];
				}
				s_l2 = sqrt(s_l2);

				for (tr = 0; tr < rg[tg].size; tr++){
					// ti: From Region
					ti = rg[tg].r_no[tr];

					// Update v_f_t[ti][tk] (group sparse solution)
					if (s_l2 <= beta){
						v_tmp = EPSILON;
					}
					// Update v_f_t[ti][tk] (= lambda[tr] - s_ast[tr] / (2.0 * mu[tr]))
					else{
						s_ast[tr] = (beta / s_l2) * s[tr];
						v_tmp = lambda[tr] - s_ast[tr] / (2.0 * mu[tr]);
						if (v_tmp < EPSILON) v_tmp = EPSILON;	// NTF
					}
					chng_avg += fabs(v_tmp - v_f_t[ti][tk]) / v_f_t[ti][tk];
					chng_num++;
					v_pre = v_f_t[ti][tk];
					v_f_t[ti][tk] = v_tmp;

					// Estimate transition counts --> a_hat
					for (tu = 0; tu < euser_num; tu++){
						// Continue if elements are unobserved
						a = ttensor[tu][ti][0];
						if (a == -1.0) continue;
						for (tj = 0; tj < r_num; tj++){
							a_hat[tu][ti][tj] += (v_f_t[ti][tk] - v_pre) * v_t_f[tj][tk];
						}
					}
				}
			}
		}

		// Update v_t_f (To Region x From Region)
		// tk: Feature dimension
		for (tk = 0; tk < k_f_t; tk++){
			// tg: group of To Region
			for (tg = 0; tg < rg_num; tg++){
				// Compute [lambda]+, mu, s, and the L2 norm of s --> lambda, mu, s, s_l2
				s_l2 = 0.0;
				for (tr = 0; tr < rg[tg].size; tr++){
					// tj: To Region
					tj = rg[tg].r_no[tr];

					numerator = denominator = 0.0;
					// tu: User
					for (tu = 0; tu < euser_num; tu++){
						// ti: From Region
						for (ti = 0; ti < r_num; ti++){
							// Continue if elements are unobserved
							a = ttensor[tu][ti][0];
							if (a == -1.0) continue;

							// Transition count --> a
							a = ttensor[tu][ti][tj];

							// Numerator --> numerator
							numerator += (a - a_hat[tu][ti][tj] + v_f_t[ti][tk] * v_t_f[tj][tk]) * v_f_t[ti][tk];
							// Denominator --> denominator
							denominator += v_f_t[ti][tk] * v_f_t[ti][tk];
						}
					}

					if (denominator == 0.0){
						mu[tr] = EPSILON;
						lambda[tr] = EPSILON;
					}
					else{
						mu[tr] = denominator;
						lambda[tr] = numerator / denominator;
						if (lambda[tr] < EPSILON) lambda[tr] = EPSILON;
					}
					s[tr] = 2.0 * mu[tr] * lambda[tr];
					s_l2 += s[tr] * s[tr];
				}
				s_l2 = sqrt(s_l2);

				for (tr = 0; tr < rg[tg].size; tr++){
					// tj: To Region
					tj = rg[tg].r_no[tr];

					// Update v_t_f[ti][tk] (group sparse solution)
					if (s_l2 <= beta){
						v_tmp = EPSILON;
					}
					// Update v_t_f[ti][tk] (= lambda[tr] - s_ast[tr] / (2.0 * mu[tr]))
					else{
						s_ast[tr] = (beta / s_l2) * s[tr];
						v_tmp = lambda[tr] - s_ast[tr] / (2.0 * mu[tr]);
						if (v_tmp < EPSILON) v_tmp = EPSILON;	// NTF
					}
					chng_avg += fabs(v_tmp - v_t_f[tj][tk]) / v_t_f[tj][tk];
					chng_num++;
					v_pre = v_t_f[tj][tk];
					v_t_f[tj][tk] = v_tmp;

					// Estimate transition counts --> a_hat
					for (tu = 0; tu < euser_num; tu++){
						for (ti = 0; ti < r_num; ti++){
							// Continue if elements are unobserved
							a = ttensor[tu][ti][0];
							if (a == -1.0) continue;
							a_hat[tu][ti][tj] += v_f_t[ti][tk] * (v_t_f[tj][tk] - v_pre);
						}
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
					// Continue if elements are unobserved
					a = ttensor[tu][ti][0];
					if (a == -1.0) continue;

					// To Region
					for (tj = 0; tj < r_num; tj++){
						// Transition count --> a
						a = ttensor[tu][ti][tj];

						// Numerator --> numerator
						numerator += (a - a_hat[tu][ti][tj] + v_u_t[tu][tk] * v_t_u[tj][tk]) * v_t_u[tj][tk];
						// Denominator --> denominator
						denominator += v_t_u[tj][tk] * v_t_u[tj][tk];
					}
				}

				// Update v_u_t[tu][tk] (= [numerator]+ / (denominator + alpha))
				if (denominator + alpha == 0.0) v_tmp = EPSILON;
				else{
					v_tmp = numerator / (denominator + alpha);
					if (v_tmp < EPSILON) v_tmp = EPSILON;	// NTF
				}
				chng_avg += fabs(v_tmp - v_u_t[tu][tk]) / v_u_t[tu][tk];
				chng_num++;
				v_pre = v_u_t[tu][tk];
				v_u_t[tu][tk] = v_tmp;

				// Estimate transition counts --> a_hat
				for (ti = 0; ti < r_num; ti++){
					// Continue if elements are unobserved
					a = ttensor[tu][ti][0];
					if (a == -1.0) continue;
					for (tj = 0; tj < r_num; tj++){
						a_hat[tu][ti][tj] += (v_u_t[tu][tk] - v_pre) * v_t_u[tj][tk];
					}
				}
			}
		}

		// Update v_t_u (To Region x User)
		// tk: Feature dimension
		for (tk = 0; tk < k_u_t; tk++){
			// tg: group of To Region
			for (tg = 0; tg < rg_num; tg++){
				// Compute [lambda]+, mu, s, and the L2 norm of s --> lambda, mu, s, s_l2
				s_l2 = 0.0;
				for (tr = 0; tr < rg[tg].size; tr++){
					// tj: To Region
					tj = rg[tg].r_no[tr];

					numerator = denominator = 0.0;
					// User
					for (tu = 0; tu < euser_num; tu++){
						// From Region
						for (ti = 0; ti < r_num; ti++){
							// Continue if elements are unobserved
							a = ttensor[tu][ti][0];
							if (a == -1.0) continue;

							// Transition count --> a
							a = ttensor[tu][ti][tj];

							// Numerator --> numerator
							numerator += (a - a_hat[tu][ti][tj] + v_u_t[tu][tk] * v_t_u[tj][tk]) * v_u_t[tu][tk];
							// Denominator --> denominator
							denominator += v_u_t[tu][tk] * v_u_t[tu][tk];
						}
					}

					if (denominator == 0.0){
						mu[tr] = EPSILON;
						lambda[tr] = EPSILON;
					}
					else{
						mu[tr] = denominator;
						lambda[tr] = numerator / denominator;
						if (lambda[tr] < EPSILON) lambda[tr] = EPSILON;
					}
					s[tr] = 2.0 * mu[tr] * lambda[tr];
					s_l2 += s[tr] * s[tr];
				}
				s_l2 = sqrt(s_l2);

				for (tr = 0; tr < rg[tg].size; tr++){
					// tj: To Region
					tj = rg[tg].r_no[tr];

					// Update v_t_u[tj][tk] (group sparse solution)
					if (s_l2 <= beta){
						v_tmp = EPSILON;
					}
					// Update v_t_u[tj][tk] (= lambda[tr] - s_ast[tr] / (2.0 * mu[tr]))
					else{
						s_ast[tr] = (beta / s_l2) * s[tr];
						v_tmp = lambda[tr] - s_ast[tr] / (2.0 * mu[tr]);
						if (v_tmp < EPSILON) v_tmp = EPSILON;	// NTF
					}
					chng_avg += fabs(v_tmp - v_t_u[tj][tk]) / v_t_u[tj][tk];
					chng_num++;
					v_pre = v_t_u[tj][tk];
					v_t_u[tj][tk] = v_tmp;

					// Estimate transition counts --> a_hat
					for (tu = 0; tu < euser_num; tu++){
						for (ti = 0; ti < r_num; ti++){
							// Continue if elements are unobserved
							a = ttensor[tu][ti][0];
							if (a == -1.0) continue;
							a_hat[tu][ti][tj] += v_u_t[tu][tk] * (v_t_u[tj][tk] - v_pre);
						}
					}
				}
			}
		}

		chng_avg /= (double)chng_num;

		j++;

		if (chng_avg < ANLS_CHNG_THR) converge = 1;
		if (converge == 1 || j == MAX_ITERATION) break;
	}
	if (j == MAX_ITERATION) printf("GSTF: Reached the maximum number of iterations (=%d).\n", MAX_ITERATION);
	else printf("GSTF: Converged at the %d-th iterations\n", j);

	// Free
	free3D(a_hat, euser_num, r_num);
	free1D(mu);
	free1D(lambda);
	free1D(s);
	free1D(s_ast);
}
