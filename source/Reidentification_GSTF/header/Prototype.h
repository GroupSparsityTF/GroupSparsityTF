/**
** @file Prototype.h
** 
** Copyright (c) 2017 National Institute of Advanced Industrial Science and Technology (AIST)
** 
** This software is released under the MIT License.
** http://opensource.org/licenses/mit-license.php
*/

#ifndef PROTOTYPE
#define PROTOTYPE

#include "Struct.h"

/* CommonFunction.cpp */
void RndPerm(int *rndperm, int x);
double getRandomValue();
int compare_double(const void *a1, const void *a2);
void Exit(int rc);

/* FileIO.cpp */
iniParameter readIniParameter();
void FileOpen(FILE **file, char *filename, char *mode);
void fileout(FILE *fp, const char * Format, ...);

/* TF.cpp */
void TF(
	double **v_u_f, double **v_f_u,
	double **v_f_t,	double **v_t_f,
	double **v_u_t,	double **v_t_u, 
	int r_num, int euser_num,
	int k_u_f, int k_f_t, int k_u_t,
	double ***ttensor, double alpha
	);

/* GSTF.cpp */
void GSTF(
	double **v_u_f, double **v_f_u,
	double **v_f_t, double **v_t_f, 
	double **v_u_t, double **v_t_u, 
	int r_num, int euser_num,
	int k_u_f, int k_f_t, int k_u_t,
	double ***ttensor, double alpha,
	double beta, struct rg_type *rg, int rg_num
	);

/* Train.cpp */
void AssignRegionID(struct trace **tetrace, int user_num, int *trace_num, int v_num, int h_num, double v_min, double v_max, double h_min, double h_max);
void MakeTrans(struct trace **tetrace, int user_num, int *ttrace_num, int *tevent_num, struct trans **ttrans, int *ttrans_num);
void Trans2Tensor(struct trans **ttrans, int *ttrans_num, int *user_rnd, int euser_num, int r_num, double ***ttensor, double **tmatrix_all);
void TFParam2Tensor(double **v_u_f, double **v_f_u, double **v_f_t, double **v_t_f, double **v_u_t, double **v_t_u, int k_u_f, int k_f_t, int k_u_t, int euser_num, int r_num, double ***ttensor);
void ML(double ***ttensor, int euser_num, int r_num);
void Tensor2Prior(double ***ttensor, int euser_num, int r_num, double **prior);
void CalcPopMat(double ***ttensor, double ***ttensor_n, int euser_num, int r_num);

/* Reidentification.cpp */
void CalcLogLikeli(double ***ttensor, double ***ttensor_n, double **prior, double **prior_n, int user_num, int euser_num, int *user_rnd, int r_num, int v_num, int h_num,struct trace **tetrace, int *ttrace_num, int *etrace_num, double ***loglikeli, double ***loglikeli_n);
void CalcFPRTPR(char *outfile, double *e_loglikeli_ratio, int e_llr_num, double *n_loglikeli_ratio, int n_llr_num, double max_sim, double min_sim, double thr_shift);

#endif