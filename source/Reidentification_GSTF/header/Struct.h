/**
** @file Struct.h
** 
** Copyright (c) 2017 National Institute of Advanced Industrial Science and Technology (AIST)
** 
** This software is released under the MIT License.
** http://opensource.org/licenses/mit-license.php
*/

#ifndef STRUCT
#define STRUCT

#include "parameter.h"

/** Trace structure */
struct trace{
	/** Vertical region no. */
	int *v_no;
	/** Horizontal region no. */
	int *h_no;
	/** Region no. */
	int *r_no;
	/** Latitude */
	double *v;
	/** Longitude */
	double *h;
	/** Number of events */
	int event_num;
};

/** Transition pattern structure */
struct trans{
	/** From Region no. */
	int fr_no;
	/** To Region no. */
	int tr_no;
	/** Transition count */
	double num;
};

/** region group structure */
struct rg_type{
	/** Number of regions in the group*/
	int size;
	/** Region no. */
	int *r_no;
};

/** INI parameter structure */
struct iniParameter{
	/** extracted location data */
	char data_extr[1024];
	/** MCL(Markov Clustering) directory */
	char mcl_dir[1024];
	/** Rate of training events */
	double tevent_rate;
	/** Rate of training events (char array) */
	char tevent_rate_s[1024];
	/** Number of target users */
	int euser_num;
	/** Total number of users (target users and non-target users) */
	int tuser_num;
	/** Number of ways to randomly divide users into target users and non-target users */
	int euser_rand_num;
	/** Number of regions in a vertical way */
	int v_num;
	/** Number of regions in a horizontal way */
	int h_num;
	/** Minimum latitude */
	double v_min;
	/** Maximum latitude */
	double v_max;
	/** Minimum longitude */
	double h_min;
	/** Maximum longitude */
	double h_max;
	/** Trainng method (0:ML, 1:TF, 2:GSTF) */
	int train_method;
	/** Whether to make the MCL files (0: do not make,  1: make (when TRAIN_METHOD = 1)) */
	int MCL_make_flag;
	/** Feature dimensionality of model parameters between USER and From Region */
	int k_u_f;
	/** Feature dimensionality of model parameters between From Region and To Region */
	int k_f_t;
	/** Feature dimensionality of model parameters between To Region and USER */
	int k_u_t;
	/** Regularization parameter alpha */
	double alpha;
	/** Regularization parameter alpha (char array) */
	char alpha_s[1024];
	/** Regularization parameter beta */
	double beta;
	/** Regularization parameter beta (char array) */
	char beta_s[1024];
};

#endif