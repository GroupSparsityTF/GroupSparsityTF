/**
** @file parameter.h
** 
** Copyright (c) 2017 National Institute of Advanced Industrial Science and Technology (AIST)
** 
** This software is released under the MIT License.
** http://opensource.org/licenses/mit-license.php
*/

#ifndef PARAMETER
#define PARAMETER

/** INI file path */
#define INI_PATH		".\\Reidentification_GSTF.ini"
/** Maximum number of users */
#define MAX_USER_NUM	400
/** Maximum number of events per trace */
#define MAX_EVENT_NUM	1000
/** Minimum number of events */
#define MIN_EVENT_NUM	100

/** Threshold for the rate of change in ANLS (used in TF/GSTF) */
#define ANLS_CHNG_THR	0.001
/** Maximum number of iterations (used in TF/GSTF) */
#define MAX_ITERATION	200
/** Threshold for the rate of change of prior probabilities */
#define PRIOR_CHNG_THR	0.00001
/** Epsilon value in NTF (Non-negative Tensor Factorization) */
#define EPSILON			1.0e-16
/** Minimum transition probability */
#define PROB_MIN		1.0e-8
/** Maximum number of regions */
#define MAX_REGION_NUM	256

#endif