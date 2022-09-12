/**
** @file CommonFunction.cpp
** 
** @brief
** Common functions.
** 
** Copyright (c) 2017 National Institute of Advanced Industrial Science and Technology (AIST)
** 
** This software is released under the MIT License.
** http://opensource.org/licenses/mit-license.php
*/

#include <stdio.h>
#include "mt19937ar.c"
#include "MemoryOperation.h"

/**
** @brief 			Randomly generate a permutation of 0, 1, 2, ..., x-1.
** @param rndperm 	Permutation of 0, 1, 2, ..., x-1
** @param x	 		Natural number
*/
void RndPerm(int *rndperm, int x){
	int rnd;
	int *ordperm;
	int i, j;

	malloc1D(&ordperm, x);
	for (i = 0; i < x; i++) ordperm[i] = i;

	for (i = 0; i < x; i++){
		rnd = genrand_int32() % (x - i);
		rndperm[i] = ordperm[rnd];
		for (j = rnd + 1; j < x - i; j++){
			ordperm[j - 1] = ordperm[j];
		}
	}

	free1D(ordperm);
}

/**
** @brief 			Get a random value.
** @return			Random value
*/
double getRandomValue()
{
	return genrand_real3();
}

/**
** @brief 				Comparator function used in qsort.
** @param a1	 		1st value
** @param a2	 		2st value
** @return				1: (*a1 > *a2), -1: (*a1 <= *a2)
*/
int compare_double(const void *a1, const void *a2){
	double *b1, *b2;
	b1 = (double *)a1;
	b2 = (double *)a2;
	if (*b1 > *b2) return 1;
	else return -1;
}

/**
** @brief 			Exit.
** @param rc	 	Return code
*/
void Exit(int rc){
	if (rc != 0)
		printf("\nERROR CODE:%d\nPress ENTER to exit.\n",rc);
	else
		printf("\nPress ENTER to exit.\n");
	getchar();
	exit(rc);
}
