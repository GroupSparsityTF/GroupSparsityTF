/**
** @file main.cpp
** 
** @brief
** Main function.
** 
** Copyright (c) 2017 National Institute of Advanced Industrial Science and Technology (AIST)
** 
** This software is released under the MIT License.
** http://opensource.org/licenses/mit-license.php
*/

#include <algorithm>
#include <Windows.h>
#include <direct.h>
#include <shlwapi.h>
#include <time.h>
#include "MemoryOperation.h"
#include "parameter.h"
#include "Prototype.h"
#include "Struct.h"

/** Traces (training traces and testing traces) */
struct trace **tetrace;

/**
** @brief 			Perform experiments in a biased setting (i.e. Section VII in [Murakami+, TIFS17]).
** @return			0: normal, -1: abnormal
*/
int main(void){
	int r_num;								// Number of regions
	int user_num = 0;						// Number of users
	int nuser_num;							// Number of non-target users

	int trace_num[MAX_USER_NUM] = { 0 };	// Number of traces
	int ttrace_num[MAX_USER_NUM]; 			// Number of training traces
	int etrace_num[MAX_USER_NUM]; 			// Number of testing traces
	int etrace_sum;							// Total number of testing traces

	int tot_event_num[MAX_USER_NUM];		// Total number of events
	int tevent_num[MAX_USER_NUM];			// Number of training events

	int trace_hist[MAX_EVENT_NUM];			// Histogram (the number of traces)
	int event_hist[MAX_EVENT_NUM];			// Histogram (the number of events)

	int **user_rnd;							// Randomized user IDs

	struct trans **ttrans;					// Transition pattern (User x Transition Pattern)
	int *ttrans_num;						// Number of transition patterns

	// Model parameters in PITF
	double **v_u_f;							// User - From Region
	double **v_f_u;							// From Region - User
	double **v_f_t;							// From Region - To Region
	double **v_t_f;							// To Region - From Region
	double **v_u_t;							// User - To Region
	double **v_t_u;							// To Region - User

	struct rg_type *rg;						// Region group
	int rg_num = 0;							// Number of region groups

	double **tmatrix_all;					// Transition count matrix for all users (From Region x To Region)

	double ***ttensor;						// Personalized transition matrices (User x From Region x To Region)
	double ***ttensor_n;					// Population transition matrices (User x From Region x To Region)
	double **prior;							// Personalized prior distribution (User x Region)
	double **prior_n;						// Population prior distribution (User x Region)

	// Log-likelihood computed using the personalized transitoin matrix (Testing User x Testing Trace x Target User)
	double ***loglikeli;
	// Log-likelihood computed using the population transition matrix (Testing User x Testing Trace x Target User)
	double ***loglikeli_n;

	// Log-likelihood ratios
	double *e_loglikeli_ratio;		// Log-likelihood ratio for each target user
	double *n_loglikeli_ratio;		// Log-likelihood ratio for each non-target user
	double *e_loglikeli_ratio_ast;	// Log-likelihood ratio for each target user (when the testing traces are combined into one)
	double *n_loglikeli_ratio_ast;	// Log-likelihood ratio for each non-target user (when the testing traces are combined into one)
	int e_llr_num = 0; 				// Number of log-likelihood ratios for target users
	int n_llr_num = 0;				// Number of log-likelihood ratios for non-target users
	int e_llr_ast_num = 0;			// Number of log-likelihood ratios for target users (when the testing traces are combined into one)
	int n_llr_ast_num = 0;			// Number of log-likelihood ratios for non-target users (when the testing traces are combined into one)

	// trace type based on the trace length (the number of events)
	int trace_type;
	int trace_interval[5] = { 5, 10, 15, 20, 25 };

	// Log-likelihood ratios for each trace type (classified based on the trace length)
	double **e_loglikeli_ratio_tr;
	double **n_loglikeli_ratio_tr;
	int e_llr_num_tr[6] = { 0 };
	int n_llr_num_tr[6] = { 0 };

	// others
	double tv_tmp[MAX_EVENT_NUM], th_tmp[MAX_EVENT_NUM];
	char file_extr[1024];
	int r_trace_num, r_event_num;
	double tv, th;
	char  *tok, *ctx;
	int tim;
	double llr_sum;
	char s[1024];
	char outfile[1024];
	FILE *fp,*fp2;
	errno_t error;
	int i, j, u;
	int en;
	int etn;
	int event_num_tmp;

	/******************************************************* Read INI File *******************************************************/
	// INI file parameters --> ini
	struct iniParameter ini = readIniParameter();

	// Total number of regions --> r_num
	r_num =ini.v_num *ini.h_num;

	if (r_num > MAX_REGION_NUM){
		printf("The total number of regions exceeds the maximum value (=%d).\n", MAX_REGION_NUM);
		Exit(-1);
	}

	/******************************************************** Read Traces ********************************************************/
	HANDLE hFind;
	WIN32_FIND_DATA fd;
	int event_total;
	int err_code;

	// Find a trace file
	sprintf_s(file_extr, "%s\\*.txt", ini.data_extr);
	hFind = FindFirstFile(file_extr, &fd);
	if (hFind == INVALID_HANDLE_VALUE) {
		printf("No user files.\nWildCard:%s\n",file_extr);
		Exit(-1);
	}

	// Malloc (tetrace)
	tetrace = (struct trace **)malloc(sizeof(struct trace *) * MAX_USER_NUM);

	// Read traces for each user (i: user ID)
	i = 0;

	FileOpen(&fp2, "read_trace_files.log", "w");

	do{
		event_total = 0;
		err_code = 0;
		sprintf_s(file_extr, "%s\\%s", ini.data_extr, fd.cFileName);

		// Open the trace file
		if (error = fopen_s(&fp, file_extr, "r") != 0){
			fileout(fp2, "not found:%s\n", fd.cFileName);
			continue;
		}

		// +++++++++++ Check whether the total number of events is more than or equal to MIN_EVENT_NUM. +++++++++++ //
		// Read each location from the trace file. 
		r_trace_num = r_event_num = 0;
		while (fgets(s, 1023, fp) != NULL){
			// Ignore the events after the number of events exceeds MAX_EVENT_NUM.
			if (r_event_num >= MAX_EVENT_NUM) continue;

			// End of the trace (line break)
			if (strcmp(s, "\n") == 0){
				// Update the total number of events --> event_total
				event_total += r_event_num;
				// If the number of events in the trace is more than one, increase the number of traces by one.
				if (r_event_num >= 2) r_trace_num++;
				// Initialize the number of events in the trace
				r_event_num = 0;
				continue;
			}

			// Latitude, latitude, time --> tv, th, tim
			if (EOF == sscanf_s(s, "%lf %lf %*d %d", &tv, &th, &tim)){
				err_code = -1;
				break;
			}

			// Do not read a location out of range
			if ((tv < ini.v_min || tv >= ini.v_max) || (th < ini.h_min || th >= ini.h_max)){
				// Update the total number of events --> event_total
				event_total += r_event_num;
				// If the number of events in the trace is more than one, increase the number of traces by one.
				if (r_event_num >= 2) r_trace_num++;
				// Initialize the number of events in the trace
				r_event_num = 0;
				continue;
			}

			// Update the number of events --> r_event_num
			r_event_num++;
		}

		// If the file format is incorrect, read the next file.
		if (err_code != 0) continue;

		// Update the total number of events --> event_total
		event_total += r_event_num;
		// If the number of events in the trace is more than one, increase the number of traces by one.
		if (r_event_num >= 2) r_trace_num++;

		// Number of traces --> trace_num[i]
		trace_num[i] = r_trace_num;

		// Malloc (tetrace)
		tetrace[i] = (struct trace *)malloc(sizeof(struct trace) * trace_num[i]);

		// Close the trace file
		fclose(fp);

		// If the number of events is smaller than MIN_EVENT_NUM, read the next file.
		if (event_total < MIN_EVENT_NUM){
			fileout(fp2, "The number of events is less than %d:%s\n", MIN_EVENT_NUM, fd.cFileName);
			continue;
		}

		// +++++++++++ If the number of events is more than or equal to MIN_EVENT_NUM, read the traces. +++++++++++ //
		fileout(fp2,"Add:%s\n", fd.cFileName);
		
		// Open the trace file
		if ((error = fopen_s(&fp, file_extr, "r")) != 0) continue;

		// Read each location from the trace file
		r_trace_num = r_event_num = 0;
		while (fgets(s, 1023, fp) != NULL){
			// Ignore the events after the number of events exceeds MAX_EVENT_NUM.
			if (r_event_num >= MAX_EVENT_NUM) continue;

			// End of the trace (line break)
			if (strcmp(s, "\n") == 0){
				// If the number of events in the trace is more than one, read the trace.
				if (r_event_num >= 2){
					// Malloc (tetrace)
					malloc1D(&(tetrace[i][r_trace_num].v_no), r_event_num);
					malloc1D(&(tetrace[i][r_trace_num].h_no), r_event_num);
					malloc1D(&(tetrace[i][r_trace_num].r_no), r_event_num);
					malloc1D(&(tetrace[i][r_trace_num].v), r_event_num);
					malloc1D(&(tetrace[i][r_trace_num].h), r_event_num);

					// Read the trace --> tetrace[i][r_trace_num]
					for (j = 0; j < r_event_num; j++){
						tetrace[i][r_trace_num].v[j] = tv_tmp[j];
						tetrace[i][r_trace_num].h[j] = th_tmp[j];
					}
					tetrace[i][r_trace_num].event_num = r_event_num;

					// Increase the number of traces by one
					r_trace_num++;
				}
				// Initialize the number of events in the trace
				r_event_num = 0;
				continue;
			}

			// Latitude, latitude, time --> tv, th, tim
			sscanf_s(s, "%lf %lf %d", &tv, &th, &tim);

			// Do not read a location out of range
			if ((tv < ini.v_min || tv >= ini.v_max) || (th < ini.h_min || th >= ini.h_max)){
				// If the number of events in the trace is more than one, increase the number of traces by one.
				if (r_event_num >= 2){
					// Malloc (tetrace)
					malloc1D(&(tetrace[i][r_trace_num].v_no), r_event_num);
					malloc1D(&(tetrace[i][r_trace_num].h_no), r_event_num);
					malloc1D(&(tetrace[i][r_trace_num].r_no), r_event_num);
					malloc1D(&(tetrace[i][r_trace_num].v), r_event_num);
					malloc1D(&(tetrace[i][r_trace_num].h), r_event_num);

					// Read the trace --> tetrace[i][r_trace_num]
					for (j = 0; j < r_event_num; j++){
						tetrace[i][r_trace_num].v[j] = tv_tmp[j];
						tetrace[i][r_trace_num].h[j] = th_tmp[j];
					}
					tetrace[i][r_trace_num].event_num = r_event_num;

					// Increase the number of traces by one
					r_trace_num++;
				}
				// Initialize the number of events in the trace
				r_event_num = 0;
				continue;
			}

			// Read the events --> tv_tmp, th_tmp
			tv_tmp[r_event_num] = tv;
			th_tmp[r_event_num] = th;

			// Update the number of events --> r_event_num
			r_event_num++;
		}

		// If the number of events in the trace is more than one, increase the number of traces by one.
		if (r_event_num >= 2){
			// Malloc (tetrace)
			malloc1D(&(tetrace[i][r_trace_num].v_no), r_event_num);
			malloc1D(&(tetrace[i][r_trace_num].h_no), r_event_num);
			malloc1D(&(tetrace[i][r_trace_num].r_no), r_event_num);
			malloc1D(&(tetrace[i][r_trace_num].v), r_event_num);
			malloc1D(&(tetrace[i][r_trace_num].h), r_event_num);

			// Read the trace --> tetrace[i][r_trace_num]
			for (j = 0; j < r_event_num; j++){
				tetrace[i][r_trace_num].v[j] = tv_tmp[j];
				tetrace[i][r_trace_num].h[j] = th_tmp[j];
			}
			tetrace[i][r_trace_num].event_num = r_event_num;

			// Increase the number of traces by one
			r_trace_num++;
		}

		fclose(fp);

		i++;
		if (i == ini.tuser_num) break;
	} while (FindNextFile(hFind, &fd) && i <= MAX_USER_NUM);
	FindClose(hFind);
	fclose(fp2);
	user_num = i;

	/********************************************** Output Statistical Information ***********************************************/
	// Total number of user i's events  --> tot_event_num[i]
	for (i = 0; i < user_num; i++){
		tot_event_num[i] = 0;
		for (j = 0; j < trace_num[i]; j++) tot_event_num[i] += tetrace[i][j].event_num;
	}

	// Number of user i's training events --> tevent_num[i]
	for (i = 0; i < user_num; i++) tevent_num[i] = (int)round(tot_event_num[i] * ini.tevent_rate);

	// Number of user i's training events --> ttrace_num[i]
	for (i = 0; i < user_num; i++){
		event_num_tmp = 0;
		for (j = 0; j < trace_num[i]; j++){
			event_num_tmp += tetrace[i][j].event_num;
			if (event_num_tmp >= tevent_num[i]) break;
		}
		ttrace_num[i] = j + 1;
	}

	// Number of user i's testing traces --> etrace_num[i]
	for (i = 0; i < user_num; i++) etrace_num[i] = trace_num[i] - ttrace_num[i];

	// Total number of testing traces --> etrace_sum
	etrace_sum = 0;
	for (u = 0; u < user_num; u++) etrace_sum += etrace_num[u];

	// Output the number of traces and events
	FileOpen(&fp, "trace_num.csv", "w");
	fileout(fp, "User ID,#traces,#traces(train),#traces(test),#events,#events(train)\n");
	for (i = 0; i < user_num; i++){
		fileout(fp, "%d,%d,%d,%d,%d,%d\n", i, trace_num[i], ttrace_num[i], etrace_num[i], tot_event_num[i], tevent_num[i]);
	}
	fclose(fp);

	// Histogram (the number of traces) --> trace_hist
	for (i = 0; i < MAX_EVENT_NUM; i++) trace_hist[i] = 0;
	for (i = 0; i < user_num; i++){
		if (trace_num[i] < MAX_EVENT_NUM) trace_hist[trace_num[i]]++;
		else trace_hist[MAX_EVENT_NUM - 1]++;
	}

	// Output the histogram (the number of traces)
	FileOpen(&fp, "trace_hist.csv", "w");
	fileout(fp, "#traces,count\n");
	for (i = 0; i < MAX_EVENT_NUM; i++) fileout(fp, "%d,%d\n", i, trace_hist[i]);
	fclose(fp);

	// Histogram (the number of events) --> event_hist
	for (i = 0; i < MAX_EVENT_NUM; i++) event_hist[i] = 0;
	for (i = 0; i < user_num; i++){
		for (j = 0; j < trace_num[i]; j++){
			event_hist[tetrace[i][j].event_num]++;
		}
	}

	// Output the histogram (the number of events)
	FileOpen(&fp, "event_hist.csv", "w");
	fileout(fp, "#events,count\n");
	for (i = 0; i < MAX_EVENT_NUM; i++) fileout(fp, "%d,%d\n", i, event_hist[i]);
	fclose(fp);

	/************************************ Divide Users into Target Users and Non-target Users ************************************/
	// Malloc (user_rnd)
	malloc2D(&user_rnd, ini.euser_rand_num, user_num);

	// Number of target users, non-target users --> ini.euser_num, nuser_num
	nuser_num = user_num - ini.euser_num;
	if (ini.euser_num <= 0){
		printf("The number of target users (=%d) is less than 0.\n", ini.euser_num);
		Exit(-1);
	}
	if (nuser_num < 0){
		printf("The number of non-target users (=%d) is less than 0.\n", nuser_num);
		Exit(-1);
	}

	printf("#target-users:%d, #non-target users:%d\n", ini.euser_num, nuser_num);

	// Randomly shuffle users (we attempt "ini.euser_rand_num" ways to randomly shuffle users) --> user_rnd
	for (i = 0; i < ini.euser_rand_num; i++) RndPerm(user_rnd[i], user_num);

	/********************************************************** Malloc ***********************************************************/
	ttrans = (struct trans **)malloc(sizeof(struct trans *) * user_num);
	for (u = 0; u < user_num; u++) ttrans[u] = (struct trans *)malloc(sizeof(struct trans) * tot_event_num[u]);
	malloc1D(&ttrans_num, user_num);
	malloc3D(&ttensor, ini.euser_num, r_num, r_num);
	malloc2D(&tmatrix_all, r_num, r_num);
	malloc1D(&rg, r_num);
	for (i = 0; i < r_num; i++) malloc1D(&(rg[i].r_no), r_num);

	malloc2D(&v_u_f, ini.euser_num, ini.k_u_f);
	malloc2D(&v_f_u, r_num, ini.k_u_f);
	malloc2D(&v_f_t, r_num, ini.k_f_t);
	malloc2D(&v_t_f, r_num, ini.k_f_t);
	malloc2D(&v_u_t, ini.euser_num, ini.k_u_t);
	malloc2D(&v_t_u, r_num, ini.k_u_t);

	malloc2D(&prior, ini.euser_num, r_num);
	malloc2D(&prior_n, ini.euser_num, r_num);
	malloc3D(&ttensor_n, ini.euser_num, r_num, r_num);

	malloc1D(&e_loglikeli_ratio, ini.euser_rand_num * etrace_sum * (1 + user_num - ini.euser_num));
	malloc1D(&n_loglikeli_ratio, ini.euser_rand_num * etrace_sum * (1 + user_num - ini.euser_num));
	malloc1D(&e_loglikeli_ratio_ast, ini.euser_rand_num * etrace_sum * (1 + user_num - ini.euser_num));
	malloc1D(&n_loglikeli_ratio_ast, ini.euser_rand_num * etrace_sum * (1 + user_num - ini.euser_num));
	malloc2D(&e_loglikeli_ratio_tr, 6, ini.euser_rand_num * etrace_sum * (1 + user_num - ini.euser_num));
	malloc2D(&n_loglikeli_ratio_tr, 6, ini.euser_rand_num * etrace_sum * (1 + user_num - ini.euser_num));

	/*********************************************** Training & Re-identification ************************************************/
	// Decide region boundaries and assign region IDs to events in traces --> tetrace
	AssignRegionID(tetrace, user_num, trace_num, ini.v_num, ini.h_num, ini.v_min, ini.v_max, ini.h_min, ini.h_max);

	// Compute a transition pattern structure from training traces--> ttrans, ttrans_num
	MakeTrans(tetrace, user_num, ttrace_num, tevent_num, ttrans, ttrans_num);

	// Conduct experiments for each way to randomly divide users into target users and non-target users
	for (en = 0; en < ini.euser_rand_num; en++){
		printf("Training & Re-identification #%d:\n", en + 1);
		// Compute a transition count tensor from the transition pattern structure --> ttensor
		// Compute a transition matrix from all users --> tmatrix_all
		Trans2Tensor(ttrans, ttrans_num, user_rnd[en], ini.euser_num, r_num, ttensor, tmatrix_all);

		// Make an MCL input file and MCL command file
		if (ini.MCL_make_flag == 1)
		{
			double ijn, jin;
			int ix, iy, jx, jy, l1;

			// If there is no directory, make it.
			sprintf_s(s, "%s", ini.mcl_dir);
			if (!PathFileExists(s)) _mkdir(s);

			// Make an MCL input file
			sprintf_s(s, "%s\\mcl_en%d.csv", s, en + 1);

			FileOpen(&fp, s, "w");
			for (i = 0; i < r_num; i++){
				for (j = i + 1; j < r_num; j++){
					ix = i %ini.h_num;
					iy = i / ini.h_num;
					jx = j %ini.h_num;
					jy = j / ini.h_num;
					l1 = abs(ix - jx) + abs(iy - jy);

					if (tmatrix_all[i][j] != -1) ijn = tmatrix_all[i][j];
					else ijn = 0.0;
					if (tmatrix_all[j][i] != -1) jin = tmatrix_all[j][i];
					else jin = 0.0;
					if (l1 <= 1) ijn += 0.01;
					fileout(fp, "%d %d %f\n", i, j, ijn + jin);
				}
			}
			fclose(fp);

			// Make an MCL command file
			sprintf_s(s, "%s\\command.txt", ini.mcl_dir, ini.euser_num);
			if(en == 0) FileOpen(&fp, s, "w");
			else FileOpen(&fp, s, "a");
			fileout(fp, "mcl mcl_en%d.csv -I 6 --abc -o group_en%d.txt\n", en + 1, en + 1);
			fclose(fp);
		}

		// Read region groups
		if (ini.train_method == 3){
			sprintf_s(s, "%s\\group_en%d.txt", ini.mcl_dir, en + 1);

			FileOpen(&fp, s, "r");

			rg_num = 0;
			while (fgets(s, 1023, fp) != NULL){
				s[strlen(s) - 1] = '\0';
				i = 0;
				if ((tok = strtok_s(s, "\t", &ctx)) != NULL){
					rg[rg_num].r_no[i++] = atoi(tok);
					while ((tok = strtok_s(NULL, "\t", &ctx)) != NULL) rg[rg_num].r_no[i++] = atoi(tok);
					rg[rg_num].size = i;
					rg_num++;
				}
			}
		}

		// +++++++++++++++++++++++++++++++++++++ Training transition matrices +++++++++++++++++++++++++++++++++++++ //
		switch (ini.train_method){
		// ML (Maximum Likelihood)
		case 1:
			// ML
			ML(ttensor, ini.euser_num, r_num);
			break;
		// TF (Tensor Factorization)
		case 2:
			// TF
			TF(v_u_f, v_f_u, v_f_t, v_t_f, v_u_t, v_t_u, r_num, ini.euser_num, ini.k_u_f, ini.k_f_t, ini.k_u_t, ttensor, ini.alpha);
			// Compute a transition probability tensor (i.e. personalized transition matrices) --> ttensor
			TFParam2Tensor(v_u_f, v_f_u, v_f_t, v_t_f, v_u_t, v_t_u, ini.k_u_f, ini.k_f_t, ini.k_u_t, ini.euser_num, r_num, ttensor);
			break;
		// GSTF(Group Sparsity Tensor Factorization)
		case 3:
			// GSTF
			GSTF(v_u_f, v_f_u, v_f_t, v_t_f, v_u_t, v_t_u, r_num, ini.euser_num, ini.k_u_f, ini.k_f_t, ini.k_u_t, ttensor, ini.alpha, ini.beta, rg, rg_num);
			// Compute a transition probability tensor (i.e. personalized transition matrices) --> ttensor
			TFParam2Tensor(v_u_f, v_f_u, v_f_t, v_t_f, v_u_t, v_t_u, ini.k_u_f, ini.k_f_t, ini.k_u_t, ini.euser_num, r_num, ttensor);

			break;
		}

		// Compute a personalized prior distribution (stationary distribution) --> prior
		Tensor2Prior(ttensor, ini.euser_num, r_num, prior);

		// Compute a population transition matrix for each target user --> ttensor_n
		CalcPopMat(ttensor, ttensor_n, ini.euser_num, r_num);

		// Compute a population prior distribution (stationary distribution) --> prior_n
		Tensor2Prior(ttensor_n, ini.euser_num, r_num, prior_n);

		// Assign a very small positive value (PROB_MIN) to an element whose value is 0 --> ttensor, prior
		for (u = 0; u < ini.euser_num; u++){
			for (i = 0; i < r_num; i++){
				for (j = 0; j < r_num; j++){
					if (ttensor[u][i][j] < PROB_MIN) ttensor[u][i][j] = PROB_MIN;
				}
				if (prior[u][i] < PROB_MIN) prior[u][i] = PROB_MIN;
			}
		}

		// Assign a very small positive value (PROB_MIN) to an element whose value is 0 --> ttensor_n, prior_n
		for (u = 0; u < ini.euser_num; u++){
			for (i = 0; i < r_num; i++){
				for (j = 0; j < r_num; j++){
					if (ttensor_n[u][i][j] < PROB_MIN) ttensor_n[u][i][j] = PROB_MIN;
				}
				if (prior_n[u][i] < PROB_MIN) prior_n[u][i] = PROB_MIN;
			}
		}

		// +++++++++++++++++++++++ Re-identification (Computation of Log-likelihood Ratios) +++++++++++++++++++++++ //
		// Malloc (loglikeli, loglikeli_n)
		loglikeli = (double ***)malloc(sizeof(double **)* user_num);
		loglikeli_n = (double ***)malloc(sizeof(double **)* user_num);
		for (u = 0; u < user_num; u++){
			etn = etrace_num[user_rnd[en][u]];
			malloc2D(&(loglikeli[u]), etn, ini.euser_num);
			malloc2D(&(loglikeli_n[u]), etn, ini.euser_num);
		}

		// Compute a log-likelihood using the personalized/population transition matrix --> loglikeli, loglikeli_n
		CalcLogLikeli(ttensor, ttensor_n, prior, prior_n, user_num, ini.euser_num, user_rnd[en], r_num,ini.v_num,ini.h_num, tetrace, ttrace_num, etrace_num, loglikeli, loglikeli_n);

		// Decide whether an anoymized trace is generated from target user e (two-class classification)
		// e: target user
		for (int e = 0; e < ini.euser_num; e++){
			// u: testing user
			for (u = 0; u < user_num; u++){
				// Number of user u's traces --> etn
				etn = etrace_num[user_rnd[en][u]];

				// Continue if user u is a target user other than user e
				if (u < ini.euser_num && u != e) continue;

				// Compute log-likelihood ratios
				// i: trace of user u
				for (i = 0; i < etn; i++){
					// Classify the trace based on the trace length (the number of events) --> trace_type
					for (j = 0; j < 5; j++){
						if (tetrace[user_rnd[en][u]][i].event_num <= trace_interval[j]) break;
					}
					trace_type = j;

					// Log-likelihood ratio for target user e --> e_loglikeli_ratio
					if (u == e){
						e_loglikeli_ratio[e_llr_num] = (loglikeli[u][i][e] - loglikeli_n[u][i][e]);
						e_loglikeli_ratio_tr[trace_type][e_llr_num_tr[trace_type]] = (loglikeli[u][i][e] - loglikeli_n[u][i][e]);

						e_llr_num++;
						e_llr_num_tr[trace_type]++;
					}
					// Log-likelihood ratio for a non-target user --> n_loglikeli_ratio
					else{
						n_loglikeli_ratio[n_llr_num] = (loglikeli[u][i][e] - loglikeli_n[u][i][e]);
						n_loglikeli_ratio_tr[trace_type][n_llr_num_tr[trace_type]] = (loglikeli[u][i][e] - loglikeli_n[u][i][e]);

						n_llr_num++;
						n_llr_num_tr[trace_type]++;
					}
				}

				// Compute log-likelihood ratios in the case when the testing traces are combined into one
				// Log-likelihood ratio for target user e --> e_loglikeli_ratio_ast
				if (u == e){
					llr_sum = 0.0;
					for (j = 0; j < etn; j++) llr_sum += (loglikeli[u][j][e] - loglikeli_n[u][j][e]);
					e_loglikeli_ratio_ast[e_llr_ast_num] = llr_sum;

					e_llr_ast_num++;
				}
				// Log-likelihood ratio for a non-target user --> n_loglikeli_ratio_ast
				else{
					llr_sum = 0.0;
					for (j = 0; j < etn; j++) llr_sum += (loglikeli[u][j][e] - loglikeli_n[u][j][e]);
					n_loglikeli_ratio_ast[n_llr_ast_num] = llr_sum;

					n_llr_ast_num++;
				}
			}
		}

		// Free (loglikeli, loglikeli_n)
		for (u = 0; u < user_num; u++){
			etn = etrace_num[user_rnd[en][u]];
			free2D(loglikeli[u], etn);
		}
		free1D(loglikeli);
		for (u = 0; u < user_num; u++){
			etn = etrace_num[user_rnd[en][u]];
			free2D(loglikeli_n[u], etn);
		}
		free1D(loglikeli_n);

		printf("done.\n");
	}

	/******************************************************** Evaluation *********************************************************/
	// Compute FPR, TPR, and AUC
	// Output-file name --> outfile
	sprintf_s(outfile, "result_tm%d.csv", ini.train_method);

	// Compute FPR, TPR, and AUC, and output them to the file
	CalcFPRTPR(outfile, e_loglikeli_ratio, e_llr_num, n_loglikeli_ratio, n_llr_num, 500.0, -1500.0, 1.0);

	// Compute FPR, TPR, and AUC in the case when the testing traces are combined into one
	// Output-file name --> outfile
	sprintf_s(outfile, "result_tm%d_ast.csv", ini.train_method);

	// Compute FPR, TPR, and AUC, and output them to the file
	CalcFPRTPR(outfile, e_loglikeli_ratio_ast, e_llr_ast_num, n_loglikeli_ratio_ast, n_llr_ast_num, 500.0, -1500.0, 1.0);

	// Compute FPR, TPR, and AUC for each trace type
	for (i = 0; i < 6; i++){
		// Output-file name --> outfile
		sprintf_s(outfile, "result_tm%d_tt%d.csv", ini.train_method, i + 1);

		// Compute FPR, TPR, and AUC, and output them to the file
		CalcFPRTPR(outfile, e_loglikeli_ratio_tr[i], e_llr_num_tr[i], n_loglikeli_ratio_tr[i], n_llr_num_tr[i], 500.0, -1500.0, 1.0);
	}

	// Free
	for (i = 0; i < user_num; i++){
		for (j = 0; j < trace_num[i]; j++){
			free1D(tetrace[i][j].v_no);
			free1D(tetrace[i][j].h_no);
			free1D(tetrace[i][j].r_no);
			free1D(tetrace[i][j].v);
			free1D(tetrace[i][j].h);
		}
		free1D(tetrace[i]);
	}
	free1D(tetrace);

	free2D(user_rnd, ini.euser_rand_num);

	for (u = 0; u < user_num; u++) free1D(ttrans[u]);
	free1D(ttrans);
	free1D(ttrans_num);
	free3D(ttensor, ini.euser_num, r_num);
	free2D(tmatrix_all, r_num);
	for (i = 0; i < r_num; i++) free1D(rg[i].r_no);
	free1D(rg);

	free2D(v_u_f, ini.euser_num);
	free2D(v_f_u, r_num);
	free2D(v_f_t, r_num);
	free2D(v_t_f, r_num);
	free2D(v_u_t, ini.euser_num);
	free2D(v_t_u, r_num);
	free2D(prior, ini.euser_num);
	free3D(ttensor_n, ini.euser_num, r_num);
	free2D(prior_n, ini.euser_num);

	free1D(e_loglikeli_ratio);
	free1D(n_loglikeli_ratio);
	free1D(e_loglikeli_ratio_ast);
	free1D(n_loglikeli_ratio_ast);

	free1D(e_loglikeli_ratio_tr);
	free1D(n_loglikeli_ratio_tr);

	printf("Proceed Success.");
	Exit(0);
}
