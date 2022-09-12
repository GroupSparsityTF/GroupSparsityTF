/**
** @file Source.cpp
** 
** @brief
** Extract location traces from the Gowalla dataset.
** 
** Copyright (c) 2017 National Institute of Advanced Industrial Science and Technology (AIST)
** 
** This software is released under the MIT License.
** http://opensource.org/licenses/mit-license.php
*/

#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <direct.h>
#include <shlwapi.h>
#include <time.h>

#define INI_PATH	".\\ExtrLocData_Gow.ini"
#define MAX_LINE	6442892

int Exit(int rc){
	printf("\nPress ENTER to exit.\n"); getchar();
	exit(rc);
}

int line_p[MAX_LINE];

int main(void){
	char infile[1024];
	char outdir[1024], allfile[1024], tracefile[1024];
	int time_int_min, time_int_max;

	FILE *fp, *fp2;
	errno_t error;

	double v_min, v_max;
	double h_min, h_max;
	double v, h;
	int user_id;
	int user_num;
	int pre_user_id;
	char s[1024];
	int i, j;
	int t1, t2, t_dif;
	struct tm tim;
	int line;

	// Read INI file
	GetPrivateProfileString("PATH", "GOWALLA_FILE", "default", infile, sizeof(infile), INI_PATH);
	GetPrivateProfileString("PATH", "GOWALLA_EXTR", "default", outdir, sizeof(outdir), INI_PATH);
	time_int_min = GetPrivateProfileInt("PARAMETER", "TIME_INTERVAL_MIN", 0, INI_PATH);
	time_int_max = GetPrivateProfileInt("PARAMETER", "TIME_INTERVAL_MAX", 0, INI_PATH);

	GetPrivateProfileString("PARAMETER", "V_MIN", "default", s, sizeof(s), INI_PATH);
	v_min = atof(s);
	GetPrivateProfileString("PARAMETER", "V_MAX", "default", s, sizeof(s), INI_PATH);
	v_max = atof(s);
	GetPrivateProfileString("PARAMETER", "H_MIN", "default", s, sizeof(s), INI_PATH);
	h_min = atof(s);
	GetPrivateProfileString("PARAMETER", "H_MAX", "default", s, sizeof(s), INI_PATH);
	h_max = atof(s);

	// If there is no directory, make it.
	if (!PathFileExists(outdir)) _mkdir(outdir);

	// All-trace file --> allfile
	sprintf_s(allfile, "%s\\\\all_trace.csv", outdir);

	// Open the original Gowalla dataset
	if (error = fopen_s(&fp, infile, "r") != 0){
		printf("cannot open %s\n", infile);
		Exit(-1);
	}

	// Open all-trace file
	if (error = fopen_s(&fp2, allfile, "w") != 0){
		printf("cannot open %s\n", allfile);
		Exit(-1);
	}
	// header
	fprintf_s(fp2, "user,id,latitude,longitude,time,interval\n");

	user_num = 0;
	pre_user_id = -1;
	t1 = 0;
	// Read each line
	while (fgets(s, 1023, fp) != NULL){
		sscanf_s(s, "%d	%04d-%02d-%02dT%02d:%02d:%02dZ	%lf %lf", &user_id, &(tim.tm_year), &(tim.tm_mon), &(tim.tm_mday), &(tim.tm_hour), &(tim.tm_min), &(tim.tm_sec), &v, &h);

		// Exclude locations out of range
		if ((v < v_min || v >= v_max) || (h < h_min || h >= h_max)){
			continue;
		}

		// UNIX time --> t2
		tim.tm_year -= 1900;
		tim.tm_mon -= 1;
		t2 = (int)mktime(&tim);

		// time interval --> t_dif
		if (pre_user_id == user_id) t_dif = t1 - t2;
		else t_dif = 0;

		// If time interval is less than time_int_min, continue
		if ((pre_user_id == user_id) && (t_dif < time_int_min)) continue;

		// If a new user appears, update the number of users --> user_num
		if (pre_user_id != user_id) user_num++;

		// Output the event to the all-trace file
		fprintf_s(fp2, "%d,%d,%lf,%lf,%d,%d\n", user_num - 1, user_id, v, h, t2, t_dif);

		pre_user_id = user_id;
		t1 = t2;
	}

	fclose(fp);
	fclose(fp2);

	printf("[Gowalla dataset]\n");
	printf("Number of users: %d\n", user_num);

	// Open the all-trace file
	if (error = fopen_s(&fp, allfile, "r") != 0){
		printf("Cannot open %s\n", allfile);
		Exit(-1);
	}

	// Skip the first line
	fgets(s, 1023, fp);

	// Read each line
	printf("Making each trace file (#: 100 users): \n");
	pre_user_id = -1;
	while (fgets(s, 1023, fp) != NULL){
		sscanf_s(s, "%d,%*d,%lf,%lf,%d,%d", &user_id, &v, &h, &t2, &t_dif);

		// If a new user appears
		if (pre_user_id != user_id){
			// Open a new trace file
			sprintf_s(tracefile, "%s\\\\new_%04d.txt", outdir, user_id);
			if (error = fopen_s(&fp2, tracefile, "w") != 0){
				printf("cannot open %s\n", tracefile);
				Exit(-1);
			}

			// Output the event to the trace file
			fprintf_s(fp2, "%lf %lf %d\n", v, h, t2);

			fclose(fp2);

			if (user_id % 100 == 99) printf("#");
		}
		else{
			// Open the trace file
			if (error = fopen_s(&fp2, tracefile, "a") != 0){
				printf("cannot open %s\n", tracefile);
				Exit(-1);
			}

			// If time interval exceeds time_int_max, we assume that the event is in another trace
			if (time_int_max != -1 && (t_dif > time_int_max)) fprintf_s(fp2, "\n");

			// Output the event to the trace file
			fprintf_s(fp2, "%lf %lf %d\n", v, h, t2);

			fclose(fp2);
		}

		pre_user_id = user_id;
	}
	printf(" done.\n");

	fclose(fp);


	// Inverse each line in the trace files
	printf("Reversing each trace file (#: 100 users): \n");
	for (i = 0; i < user_num; i++){
		sprintf_s(tracefile, "%s\\\\new_%04d.txt", outdir, i);
		if (error = fopen_s(&fp, tracefile, "r") != 0){
			printf("cannot open %s\n", tracefile);
			Exit(-1);
		}

		line = 0;
		line_p[0] = 0;
		while (fgets(s, 1023, fp) != NULL){
			line++;
			line_p[line] = line_p[line - 1] + strlen(s) + 1;
		}

		if (error = fopen_s(&fp2, "tmp.txt", "w") != 0){
			printf("cannot open tmp.txt\n");
			Exit(-1);
		}

		for (j = 0; j < line; j++){
			fseek(fp, line_p[line - j - 1], SEEK_SET);
			fgets(s, 1023, fp);
			fputs(s, fp2);
		}

		fclose(fp);
		fclose(fp2);

		if (remove(tracefile) != 0){
			printf("Error: failed to remove %s.\n", tracefile);
		}

		if (rename("tmp.txt", tracefile) != 0){
			printf("Error: failed to move tmp.txt to %s.\n", tracefile);
		}

		if (i % 100 == 99) printf("#");
	}
	printf(" done.\n");

	Exit(0);
}
