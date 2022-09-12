/**
** @file FileIO.cpp
** 
** @brief
** File I/O functions.
** 
** Copyright (c) 2017 National Institute of Advanced Industrial Science and Technology (AIST)
** 
** This software is released under the MIT License.
** http://opensource.org/licenses/mit-license.php
*/

#include <stdio.h>
#include <shlwapi.h>
#include <Windows.h>
#include "MemoryOperation.h"
#include "parameter.h"
#include "Prototype.h"
#include "Struct.h"

/**
** @brief 			Read parameters in INI file.
** @return			parameters in INI file
*/
iniParameter readIniParameter()
{
	char s[1024];
	iniParameter ret;
	
	if (!PathFileExists(INI_PATH))
	{
		GetFullPathName(INI_PATH, sizeof(s), s, NULL);
		printf("readIniParameter: ini file not exists. %s", s);
		Exit(-1);
	}

	// Read INI file
	GetPrivateProfileString("PATH", "DATA_EXTR", "default", s, sizeof(s), INI_PATH);
	GetFullPathName(s, sizeof(ret.data_extr), ret.data_extr, NULL);
	GetPrivateProfileString("PATH", "MCL_DIR", "default", s, sizeof(s), INI_PATH);
	GetFullPathName(s, sizeof(ret.mcl_dir), ret.mcl_dir, NULL);

	GetPrivateProfileString("PARAMETER", "TEVENT_RATE", "default", ret.tevent_rate_s, sizeof(ret.tevent_rate_s), INI_PATH);
	ret.tevent_rate = atof(ret.tevent_rate_s);

	ret.euser_num = GetPrivateProfileInt("PARAMETER", "EUSER_NUM", 0, INI_PATH);
	ret.tuser_num = GetPrivateProfileInt("PARAMETER", "TUSER_NUM", 0, INI_PATH);
	ret.euser_rand_num = GetPrivateProfileInt("PARAMETER", "EUSER_RAND_NUM", 0, INI_PATH);

	ret.v_num = GetPrivateProfileInt("PARAMETER", "VREGION_NUM", 0, INI_PATH);
	ret.h_num = GetPrivateProfileInt("PARAMETER", "HREGION_NUM", 0, INI_PATH);

	GetPrivateProfileString("PARAMETER", "V_MIN", "default", s, sizeof(s), INI_PATH);
	ret.v_min = atof(s);
	GetPrivateProfileString("PARAMETER", "V_MAX", "default", s, sizeof(s), INI_PATH);
	ret.v_max = atof(s);
	GetPrivateProfileString("PARAMETER", "H_MIN", "default", s, sizeof(s), INI_PATH);
	ret.h_min = atof(s);
	GetPrivateProfileString("PARAMETER", "H_MAX", "default", s, sizeof(s), INI_PATH);
	ret.h_max = atof(s);

	ret.train_method = GetPrivateProfileInt("PARAMETER", "TRAIN_METHOD", 1, INI_PATH);

	ret.k_u_f = GetPrivateProfileInt("PARAMETER", "K_U_F", 0, INI_PATH);
	ret.k_f_t = GetPrivateProfileInt("PARAMETER", "K_F_T", 0, INI_PATH);
	ret.k_u_t = GetPrivateProfileInt("PARAMETER", "K_U_T", 0, INI_PATH);

	GetPrivateProfileString("PARAMETER", "ALPHA", "default", ret.alpha_s, sizeof(ret.alpha_s), INI_PATH);
	ret.alpha = atof(ret.alpha_s);
	GetPrivateProfileString("PARAMETER", "BETA", "default", ret.beta_s, sizeof(ret.beta_s), INI_PATH);
	ret.beta = atof(ret.beta_s);

	ret.MCL_make_flag = GetPrivateProfileInt("PARAMETER", "MCL_MAKE", 0, INI_PATH);

	if (ret.train_method != 1)
		ret.MCL_make_flag = 0;
	return ret;
}

/**
** @brief 			Open a file.
** @param fp	 	File pointer
** @param filename 	File name
** @param mode	 	Mode
*/
void FileOpen(FILE **fp, char *filename, char *mode){
	errno_t error;

	if (error = fopen_s(fp, filename, mode) != 0){
		printf("FileOpen: Canot open %s\n", filename);
		Exit(-1);
	}
}

/**
** @brief 			Output formatted data to a file.
** @param fp	 	File pointer
** @param Format	C string
*/
void fileout(FILE *fp, const char * Format, ...)
{
	va_list st;
	va_start(st, Format);
	if (fp != NULL) vfprintf(fp, Format, st);
}
