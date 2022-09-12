# GroupSparsityTensorFactorization

This is a source code of re-identification attacks in an open scenario based on group sparsity tensor factorization (GSTF), which is described in the following article: 

[Murakami+, TIFS17] Takao Murakami, Atsunori Kanemura, Hideitsu Hino, "Group Sparsity Tensor Factorization for Re-identification of Open Mobility Traces," 
IEEE Transactions on Information Forensics and Security, Vol.12, No.3, pp.689-704, 2017.

This source code is written in c/c++. It is made public to show how the proposed training method in [Murakami+, TIFS17] (i.e., GSTF) is implemented.
The source code can also be used to conduct experiments in a biased setting using the Gowalla dataset (which is described in Section VII of [Murakami+, TIFS17]).
For the contents of the source code, please see comments in the code.
A supplemental material (supplemental_material.pdf) is also provided to describe how to solve the optimization problems (19), (20), (22), and (24) in GSTF in details.

Please cite [Murakami+, TIFS17] when you publish a paper using this code.

# Execution Environment
Windows 10, Visual Studio 2013

# How to Compile the Source Codes

(Steps 1 and 2 are required to build the Reidentification_GSTF solution file.)

1. Download the Mersenne Twister source code (mt19937ar.sep.tgz) from the following website:
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html

2. Move mt19937ar.c and mt19937ar.h as follows:
   - "mt19937ar.c" --> "source\Reidentification_GSTF\cpp\mt19937ar.c"
   - "mt19937ar.h" --> "source\Reidentification_GSTF\header\mt19937ar.h"

3. Build the solution files.

4. Copy "ExtrLocData_Gow.exe" and "Reidentification_GSTF.exe" to "release" directory.

# How to Execute the Programs

(Steps 6 and 7 can be skipped.)

1. Download the Gowalla dataset (loc-gowalla_totalCheckins.txt.gz): https://snap.stanford.edu/data/loc-gowalla.html

2. Unzip the downloaded dataset and move "Gowalla_totalCheckins.txt" to "dataset" directory.

3. Go to "release" directory and run "ExtrLocData_Gow.exe".

   Location traces are extracted, and are output in "dataset\Data_30min" directory. 
   Each trace is separated by a line break, and the 1st, 2nd, and 3rd columns are [latitude], [longitude], and [UNIX time], respectively.

4. Copy "Reidentification_GSTF_tm1.ini" to "Reidentification_GSTF.ini", and run "Reidentification_GSTF.exe".

   Experimental results of ML (and ML*) are output in "release" directory as follows: 
   - event_hist.csv: #events (locations) per trace v.s. frequency
   - result_tm1.csv: FPR, TPR, and AUC of ML
   - result_tm1_ast.csv: FPR, TPR, and AUC of ML*
   - result_tm1_tt[1|2|3|4|5|6].csv: FPR, TPR, and AUC of ML for the testing trace of the size [2-5|6-10|11-15|16-20|21-25|26-].
   - trace_hist.csv: #traces per user v.s. frequency
   - trace_num.csv: #traces and #events for each user

5. Copy "Reidentification_GSTF_tm2.ini" to "Reidentification_GSTF.ini", and run "Reidentification_GSTF.exe".

   Experimental results of TF (and TF*) are output in "release" directory.

6. Download the MCL tool (mcl-14-137.tar) from http://www.micans.org/mcl/ and install the MCL tool. 
   For example, if you have installed Cygwin (with automake-1.15, gcc, and make packages), move to "mcl-14-137" directory and run the following commands:
   $ cp /usr/share/automake-1.15/config.guess autofoo/config.guess
   $ ./configure
   $ make
   $ make install

7. Go to "release\mcl" directory, and run the MCL commands written in "command.txt".
   Region group files (group_en[1-10].txt) are output in "release\mcl" directory. Each line represents a region group with region IDs separated by TAB.

8. Copy "Reidentification_GSTF_tm3.ini" to "Reidentification_GSTF.ini", and run "Reidentification_GSTF.exe".

   Experimental results of GSTF (and GSTF*) are output in "release" directory.
   The experimental results in Section VII of [Murakami+, TIFS17] (see "result" directory) can be obtained from the files output in Steps 4, 5, and 8.

This software is released under the MIT License, see LICENSE.txt.
