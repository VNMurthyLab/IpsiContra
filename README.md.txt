#Bilateral Alignment of Receptive Fields in the Olfactory Cortex Points to Non-Random Connectivity

## OVERVIEW

This folder contains all the Matlab data and code necessary to replicate the results of the following article:

* **Authors**: Grimaud J, Dorrell W, Pehlevan C, Murthy V.
* **Title**: Bilateral Alignment of Receptive Fields in the Olfactory Cortex Points to Non-Random Connectivity.
* **Published in** bioRxiv on [INSERT DATE].

Details on how each variable was calculated can be found in the Methods section of the article.
If you use this code or data please cite the article mentioned above.

Below is a description of all the Matlab code and data files.
Matlab version: R2016a.

## 1- CODE

### Figures 2, 3, 4, S1, S3, and S4

There is one Matlab script per figure.
The name of each Matlab code is self-explanatory.
For example: code_for_figure_2.mat re-creates the plots from Figure 2.
To run this code, keep the code files in the same folder as the data files.

### Figures 5, 6, S5, and S7

The code corresponding to the decoding and modeling figures can be found in the folder called "decoding_and_modeling".
To run this code, keep all the files from "decoding_and_modeling" in the same folder. Within this one folder, create a subfolder called "Data_and_Notes", in which you will need to add the following data files:
* SniffTime
* SniffTimeOB
* tetrodeRecordings_OC_2s
* tetrode_recordings_formatted_for_alex_v2_ordered
* tetrode_recordings_formatted_for_alex_v2_ob_ordered

For details about the data files, see below.
For more details and the decoding and modeling code, see the Readme.txt file inside the "decoding_and_modeling" folder.

### Other Figures

Figures 1 and S2 do not require any script.
Figure S6 is made of various plots from other figures.

## 2- DATA FILES RELATED TO OC TETRODE RECORDINGS

### tetrodeRecordings_OC_2s

This file contains 4 variables: A, B, BasalFR, and C.

A is a 1x10 cell.
For mouse m, M = A{m}, M is a 4-D matrix.
M(n,o,s,t) is the response of neuron n to odor o presented on the side s of the mask for the t-th time.
The response was calculated as described in the Methods section.
side 1 is ipsi, side 2 is contra.

B is a 1x10 cell.
For mouse m, M = B{m}, M is a 1-D matrix.
M(n) is the recording session identity of neuron n.
For example, all neurons with M(n)==2 were recorded during the 2nd recording session.

BasalFR is a 1x10 cell.
For mouse m, M = BasalFR{m}, M is a 1-D matrix.
M(n) is the basal firing rate of neuron n (in Hz).

C is a 1x10 cell.
For mouse m, M = C{m}, M is a 3-D matrix.
M(n,o,s) is equal to 1 if neuron n responded significantly to odor o presented on the side s of the mask, 0 otherwise.
Significance was calculated as described in the Methods section.
side 1 is ipsi, side 2 is contra.

### tetrodeRecordings_OC_2s_notBaselineSubstracted

This file is similar to tetrodeRecordings_OC_2s.
The only difference is:
For A, baseline activities were not substracted from odor-induced activities.

### tetrodeRecordings_OC_waveform

This file contains 2 variables: SpikeTrace and SpikeWidth

SpikeTrace is a 1x10 cell.
For mouse m, M = SpikeTrace{m}, M is a 2-D matrix.
M(n,:) is the average spike trace of neuron n (acquisition rate: 30kHz).

SpikeWidth is a 1x10 cell.
For mouse m, M = SpikeWidth{m}, M is a 1-D matrix.
M(n) is the width at half max amplitude of the average spike trace of neuron n (in microseconds).

### tetrode_recordings_formatted_for_alex_v2_ordered_aon

This file contains 3 variables. All variables are as follows: Array[XXXX].
[XXXX] is the name of a mouse recorded in the AON.

For mouse [XXXX], neuron n, odor o, side s, M = Array[XXXX]{n,o,s}, M is a 2-D matrix.
M is the PSTH of neuron n for odor o presented on side s (7 repeats).
The bin size of the PSTH is 10ms. Each bin contains a spike count.
For each repeat, the first 599 bins correspond to 6s preceeding odor onset.
The odor valve opens at time bin #600 and closes at bin #799 (total: 2s of odor delivery).
The rest of the PSTH (bins #800 to #1100) corresponds to the activity during the 4s following the end of odor onset.
side 3 is ipsi, side 2 is contra. Note: s==1 is unused in this article.

### tetrode_recordings_formatted_for_alex_v2_ordered_apc

This is similar to tetrode_recordings_formatted_for_alex_v2_ordered_aon, for mice recorded in the APC

### tetrode_recordings_formatted_for_alex_v2_ordered_ppc

This is similar to tetrode_recordings_formatted_for_alex_v2_ordered_aon, for mice recorded in the PPC

### SniffTime

This file contains one variable: SniffTime.

SniffTime is a 4-D cell.
For mouse m, neuron n, odor o, side s, M = SniffTime{m,n,o,s}, M is a 1-D matrix.
For trial t, M(t) is the time from odor onset to the first sniff after odor onset (in seconds).
If no sniff was detected between odor onset and 0.5s after odor onset, then M(t) = 0.
side 3 is ipsi, side 2 is contra. Note: s==1 is unused in this article.

## 3- DATA FILES RELATED TO OB TETRODE RECORDINGS

### tetrodeRecordings_OB_new

This file contains 2 variables: A and C.

A is a 1x1 cell.
M = A{1}, M is a 3-D matrix.
M(n,o,s) is the average response of neuron n to odor o presented on the side s.
side 1 is ipsi, side 2 is contra.

C is a 1x1 cell.
M = C{1}, M is a 3-D matrix.
M(n,o,s) is equal to 1 if neuron n responded significantly to odor o presented on the side s of the mask, 0 otherwise.
side 1 is ipsi, side 2 is contra.

### tetrode_recordings_formatted_for_alex_v2_ob_ordered

This file is similar to tetrode_recordings_formatted_for_alex_v2_ordered, except it contains OB instead of OC recordings.

### SniffTimeOB

This file is similar to SniffTime, except (1) it reports sniff times for OB instead of OC recordings, and (2) the neurons from all mice have been pooled together (as only a handful of neuron was recorded in each mouse).

## 4- DATA FILES RELATED THE PID RECORDINGS

### pid_test_facialMask

This file contains 3 variables: odorlist, PIDtracesL, and PIDtracesR.

odorlist is a 2-D matrix.
The first dimension is the number of trials on the PID (i.e. the number of odor pulses delivered to the PID nozzle).
For trial t, M = odorlist(t,:), M contains two elements.
M(1) is the identity of the odor delivered during trial t.
M(2) is the side of the mask through which the odor was delivered.
side 1 is right, side 2 is left.

PIDtracesL is a 2-D matrix.
For trial t, M = PIDtracesL(:,t) is the PID trace recorded with the PID nozzle placed in the left part of the mask .
M shows 4 seconds of recording - 1s before odor delivery, 2s of odor delivery, and 1s after odor delivery (acquisition rate: 100Hz).

PIDtracesR is identical to PIDtracesL, except the PID nozzle was placed in the right part of the mask.

## 5- DATA FILES RELATED TO CALCIUM IMAGING

### gcamp3_mouse1_intensityTraces_glomerularMasks

This file contains 4 variables: Glomeruli, intensityTracesL, intensityTracesR, and mipbaseline.
They all correspond to the first of three mice imaged above their olfactory bulbs for glomerular calcium signals.

Glomeruli is a 3-D matrix.
For glomerulus g, Glomeruli(:,:,g) is the mask of g (a binary image where pixel value is 0 on top of the glomerulus, 1 otherwise).

intensityTracesL is a 3-D matrix.
For glomerulus g, odor o, M = intensityTracesL(o:16:16*7,:,g), M is a 2-D matrix.
M contains the raw calcium traces of g responding to o presented on the left side of the mask.
There are 7 repeats of each odor. Acquisition rate: 4Hz.
Total duration of each trace: 9.25 seconds, including 4s pre-odor delivery and 2s of odor delivery.

intensityTracesR is identical to intensityTracesL, except odors were delivered on the right side of the mask.

mipbaseline is a 2-D matrix.
It is the image of the craniotomy before any odor was applied.
Glomerulus masks from Glomeruli can be overimpose on mipbaseline.

### gcamp3_allmice_glomResponsesIpsiSignif

This file contains 12 variables. All variables are as follows: dFFOneSide[X]_m[Y]_b[Z]
X is empty for ipsilateral odor deliveries, or X==2 for contralateral deliveries.
Y is the mouse number (note: we imaged 3 mice).
Z is the olfactory bulb number (note: we image both bulbs in each mouse).
For example: dFFOneSide_m2_b1 corresponds to mouse #2, bulb #1, with odors delivered ipsilaterally.

dFFOneSide[X]_m[Y]_b[Z] is a 1-D matrix.
It contains the dF/F values of all glomerulus-odor pairs significantly responding to ipsilateral presentations.
Note that dFFOneSide_m[Y]_b[Z] and dFFOneSide2_m[Y]_b[Z] contain the same glomerulus-odor pairs.
See Methods section of the manuscript for more details.

## 6- DATA FILES RELATED TO BREATHING MONITORING TESTS

### airFlowSensor_versus_cannula_test

This file contains 6 variables: BR_exh_airflow, BR_exh_cannula, BR_inh_airflow, BR_inh_cannula, Shift_exh, and Shift_inh.

BR_exh_airflow is a 1-D matrix.
It contains all the instantaneous breathing rates (in Hz) calculated from the exhalation peaks detected with the airflow sensor (face mask).

BR_inh_airflow is identical to BR_exh_airflow, except we used the inhalation peaks.

BR_exh_cannula is a 1-D matrix.
It contains all the instantaneous breathing rates (in Hz) calculated from the exhalation peaks detected with the pressure sensor (chronic intranasal cannula).

BR_inh_cannula is identical to BR_exh_cannula, except we used the inhalation peaks.

Shift_exh is a 1-D matrix.
It contains, for all exhalation peaks detected with the cannula, the delay (in seconds) from each of these peaks to the equivalent peak on the airflow signal.

Shift_inh is identical to Shift_exh, except we used the inhalation peaks.