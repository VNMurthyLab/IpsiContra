Code for decoding and modelling in:

"Matched Bilateral Odor Responses in the Olfactory Cortex Point to Non-Random
Connectivity", Grimaud et al.

There are three files that you should run in order to generate the modelling 
and decoding analysis figures (i.e. Figures 5, 6, S5 and S7)

1) DecodingAnalysis.m analyses the data and producing the decoding curves as 
	seen in fig 5 B-D and fig S5.
2) DataExtractor.m pulls from the data the information that we will later use
	to validate the model. Similar to other pieces of code used to
	generate figures 3E-G.
3) PiriformModelling.m performs the modelling simulations presented in fig 6
	and fig S7.

Prior to running the code, please load all variables from the following data files:
   - tetrode_recordings_formatted_for_alex_v2_ordered_aon
   - tetrode_recordings_formatted_for_alex_v2_ordered_apc
   - tetrode_recordings_formatted_for_alex_v2_ordered_ppc
And save them together into a data file named "tetrode_recordings_formatted_for_alex_v2_ordered".
