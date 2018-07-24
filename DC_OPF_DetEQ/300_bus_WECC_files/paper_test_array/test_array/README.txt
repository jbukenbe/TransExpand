Poorly maintained right now as experiments are ongoing.

Currently setup to run by calling LF_array_call which was designed to be run in an HPC array
I think this only calls run_plan_in_gams, which is the other function that is useful here.
Note some sloppy indexing that was useful for debugging runs that got stuck for some reason on the HPC,
some indexing is done in LF_array_call and some is done in run_plan_in_gams, and saved periodically so partial results can be debugged.

Data files are not clear to read, as they are stored compactly in a way the lets them be multiplied together quicky with matlab vector notation. 
 
OPF_file is the opf GAMS file that is used.