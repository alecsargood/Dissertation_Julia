# Dissertation_Julia
Julia code developed for the dissertation titled "Gene-Expression Time Delays in Reaction Diffusion Systems"

Numerical simulations for stiff Reaction-Diffusion systems, with Schnakenberg or Gierer-Meinhardt kinetics, can be run from the main.jl file. The script allows
incorporation of gene-expression time delays, in the form of a fixed delay, or distributed delay. Distributions considered include the skew and symmetric truncated 
Gaussian pdfs. 

main.jl: The main file where numerical solutions are generated.
func_mod.jl: Script where parameter structure is initialised, and relevent functions are defined.
Remaining scripts are auxillary files where problem dependent information is held (e.g. RHS of system with fixed delay, distributed delay)
