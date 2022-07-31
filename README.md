# Dissertation_Julia
Julia and MATLAB code developed for the paper titled "Fixed and Distributed Gene Expression Time Delays in Reaction-Diffusion Systems"

The code allows for both analtyical and numerical exploration of reaction-diffusion models with delay. Numerical simulations for stiff Reaction-Diffusion systems, with Schnakenberg or Gierer-Meinhardt kinetics, can be run from the 'main.jl' file. The script allows incorporation of gene-expression time delays, in the form of a fixed delay, or distributed delay. Distributions considered include the skew and symmetric truncated 
Gaussian pdfs. 

Below is an overview of the code used to generate the figures in the given paper:

Figure 1: Produced uisng the 'Distributions.m' script - by setting the relevant parameter values and running the script (sections).

Figure 2: Produced using the 'TuringSpace.m' script - running either section 1 or 2 of the code for the relevant model kinetics, then running section 3 to  
          compute and plot the Turing space.
          
Figures 3 and 4: Produced using the 'fixed.m' script - running the first section to pick the model, followed by the second section titled 'Bifurcation  
                 diagram' with the selected parameter values. 
                 
Figure 5: The Fixed delay curve produced analogously to Figure 6 (see below). The skewed curves can be computed by running the
          'dist.m' script (with the given parameter values). The curve is plotted by plotting 'tauvec' against 'results_vec'.

Figure 6: Produced using the 'fixed.m' script - running the first section to pick the model, followed by the third section title 'Dispersion relation' with
          the selected parameter values.
          
Figures 7,8,9 - The numerical results are produced by running the 'main.jl' script (for fixed delay, with set threshold and initial conditions), and then
                then plotting 'tauvec' against 'ttakenvec'. The 'analytical' results in Figure are produced by using eqn 24. The values for $A_k(T^*)$ are
                stored in 'ATvec', and the value of $A_k$(0) stored as 'A0'. We note that the dispersion analysis (as produced in Figure 6) to first be
                conducted to find the dominant k and $\lambda_k$ values. The relevant model kinetics and steady-states must be ammended in the 
                'func_mod.jl' script (e.g. LI or GM kinetics).
                
Figures 10,11 - Produced by running the 'main.jl' script - by setting the relevent problem type and parameters for numerical simulation. The solution can
                be plotted by using the 'turing_plot()' function call - commented out on line 70. 

Figures 12,13 - Plots (a) produced analogously to Figure 2. Plots (b), (c) and (d) produced analogously to Figures 10 and 11.

Figure 14 - Taken from cited paper.

Figures 15,16 - Produced analogously to Figures 10,11 with different initial conditions set in 'func_mod.jl' script.

Figure 17 - Produced analogously to Figures 10,11 with different boundary conditions set in 'main.jl' script as 'Dirichlet'. 
              
          
