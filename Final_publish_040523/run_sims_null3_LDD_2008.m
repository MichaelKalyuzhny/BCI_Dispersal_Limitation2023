%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run simulation of null model 3 (increased long-distance dispesal) for every species
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Michael Kalyuzhny
%
% Date created: 23/08/2020
% Date last modified: 06/01/2021
%
% Load the species table and run the simulations for the parameters of each species. The simulations are run on a
% lanscape of 1200X1200 meters but with the same density as in BCI, which ia assounted for in the calculation of J.
% Immigration is small, only to maintain some diversity.
% In this version of the null model, a proportion of seeds inp.prop_random is dispersed at random
%% set up:

load('species_data1.mat','sp_dat')

cd('Results') %generate all results in this folder

species_run = find(~isnan(sp_dat{:,'HML2008AlphaFitted'})); %these are the species to run the simulation for

disp('Running Null3 (Long Distance Dispersal) with H2008 distances: ')
tic
parfor spp = 1:length(species_run)
    sp = species_run(spp)
    inp = struct;
    inp.J = round(1200*1200/(1000*500) * sp_dat{sp,'comparable_trees'}); %number of individuals in the forest
    inp.L = 1200; %landscape edge (meters)
    inp.S_reg = 300; %number of species in regional pool
    inp.b = 2; %'b' parameter of the 2DT kernel (sensusp= Nathan 2012), = p+1 (sensu Clark 1999), representing heavy tailedness
    inp.a = sqrt(exp(sp_dat{sp,'HML2008AlphaFitted'})); %second parameter of the 2DT
    inp.imm_prob = 10/inp.J; %number of immigrating seeds to add
    inp.prop_random = 0.1; %number of time steps a tree "waits" in the "stack".
    
    % Time and sampling parameters (time in sweeps):
    inp.samps_tot = 1000; %how many samples to take?
    inp.samp_freq = 10; %how many sweeps between samples?
    inp.first_samp = 5000; %after how many sweeps to take first sample?
    inp.print_freq = 1000;
    inp.output_file = ['null3_2008_' sp_dat{sp,'SpeciesCode'}{1}];
    
    [~, ~, ~] = Null3_LDD(inp)
end
toc

cd('..') %Go back to parent folder