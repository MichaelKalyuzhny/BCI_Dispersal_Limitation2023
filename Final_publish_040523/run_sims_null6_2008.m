%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run simulation of null model 6 for every species
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Michael Kalyuzhny
%
% Date created: 23/08/2020
% Date last modified: 02/03/2021
%
% Load the species table and run the simulations for the parameters of each species. The simulations are run on a
% lanscape of 500X1000 meters, with EXACTLY the same tree location then in 2015 in BCI.
% Immigration is small, only to maintain some diversity.
%% set up:  

load('species_data1.mat','sp_dat')
load('Results/spatial_distribution_all_trees.mat');

cd('Results') %generate all results in this folder

species_run = find(~isnan(sp_dat{:,'HML2008AlphaFitted'})); %these are the species to run the simulation for

disp('Running Null6 with H2008 distances: ')
tic
parfor spp = 1:length(species_run)
    sp = species_run(spp)
    inp = struct;
    inp.J = sum(point_all_by_species(:, sp)); %number of individuals in the forest
    inp.Lx = 1000; %landscape edge (meters)
    inp.Ly = 500; %landscape edge (meters)
    inp.S_reg = 300; %number of species in regional pool
    inp.b = 2; %'b' parameter of the 2DT kernel (sensusp= Nathan 2012), = p+1 (sensu Clark 1999), representing heavy tailedness
    inp.a = sqrt(exp(sp_dat{sp,'HML2008AlphaFitted'})); %second parameter of the 2DT
    inp.imm_prob = 2/inp.J; %number of immigrating seeds to add
    inp.trees = points_all(point_all_by_species(:, sp),:);
    inp.D = 50; %the size of quadrats. If a smaller number is used, for some trees some quadrats will be empty
    
    % Time and sampling parameters (time in sweeps):
    inp.samps_tot = 1000; %how many samples to take?
    inp.samp_freq = 10; %how many sweeps between samples?
    inp.first_samp = 5000; %after how many sweeps to take first sample?
    inp.print_freq = 1000;
    inp.output_file = ['null6_Fiexd_Trees_2008_' sp_dat{sp,'SpeciesCode'}{1}];
    
    [~, ~, ~] = Null6_Fixed_Trees(inp)
end
toc

cd('..') %Go back to parent folder