%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forest simulator - Neutral, specieation and no boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Michael Kalyuzhny
%
% First written: 21/08/2020
% Last Modified: 28/12/2020
%
%% Description
%  Main community matrix - 3 columns: x,y,species ID. Seeds are generated
% (with same attributes) and dispersed with torus BC. with some probability species will be a new one All time is in sweeps. Fixed number of samples is taken, after equilibration. 
% Samples are saved in a preset folder and then merged
%
% Input: 'inp' struct with various fields
% output - spatial community 3D array: first dimension is the invividual, second dimension has three elements:
% x,y,species ID, third dimension is the sample. So for example, com_samp(100,:,20) will return the attributes of the
% 100th individual in the 20th sample: [x y species_ID].

function [com_samp] = sim2_N_spec_nobc_continue(inp, next_samp)

% %% Input patameters - example
% % This can be uncommented when transforming this function into a script.
% 
% inp.J = 1600; %number of individuals in the forest
% inp.L = 200; %landscape edge (meters)
% inp.S_init = 10; %number of species in regional pool
% inp.b = 2; %'b' parameter of the 2DT kernel (sensu Nathan 2012), = p+1 (sensu Clark 1999), representing heavy tailedness
% inp.a = 4; %second parameter of the 2DT
% inp.mu = 0.001; %number of immigrating seeds to add
% 
% % Time and sampling parameters (time in sweeps):
% inp.samps_tot = 100; %how many samples to take?
% inp.samp_freq = 10; %how many sweeps between samples?
% inp.first_samp = 0; %after how many sweeps to take first sample?
% inp.print_freq = 10;
% inp.output_file = 'sim2x_1.mat';

%% Initialization

%all time is in sweeps (=generations, J events)
t_now = 0;
next_samp_time = inp.samp_freq;
stop_samp_ind = inp.samps_tot + 1; %when this is the next sample - stop the run
next_samp_ind = next_samp; %index of next sample
next_print = inp.print_freq;
time_add = 1/inp.J; %the amount of time (in sweeps) that = a time step

% Import initial community:

l = load([inp.output_file '/s_' num2str(next_samp_ind-1) '_' inp.output_file '.mat']);
com = l.com;

next_sp = max(com(:,3)) + 100;

%% Main Loop:
tic

while next_samp_ind < stop_samp_ind

    seed = com(randi(inp.J,'single'),:); %choose a random parent
    u = rand; %temporary variable, to be transformed into distance:
    r = sqrt(((1-u)./(inp.a.^(2*inp.b-2))).^(1./(1-inp.b)) - inp.a^2 ); %transform to 2DT distance distribution
    theta = rand*2*pi; %direction
    seed(1) = seed(1) + r*cos(theta); %perform seed displacement - x coord.
    seed(2) = seed(2) + r*sin(theta); %perform seed displacement - y coord.

    % If seed is outside - repeat:
    while any(seed(1:2) > inp.L) || any(seed(1:2) < 0)
        seed = com(randi(inp.J,'single'),:); %choose a random parent
        u = rand; %temporary variable, to be transformed into distance:
        r = sqrt(((1-u)./(inp.a.^(2*inp.b-2))).^(1./(1-inp.b)) - inp.a^2 ); %transform to 2DT distance distribution
        theta = rand*2*pi; %direction
        seed(1) = seed(1) + r*cos(theta); %perform seed displacement - x coord.
        seed(2) = seed(2) + r*sin(theta); %perform seed displacement - y coord.
    end
    
    if rand < inp.mu %speciation
        seed(3) = next_sp;
        next_sp = next_sp + 1;
    end

    com(randi(inp.J,'single'),:) = seed; % replace a random tree with the seed
    
    %print time:
    if t_now >= next_print
        disp(['Runtime: ' sprintf('%0.8g',toc) ' sec., Time: ' sprintf('%0.8g',t_now) ' gen.'])
        next_print = next_print + inp.print_freq;
    end
  
    if t_now >= next_samp_time %if it's time to sample
        time_samp(next_samp_ind) = t_now; %record time
        
        save([inp.output_file '/s_' num2str(next_samp_ind) '_' inp.output_file '.mat'],'com');
        
        %next sampling:
        next_samp_time = next_samp_time + inp.samp_freq;
        next_samp_ind = next_samp_ind + 1;
    end
    
    t_now = t_now + time_add; %progress time
end

% Compose the pieces and save:

com_samp = NaN(inp.J,3,inp.samps_tot,'single');

for ss = 1:inp.samps_tot
    l = load([inp.output_file '/s_' num2str(ss) '_' inp.output_file '.mat'],'com');
    com_samp(:,:,ss) = l.com;
end

save([inp.output_file '.mat'], '-v7.3')
rmdir(inp.output_file, 's') % delete temporary directory for the samples
pause(3); %wait until delete is done, just in case

disp(['Finished! Runtime: ' sprintf('%0.8g',toc) ' sec., Time: ' sprintf('%0.8g',t_now) ' gen.'])

end