%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% what is the effect of a combination of CNDD and HNDD?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Michael Kalyuzhny
%
% Date created: 17/04/2023
% Date last modified: 27/04/2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% Many of the parameters below, aimed at creating differences beween CNDD and HNDD
% are not used.
% Evaluates statistics witn and without edge effects
% the third regime is not used as it gives similar results.
%% Set simulations to analyze

inp_all2 = struct;

inp_all2.J = 22000; %number of individuals in the forest
inp_all2.L = 1200; %landscape edge (meters)
inp_all2.S_reg = 300; %number of species in regional pool
inp_all2.b = 2; %'b' parameter of the 2DT kernel (sensu Nathan 2012), = p+1 (sensu Clark 1999), representing heavy tailedness
inp_all2.a = 25; %second parameter of the 2DT
inp_all2.Q_HNDD = 10; %magnitude of HNDD
inp_all2.Q_CNDD_min = 10; %minimal value of CNDD
inp_all2.Q_CNDD_max = 10;
inp_all2.CNDD_dist = 7; %spatial scale at which CNDD drops
inp_all2.CNDD_steep = 6; %steepness of CNDD drop with distance
inp_all2.H_dist = 7;
inp_all2.H_steep = 6;    
inp_all2.D = 20;
inp_all2.alpha = 1; %total competition
inp_all2.imm_prob = 0.0005; %number of immigrating seeds to add

Lx = 1000;
Ly = 500;

% Time and sampling parameters (time in sweeps):
inp_all2.samps_tot = 1000; %how many samples to take?
inp_all2.samp_freq = 10; %how many sweeps between samples?
inp_all2.first_samp = 1000; %after how many sweeps to take first sample?
inp_all2.print_freq = 1000;
inp_all2.output_file = 'tmp';

magnitude_vals = {[0 0],[3 3], [6 0],[12 0]}; %first number is CNDD, second is HNDD
dist_vals = [7 20];

inp_all2 = repmat(inp_all2,length(magnitude_vals),length(dist_vals));
%Dimensions: 1. DD regime; 2. dispersal distance

for nn=1:length(magnitude_vals)
        for dd = 1:length(dist_vals)
            inp_all2(nn,dd).Q_HNDD = magnitude_vals{nn}(2);
            inp_all2(nn,dd).Q_CNDD_min = magnitude_vals{nn}(1);
            inp_all2(nn,dd).Q_CNDD_max = magnitude_vals{nn}(1);
            inp_all2(nn,dd).a = dist_vals(dd)*2*gamma(1)/(sqrt(pi)*gamma(0.5)); %set disperdal distance
            
            if nn == 1 %neutral case
                inp_all2(nn,dd).samps_tot = 5000;
            end
                      
            inp_all2(nn,dd).output_file = ['Edge2_QC_'  sprintf('%0.4g',magnitude_vals{nn}(1)) '_QH_' sprintf('%0.4g',magnitude_vals{nn}(2)) '_D_' sprintf('%0.4g',dist_vals(dd))];
        end
end

save('Edge2_21.mat','inp_all2', 'magnitude_vals', 'dist_vals','Lx','Ly')

%% Run:
parfor nn=1:numel(inp_all2)
    inp = inp_all2(nn);
    tmp = rand(1,nn); %create differences between simulations with similar seed

    [~, ~, ~] = sim3(inp);
end


%% Analyze data both usual ways:

% Set analysis parameters:
abd_bin_edges = [1 5:2:51 55:5:100 110:10:200 250:50:600 700:100:1000 1250 1500 5500];
dist_bin_edges = 0:10:200;

load('Edge2_21.mat')

parfor nn = 1:numel(inp_all2)
    if isfile([inp_all2(nn).output_file '.mat'])
        l = load(inp_all2(nn).output_file);
        com_samp = l.com_samp;

        Stats_regime_meta(com_samp, abd_bin_edges, dist_bin_edges, 5, 1500, inp_all2(nn).L, Lx, Ly, [inp_all2(nn).output_file '_regular_analysis_results.mat']); %run analysis for this regime and save result separately
        Stats_regime_unobs(com_samp, abd_bin_edges, dist_bin_edges, 5, 1500, inp_all2(nn).L, Lx, Ly, [inp_all2(nn).output_file '_unobs_analysis_results.mat']); %run analysis for this regime and save result separately

    end
end

%% Load all data:

load('Edge2_21.mat')

res_torus = struct;
res_unobs = struct;
res_torus = repmat(res_torus, size(inp_all2));
res_unobs = repmat(res_unobs, size(inp_all2));

for nn=1:length(magnitude_vals)
    for dd = 1:length(dist_vals)
        
        % Regular analysis:
        fn = [inp_all2(nn,dd).output_file '_regular_analysis_results.mat'];
        if isfile(fn)
            if (nn == 1) && (dd == 1)
                res_torus = load(fn);
            else
                res_torus(nn,dd) = load(fn);
            end
        end
        
        % Analysis with unobserved distances:
        fn = [inp_all2(nn,dd).output_file '_unobs_analysis_results.mat'];
        if isfile(fn)
            if isfile(fn)
                if (nn == 1) && (dd == 1)
                    res_unobs = load(fn);
                else
                    res_unobs(nn,dd) = load(fn);
                end
            end
            
        end
    end
end


%% Compute EA and ED:

% The assumption is the first level of CNDD is the neutral, and all is
% compared to it.

abd_thresh = 20;

% Compute levels:
num_abd_bins = length(res_torus(1,1).samps_per_bin);
len_cndd = length(magnitude_vals)-1;
len_dist = length(dist_vals);

% Preallocate:
EAs_torus = nan(num_abd_bins, len_dist, len_cndd);
EAs_unobs = nan(num_abd_bins, len_dist, len_cndd);
EDs_torus = nan(num_abd_bins, len_dist, len_cndd);
EDs_unobs = nan(num_abd_bins, len_dist, len_cndd);

for dd = 1:len_dist %dispersal distance
    for cc = 1:len_cndd %CNDD regimes
        EAs_torus(:, dd, cc) = res_torus(cc + 1,dd).log_N20 - res_torus(1,dd).log_N20; %the 1st cndd level is neutral
        EAs_unobs(:, dd, cc) = res_unobs(cc + 1,dd).log_N20 - res_unobs(1,dd).log_N20;
        
        EDs_unobs(:, dd, cc) = res_unobs(cc + 1,dd).mean_log_dist - res_unobs(1,dd).mean_log_dist;
        EDs_torus(:, dd, cc) = res_torus(cc + 1,dd).mean_log_dist - res_torus(1,dd).mean_log_dist;
        
        % remove cases with too few observations:
        to_rem_torus = (res_torus(cc + 1,dd).samps_per_bin < abd_thresh) | (res_torus(1,dd).samps_per_bin < abd_thresh);
        to_rem_unobs = (res_unobs(cc + 1,dd).samps_per_bin < abd_thresh) | (res_unobs(1,dd).samps_per_bin < abd_thresh);
        
        EAs_torus(to_rem_torus, dd, cc) = NaN;
        EDs_torus(to_rem_torus, dd, cc) = NaN;
        
        EAs_unobs(to_rem_unobs, dd, cc) = NaN;
        EDs_unobs(to_rem_unobs, dd, cc) = NaN;

    end
end

abds = res_torus.abd_bin_centers; %assuming the abundance bins are the same across the board

save('results_edge.mat','EDs_unobs','EAs_unobs','EDs_torus','EAs_torus','dist_vals','abds')

%% Generate plot:

load('results_edge.mat')

Regime_titles = {'{\itQ_C} = {\itQ_H} = 3', '{\itQ_C} = 6, {\itQ_H} = 0', '{\itQ_C} = 12, {\itQ_H} = 0'};

figure()

% blue = torus, red = no edge
% dash = D = 7, solid = D = 20


s1 = subplot(2,2,1);
hold on
set(s1,'FontSize',16,'FontWeight','bold','LineWidth',1.5,'XScale','log');
ylabel(Regime_titles{1})
plot(abds,EAs_torus(:,1,1),'b--', abds,EAs_torus(:,2,1),'b-', ...
    abds,EAs_unobs(:,1,1), 'r--', abds,EAs_unobs(:,2,1), 'r-', 'LineWidth',1)
title('Excess Abundance ({\itEA}(20))')
text(100,-0.7,'A','FontSize',16,'FontWeight','bold')
xticks([10 100 1000])

s1 = subplot(2,2,2);
hold on
set(s1,'FontSize',16,'FontWeight','bold','LineWidth',1.5,'XScale','log');
plot(abds,EDs_torus(:,1,1),'b--', abds,EDs_torus(:,2,1),'b-', ...
    abds,EDs_unobs(:,1,1), 'r--', abds,EDs_unobs(:,2,1), 'r-', 'LineWidth',1)
title('Excess Distance ({\itED})')
text(100,0.3,'B','FontSize',16,'FontWeight','bold')
xticks([10 100 1000])

% s1 = subplot(3,2,3);
% hold on
% set(s1,'FontSize',16,'FontWeight','bold','LineWidth',1.5,'XScale','log');
% ylabel(Regime_titles{2})
% plot(abds,EAs_torus(:,1,2),'b--', abds,EAs_torus(:,2,2),'b-', ...
%     abds,EAs_unobs(:,1,2), 'r--', abds,EAs_unobs(:,2,2), 'r-', 'LineWidth',1)
% text(100,-0.7,'C','FontSize',16,'FontWeight','bold')
% xlabel('Abundance')
% 
% s1 = subplot(3,2,4);
% hold on
% set(s1,'FontSize',16,'FontWeight','bold','LineWidth',1.5,'XScale','log');
% plot(abds,EDs_torus(:,1,2),'b--', abds,EDs_torus(:,2,2),'b-', ...
%     abds,EDs_unobs(:,1,2), 'r--', abds,EDs_unobs(:,2,2), 'r-', 'LineWidth',1)
% text(100,0.7,'D','FontSize',16,'FontWeight','bold')
% legend({'Torus, {\itD} = 7 m.','Torus, {\itD} = 20 m.', 'No edge, {\itD} = 7 m.', 'No edge, {\itD} = 20 m.'})


s1 = subplot(2,2,3);
hold on
set(s1,'FontSize',16,'FontWeight','bold','LineWidth',1.5,'XScale','log');
ylabel(Regime_titles{3})
plot(abds,EAs_torus(:,1,3),'b--', abds,EAs_torus(:,2,3),'b-', ...
    abds,EAs_unobs(:,1,3), 'r--', abds,EAs_unobs(:,2,3), 'r-', 'LineWidth',1)
text(100,-0.7,'C','FontSize',16,'FontWeight','bold')
xlabel('Abundance')
xticks([10 100])

s1 = subplot(2,2,4);
hold on
set(s1,'FontSize',16,'FontWeight','bold','LineWidth',1.5,'XScale','log');
plot(abds,EDs_torus(:,1,3),'b--', abds,EDs_torus(:,2,3),'b-', ...
    abds,EDs_unobs(:,1,3), 'r--', abds,EDs_unobs(:,2,3), 'r-', 'LineWidth',1)
text(100,0.7,'D','FontSize',16,'FontWeight','bold')
legend({'Torus, {\itD} = 7 m.','Torus, {\itD} = 20 m.', 'No edge, {\itD} = 7 m.', 'No edge, {\itD} = 20 m.'})
xticks([10 100])

