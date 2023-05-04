%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% What is the effect of Dioecious
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Michael Kalyuzhny
%
% Date created: 07/04/2020
% Date last modified: 27/04/2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compare statistics between Dioecious and Monoecious species - only for null model
%
%% Set standard null simulations

inp_null = struct;

inp_null.J = 5500; %number of individuals in the forest
inp_null.L = 600; %landscape edge (meters)
inp_null.S_reg = 300; %number of species in regional pool
inp_null.b = 2; %'b' parameter of the 2DT kernel (sensu Nathan 2012), = p+1 (sensu Clark 1999), representing heavy tailedness
inp_null.a = 25; %second parameter of the 2DT
% inp_null.Q_HNDD = 10; %magnitude of HNDD
% inp_null.Q_CNDD_min = 10; %minimal value of CNDD
% inp_null.Q_CNDD_max = 10;
% inp_null.CNDD_dist = 7; %spatial scale at which CNDD drops
% inp_null.CNDD_steep = 6; %steepness of CNDD drop with distance
% inp_null.H_dist = 7;
% inp_null.H_steep = 6;    
% inp_null.D = 20;
% inp_null.alpha = 1; %total competition
inp_null.imm_prob = 0.001; %number of immigrating seeds to add

% Time and sampling parameters (time in sweeps):
inp_null.samps_tot = 5000; %how many samples to take?
inp_null.samp_freq = 10; %how many sweeps between samples?
inp_null.first_samp = 5000; %after how many sweeps to take first sample?
inp_null.print_freq = 1000;
inp_null.output_file = 'tmp';

dist_vals = [7 20 60];
J_vals = [5500];

inp_null = repmat(inp_null,length(dist_vals),length(J_vals));
%Dimensions: 1. dispersal distance 2. J

for jj=1:length(J_vals)
    for dd = 1:length(dist_vals)
        inp_null(dd,jj).a = dist_vals(dd)*2*gamma(1)/(sqrt(pi)*gamma(0.5)); %set disperdal distance
        inp_null(dd,jj).J = J_vals(jj);
        inp_null(dd,jj).output_file = ['Standard_D_'  sprintf('%0.4g',dist_vals(dd)) '_J_' sprintf('%0.4g', J_vals(jj))];
    end
end

save('inp_standard_null.mat','inp_null', 'J_vals', 'dist_vals')

%% Run standard null:
load('inp_standard_null.mat','inp_null')

parfor cc=1:numel(inp_null)
    inp = inp_null(cc);
    tmp = rand(1,cc); %create differences between simulations with similar seed

    [~, ~] = sim2_N(inp);
end

%% set Dioecious null:

% Set input (change filenames, leave everything else):

load('inp_standard_null.mat','inp_null', 'J_vals', 'dist_vals')
inp_Die = inp_null;

for jj=1:length(J_vals)
    for dd = 1:length(dist_vals)

        inp_Die(dd,jj).output_file = ['Die_D_'  sprintf('%0.4g',dist_vals(dd)) '_J_' sprintf('%0.4g', J_vals(jj))];
    end
end

save('inp_Die.mat','inp_Die', 'J_vals', 'dist_vals')

%% Run di null:

load('inp_Die.mat','inp_Die', 'J_vals', 'dist_vals')

parfor cc=1:numel(inp_Die)
    inp = inp_Die(cc);
    tmp = rand(1,cc); %create differences between simulations with similar seed

    [~, ~] = sim2_Die(inp);
end

%% Analyze standard null results:

% Set analysis parameters:
abd_bin_edges = [1 5:2:51 55:5:100 110:10:200 250:50:600 700:100:1000 1250 1500 5500];
dist_bin_edges = 0:10:200;

load('inp_standard_null.mat','inp_null', 'J_vals', 'dist_vals')

parfor cc = 1:numel(inp_null)
    f = [inp_null(cc).output_file '.mat']; %filename to load
    if isfile(f)
        l = load(f);
        com_samp = l.com_samp;

        Stats_regime3(com_samp, abd_bin_edges, dist_bin_edges, 0.25, 5, inp_null(cc).L, [inp_null(cc).output_file '_analysis_results.mat']); %run analysis for this regime and save result separately
    end
end

%% Analyze Dioecious simulation results:

% Set analysis parameters:
abd_bin_edges = [1 5:2:51 55:5:100 110:10:200 250:50:600 700:100:1000 1250 1500 5500];
dist_bin_edges = 0:10:200;

load('inp_Die.mat','inp_Die', 'J_vals', 'dist_vals')

parfor cc = 1:numel(inp_Die)
    f = [inp_Die(cc).output_file '.mat']; %filename to load
    if isfile(f)
        l = load(f);
        com_samp = l.com_samp;

        Stats_regime3(com_samp, abd_bin_edges, dist_bin_edges, 0.25, 5, inp_Die(cc).L, [inp_Die(cc).output_file '_analysis_results.mat']); %run analysis for this regime and save result separately
    end
end

%% Load data:

% Load standard:
load('inp_standard_null.mat')

res_null = cell(1);
res_null = repmat(res_null,size(inp_null));

for ii = 1:size(inp_null,1)
    for jj = 1:size(inp_null,2)
        fn = [inp_null(ii,jj).output_file '_analysis_results.mat']; %filename to try loading
        if isfile(fn)
            res_null{ii,jj} = load([inp_null(ii,jj).output_file '_analysis_results.mat']);
        end
    end
end

% Load Dioecious:
load('inp_Die.mat')

res_Die = cell(1);
res_Die = repmat(res_Die,size(inp_Die));

for ii = 1:size(inp_Die,1)
    for jj = 1:size(inp_Die,2)
        fn = [inp_Die(ii,jj).output_file '_analysis_results.mat']; %filename to try loading
        if isfile(fn)
            res_Die{ii,jj} = load([inp_Die(ii,jj).output_file '_analysis_results.mat']);
        end
    end
end

save('results_di_vs_mono.mat','res_Die','res_null','dist_vals')

%% Plot:

load('results_di_vs_mono.mat')

figure()

for dd = 1:3
    s1 = subplot(2,3,dd);
    hold on
    set(s1,'FontSize',16,'FontWeight','bold','LineWidth',1.5,'XScale','log');
    
    plot(res_null{dd,1}.abd_bin_centers(1:end-1), res_Die{dd,1}.log_N20(1:end-1) - res_null{dd,1}.log_N20(1:end-1), '-', 'LineWidth',1)
    
    title(['{\itD} = ' num2str(dist_vals(dd)) ' m.']);
    %line([res_null{dd,1}.abd_bin_centers(2) res_null{dd,1}.abd_bin_centers(end)],[0 0], 'LineStyle', '--','color',[.5 .5 .5], 'LineWidth',1.5)
    line([res_null{dd,1}.abd_bin_centers(2) s1.XLim(2)],[0 0], 'LineStyle', '--','color',[.5 .5 .5], 'LineWidth',1.5)
    s1.YLim(1) = -.01;
    xlim([5 1500])
    xticks([10 100 1000])
    
    if dd == 1
        ylabel('{\itEA}(20) of Dioecious')
        annotation('arrow',[.95 .95],[.55 .65],'LineWidth',2.5)
        t = text(1000,.2,'Overdispersion', 'FontSize', 14', 'FontWeight', 'bold');
        set(t,'Rotation',90)
    end
    
end

for dd = 1:3
    s1 = subplot(2,3,3+dd);
    hold on
    set(s1,'FontSize',16,'FontWeight','bold','LineWidth',1.5,'XScale','log');
    
    plot(res_null{dd,1}.abd_bin_centers(1:end-1), res_Die{dd,1}.mean_log_dist(1:end-1) - res_null{dd,1}.mean_log_dist(1:end-1), '-', 'LineWidth',1)
    line([res_null{dd,1}.abd_bin_centers(2) s1.XLim(2)],[0 0], 'LineStyle', '--','color',[.5 .5 .5], 'LineWidth',1.5)
    xlim([5 1500])
    xticks([10 100 1000])
    
    if dd == 1
        ylabel('{\itED} of Dioecious')
        annotation('arrow',[.95 .95],[.45 .35],'LineWidth',2.5)
        t = text(1000,.2,'Overdispersion', 'FontSize', 14', 'FontWeight', 'bold');
        set(t,'Rotation',90)
    end
    
    if dd == 2
        xlabel('Abundance')
    end
end
