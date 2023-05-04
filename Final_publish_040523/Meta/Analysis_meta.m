%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Does the torus BC with 'm' approximate well a metacommunity?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Michael Kalyuzhny
%
% Date created: 13/04/2023
% Date last modified: 19/04/2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compare statistics between Dioecious and Monoecious species
%
%% Set standard null simulations

inp_meta1 = struct;

inp_meta1.J = 550000; %number of individuals in the forest
inp_meta1.L = 6000; %landscape edge (meters)
inp_meta1.S_init = 200; %number of species in regional pool
inp_meta1.b = 2; %'b' parameter of the 2DT kernel (sensu Nathan 2012), = p+1 (sensu Clark 1999), representing heavy tailedness
inp_meta1.a = 25; %second parameter of the 2DT
inp_meta1.mu = 0.00005; %number of immigrating seeds to add

% Time and sampling parameters (time in sweeps):
inp_meta1.samps_tot = 2000; %how many samples to take?
inp_meta1.samp_freq = 10; %how many sweeps between samples?
inp_meta1.first_samp = 10000; %after how many sweeps to take first sample?
inp_meta1.print_freq = 100;
inp_meta1.output_file = 'tmp';

dist_vals = [7 20 60];

inp_meta1 = repmat(inp_meta1,length(dist_vals),1);
%Dimensions: 1. dispersal distance 2. J


    for dd = 1:length(dist_vals)
        inp_meta1(dd).a = dist_vals(dd)*2*gamma(1)/(sqrt(pi)*gamma(0.5)); %set disperdal distance
        inp_meta1(dd).output_file = ['meta_nobc2_D_'  sprintf('%0.4g',dist_vals(dd))];
    end

save('inp_meta1.mat','inp_meta1', 'dist_vals')

%% Run meta null:
load('inp_meta1.mat')

parfor cc=1:numel(inp_meta1)
    inp = inp_meta1(cc);
    tmp = rand(1,cc); %create differences between simulations with similar seed

    [~, ~] = sim2_N_spec_nobc(inp);
end

%% Analyze meta:

abd_bin_edges = [1 5:2:51 55:5:100 110:10:200 250:50:600 700:100:1000 1250 1500 5500];
dist_bin_edges = 0:10:200;
Lx = 1000;
Ly = 500;

load('inp_meta1.mat')

parfor nn = 1:numel(inp_meta1)
    if isfile([inp_meta1(nn).output_file '.mat'])
        l = load(inp_meta1(nn).output_file);
        com_samp = l.com_samp;

        Stats_regime_meta(com_samp, abd_bin_edges, dist_bin_edges, 5, 1500, inp_meta1(nn).L, Lx, Ly, [inp_meta1(nn).output_file '_regular_analysis_results.mat']); %run analysis for this regime and save result separately        Stats_regime_meta(com_samp, abd_bin_edges, dist_bin_edges, 5, 1500, inp_meta1(nn).L, Lx, Ly, [inp_meta1(nn).output_file '_regular_analysis_results.mat']); %run analysis for this regime and save result separately

    end
end

%% Set parameters for smaller landscapes

inp_meta_s1 = struct;

inp_meta_s1.J = 22000; %number of individuals in the forest
inp_meta_s1.L = 1200; %landscape edge (meters)
inp_meta_s1.S_reg = 300; %number of species in regional pool
inp_meta_s1.b = 2; %'b' parameter of the 2DT kernel (sensu Nathan 2012), = p+1 (sensu Clark 1999), representing heavy tailedness
inp_meta_s1.a = 25; %second parameter of the 2DT
inp_meta_s1.imm_prob = 0.00005; %number of immigrating seeds to add

% Time and sampling parameters (time in sweeps):
inp_meta_s1.samps_tot = 8000; %how many samples to take?
inp_meta_s1.samp_freq = 10; %how many sweeps between samples?
inp_meta_s1.first_samp = 1000; %after how many sweeps to take first sample?
inp_meta_s1.print_freq = 1000;
inp_meta_s1.output_file = 'tmp';

dist_vals = [7 20 60];
J = inp_meta_s1.J; 
m_vals = [0.2/J 2/J 10/J 100/J];

inp_meta_s1 = repmat(inp_meta_s1,length(dist_vals),length(m_vals));
%Dimensions: 1. dispersal distance 2. m

for mm = 1:length(m_vals)
    for dd = 1:length(dist_vals)
        inp_meta_s1(dd,mm).a = dist_vals(dd)*2*gamma(1)/(sqrt(pi)*gamma(0.5)); %set disperdal distance
        inp_meta_s1(dd,mm).imm_prob = m_vals(mm);
        inp_meta_s1(dd,mm).output_file = ['meta_small1_D_'  sprintf('%0.4g',dist_vals(dd)) '_m_'  sprintf('%0.4g',m_vals(mm))];
    end
end

save('inp_meta_s1.mat','inp_meta_s1', 'dist_vals','m_vals')

%% Run smaller landscapes with m:

load('inp_meta_s1.mat')

parfor cc=1:numel(inp_meta_s1)
    inp = inp_meta_s1(cc);
    tmp = rand(1,cc); %create differences between simulations with similar seed

    [~, ~] = sim2_N(inp);
end

%% Analyze smaller landscapes:

abd_bin_edges = [1 5:2:51 55:5:100 110:10:200 250:50:600 700:100:1000 1250 1500 5500];
dist_bin_edges = 0:10:200;
Lx = 1000;
Ly = 500;

load('inp_meta_s1.mat')

parfor nn = 1:numel(inp_meta_s1)
    f = [inp_meta_s1(nn).output_file '.mat'];
    if isfile(f)
        l = load(f);
        com_samp = l.com_samp;

        Stats_regime_meta(com_samp, abd_bin_edges, dist_bin_edges, 5, 1500, inp_meta_s1(nn).L, Lx, Ly, [inp_meta_s1(nn).output_file '_regular_analysis_results.mat']); %run analysis for this regime and save result separately        Stats_regime_meta(com_samp, abd_bin_edges, dist_bin_edges, 5, 1500, inp_meta1(nn).L, Lx, Ly, [inp_meta1(nn).output_file '_regular_analysis_results.mat']); %run analysis for this regime and save result separately

    end
end

%% Import all data:
load('inp_meta1.mat')

% Import large metacommunity data:

for dd = 1:length(dist_vals)

    if dd == 1
        large = load([inp_meta1(dd).output_file '_regular_analysis_results.mat']);
        large = repmat(large,length(dist_vals),1);
    else
        large(dd) = load([inp_meta1(dd).output_file '_regular_analysis_results.mat']);
    end
end

% Import small local community data:

load('inp_meta_s1.mat')

for dd = 1:length(dist_vals)
    for mm = 1:length(m_vals)

        if (dd == 1) && (mm == 1)
            small = load([inp_meta_s1(dd,mm).output_file '_regular_analysis_results.mat']);
            small = repmat(small, length(dist_vals), length(m_vals));
        else
            small(dd,mm) = load([inp_meta_s1(dd,mm).output_file '_regular_analysis_results.mat']);
        end
    end
end

save('meta_results.mat','small','large', 'm_vals','dist_vals')

%% Generate figures

load('meta_results.mat')

figure()

abds = small(1).abd_bin_centers; %assuming the abundance bins are the same across the board

for dd = 1:3
    
    % Densities:
    s1 = subplot(2,3,dd);
    hold on
    set(s1,'FontSize',14,'FontWeight','bold','LineWidth',1.5,'XScale','log');
    
    plot(abds(1:end-1),large(dd).N20(1:end-1),'b--','LineWidth',1)
    for mm = 1:length(m_vals)
        plot(abds(1:end-1),small(dd,mm).N20(1:end-1),'LineWidth',1)
    end
    xticks([10 100 1000])
    title(['{\itD} = ' num2str(dist_vals(dd)) ' m.']); 
    if dd == 1
        ylabel('Neighbors within 20 m.')
    end
    
    if dd == 3
        legend({'Metacomm.','{\itk} = 0.2/{\itJ}', '{\itk} = 2/{\itJ}', '{\itk} = 10/{\itJ}', '{\itk} = 100/{\itJ}'});
    end
    
    % Distances:
    s1 = subplot(2,3, 3 + dd);
    hold on
    set(s1,'FontSize',14,'FontWeight','bold','LineWidth',1.5,'XScale','log');
    
    plot(abds(1:end-1),large(dd).log_mean_dist(1:end-1),'b--','LineWidth',1)
    
    for mm = 1:length(m_vals)
        plot(abds(1:end-1),small(dd,mm).log_mean_dist(1:end-1),'LineWidth',1)
    end
    xticks([10 100 1000])
    
    if dd == 1
        ylabel('Nearest neighbor distance (m.)')
    end
    
    if (dd == 2)
        xlabel('Abundance')
    end
    
end
