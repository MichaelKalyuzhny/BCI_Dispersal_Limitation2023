%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Power analysis 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Michael Kalyuzhny
%
% Date created: 17/04/2023
% Date last modified: 17/04/2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% Many of the parameters below, aimed at creating differences beween CNDD and HNDD
% are not used.
% Evaluates statistics witn and without edge effects

%% Analyze to obtain raw, sample-level data:

% Set analysis parameters:
abd_bin_edges = [1 5:2:51 55:5:100 110:10:200 250:50:600 700:100:1000 1250 1500 5500];
min_abd = 5;
max_abd = 1500;
prealloc = 500000;
Tot_L = 1200;

%% Run analysis:

load('Edge2_21.mat')

parfor nn = 1:numel(inp_all2)
    if isfile([inp_all2(nn).output_file '.mat'])
        l = load(inp_all2(nn).output_file);
        com_samp = l.com_samp;

        Stats_regime_meta_export(com_samp, abd_bin_edges, min_abd, max_abd, prealloc, Tot_L, Lx, Ly, [inp_all2(nn).output_file '_raw_stats.mat']);

    end
end

%% Create intervals between low and high percentiles:

load('Edge2_21.mat')

thresh_null = 50;
pctiles = [2.5 97.5];

num_abd_bins = length(abd_bin_edges) - 1;
len_cndd = length(magnitude_vals)-1;
len_dist = length(dist_vals);

n20_envelopes = nan(2,num_abd_bins,len_dist); %first: upper/lower, 2nd: bins
nnd_envelopes = nan(2,num_abd_bins,len_dist); %first: upper/lower, 2nd: bins

for dd = 1:len_dist
    l = load([inp_all2(1,dd).output_file '_raw_stats.mat']); %NEUTRAL
    
    n20_envelopes(:,:,dd) = prctile(l.log_N20s,[pctiles(2) pctiles(1)],2)';
    nnd_envelopes(:,:,dd) = prctile(l.log_NNDs,[pctiles(2) pctiles(1)],2)';

    
    to_rem = sum(~isnan(l.log_N20s),2) < thresh_null;
    n20_envelopes(:,to_rem',dd) = NaN;
    nnd_envelopes(:,to_rem',dd) = NaN;
    
end


%% Load data and compute power:

load('Edge2_21.mat')

min_compute_power = 20;

% Compute levels:
num_abd_bins = length(abd_bin_edges) - 1;
len_cndd = length(magnitude_vals)-1;
len_dist = length(dist_vals);

% Preallocate:
EA20_power = nan(num_abd_bins, len_cndd, len_dist);
ED_power = nan(num_abd_bins, len_cndd, len_dist);

% Import raw data:
for nn=1:len_cndd
    for dd = 1:len_dist
        
        fn = [inp_all2(nn+1,dd).output_file '_raw_stats.mat'];
        if isfile(fn)
            if (nn == 1) && (dd == 1)
                res = load(fn);
                res = repmat(res,len_cndd, len_dist);
            else
                res(nn,dd) = load(fn);
            end
            
            % Compute power:
            for bb = 1:num_abd_bins
                
                % Power of EA20:
                vals = res(nn,dd).log_N20s(bb,:);
                vals = vals(~isnan(vals));
                
                if length(vals) >= min_compute_power
                    EA20_power(bb,nn,dd) = sum((vals < n20_envelopes(2,bb,dd)) | (vals > n20_envelopes(1,bb,dd)))/length(vals);
                end

                % Power of ED:
                vals = res(nn,dd).log_NNDs(bb,:);
                vals = vals(~isnan(vals));
                
                if length(vals) >= min_compute_power
                    ED_power(bb,nn,dd) = sum((vals < nnd_envelopes(2,bb,dd)) | (vals > nnd_envelopes(1,bb,dd)))/length(vals);
                end
            end
        end
        
    end
end

abds = abd_bin_edges(1:end-1) + diff(abd_bin_edges)/2; %assuming the abundance bins are the same across the board


save('power_results.mat', 'abds', 'EA20_power', 'ED_power')

%% Generate plot - plot only two regimes:

load('power_results.mat')

figure()

% blue = CNDD regime 1, red = CNDD regime 3
% dash = D = 7, solid = D = 20

s1 = subplot(1,3,1);
hold on
set(s1,'FontSize',14,'FontWeight','bold','LineWidth',1.5,'XScale','log');
ylabel('Power')
plot(abds, EA20_power(:,3,1), 'r--', abds, EA20_power(:,3,2), 'r-', ...
    abds, EA20_power(:,1,1), 'b--', abds, EA20_power(:,1,2), 'b-' , 'LineWidth',1)

title('Excess Abundance ({\itEA}(20))')
text(10,.9,'A','FontSize',16,'FontWeight','bold')
legend({'{\itQ_C} = 12, {\itQ_H} = 0, {\itD} = 7 m.', '{\itQ_C} = 12, {\itQ_H} = 0, {\itD} = 20 m.', ...
    '{\itQ_C} = {\itQ_H} = 3, {\itD} = 7 m.', '{\itQ_C} = {\itQ_H} = 3, {\itD} = 20 m.'})
xticks([10 100 1000])

s1 = subplot(1,3,2);
hold on
set(s1,'FontSize',14,'FontWeight','bold','LineWidth',1.5,'XScale','log');
plot(abds, ED_power(:,3,2), 'r-', abds, ED_power(:,3,1), 'r--', ...
    abds, ED_power(:,1,2), 'b-', abds, ED_power(:,1,1), 'b--', 'LineWidth',1)

title('Excess Distance ({\itED})')
text(10,.9,'B','FontSize',16,'FontWeight','bold')
xticks([10 100 1000])
xlabel('Abudnance')

% %% Generate plot - original:
% 
% figure()
% 
% % blue = CNDD regime 1, red = CNDD regime 2, magenta = CNDD regime 3
% % dash = D = 7, solid = D = 20
% 
% abds = abd_bin_edges(1:end-1) + diff(abd_bin_edges)/2; %assuming the abundance bins are the same across the board
% 
% s1 = subplot(1,3,1);
% hold on
% set(s1,'FontSize',14,'FontWeight','bold','LineWidth',1.5,'XScale','log');
% ylabel('Power')
% plot(abds, EA20_power(:,3,1), 'm--', abds, EA20_power(:,3,2), 'm-', ...
%     abds, EA20_power(:,2,1), 'r--', abds, EA20_power(:,2,2), 'r-', ...
%     abds, EA20_power(:,1,1), 'b--', abds, EA20_power(:,1,2), 'b-' , 'LineWidth',1)
% 
% title('Excess Abundance ({\itEA}(20))')
% text(10,.9,'A','FontSize',16,'FontWeight','bold')
% legend({'{\itQ_C} = 12, {\itQ_H} = 0, {\itD} = 7 m.', '{\itQ_C} = 12, {\itQ_H} = 0, {\itD} = 20 m.', ...
%     '{\itQ_C} = 6, {\itQ_H} = 0, {\itD} = 7 m.', '{\itQ_C} = 6, {\itQ_H} = 0, {\itD} = 20 m.', ...
%     '{\itQ_C} = {\itQ_H} = 3, {\itD} = 7 m.', '{\itQ_C} = {\itQ_H} = 3, {\itD} = 20 m.'})
% xticks([10 100])
% 
% s1 = subplot(1,3,2);
% hold on
% set(s1,'FontSize',14,'FontWeight','bold','LineWidth',1.5,'XScale','log');
% plot(abds, ED_power(:,3,2), 'm-', abds, ED_power(:,3,1), 'm--', ...
%     abds, ED_power(:,2,2), 'r-', abds, ED_power(:,2,1), 'r--', ...
%     abds, ED_power(:,1,2), 'b-', abds, ED_power(:,1,1), 'b--', 'LineWidth',1)
% 
% title('Excess Distance ({\itED})')
% text(10,.9,'B','FontSize',16,'FontWeight','bold')
% xticks([10 100])
% xlabel('Abudnance')
