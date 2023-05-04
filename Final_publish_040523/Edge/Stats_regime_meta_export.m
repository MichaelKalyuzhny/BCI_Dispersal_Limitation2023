%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute multiple statistics for the samples of a parameter regime
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Michael Kalyuzhny
% Date created: 29/12/2020
% Date last modified: 17/04/2023
%
% Input: 
%       1. 'com_samp' array, where first dimension is the invividual, second dimension has three elements: x,y,species ID, third dimension is the sample.
%       2. 'abd_bin_edges' array - the left edge of each abundance bin. first entry should be one, last should be
%       community size
%       3. 'dist_bin_edges' array of the left edge of every bin. last value is the final edge.
%       4. 'prealloc_size' - fraction of total observations to preallocate bin
%       5. 'min_abd' - the minimal abundance to consider
%       6. 'L' - the edge length of the landscape
%       7. 'filename' - if not NaN, the output will be saved to this location.
%c
% Output:
%       1. 'log_mean_dist' - the mean of the log conspecific nearest neighbor distance in each abundance bin
%       2. 'mean_log_dist - the log of the mean conspecific nearest neighbor distance in each abundance bin
%       3. 'CE' - the mean Clark-Evans statistic in each abundance bin
%       4. 'N20' - the number of neighbors within 20 meters, by abundance bins
%       5. log_N20 the log number of individuals 
%       6. 'ND' - the number of neighbors per distance bin (neighbor distance)
%       7. log_ND - abundnace bins by distance
%       8. 'samps_per_bin' - how many samples were obtained for each abundance bin?
%       9. 'abd_bin_centers' - the centers of abundance bins
%       10. 'dist_bin_centers' - the centers of distance bins
%       11. 'PCF' - the pair correlation function averaged for all species with abundnace above min_abd
%       12. RND20 - relative neighborhood density at 20 meters
%
% Procedure:
%       1. Compute the composition througth time matrix - 'composition_ts', where species are columns and samples (time)
%       are rows
%       2. For every entry that has sufficient abundance, calculate log_dist, CE, N20 and ND. Comments:
%           - ND is a matrix of samples by distance bins
%           - for all others - create a matrix of samples by species and then allocate to bins
%
% This version performs the computation for only the central part of the metacommunity landscape, with an a-priori
% unknown number of species. It  exports the raw values of the
% statistics - log_N20 and mean_log_dist
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [log_N20s, log_NNDs, abd_bin_centers, samps_per_bin] = Stats_regime_meta_export(com_samp, abd_bin_edges, min_abd, max_abd, prealloc, Tot_L, Lx, Ly, filename)
% This function computes the statistics for a parameter regime

% % Example input:
% abd_bin_edges = [1 5:2:51 55:5:100 110:10:200 250:50:600 700:100:1000 1250 1500 5500];
% min_abd = 5;
% max_abd = 1500;
% Tot_L = 600;
% Lx = 400;
% Ly = 200;
% prealloc = 10000;
% filename = 'Stats_regime2_test.mat';

% Precalculations:
abd_bin_centers = abd_bin_edges(1:end-1) + (diff(abd_bin_edges)-1)/2; %centers of abundance bins
abd_bin_edges_use = abd_bin_edges - 0.01;
abd_bins = length(abd_bin_edges) - 1;

% Preallocate: 
log_N20s = nan(abd_bins, prealloc); %preallocate neighborhood density at one short range
log_NNDs = nan(abd_bins, prealloc); %preallocate log nearest neighbor distance

counters = ones(abd_bins,1); %how many observations in each bin?

ly = size(com_samp,3); %samples (years) in the data
ly100 = round(ly/100); %one percent of the total calculations to do. Calculated to show progress.

%run over all species-year combinations:
disp(['total samples (years): ' num2str(ly)])
for ii = 1:ly
    
    %Truncate the data:
    
    samp_now = com_samp(:,:,ii);
    outside = (samp_now(:,1) < (Tot_L-Lx)/2) | (samp_now(:,1) > (Tot_L+Lx)/2) | (samp_now(:,2) < (Tot_L-Ly)/2) | (samp_now(:,2) > (Tot_L+Ly)/2);
    samp_now(outside,:) = [];
    
    % Find all unique species in focal area:
    u = unique(samp_now(:,3));
    
    % Some initial calculations for this species:
    for ss = 1:length(u)
        
        % check species
        sp_inds = samp_now(:,3) == u(ss);
        this_abd = sum(sp_inds);
        
        if (this_abd >= min_abd) && (this_abd <= max_abd) %if it's within abundance range - analyze it
            points = samp_now(sp_inds,1:2);
            this_bin = find(abd_bin_edges_use-this_abd>0,1)-1;

            % Compute nearest neighbor distances:

            xdists = abs(points(:,1)-points(:,1)'); %distance between every pair of points on a torus along x axis
            xdists = min(xdists,Lx-xdists); %torus!
            ydists = abs(points(:,2)-points(:,2)'); %distance between every pair of points on a torus along y axis
            ydists = min(ydists,Ly-ydists); %torus!
            dd = sqrt(xdists.^2+ydists.^2); %distances of first 'num_samples' individuals to all 

            %replace diagonal with Inf (so that your distance to yourseld does not get counted in any distance bin):
            dd(1:(this_abd+1):(this_abd^2)) = Inf;

            %find minimum in every row:
            min_dists = min(dd,[],2);

            %assign to preallocated arrays:
            log_NNDs(this_bin, counters(this_bin)) = log(mean(min_dists));
            tmp = mean(sum(dd<20));
            
            if tmp == 0 %Correct infinite values when no neighbors are observed:
                log_N20s(this_bin, counters(this_bin)) = log(1/this_abd);
            else
                log_N20s(this_bin, counters(this_bin)) = log(tmp);
            end
            
            counters(this_bin) = counters(this_bin) + 1;
        end
    end
    
    if rem(ii,ly100) == 0
        disp(['total samples: ' num2str(ly) ', done: ' num2str(ii)])
    end
end

samps_per_bin = counters - 1;

log_N20s(:,max(counters):end) = [];
log_NNDs(:,max(counters):end) = [];


if ~isnan(filename) %save if required to do so
    save(filename, 'log_N20s', 'log_NNDs', 'abd_bin_centers','samps_per_bin')
end

end
