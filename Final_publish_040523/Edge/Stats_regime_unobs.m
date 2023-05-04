%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute multiple statistics for the samples of a parameter regime
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Michael Kalyuzhny
% Date created: 29/12/2020
% Date last modified: 03/05/2023
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
% unknown number of species. The distances and densities are calculated with respect to ALL points. Here distances and
% densities do not use a torus BC
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mean_log_dist, log_mean_dist, CE, log_N20, N20, RND20, log_ND, ND, PCF, abd_bin_centers, dist_bin_centers] = Stats_regime_unobs(com_samp, abd_bin_edges, dist_bin_edges, min_abd, max_abd, Tot_L, Lx, Ly, filename)
% This function computes the statistics for a parameter regime

% Example input:
% abd_bin_edges = [1 5:2:51 55:5:100 110:10:200 250:50:600 700:100:1000 1250 1500 5500];
% dist_bin_edges = 0:10:200;
% min_abd = 5;
% max_abd = 1500;
% Tot_L = 1200;
% Lx = 1000;
% Ly = 500;
% filename = 'Stats_regime1_test.mat';

% Precalculations:
dist_bin_centers = dist_bin_edges(1:end-1) + diff(dist_bin_edges)/2; %centers of distance bins
abd_bin_centers = abd_bin_edges(1:end-1) + (diff(abd_bin_edges)-1)/2; %centers of abundance bins
abd_bin_edges_use = abd_bin_edges - 0.01;
abd_bins = length(abd_bin_edges) - 1;
dist_bins = length(dist_bin_edges) - 1;
con_area = pi*dist_bin_edges.^2; %the area of annuli with radii 'dist_bin_edges'
one_over_prop_area = (Lx*Ly)./diff(con_area); %what is the proportion of area (out of the total) in each annulus?
one_over_prop_area20 = (Lx*Ly)./(pi*400);

% Preallocate: 
ND_prealloc = zeros(abd_bins, dist_bins); %preallocate neighborhood density at a longer scale
N20_prealloc = zeros(abd_bins, 1); %preallocate neighborhood density at one short range
RND20_prealloc = zeros(abd_bins, 1); %preallocate rel. neighborhood density at one short range
log_ND_prealloc = zeros(abd_bins, dist_bins); %preallocate neighborhood density at a longer scale
log_N20_prealloc = zeros(abd_bins, 1); %preallocate neighborhood density at one short range
nnd_mean_prealloc = zeros(abd_bins, 1); %preallocate nearest neighbor distance
nnld_mean_prealloc = zeros(abd_bins, 1); %preallocate log nearest neighbor distance

PCF_prealloc = zeros(abd_bins, dist_bins);

counters = zeros(abd_bins,1); %how many observations in each bin?

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
            points_inside = samp_now(sp_inds,1:2);
            this_bin = find(abd_bin_edges_use-this_abd>0,1)-1;
            
            points_all_inds = com_samp(:,3,ii) == u(ss);
            points_all = com_samp(points_all_inds,:,ii);

            % Compute nearest neighbor distances:
% 
%             xdists = abs(points_inside(:,1)-points_inside(:,1)'); %distance between every pair of points on a torus along x axis
%             xdists = min(xdists,Lx-xdists); %torus!
%             ydists = abs(points_inside(:,2)-points_inside(:,2)'); %distance between every pair of points on a torus along y axis
%             ydists = min(ydists,Ly-ydists); %torus!
%             dd = sqrt(xdists.^2+ydists.^2); %distances of first 'num_samples' individuals to all 
%             
            dd = sqrt((points_inside(:,1) - points_all(:,1)').^2 + (points_inside(:,2) - points_all(:,2)').^2); %compute distances the simple way, without torus

            %replace zero distances with Inf (so that your distance to yourseld does not get counted in any distance bin):
            dd(dd == 0) = Inf;

            %find minimum in every row:
            min_dists = min(dd,[],2);

            % Calculate neighborhood density:
            dists_binned = nan(this_abd,dist_bins); %the thing to average on

            for pp = 1:this_abd
                [dists_binned(pp,:) , ~] = histcounts(dd((pp),:),dist_bin_edges); %for every indivisual, count neighbors within bins of distance
            end

            ND_now = mean(dists_binned,1); %mean neighbors in rings (over individuals)
            log_ND_now = log(ND_now);
            
            %assign to preallocated arrays:
            nnd_mean_prealloc(this_bin) = nnd_mean_prealloc(this_bin) + mean(min_dists);
            nnld_mean_prealloc(this_bin) = nnld_mean_prealloc(this_bin) + log(mean(min_dists));
            tmp = mean(sum(dd<20,2));
            N20_prealloc(this_bin) = N20_prealloc(this_bin) + tmp;
            RND20_prealloc(this_bin) = RND20_prealloc(this_bin) + tmp*one_over_prop_area20/this_abd;
            ND_prealloc(this_bin,:) = ND_prealloc(this_bin,:) + ND_now;
            PCF_prealloc(this_bin,:) = ND_now.*one_over_prop_area/this_abd;
            
            if tmp == 0 %Correct infinite values when no neighbors are observed:
                log_N20_prealloc(this_bin) = log_N20_prealloc(this_bin) + log(1/this_abd);
            else
                log_N20_prealloc(this_bin) = log_N20_prealloc(this_bin) + log(tmp);
            end
            
            log_ND_now(~isfinite(log_ND_now)) = log(1/this_abd);
            log_ND_prealloc(this_bin,:) = log_ND_prealloc(this_bin,:) + log_ND_now;

            counters(this_bin) = counters(this_bin) + 1;
        end
    end
    
    if rem(ii,ly100) == 0
        disp(['total samples: ' num2str(ly) ', done: ' num2str(ii)])
    end
end

%Final calculations:
mean_log_dist = nnld_mean_prealloc./counters;
log_mean_dist = nnd_mean_prealloc./counters;
CE = (nnd_mean_prealloc./counters).*2.*sqrt(abd_bin_centers'/(Lx*Ly^2));
N20 = (N20_prealloc)./counters;
ND = ND_prealloc./(counters);
log_N20 = log_N20_prealloc./counters;
log_ND = log_ND_prealloc./counters;
PCF = (PCF_prealloc)./counters;
RND20 = (RND20_prealloc)./counters;

samps_per_bin = counters;

if ~isnan(filename) %save if required to do so
    save(filename, 'mean_log_dist', 'PCF', 'log_mean_dist', 'CE', 'N20', 'RND20', 'log_N20', 'ND', 'log_ND', 'abd_bin_centers', 'dist_bin_centers','samps_per_bin')
end

end
