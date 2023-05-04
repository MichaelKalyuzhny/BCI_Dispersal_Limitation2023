  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute statistics for a null model for a set of species
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Michael Kalyuzhny
% Date created: 17/01/2021
% Date last modified: 05/04/2023
%
% Input: 
%   1. Lx, Ly - sizes of landscape
%   2. dist_bind_edges - the left edge of the bins of the Relative Neighborhood density R(r). The last number is the
%   right edge
%   3. species_analyze - the species to run the analysis for
%   4. filenames_nulls - the filenames of the null model runs 
%   5. fileneme_save - filename to save the results
%   6. truncate_space - should space be truncated to the size of BCI?
%
% Output:
%    R(20), R75_125, ED_mean, ED_median and their p values; R(r) and it's 95% simulation envelope; the centers of distance bins
% 
% Procedure:
%   1. Compute the set of statistics for every one of the observed species (each represented by a set of points)
%   2. Import the Null data for that species and compute the distributions of these statistics (including simulation
%   envelopes) for a species with a similar binning
%   3. Compute the normalized statistics

% Version 4 computes global envelopes by outputing the simulation results to an R file, which employs the GET package.

% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R20l, R20l_pval, R75_125l, R75_125l_pval, R, R_sim_envelope, ED_mean, ED_mean_pval, ED_median, ED_median_pval, comparable_samps_null, dist_bin_centers, R_global_envelope] = Analyze_nulls4(species_analyze, filenames_nulls, Lx, Ly, dist_bin_edges, fileneme_save, truncate_space)

%% Precalculations:

% Temporart version of the input:
load('species_data1.mat','sp_dat')

% Example input:
% species_analyze = 66;%find(~isnan(sp_dat{:,'HML2008AlphaFitted'})); %these are the species to run the simulation for
% filenames_nulls = [ repmat('../Results/null1_2008_', length(species_analyze), 1), char(sp_dat{species_analyze,'SpeciesCode'}) repmat('.mat', length(species_analyze), 1) ];
% Lx = 1000; %length (in meters) of the X dimension
% Ly = 500;
% dist_bin_edges = 10:10:200; %'dist_bin_edges' array of the left edge of every bin. last value is the final edge.

num_sp = length(species_analyze);
num_dist_bins = length(dist_bin_edges) - 1;
dist_bin_centers = dist_bin_edges(1:end-1) + diff(dist_bin_edges)/2; %centers of distance bins

%% Compute observed (raw) statistics:

load('Results/spatial_distributions.mat', 'point_2015')

% Empty arrays of results:
R20_raw = nan(num_sp,1);
R75_125_raw = nan(num_sp,1);
R_raw = nan(num_sp, num_dist_bins);
NND_mean_raw = nan(num_sp,1);
NND_median_raw = nan(num_sp,1);

for sp = 1:num_sp
    % data:
    points = point_2015{species_analyze(sp)};
    this_abd = size(points,1);
    
    % generate distance matrix:
    xdists = abs(points(:,1)-points(:,1)'); %distance between every pair of points on a torus along x axis
    xdists = min(xdists,Lx-xdists); %torus!
    ydists = abs(points(:,2)-points(:,2)'); %distance between every pair of points on a torus along y axis
    ydists = min(ydists,Ly-ydists); %torus!
    dd = sqrt(xdists.^2+ydists.^2); %distance matrix

    %replace diagonal with Inf (so that your distance to yourseld does not get counted in any distance bin):
    dd(1:(this_abd+1):(this_abd^2)) = Inf;

    %find minimum in every row:
    min_dists = min(dd,[],2); %This is the vector of nearest neighbor distances
    
    % Calculate neighborhood density:
    dists_binned = nan(this_abd,num_dist_bins); %the thing to average on

    for pp = 1:this_abd
        [dists_binned(pp,:) , ~] = histcounts(dd((pp),:),dist_bin_edges); %for every indivisual, count neighbors within bins of distance
    end
    
    % Compute raw statistics: 
    R20_raw(sp) = log(mean(sum(dd<20))); %number of neighbors within 20 meters
    R75_125_raw(sp) = log(mean(sum( (dd>75)&(dd<125) )));
    R_raw(sp,:) = log(mean(dists_binned,1)); % neighborhood density ~ distance
    NND_mean_raw(sp) = log(mean(min_dists)); % mean nearest neighbor distance
    NND_median_raw(sp) = log(median(min_dists)); % median nearest neighbor distance
    
    % Correct Rs that are zeros: assume 1 individual had a neighbor (in real data minimum is 2):
    if R20_raw(sp)==-Inf
        R20_raw(sp) = log(1/this_abd);
    end
    R_raw(sp,~isfinite(R_raw(sp,:))) = log(1/this_abd);
    
    if R75_125_raw(sp)==-Inf
        R75_125_raw(sp) = log(1/this_abd);
    end
end

%% For every species: analyze it's Null and then compute relativized statistics:

% Empty arrays of results:
R20l = nan(num_sp,1);
R20l_pval = nan(num_sp,1);
R75_125l = nan(num_sp,1);
R75_125l_pval = nan(num_sp,1);
R = nan(num_sp, num_dist_bins);
R_sim_envelope = nan(2, num_dist_bins, num_sp); %first dimension: lower/upper, second: distance bin, 3rd: species index
R_global_envelope = nan(2, num_dist_bins, num_sp); %first dimension: lower/upper, second: distance bin, 3rd: species index
ED_mean = nan(num_sp,1);
ED_mean_pval = nan(num_sp,1);
ED_median = nan(num_sp,1);
ED_median_pval = nan(num_sp,1);
comparable_samps_null = nan(num_sp,1);
R_ress = cell(num_sp,1);

parfor sp = 1:num_sp %run over all species
    sp
    this_abd = sp_dat{species_analyze(sp),'abd_adults_2015'}; %the abundance of the species to analyze
    
    % Load null and remove edges:
    l = load(filenames_nulls(sp,:)); %load simulation results for this species
    com_samp = l.com_samp;
    
    if truncate_space
        %transform trees outside of bounds to NAN:
        for ss=1:size(com_samp,3)
            outside = (com_samp(:,1,ss) < 100) | (com_samp(:,1,ss) > 1100) | (com_samp(:,2,ss) < 350) | (com_samp(:,2,ss) > 850); 
            com_samp(outside,:,ss) = nan;
        end
    end
    
    % Find comparable species:
    bin_width = 0.3*(this_abd)^0.75;

    composition_ts = sum_pop(com_samp); %create an abundance-time matrix

    %find when species fall within the range:
    [year,spp] = find((composition_ts > round(this_abd - bin_width)) & (composition_ts < round(this_abd + bin_width))); %pairs of year, sp where there are observations
    comparable_samps_null(sp) = length(year); 
    
    % Empty arrays of resamples:
    R20_res = nan(comparable_samps_null(sp),1);
    R75_125_res = nan(comparable_samps_null(sp),1);
    R_res = nan(comparable_samps_null(sp), num_dist_bins);
    NND_mean_res = nan(comparable_samps_null(sp),1);
    NND_median_res = nan(comparable_samps_null(sp),1);
    
    %run over comparable samples in the null model:
    for ii = 1:comparable_samps_null(sp)
        sp_inds = com_samp(:,3,year(ii)) == spp(ii); %the indexes of the sample in the null data
        points = com_samp(sp_inds,1:2,year(ii));
        this_abd = size(points,1); %since there is some variability in the abundane, it should be 

        % generate distance matrix:
        xdists = abs(points(:,1)-points(:,1)'); %distance between every pair of points on a torus along x axis
        xdists = min(xdists,Lx-xdists); %torus!
        ydists = abs(points(:,2)-points(:,2)'); %distance between every pair of points on a torus along y axis
        ydists = min(ydists,Ly-ydists); %torus!
        dd = sqrt(xdists.^2+ydists.^2); %distance matrix

        %replace diagonal with Inf (so that your distance to yourseld does not get counted in any distance bin):
        dd(1:(this_abd+1):(this_abd^2)) = Inf;

        %find minimum in every row:
        min_dists = min(dd,[],2); %This is the vector of nearest neighbor distances

        % Calculate neighborhood density:
        dists_binned = nan(this_abd,num_dist_bins); %the thing to average on

        for pp = 1:this_abd
            [dists_binned(pp,:) , ~] = histcounts(dd((pp),:),dist_bin_edges); %for every indivisual, count neighbors within bins of distance
        end
        
        R20_res(ii) = log(mean(sum(dd<20))); %number of neighbors within 20 meters
        R75_125_res(ii) = log(mean(sum( (dd>75)&(dd<125) )));
        R_res(ii,:) = log(mean(dists_binned,1)); %
        NND_mean_res(ii) = log(mean(min_dists)); % mean nearest neighbor distance
        NND_median_res(ii) = log(median(min_dists)); % median nearest neighbor distance
        
        if R20_res(ii)==-Inf
            R20_res(ii) = log(1/this_abd);
        end
        R_res(ii,~isfinite(R_res(ii,:))) = log(1/this_abd);
        
        if R75_125_res(ii)==-Inf
            R75_125_res(ii) = log(1/this_abd);
        end
    end
    
    % Finally, compute all statistics for this species:
    R75_125l(sp) = R75_125_raw(sp) - mean(R75_125_res);
    R75_125l_pval(sp) = dist_pval(R75_125_res, R75_125_raw(sp));
    R20l(sp) = R20_raw(sp) - mean(R20_res);
    R20l_pval(sp) = dist_pval(R20_res, R20_raw(sp));
    R(sp,:) = R_raw(sp,:) - mean(R_res);
    ED_mean(sp) = NND_mean_raw(sp) - mean(NND_mean_res);
    ED_mean_pval(sp) = dist_pval(NND_mean_res, NND_mean_raw(sp));
    ED_median(sp) = NND_median_raw(sp) - median(NND_median_res);
    ED_median_pval(sp) = dist_pval(NND_median_res, NND_median_raw(sp));

    %Compute dimulation envelopes:
    s = sort(R_res);
    bounds = [round(0.025*comparable_samps_null(sp)) round(0.975*comparable_samps_null(sp))]; %the 2.5 and 97.5 percentiles of the resampled densitites
    R_sim_envelope(:,:,sp) = [s(bounds(1),:) - mean(R_res) ; s(bounds(2),:) - mean(R_res)]; %the simulation envelope
    
    % Output centralized simulation results (to be used for constructing global envelopes)
    R_ress{sp} = R_res - mean(R_res);
end

% Compute global envelopes:

for sp = 1:num_sp %run over all species
    
    % Export:
    R_res2 = R_ress{sp};
    R_export = R(sp,:);
    writematrix(R_res2,'r_res.csv')
    writematrix(R_export,'r.csv')

    % Run R script:

    % First, make sure that R can be called from the command line. In windoes, this includes ensuring that the folder where R is installed is in the 'Path', so it can be ran in command line. Here is a
    % good link for how to do that:
    % https://www.java.com/en/download/help/path.html
    % Also, check the R script and install the GET package.
    
    !Rscript Export_ERL_env.R 

    % Import:
    
    inp = readmatrix('export_rel.csv');
    R_global_envelope(:,:,sp) = inp(:,2:3)'; %order of upper and lower is inverted here
end

% Save results:

if ~isempty(fileneme_save)
    save(fileneme_save, 'R20l', 'R20l_pval', 'R75_125l', 'R75_125l_pval', 'R', 'R_sim_envelope', 'ED_mean', 'ED_mean_pval', 'ED_median', 'ED_median_pval', 'comparable_samps_null', 'R_global_envelope')
end

end