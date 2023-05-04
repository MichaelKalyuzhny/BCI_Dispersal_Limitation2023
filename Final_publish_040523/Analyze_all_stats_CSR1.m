%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform the analysis of statistics with the Random null
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Michael Kalyuzhny
% Date created: 10/04/2023
% Date last modified: 27/04/2023
%
% Computes the statistics with the Random null (CSR) as reference and adds them to the table
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data and prepare results table:

ittr = 2000;

load('species_data3.mat','sp_dat')

all_stats_names = {'R20l_CSR', 'R20l_CSR_pval', 'ED_mean_CSR', 'ED_mean_CSR_pval'}; %names of statistics to add
stats = array2table(nan(size(sp_dat,1),length(all_stats_names)),"VariableNames",all_stats_names);

%% Compute raw stats:

load('Results/spatial_distributions.mat', 'point_2015')

num_sp = size(sp_dat,1);
Lx = 1000; %size of BCI
Ly = 500;

% Empty arrays of results:
R20_raw = nan(num_sp,1);
NND_mean_raw = nan(num_sp,1);

for sp = 1:num_sp
    % data:
    points = point_2015{sp};
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

    % Compute raw statistics: 
    R20_raw(sp) = log(mean(sum(dd<20))); %number of neighbors within 20 meters
    NND_mean_raw(sp) = log(mean(min_dists)); % mean nearest neighbor distance
    
    % Correct Rs that are zeros: assume 1 individual had a neighbor (in real data minimum is 2):
    if R20_raw(sp)==-Inf
        R20_raw(sp) = log(1/this_abd);
    end

end

%% Generate CSR null and compute final stats:

for sp = 1:num_sp
    
    sp
    
    R20_res = nan(1,ittr);
    NND_mean_res = nan(1,ittr);
    
    for ii = 1:ittr
        [sp ii]
        
        pop_generate = max(2,poissrnd(sp_dat.abd_adults_2015(sp)));
        
        % data:
        points = rand(pop_generate,2).*[Lx Ly];
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

        % Compute raw statistics: 
        R20_res(ii) = log(mean(sum(dd<20))); %number of neighbors within 20 meters
        NND_mean_res(ii) = log(mean(min_dists)); % mean nearest neighbor distance

        % Correct Rs that are zeros: assume 1 individual had a neighbor (in real data minimum is 2):
        if R20_res(ii)==-Inf
            R20_res(ii) = log(1/this_abd);
        end
    end
    
    stats.ED_mean_CSR(sp) = NND_mean_raw(sp) - mean(NND_mean_res);
    stats.ED_mean_CSR_pval(sp) = dist_pval(NND_mean_res, NND_mean_raw(sp));
    stats.R20l_CSR(sp) = R20_raw(sp) - mean(R20_res);
    stats.R20l_CSR_pval(sp) = dist_pval(R20_res, R20_raw(sp));

end

%% export extended dataset:

sp_dat = [sp_dat stats];      
save('species_data4.mat','sp_dat') %THIS INTERMEDIATE FILE IS PROVIDED