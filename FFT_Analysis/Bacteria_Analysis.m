%% WARNING: This was designed using functions introduced in MATLAB 2016B. 
%  Results may vary if used with earlier MATLAB releases.

%% Initialize Environment
% clear all; % Do not clear all if reusing metadata table
close all;
clc;

%% Function Setup

% Day range for analysis 
range = 1:70;

% Target Directory for Output Plots
output_dir = '';
if((output_dir ~= '') & (~exist(output_dir,'dir')))
    mkdir(output_dir);
end

% Saving Plots Enabled (1) or Disabled (0)
save_plots = 0;

% Load CSV Data and Taxonomy Files
if ~exist('SRA_metadata','var')
    SRA_metadata = readtable('DataPreparation/MergedDataSet.csv');
end
if ~exist('taxonomy_table','var')
    taxonomy_table = readtable('taxatable.csv');
end

% For ease of analysis, "description_s" data was replaced with the following
% Person A = -100, Person B = -200, Person A stool = -101
StoolA_index = find(SRA_metadata.description_s == -100);
SalivaA_index = find(SRA_metadata.description_s == -101);

%% Loop with FFT Function Call

% Choose Sequence Numbers for Analysis
for i = 1:5
    % Initialize Data Vectors with NaN for Future Processing
    StoolA = ones(365,1)*NaN;
    SalivaA = ones(365,1)*NaN;
    
    % Grab Bacteria Data from Appropriate Table Column
    % (sequence data starts at column 4)
    bact = [SRA_metadata.collection_day_s SRA_metadata.(i+3)];
    % If sequence data exists for a given day, 
    % replace NaN with proper value
    for day = 1:365
        A = find(bact(StoolA_index,1) == day);
        if(A)
            StoolA(day) = bact(StoolA_index(A(1)),2);
        end
        A_s = find(bact(SalivaA_index,1) == day);
        if(A_s) 
            SalivaA(day) = bact(SalivaA_index(A_s(1)),2);
        end
    end
    
    % Select Out the Data for Proper Time Period
    StoolA = StoolA(range);
    SalivaA = SalivaA(range);

    % Define Output File Names Based on Taxonomy Data
    ouput_name = [num2str(i),'_',taxonomy_table.Family{i},'_',taxonomy_table.Genus{i}];
    Bacteria_FFT(ouput_name,StoolA,SalivaA,range,output_dir,save_plots);
end