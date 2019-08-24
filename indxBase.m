% This should be my new base script
fetch_allSubj;

try
% set paths
addpath(genpath('Z:\Luca\functions\wave_clus-master')); % waveclus 3.0
addpath('Z:\Luca\functions'); % my functions
addpath('Z:\Luca\functions\Neuralynx_19012019'); % Neuralynx (the commons folder function doesnt work on my PC)
addpath(genpath('Z:\Luca\TREBER\Scripts'))

catch
% for server
addpath(genpath('/media/ldk898/rds-share/Luca/functions/wave_clus-master')); % waveclus 3.0
addpath('/media/ldk89h8/rds-share/Luca/functions'); % my functions
addpath('/media/ldk898/rds-share/Luca/functions/Neuralynx_19012019'); % Neuralynx (the commons folder function doesnt work on my PC)
addpath(genpath('/media/ldk898/rds-share/Luca/TREBER/Scripts'))
end

%% detect and cluster spikes (manually and automatically) using fndIndx3.m

%  Here I use cross-correlation to get rid of some artefacts
subjID = 'P04_S3';
subjID = 'P0#_S#';
disp(subjID);

xc_postClust(subjID); % does the crosscorrelation
quickfix = xc_postClust2(subjID); % rejects cluster based on crosscorrelation
if ~isempty(quickfix); open quickfix; end
rename_clNames(subjID); % fixes possible label errors

[encSpiketimes_cueLocked, encSpiketimes_respLocked, encSpiketimes_preTrial, encSpiketimes_stimLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID); % put spiketimes into the table template, output to visually check
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, encSpiketimes_preTrial, encSpiketimes_stimLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID);

% extract results for each participant in a way that works for SPSS
extractMeanRSA

% generate a mean RSA for all comparisons for hits as well as for misses that can be visualized
Average4RSA(allSubj)

% visualize all single units
visuSU

% SME
% script that runs function "SME" for all subjects (allSubj)
% then it runs "extractSME". This function summarizes the output from "SME" for all subjects (allSubj)
runSME

%% do do
% fix edge problem with retrieval in visuSU
% run extractMeanRSA for SPSS
% run visuSU to visualize all single units


%%
clear;clc
fetch_allSubj;
for sbj = 1 : size(allSubj, 1)
subjID = allSubj{sbj};
disp(subjID);
[encSpiketimes_cueLocked, encSpiketimes_respLocked, encSpiketimes_preTrial, encSpiketimes_stimLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD] = RSA(subjID); % put spiketimes into the table template, output to visually check
RSA2(encSpiketimes_cueLocked, encSpiketimes_respLocked, encSpiketimes_preTrial, encSpiketimes_stimLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, clNames, mCD, subjID);
end


%% kapa640 (desktop folder)
blockwise_diffScore_cond
drawBoxplot