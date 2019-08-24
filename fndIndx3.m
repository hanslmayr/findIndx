% ideas:
% segmentierung von spike detection sollte auf varianz basieren!

try
p2d = cd;
p2d(end+1)='\';
[encSpiketimes_cueLocked, encSpiketimes_respLocked, retSpiketimes_cueLocked, retSpiketimes_respLocked, hitsIdx, missIdx, allTrials, sessionNum, retTrigger, encTrigger]=loadLogs(p2d);
dts.start = encTrigger(1);
dts.end = max(max(retTrigger))+1; % remember that retTrigger is not linearly increasing
dts.dt = 0.001:0.001:dts.end; % create a linerally spaced vector from the beginning to the end of the experiment in steps of 1ms

catch
    disp('Using default dt values');
    dts.start = 1;
    dts.end = 4000;
    dts.dt = 3;
end

% positive detection
par = set_parameters();
par.detection = 'pos';
Get_spikes('all', 'parallel', true, 'par', par);
movefile('*.mat', 'posDetect'); % moves all .mat files that include the just detected spikes into the folder "posDetect"

% % negative detection
par.detection = 'neg';
Get_spikes('all', 'parallel', true, 'par', par);
movefile('*.mat', 'negDetect'); % moves all .mat files that include the just detected spikes into the folder "negDetect"

% positive clustering
cd posDetect
Do_clustering('all', 'parallel', true, 'make_plots', true);
% Do_clustering('CSC_antHippR1_spikes.mat', 'parallel', false,'make_plots', true);

% % negative clustering
cd ..\negDetect
Do_clustering('all', 'parallel', true, 'make_plots', true);

clear it
% positive manual inspection
cd ..\posDetect
if ~exist('it'); it=1; end
manuSorting(dts, it);
clear it

%     checkchannelLDK(it)
%     mkSpiketimes(dts ,it);

%     pnClus = dir('times_*.mat');
%     wave_clus(pnClus(it).name);
%     load('it.mat'); disp(it);

% negative manual inspection
clear it; cd ..\negDetect;
if ~exist('it'); it=1; end
manuSorting(dts, it);
clear it