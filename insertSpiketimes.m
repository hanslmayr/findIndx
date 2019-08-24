%% insertSpiketimes
% trigger: encTrigger
% spiketimescluster=spiketimescluster
% locking: cuelocked = 1; stimlocked = 2; resplocked = 3
% timeWindow: Timewindow before and after the trigger that is considered

function spiketimes_clusterTrial=insertSpiketimes(trigger, spikeTimesCluster, locking, timeWindow)
spiketimes_clusterTrial={};
for ix=1:size(trigger,1) % number of trials  
    spiketimes_clusterTrial{1,ix} = spikeTimesCluster(find(spikeTimesCluster >= trigger(ix, locking) + timeWindow(1) ...
        & spikeTimesCluster <= trigger(ix,locking) + timeWindow(2) )) - trigger(ix,locking);    
end
end