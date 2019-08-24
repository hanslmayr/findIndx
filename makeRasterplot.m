%% Function to make a Rasterplot
% allTrials: Number of Trials in that Session
% timeWindow: Timewindow before and after the trigger that is considered
% spiketimes: encSpiketimes_cueLocked for encCueLocked
% position: [2:3:11] for encCueLocked

function [n, dt]=makeRasterplot(timeWindow, allTrials, spiketimes, position)
subplot(13,3,position,'align');
hold on
dt = linspace(timeWindow(1),timeWindow(2),41);
n = [];
for ix = 1:allTrials
    x=spiketimes{1,ix};
    if isempty(x{1}) && ix~=allTrials
        continue % if that cluster has no spikes in that trial, continue with the next trial
    elseif isempty(x{1}) && ix==allTrials
        n(ix,:)=0; % otherwise n has too few rows and fr is wrong
    else
        x=cell2mat(x);
        xd = [x;x];
        y = ix*ones(1,length(x));
        y = [y-.5;y+.5];
        line(xd,y,'Color','r');
        [n(ix,:),~] = hist(x,dt);
    end
end
axis tight
xlim([timeWindow(1)+0.5 timeWindow(2)-0.5])
ylim([0 allTrials])
ylabel('Trial Number');
box on
end
