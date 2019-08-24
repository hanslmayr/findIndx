RSAmean_h_all = [];
RSAmean_c_all = [];
RSAmean_o_all = [];

for i = 1 : size(allSubj,1)
    subjID = allSubj{i};
    clear RSAmean_h RSAmean_c RSAmean_o
    
    try
        cd Z:/Luca/data
    catch
        cd /media/ldk898/rds-share/Luca/data
    end
    mSubject = subjID(1:end-3);
    mSession = subjID(end-1:end);
    
    % if the session name is called 1b then this line prevents an error during cd
    mSubject(regexp(mSubject,'_')) = [];
    if isempty(regexp(mSession,'S', 'ONCE'))
        mSession = ['S', mSession];
    end
    
    cd(mSubject)
    cd(mSession)
    abc = dir;
    cd(abc(3).name)
    cd advancedAnalysis\
    cd RSA
    abc = dir('RSAmean_*');
    load(abc.name)
    
    %%
    if ~isempty(RSAmean_h) % in case there are not enough hippocampal single units
        if size(RSAmean_h,1) == 3 % in case there are not enough misses
            RSAmean_h(4:6, 1:8) = nan;
            RSAmean_h(4:6, 9  ) = [1;2;3];
            RSAmean_h(4:6, 10 ) = [0;0;0]; % redundant
        end
        RSAmean_h(1:size(RSAmean_h,1),end+1) = i;
        RSAmean_h_all = [RSAmean_h_all; RSAmean_h];
    end
    
    if ~isempty(RSAmean_c)
        if size(RSAmean_c,2) == 3
            RSAmean_c(4:6, 1:8) = nan;
            RSAmean_c(4:6, 9  ) = [1;2;3];
            RSAmean_c(4:6, 10 ) = [0;0;0]; % redundant
        end
        RSAmean_c(1:size(RSAmean_c,1),end+1) = i;
        RSAmean_c_all = [RSAmean_c_all; RSAmean_c];
    end
    
    if ~isempty(RSAmean_o)
        if size(RSAmean_o,2) == 3
            RSAmean_o(4:6, 1:8) = nan;
            RSAmean_o(4:6, 9  ) = [1;2;3];
            RSAmean_o(4:6, 10 ) = [0;0;0]; % redundant
        end
        RSAmean_o(1:size(RSAmean_o,1),end+1) = i;
        RSAmean_o_all = [RSAmean_o_all; RSAmean_o];
    end
end

clearvars -except RSAmean_h_all RSAmean_c_all RSAmean_o_all