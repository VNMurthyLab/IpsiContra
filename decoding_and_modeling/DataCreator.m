%% DataCreator.m
% Function to take in the data in raw form and output one specific
% processed variety. Flags outlined in NowClassify.m tell you which type
% you are accessing, condensed info here:
%
% INPUTS
% Mice - List of Mice to use: 1 OB, 2 - 4 AON, 5 - 8 APC, 9 - 11 PPC.
% Side - Which odour presentation side. 1 - Both, 2 - Contra, 3 - Ipsi.
% Chosen_Odours - which odours to include. 1 - 15
% Bin_Size - The size of the window.
% Sniff_Align: Wether (1) or not (0) to shift 0 to first sniff.
% All the Data (so you only load it all in once)
% 
% OUTPUT
% The nicely processed data and labels.
% X - Neurons x Odours x Trials summed spike responds
% Y - Odours labels

function [X, Y] = DataCreator(Mice, Side, Chosen_Odours,Bin_Size, ...
    Sniff_Align, DataSets, SniffTime)

% Number of odours to include
NumOdours = 15;

% Extract data and sniff times for requested mice
TimeSeries = [];
Sniffing = [];
for i = 1:length(Mice)
    TimeSeries = [TimeSeries; DataSets{Mice(i)}];
    Sniffing = [Sniffing; squeeze(SniffTime(Mice(i), 1:size(DataSets{Mice(i)}, 1), :,Side))];
end

% Get only the requested side
TimeSeries = TimeSeries(:,:,Side);
Neurons = size(TimeSeries,1);       % tot num of neurons

% Here we jiggle the data into a more manageable one layer array (rather than cell and array)
TimeData = zeros(Neurons, NumOdours+1, 7, 1100);
for n = 1:Neurons
    for o = 1:(NumOdours+1)
        TimeData(n, o, :,:) = TimeSeries{n,o};
    end
end

% Now, starting from the time specified by Offset_Time and running
% for the time specified by Bin_Size create a single summed
% spike count
Summed_Initial = zeros(Neurons, NumOdours+1, 7);
start = 601;
for o = 1:NumOdours+1
    for n = 1:Neurons
        if Sniff_Align
            Sniff_Vector = Sniffing{n, o};
        else
            Sniff_Vector = zeros(1, 7);
        end
        for t = 1:7
            starter = round(start + Sniff_Vector(t)*100);
            if starter > 1100
                starter = 1100;
            end
            ender = starter+Bin_Size*100;
            if ender > 1100
                ender = 1100;
            end 
            Summed_Initial(n, o, t) = sum(TimeData(n,o,t,starter:ender));
        end
    end
end
X = Summed_Initial(:,1:NumOdours,:) - mean(Summed_Initial(:,16,:),3);
%X = permute(X(:,1:NumOdours,:), [2,3,1]);

Y = sort(Chosen_Odours);
% Now limit data to only trials corresponding to chosen odours
X = X(:,Y, :);
end
