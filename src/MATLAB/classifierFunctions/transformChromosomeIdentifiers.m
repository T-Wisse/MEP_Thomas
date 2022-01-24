function [data] = transformChromosomeIdentifiers(data)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    chromosomeIdentifiers = unique(data(:,1)); %Get the 17 unique chromosome identifiers used in the data
    chromosomeIdentifiers = chromosomeIdentifiers - 1132;
    data(:,1) = data(:,1)-1132;
    data(data(:,1)==92,1) = data((data(:,1)==92),1)-75;
    if unique(data(:,1)) ~= [1:17]
        disp('An issue occured with transformation of chromosome identifiers')
    end
end

% chromosomeIdentifiers = unique(data(:,1)); %Get the 17 unique chromosome identifiers used in the data
% chromosomeIdentifiers = chromosomeIdentifiers - 1132;
% data(:,1) = data(:,1)-1132;
% data(data(:,1)==92,1) = data((data(:,1)==92),1)-75;
% if unique(data(:,1)) ~= [1:17]
%     disp('An issue occured with transformation of chromosome identifiers')
% end