function [promotorData] = promotorData(posConcat,geneCoordinatesConcat,geneCount,geneSenseBinary)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
promotorData = zeros(geneCount,1);
    for kk=1:geneCount %For each gene, find the number of reads and transposons in the 100 bp upstream (promotor region)
        if geneSenseBinary(kk)==1    
            promotorData(kk)=sum(posConcat(geneCoordinatesConcat(kk,1)-100:geneCoordinatesConcat(kk,1)));
        else
            promotorData(kk)=sum(posConcat(geneCoordinatesConcat(kk,2):geneCoordinatesConcat(kk,2)+100));
        end
    end
end