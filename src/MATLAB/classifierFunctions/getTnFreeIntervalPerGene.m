function [tnFreeInterval] = getTnFreeIntervalPerGene(tnCoordinatesConcat,geneStartCoordinatesConcat,geneEndCoordinatesConcat,geneCount)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    for ii=1:geneCount %Loop trough each gene. For each, find the longest interval between transposons in bp
        ww = tnCoordinatesConcat>=geneStartCoordinatesConcat(ii)&tnCoordinatesConcat<=geneEndCoordinatesConcat(ii);
        ww = tnCoordinatesConcat(ww);
        ww = [geneStartCoordinatesConcat(ii) ; ww ; geneEndCoordinatesConcat(ii)]; %Concat the transposon locations of gene ii with the start and end of the gene
        tnFreeInterval(ii,1) = max(diff(ww)); %determine the largest difference bertween previously concatenated locations to find longest region without a transposon in the gene 
    end

end