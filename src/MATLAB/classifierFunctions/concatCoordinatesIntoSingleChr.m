function [tnCoordinatesConcat,geneStartCoordinatesConcat,geneEndCoordinatesConcat] = concatCoordinatesIntoSingleChr(T,tnCoordinates,geneCount,lengthPreviousChr)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    
    tnCoordinatesConcat = [];
    
    geneChromosomeLoc = zeros(geneCount,1);
    
    for ii = 1:geneCount
        geneChromosomeLoc(ii) = roman2num(T.Chromosome{ii});
    end
    geneChromosomeLoc(isnan(geneChromosomeLoc)) = 17; %Mitochondrial DNA (denoted as 'mito') is returned as NaN by roman2num. All this DNA is counted as chromosome 17
    geneStartCoordinatesConcat = T.StartLocation;
    geneEndCoordinatesConcat = T.EndLocation;
    
    for ii=1:17
        %NOTE: tnCoordinatesConcat includes double locations. Presumably where
        %transposons hit either direction
        tnCoordinatesConcat= [tnCoordinatesConcat; tnCoordinates(tnCoordinates(:,1)==ii,2)+lengthPreviousChr(ii)]; %Translate all transposon coordinates to a total bp. 
        geneStartCoordinatesConcat(geneChromosomeLoc==ii) = geneStartCoordinatesConcat(geneChromosomeLoc==ii)+lengthPreviousChr(ii); %Translate all gene start coordinates to a total bp. 
        geneEndCoordinatesConcat(geneChromosomeLoc==ii) = geneEndCoordinatesConcat(geneChromosomeLoc==ii)+lengthPreviousChr(ii); %Translate all gene end coordinates to a total bp. 
    end

end