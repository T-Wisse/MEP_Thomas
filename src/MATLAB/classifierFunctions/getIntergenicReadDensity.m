function [intergenicRead20kbPerGene,intergenicReadDensity20kb] = getIntergenicReadDensity(geneLength,chromosomeEndPos,readDataChrCoords,coordinatesGenes,chromosomePerGene)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    wl = 20000; % define window length
    
    tt=0;
    uu=0;
    intergenicLength20kbPerGene = zeros(length(geneLength),2);
    intergenicRead20kbPerGene = zeros(length(geneLength),2);
    
    for kk=1:length(chromosomeEndPos) %For each chromosome
        chr = strcat('ch',num2str(kk)); % create variable with all bp positions for every chromosome
        genePos.(chr) = zeros(chromosomeEndPos(kk),1); % initialize field as number of bp in the chromosome
        intergenicReadPos.(chr)= zeros(chromosomeEndPos(kk),1); % initialize field as number of bp in the chromosome
        readCountPos.(chr) = zeros(chromosomeEndPos(kk),1);
    
        for ii=tt+1:tt+sum(chromosomePerGene==kk) % For each gene in the corresponding chromosome
            genePos.(chr)(coordinatesGenes(ii,1):coordinatesGenes(ii,2)) = 1; %put 1 at the region (bp's) in the chromosome that corresponds to a gene
            tt = ii;
        end
        
        temp = readDataChrCoords(readDataChrCoords(:,3) == kk,1); % Select the read coordinates for a gene. 
        readCountPos.(chr)(temp) = readDataChrCoords(readDataChrCoords(:,3)==kk,2); %gives the read count for every bp position per chromosome
        
        genePos.(chr)=logical(genePos.(chr));
        intergenicReadPos.(chr)=readCountPos.(chr);
        intergenicReadPos.(chr)(genePos.(chr)) = 0; %why?
        
        
        for jj=uu+1:uu+sum(chromosomePerGene==kk) % every gene of each chromosome
            %Left hand side window?
            if coordinatesGenes(jj,1)-wl<1 %If the window extends outside of the chromosome at the start
                intergenicLength20kbPerGene(jj,1) = sum(1-(genePos.(chr)(1:coordinatesGenes(jj,1)))); % calculate intergenic length before gen in 20kb window
                intergenicRead20kbPerGene(jj,1) = sum(intergenicReadPos.(chr)(1:coordinatesGenes(jj,1))); % calculates number of reads in this region
            else
                intergenicLength20kbPerGene(jj,1) = sum(1-(genePos.(chr)(coordinatesGenes(jj,1)-wl:coordinatesGenes(jj,1)))); % calculate intergenic length before gen in 20kb window
                intergenicRead20kbPerGene(jj,1) = sum(intergenicReadPos.(chr)(coordinatesGenes(jj,1)-wl:coordinatesGenes(jj,1)));
            end
            
            %Right hand side window?
            if coordinatesGenes(jj,2)+wl>chromosomeEndPos(kk) %If the window extends outside of the chromosome at the end
                 intergenicLength20kbPerGene(jj,2) = sum(1-(genePos.(chr)((coordinatesGenes(jj,2):chromosomeEndPos(kk))))); % calculate intergenic length before gen in 20kb window
                 intergenicRead20kbPerGene(jj,2) = sum(intergenicReadPos.(chr)(coordinatesGenes(jj,2):chromosomeEndPos(kk)));
            else
                intergenicLength20kbPerGene(jj,2) = sum(1-(genePos.(chr)(coordinatesGenes(jj,2):coordinatesGenes(jj,2)+wl))); % calculate intergenic length before gen in 20kb window
                intergenicRead20kbPerGene(jj,2) = sum(intergenicReadPos.(chr)(coordinatesGenes(jj,2):coordinatesGenes(jj,2)+wl));
            end
            uu=jj;
        end
    end
    
    intergenicReadDensity20kb = intergenicRead20kbPerGene./(1+intergenicLength20kbPerGene); %1st column = left hand side, 2nd column = right hand side? % ADDED 1 SO THAT IT CANT GIVE NAN. WHY IT IS 0 IS UNCLEAR TO ME....

end