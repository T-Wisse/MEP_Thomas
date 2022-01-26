function [intergenicLength20kbPerGene,intergenicTn20kbPerGene,intergenic_tndensity20kb] = getIntergenicTnDensity(geneLength,chromosomeEndPos,tnCoordinates,coordinates,chromosomePerGene)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    wl = 20000; % define window length
    
    tt=0;
    uu=0;
    intergenicLength20kbPerGene = zeros(length(geneLength),2);
    
    for kk=1:length(chromosomeEndPos) %For each chromosome
        myField = strcat('ch',num2str(kk)); % create variable with all bp positions for every chromosome
        genePos.(myField) = zeros(chromosomeEndPos(kk),1); % initialize field as number of bp in the chromosome
        intergenicTnPos.(myField)= zeros(chromosomeEndPos(kk),1); % initialize field as number of bp in the chromosome
        tnPos.(myField) = zeros(chromosomeEndPos(kk),1);
    
        for ii=tt+1:tt+sum(chromosomePerGene==kk) % For each gene in the corresponding chromosome
            genePos.(myField)(coordinates(ii,1):coordinates(ii,2)) = 1; %put 1 at the region (bp's) in the chromosome that corresponds to a gene
            tt = ii;
        end
        
        temp = tnCoordinates(tnCoordinates(:,1) == kk,2); % Select the transposon coordinates for a gene. 
        tnPos.(myField)(temp) = 1; %gives a 1 for every bp position that has a transposon insertion
            
        intergenicTnPos.(myField)=tnPos.(myField);
        genePos.(myField)=logical(genePos.(myField));
        intergenicTnPos.(myField)(genePos.(myField)) = 0; %why?
        
        
        for jj=uu+1:uu+sum(chromosomePerGene==kk) % every gene of each chromosome
            if coordinates(jj,1)-wl<1
                intergenicLength20kbPerGene(jj,1) = sum(1-(genePos.(myField)(1:coordinates(jj,1)))); % calculate intergenic length before gen in 20kb window
                intergenicTn20kbPerGene(jj,1) = sum(intergenicTnPos.(myField)(1:coordinates(jj,1))); % calculates number of tn in this region
            else
                intergenicLength20kbPerGene(jj,1) = sum(1-(genePos.(myField)(coordinates(jj,1)-wl:coordinates(jj,1)))); % calculate intergenic length before gen in 20kb window
                intergenicTn20kbPerGene(jj,1) = sum(intergenicTnPos.(myField)(coordinates(jj,1)-wl:coordinates(jj,1)));
            end
            
            if coordinates(jj,2)+wl>chromosomeEndPos(kk)
                 intergenicLength20kbPerGene(jj,2) = sum(1-(genePos.(myField)((coordinates(jj,2):chromosomeEndPos(kk))))); % calculate intergenic length before gen in 20kb window
                 intergenicTn20kbPerGene(jj,2) = sum(intergenicTnPos.(myField)(coordinates(jj,2):chromosomeEndPos(kk)));
            else
                intergenicLength20kbPerGene(jj,2) = sum(1-(genePos.(myField)(coordinates(jj,2):coordinates(jj,2)+wl))); % calculate intergenic length before gen in 20kb window
                intergenicTn20kbPerGene(jj,2) = sum(intergenicTnPos.(myField)(coordinates(jj,2):coordinates(jj,2)+wl));
            end
            uu=jj;
        end
    end
    
    intergenic_tndensity20kb = intergenicTn20kbPerGene./intergenicLength20kbPerGene;

end