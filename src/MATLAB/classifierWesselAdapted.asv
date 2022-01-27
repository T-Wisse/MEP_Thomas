%https://sites.google.com/site/satayusers/

clear all
close all

%ORDER OF FILES LOADED: BEDFILE, perInsertionsFile, perGeneFile, essentialGeneFile
% File locations: C:\Users\Thomas\Desktop\Uni\SATAYDATA_211210\SATAY_Data\SATAY_Analysis
% filebase = 'yTW001_4_merged_cleaned_forward_reads_trimmed.sorted.bam';
filebase = 'yLIC137_7_merged_cleaned_forward_reads_trimmed.sorted.bam';
% filebase = 'WT_merged-techrep-a_techrep-b_trimmed.sorted.bam';
%%renamed the WT .wig file to fit this naming scheme. should check if it is
%%from the exact same run. Also renamed yLIC bed file


%% Load satay data from transposonmapper & annotation file

% Should be able to do this with path and regexp...
perInsertionsFile = append(filebase,'_pergene_insertions.txt');
perGeneFile = append(filebase,'_pergene.txt');
essentialGeneFile = append(filebase, '_peressential.txt');
wigFile = append(filebase, '.wig');
bedFile = append(filebase,'.bed.txt');
[coordinates, chromosomePerGene, chromosomeEndPos, geneID, standardGeneID, geneSenseBinary] = readFileGFF(); %Read in the gff annotation file

% [bedFile, path]=uigetfile('Pick a file'); %<-point to the adjusted .bed
% file ALTERNATIVE METHOD


%pergene_insertions.txt: Gene name	Chromosome	Start location	End location  Insertion locations	Reads per insertion location

Table1 = readtable(perInsertionsFile); %load .insertions file
Table2 = readtable(perGeneFile); %load .pergene file
Table3 = readtable(essentialGeneFile); %load peressential file
data = load(bedFile, '-ascii'); %.bed file: Chromosome name start location end location . score. Note: file must be changed manually to be numerical only
readData = loadWigfile(wigFile); %Takes a long time, can we shorten this somehow?
essentialGenes = readtable('Cerevisiae_AllEssentialGenes_List.txt'); %Read in the essential genes file

geneCount = height(Table1);
essentialGeneList = zeros(length(standardGeneID),1);
for ii = 1:height(standardGeneID)
    essentialGeneList(ii) = sum(strcmp(standardGeneID{ii},essentialGenes.AllEssentialGenesFoundInLists___Cerevisiae_EssentialGenes_List_)); %For each gene check if it is listed as essential
end
data = data(:,[1,2]); %We are only interested in the first 2 entries for the .bed file.
data(:,1) = transformChromosomeIdentifiers(data(:,1)); % transformation of used chromosome identifiers to 1:17. 

%% put all features and coordinates on one single huge concatenated chromosome
% Why do we have non-unique entries to tn locations? probably they go in reverse directions?

lengthPreviousChr = [0; cumsum(chromosomeEndPos)]; %get the length in bp of each chromosome. The first chromosome starts at 0
lengthPreviousChr = double(lengthPreviousChr(1:end-1)); % Remove the last entry as it does not act as a starting point to another chromosome
%force type double to allow for addition with other double type data. May
%be useful to force double immediately upon reading in, also for other data
tnCoordinates = data;

[tnCoordinatesConcat,geneStartCoordinatesConcat,geneEndCoordinatesConcat] = concatCoordinatesIntoSingleChr(Table1,tnCoordinates,geneCount,lengthPreviousChr);


%% reading in wig file

chrEntries = find(isnan(readData(:,1)));
chrEntries = [chrEntries; length(readData)];

for ii = 1:17 %ARE THEY IN THE SAME ORDER? OTHERWISE THIS DOES NOT MAKE SENSE. Readdata counts up 1133 - 1148. Mito last
    readData(chrEntries(ii):chrEntries(ii+1),1) = readData(chrEntries(ii):chrEntries(ii+1),1) + lengthPreviousChr(ii); 
end
readData(isnan(readData(:,1)),:) = []; %Remove rows which contain a nan value


%% Get number of reads & transposon per gene
tnPerGene = Table2.NumberOfTransposonsPerGene;
readPerGene = Table2.NumberOfReadsPerGene;

geneLength = Table1.EndLocation - Table1.StartLocation; %CAN USE GFF DATA FOR THIS
tnDensity = tnPerGene./geneLength; %define transposon density per gene

% Count number of transposon per gene minus end and beginning
[tnPerGeneMinTen,readPerGeneMinTen,tnDensityMinTen] = getTransposonsPerGeneMinTen(Table1,geneLength,geneCount); %the read count from table1 does not correspond with read count from table2. 

% Find the longest transposon free interval within each gene
tnFreeInterval = getTnFreeIntervalPerGene(tnCoordinatesConcat,geneStartCoordinatesConcat,geneEndCoordinatesConcat,geneCount);

%% intergenic regions of 20kb 
% RUNS, BUT DOES IT actuaLLY WORK AS INTENDED?

[intergenicLength20kbPerGene,intergenicTn20kbPerGene,intergenicTnDensity20kb] = getIntergenicTnDensity(geneLength,chromosomeEndPos,tnCoordinates,coordinates,chromosomePerGene);


%% creating a list with the number of transposons and reads in the 100bp 5' of each gene (promotor region)

geneCoordinatesConcat(:,1) = geneStartCoordinatesConcat;
geneCoordinatesConcat(:,2) = geneEndCoordinatesConcat;

% creating a concatinated verson of all transposons and reads.
readPosConcat = zeros(max(max(geneCoordinatesConcat)),1);
readPosConcat(readData(:,1)) = readData(:,2);
tnPosConcat = zeros(max(max(geneCoordinatesConcat)),1);
tnPosConcat(tnCoordinatesConcat) = 1;

[promotorRead] = promotorData(readPosConcat,geneCoordinatesConcat,geneCount,geneSenseBinary);
[promotorTn] = promotorData(tnPosConcat,geneCoordinatesConcat,geneCount,geneSenseBinary);

promotorReadPerTn=promotorRead./promotorTn;

%% Combine training/testing data into a single table
readnumb = 1; %TEMPORARY SO I CAN CONTINUE. VERY VERY WRONG! WHAT IS READNUMB? %READNUMB IS PROBABLY THE READ COUNT PER BP? NO PROB NOT, MAYBE TOTAL READ COUNT

readPerGenePerTn = readPerGene./(tnPerGene.*sum(readnumb)); %unused?
readPerGenePerTnMinTen = readPerGeneMinTen./(tnPerGeneMinTen.*sum(readnumb)); %unused?
tnFreeIntervalPerGeneLength = (double(tnFreeInterval)./geneLength);
intergenicTn20kbUp = intergenicTn20kbPerGene(:,1);
intergenicTn20kbDown = intergenicTn20kbPerGene(:,2);
NI20kb = tnDensity./((intergenicTnDensity20kb(:,1)+intergenicTnDensity20kb(:,2))/2); %IntergenicTnDensity20kb can be 0... --> inf values
%NI20kb is compares the density of transposons in the local region to the
%transposons in the gene

dataTableOfFeatures = table(tnPerGeneMinTen, readPerGeneMinTen, geneLength, tnDensity, tnFreeIntervalPerGeneLength, promotorTn, intergenicTn20kbUp, intergenicTn20kbDown, NI20kb, essentialGeneList);% 
dataTableOfFeatures_WT = dataTableOfFeatures;

%% Train classifier 
% [A, B, C] = trainClassifier(dataTableOfFeatures) USE BAGGED TREES !!!!
%%
% yfit = trainedModelWTStripped.predictFcn(dataTableOfFeatures);

classifier = baggedTreesModel_220125;

%% Normalization for selected datasets to combine them into a single set

%tntot1=410169; %WT 1
%readtot1=31794831; %WT 1
%tntot2=356875; %WT 2
%readtot2=1530285; %WT 2

% tnCountA = 100000; %Total number of transposons found in the experiment A
% readCountA = 1000000; %Total number of reads found in the experiment A
% 
% tableNormA = TestTable_DplI; %Select table A for normalization
% tableNormA{:,[1 4:5]} = tableNormA{:,[1, 4:5]}/tnCountA; %divide all entries based on transposon count (name includes tn) by the total number of transposons
% tableNormA{:,2} = tableNormA{:,2}/readCountA; %divide all entries based on read count by the total number of reads
% 
% 
% tnCountB = 100000; %Total number of transposons found in experiment B
% readCountB = 1000000; %Total number of reads found in experiment B
% 
% tableNormB = TestTable_DplI; %Select table B for normalization
% tableNormB{:,[1 4:5]} = tableNormB{:,[1, 4:5]}/tnCountB; %divide all entries based on transposon count (name includes tn) by the total number of transposons
% tableNormB{:,2} = tableNormB{:,2}/readCountB; %divide all entries based on read count by the total number of reads
% 
% EssentialityTable_A_B = [tableNormA ; tableNormB];

% Use Apps:  Classification Learner to train the classifier
% 
% %% Create training data set with equal number of ess and non-ess genes
% 
% A1=rand(6603,1)<=0.82;
% A2=logical(A1-essentialyesno);
% essentialyesnoreduced=essentialyesno;
% essentialyesnoreduced(A2)=[];
% tnpergenemintenred=tnPerGeneMinTen;
% tnpergenemintenred(A2)=[];
% readpergenemintenred=readPerGeneMinTen;
% readpergenemintenred(A2)=[];
% genelengthred=geneLength;
% genelengthred(A2)=[];
% tndensityred=tnDensity;
% tndensityred(A2)=[];
% tndensitymintenred=tnDensityMinTen;
% tndensitymintenred(A2)=[];
% %tndensityupred=tndensityup;
% %tndensityupred(A1)=[];
% %tndensitydownred=tndensitydown;
% %tndensitydownred(A1)=[];
% tnfreeintervalpergenelengthred=tnFreeIntervalPerGeneLength;
% tnfreeintervalpergenelengthred(A2)=[];
% promotortnred=promotortn;
% promotortnred(A2)=[];
% intergenictn20kbupred=intergenictn20kbup;
% intergenictn20kbupred(A2)=[];
% intergenictn20kbdownred=intergenictn20kbdown;
% intergenictn20kbdownred(A2)=[];
% NI20kbred=NI20kb;
% NI20kbred(A2)=[];
% 
% 
% %EssentialityTableGred = table(tnpergenemintenred', readpergenemintenred', genelengthred, tndensityred', tndensityupred, tndensitydownred, tnfreeintervalpergenelengthred', essentialyesnoreduced); %
% TestTable_DplI_red = table(tnpergenemintenred', readpergenemintenred', genelengthred, tndensityred', tnfreeintervalpergenelengthred', promotortnred', intergenictn20kbupred, intergenictn20kbdownred, NI20kbred);%
% %EssentialityTable_WTI_red = table(tnpergenemintenred', readpergenemintenred', genelengthred, tndensitymintenred', tnfreeintervalpergenelengthred', promotortnred', intergenictn20kbupred, intergenictn20kbdownred, NI20kbred, essentialyesnoreduced);
% 
% 
%% Self prediction

%  Predicted_WT2_WT1WT2training = classifier.predictFcn(EssentialityTableK_WT2);
%  
%  Predicted_WT1_WT2training_bagged_10learners = SATAY_classifier_WTII_bagged_10learners.predictFcn(EssentialityTableK_WT1_norm);
%  sum(Predicted_WT1_WT2training_bagged_10learners==essentialyesno)

%% Test DplI data set

%NOMALIZATION STEP?
% tntotDpl=590079;
% readtotDpl=15077158;
% 
% TestTable_DplI_norm = TestTable_DplI;
% 
% TestTable_DplI_norm{:,1} = TestTable_DplI{:,1}/tntotDpl;
% TestTable_DplI_norm{:,2} = TestTable_DplI{:,2}/readtotDpl;
% TestTable_DplI_norm{:,4:9} = TestTable_DplI{:,4:9}/tntotDpl;
% dataTableOfFeaturesyTW001 = dataTableOfFeatures;
% 
predicted_WT = classifier.predictFcn(dataTableOfFeatures_WT);
predicted_yLIC = classifier.predictFcn(dataTableOfFeatures_yLIC);
predicted_yTW = classifier.predictFcn(dataTableOfFeatures_yTW);

essdifference_yLIC_WT = ~(predicted_WT == predicted_yLIC);
essdifference_yLIC_yTW = ~(predicted_yTW == predicted_yLIC);


%% classification score
%USES fitcensemble ?
% temp = trainedModelWTStripped.ClassificationKNN;

kNNMdl = classifier.ClassificationEnsemble;
[labels,Class_score_yLIC137_7] = predict(kNNMdl,dataTableOfFeatures_yLIC);
[labels,Class_score_yTW001_4] = predict(kNNMdl,dataTableOfFeatures_yTW);


%% check different essential genes in DplI and DplIPsdII with >95% accuracy

%sum(Predicted_DplIPsdII_WT1WT2training==Predicted_DplI_WT1WT2training)
%Predicted_ess_in_DplI_not_in_DplIPsdII = find((Predicted_DplIPsdII_WT1WT2training-Predicted_DplI_WT1WT2training)==-1);

Predicted_Essential_yLIC_95=(Class_score_yLIC137_7(:,2)>=0.95);
Predicted_Essential_yTW_95=(Class_score_yTW001_4(:,2)>=0.95);

Predicted_ess95_in_yLIC_not_in_yTW = find((Predicted_Essential_yTW_95-Predicted_Essential_yLIC_95)==-1);
Predicted_ess95_in_yTW_not_in_yLIC = find((Predicted_Essential_yTW_95-Predicted_Essential_yLIC_95)==1);

% Histogram of classification scores

figure(1)
hist(Class_score_yLIC137_7(:,1),100)
xlabel('Classification score')
ylabel('#genes')
title('Distribution of Classification scores in yLIC137\_7')
set (gca, 'Fontsize', 20)

figure(2)
hist(Class_score_yTW001_4(:,1),100)
xlabel('Classification score')
ylabel('#genes')
title('Distribution of Classification scores in yTW001\_4')
set (gca, 'Fontsize', 20)


%% Create a table (matrix) containing all the essential genes in different genetic backgrounds
essentialGenesInDifferentBackgrounds = [dataTableOfFeatures_WT.essentiality, predicted_WT, predicted_yLIC, predicted_yTW];
essentialAtLeastOnce = (sum(essentialGenesInDifferentBackgrounds,2)>0);
onlyEssentialGenesInDifferentBackgrounds = [geneList(essentialAtLeastOnce), num2cell([dataTableOfFeatures_WT.essentiality(essentialAtLeastOnce), predicted_WT(essentialAtLeastOnce), predicted_yLIC(essentialAtLeastOnce), predicted_yTW(essentialAtLeastOnce)])];

%%

tableVars = {'annotatedEssentialWT', 'predictedEssential_WT', 'predictedEssential_yLIC', 'predictedEssential_yTW'};

A = geneList(essentialAtLeastOnce);
tempTable = array2table(essentialGenesInDifferentBackgrounds(essentialAtLeastOnce,:));
tempTable.('GeneID') = A;
essentialGenesInDifferentBackgrounds_Table = [tempTable(:,5),tempTable(:,1:4)];

for ii = 1:4
    essentialGenesInDifferentBackgrounds_Table.Properties.VariableNames{append('Var',num2str(ii))} = tableVars{ii};
end

% essentialGenesInDifferentBackgrounds_Table = join(cell2table(geneList(essentialAtLeastOnce)) , array2table(essentialGenesInDifferentBackgrounds(essentialAtLeastOnce,:)));

% Write table to excel file
% cd('D:\Users\Thomas\Studie\MEP\MEP_Thomas\src\MATLAB\classifierWorkspaces')
% writetable(essentialGenesInDifferentBackgrounds_Table,'essentialGenesInDifferentBackgrounds_Table.xlsx','Sheet',1);

%%
tempTableA = array2table(classScores);
tempTableA.('GeneID') = geneList;
classScoreTable = [tempTableA(:,3),tempTableA(:,1),tempTableA(:,2)];
% writetable(classScoreTable,'classScoreTable.xlsx','Sheet',1);

