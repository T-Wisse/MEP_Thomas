%https://sites.google.com/site/satayusers/

clear all 
close all

%ORDER OF FILES LOADED: BEDFILE, perInsertionsFile, perGeneFile, essentialGeneFile

%% Load insertion data (.bed file)

[bedfile, path]=uigetfile('Pick a file'); %<-point to the .bed file
cd(path)


%pergene_insertions.txt: Gene name	Chromosome	Start location	End location	Insertion locations	Reads per insertion location
%.bed file: Chromosome name start location end location . score

%We are only interested in the first 2 entries for the .bed file. The end location is
%just start loc + 1

data = load(bedfile, '-ascii'); %Load only works if file has no text. Note: automatically removes zeros (0) at the start of a number
data = data(:,[1,2]); %Select columns 1 and 2. It should be a score??
load('yeastGFF.mat') %<- both yeastGFF.mat should be in the same folder as the BAM file. Otherwise adapt here



%% Load insertion data per gene (.txt file)

filebase = 'yTW001_4_merged_cleaned_forward_reads_trimmed.sorted.bam';
%filebase = yLIC137_7_merged_cleaned_forward_reads_trimmed.sorted.bam
%%renamed the .wig file to fit this naming scheme
%filebase = WT_merged-techrep-a_techrep-b_trimmed.sorted.bam 

perInsertionsFile = append(filebase,'_pergene_insertions.txt');
perGeneFile = append(filebase,'_pergene.txt');
essentialGeneFile = append(filebase, '_peressential.txt');
wigFile = append(filebase, '.wig');

%Should be able to do this with path and regexp... but for now we cheat

% [perInsertionsFile, path]=uigetfile('Pick a file'); %<-point to the pergene_insertions.txt file
T = readtable(perInsertionsFile); %alternative to load, works with text.

% [perGeneFile, ~]=uigetfile('Pick a file'); %<-point to the pergene_insertions.txt file
Table2 = readtable(perGeneFile); %alternative to load, works with text.

% [essentialGeneFile, path]=uigetfile('Pick a file'); %<-point to the peressential.txt file
Table3 = readtable(essentialGeneFile); %alternative to load, works with text.
essentialGeneList = Table3.GeneName;

readData = loadWigfile(wigFile);

geneCount = height(T);


%% transformation of used chromosome identifiers to 1:17. Works but it's an ugly fix

% chromosomeIdentifiers = unique(data(:,1)); %Get the 17 unique chromosome identifiers used in the data
% chromosomeIdentifiers = chromosomeIdentifiers - 1132;
% data(:,1) = data(:,1)-1132;
% data(data(:,1)==92,1) = data((data(:,1)==92),1)-75;
% if unique(data(:,1)) ~= [1:17]
%     disp('An issue occured with transformation of chromosome identifiers')
% end

data(:,1) = transformChromosomeIdentifiers(data(:,1));




%% NOTE: RUNS. BUT I DON't USE THIS FILE SINCE THE RUN WAS DONE WITH A DIFFERENT FILE...?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%prepare chromosomal features from gff file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[genes,features,essential] = getFeaturesFromGff(gff); %

% chromosomeEndPos = cell2mat(gff(find(strcmp(gff(:,3),'omosome')),5)); %
% find chromosome lengths %DOES NOT SEEM TO WORK FOR ME AS INTENDED.
% LOCATIONS DO NOT MATCH CHROMOSOMES. DIFFEERT GFF FILE MAYBE?

temp = genes.coordinates(2:end,1)-genes.coordinates(1:end-1,1); %fun idea, doesnt work properly because of stupid reasons
chromosomeEndPos = genes.coordinates((temp<0),2);
chromosomeEndPos = [chromosomeEndPos; genes.coordinates(end,2)];
chromosomeEndPos = chromosomeEndPos + 2000; %HCACK TO BYPASS ISSUE. REALLY SHOULD NOT DO THIS...

%% OBSOLETE?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create array for all bp positions per chromosome
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% res = 200;
% 
% 
% 
% for k=1:length(chromosomeEndPos)
%     my_field = strcat('ch',num2str(k));     % create variable with all bp positions for every chromosome
%     tnpos.(my_field) = zeros(chromosomeEndPos(k),1);
%     readpos.(my_field) = zeros(chromosomeEndPos(k),1);
%     
%     x=tnCoordinates(:,1)==k; %find all positions in tncoordinates that give tn locations for specific (k) chromosome
%     y=tnCoordinates(x,2); %takes coordinates of tns on chromosome k
%     tnpos.(my_field)(y) = 1; % gives a 1 for every bp position that has a transposon insertion
%     readpos.(my_field)(y) = readnumb(x); % puts number of reads on every bp position
%     
% end

 
%% Probably works now. Why do we have non-unique entries to tn locations? probably they go in reverse directions?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%put all features and coordinates on one single huge concatenated chromosome
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lengthPreviousChr = [0; cumsum(chromosomeEndPos)]; %get the length in bp of each chromosome. The first chromosome starts at 0
lengthPreviousChr = lengthPreviousChr(1:end-1); % Remove the last entry as it does not act as a starting point to another chromosome

tnCoordinates = data;

[tnCoordinatesConcat,geneStartCoordinatesConcat,geneEndCoordinatesConcat] = concatCoordinatesIntoSingleChr(T,tnCoordinates,geneCount,lengthPreviousChr);


%% reading in wig file

chrEntries = find(isnan(readData(:,1)));
chrEntries = [chrEntries; length(readData)];

for ii = 1:17 %ARE THEY IN THE SAME ORDER? OTHERWISE THIS DOES NOT MAKE SENSE. Readdata counts up 1133 - 1148. Mito last
    readData(chrEntries(ii):chrEntries(ii+1),1) = readData(chrEntries(ii):chrEntries(ii+1),1) + lengthPreviousChr(ii); 
end
readData(isnan(readData(:,1)),:) = [];

%% OBSOLETE? Testing: Readpos, tnpos 
% 
% readPosConcat = zeros(sum(chromosomeEndPos),1); %Array with the length of all chromosome combined
% 
% tnPosConcat = zeros(sum(chromosomeEndPos),1); %Array with the length of all chromosome combined
% tnPosConcat(tnCoordinatesConcat) = 1; %Set locations of transposons in the geneome to 1
% 
% readCountConcat = [];
% for ii = 1:geneCount
%     readCountConcat = [readCountConcat str2num(T.ReadsPerInsertionLocation{ii})];
% end
% 
% readPosConcat(logical(tnPosConcat)) =  readCountConcat';

%% OBSOLETE? brute forcing readpos?

% 
% kk = 0;
% for ii = 1:geneCount
%     myField = T.Chromosome{ii}; 
%     if ~exist(strcat('readCount.(',myField,')'))
%         readCount.(myField) = [];
%     end
%     for jj = 1:length(str2num(A{ii}))
%         kk = kk+1;
%         readCount.(myField) = [readCount.(myField) str2num(T.ReadsPerInsertionLocation{ii})];
%     end
% end
%% Get number of reads & transposon per gene
tnPerGene = Table2.NumberOfTransposonsPerGene;
readPerGene = Table2.NumberOfReadsPerGene;

geneLength = T.EndLocation - T.StartLocation;
% geneLength = genes.coordinates(:,2)-genes.coordinates(:,1); %define gene length USES a DIFFERENT LIST OF GENES THAN READ DATA HAS!
tnDensity = tnPerGene./geneLength; %define transposon density per gene


%% Count number of transposon per gene minus end and beginning


% tnPerGeneMinTen = zeros(geneCount,1);
% readPerGeneMinTen = zeros(geneCount,1);
% 
% for ii=1:geneCount %Loop through each gene. May have an issue with genes on seperate chromosomes here mapping to the same location...
%     insertionLocations = str2num(T.InsertionLocations{ii});
%     readsPerInsertionLocations = str2num(T.ReadsPerInsertionLocation{ii});
%     %Find all reads (tncoordinates_copy) that start at a basepair number between the start and end of a gene (start_coor = gene.coordinates(ii,1) and end_coor = gene.coordinates(ii,2))
%     xx=insertionLocations>=T.StartLocation(ii)+geneLength(ii)*0.1 & insertionLocations<=T.EndLocation(ii)-geneLength(ii)*0.1;
%     %determine how many reads there are per gene (index of each read is stored in xx). 
%     tnPerGeneMinTen(ii) = sum(xx);
%     readPerGeneMinTen(ii) = sum(readsPerInsertionLocations(xx));
% end
% 
% tnDensityMinTen = tnPerGeneMinTen./geneLength;

[tnPerGeneMinTen,readPerGeneMinTen,tnDensityMinTen] = getTransposonsPerGeneMinTen(T,geneLength,geneCount);

%% intergenic regions of 20kb % TO DO

wl = 20000; % define window length

roman = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "mt"];
t=0;
u=0;
intergeniclength20kbpergene = zeros(length(geneLength),2);


for k=1:length(chromosomeEndPos)
    my_field = strcat('ch',num2str(k));     % create variable with all bp positions for every chromosome
    genepos.(my_field) = zeros(chromosomeEndPos(k),1);
    intergenictnpos.(my_field)= zeros(chromosomeEndPos(k),1);
    
    for i=t+1:t+sum(strcmp(roman(k),genes.chr)) % create binary with 1 for every bp that is part of a gene for every chromosome
    genepos.(my_field)(genes.coordinates(i,1):genes.coordinates(i,2)) = 1; 
    t = i;
    end
    
    tnPos.(my_field)=tnCoordinates(tnCoordinates(:,1) == k,2);
%     intergenictnpos.(my_field) = intergenictnpos; DOES NOT MAKE ANY
%     SENSE?
    genepos.(my_field)=logical(genepos.(my_field));
    intergenictnpos.(my_field)(genepos.(my_field)) = 0;
    
    for j=u+1:u+sum(strcmp(roman(k),genes.chr)) % every gene of each chromosome
        if genes.coordinates(j,1)-wl<1
            intergeniclength20kbpergene(j,1) = sum(1-(genepos.(my_field)(1:genes.coordinates(j,1)))); % calculate intergenic length before gen in 20kb window
            intergenictn20kbpergene(j,1) = sum(intergenictnpos.(my_field)(1:genes.coordinates(j,1))); % calculates number of tn in this region
        else
            intergeniclength20kbpergene(j,1) = sum(1-(genepos.(my_field)(genes.coordinates(j,1)-wl:genes.coordinates(j,1)))); % calculate intergenic length before gen in 20kb window
            intergenictn20kbpergene(j,1) = sum(intergenictnpos.(my_field)(genes.coordinates(j,1)-wl:genes.coordinates(j,1)));
        end
        
        if genes.coordinates(j,2)+wl>chromosomeEndPos(k)
             intergeniclength20kbpergene(j,2) = sum(1-(genepos.(my_field)((genes.coordinates(j,2):chromosomeEndPos(k))))); % calculate intergenic length before gen in 20kb window
             intergenictn20kbpergene(j,2) = sum(intergenictnpos.(my_field)(genes.coordinates(j,2):chromosomeEndPos(k)));
        else
            intergeniclength20kbpergene(j,2) = sum(1-(genepos.(my_field)(genes.coordinates(j,2):genes.coordinates(j,2)+wl))); % calculate intergenic length before gen in 20kb window
            intergenictn20kbpergene(j,2) = sum(intergenictnpos.(my_field)(genes.coordinates(j,2):genes.coordinates(j,2)+wl));
        end
        u=j;
    end
    
end

intergenic_tndensity20kb = intergenictn20kbpergene./intergeniclength20kbpergene;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Longest transposon free interval within a gene
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=1:geneCount %Loop trough each gene. For each, find the longest interval between transposons in bp
    ww = tnCoordinatesConcat>=geneStartCoordinatesConcat(ii)&tnCoordinatesConcat<=geneEndCoordinatesConcat(ii);
    ww = tnCoordinatesConcat(ww);
    ww = [geneStartCoordinatesConcat(ii) ; ww ; geneEndCoordinatesConcat(ii)]; %Concat the transposon locations of gene ii with the start and end of the gene
    tnFreeInterval(ii,1) = max(diff(ww)); %determine the largest difference bertween previously concatenated locations to find longest region without a transposon in the gene 
end


% %% creating a list with the number of transposons and reads in the 100bp 5' of each gene (promotor region)
% TO DO ... This I need to use the genes file for...

genesense = cell2mat(gff(features.genes,7));
genesensebinary = genesense=='+';

% creating a concatinated verson of all transposons and reads.
% tnPos_concat = [tnPos.ch1 ; tnPos.ch2 ; tnPos.ch3 ; tnPos.ch4 ; tnPos.ch5 ; tnPos.ch6 ; tnPos.ch7 ; tnPos.ch8 ; tnPos.ch9 ; tnPos.ch10 ; tnPos.ch11 ; tnPos.ch12 ; tnPos.ch13 ; tnPos.ch14 ; tnPos.ch15 ; tnPos.ch16 ; tnPos.ch17];
% readpos_concat = [readpos.ch1 ; readpos.ch2 ; readpos.ch3 ; readpos.ch4 ; readpos.ch5 ; readpos.ch6 ; readpos.ch7 ; readpos.ch8 ; readpos.ch9 ; readpos.ch10 ; readpos.ch11 ; readpos.ch12 ; readpos.ch13 ; readpos.ch14 ; readpos.ch15 ; readpos.ch16 ; readpos.ch17];
readPos_concat = readData(:,1);
genes.coordinates_concat(:,1) = geneStartCoordinatesConcat;
genes.coordinates_concat(:,2) = geneEndCoordinatesConcat;

tnPos_concat = tnCoordinatesConcat; %maybe?


for k=1:length(genes.coordinates_concat)
    if genesensebinary(k)==1    
        promotortn(k)=sum(tnPos_concat(genes.coordinates_concat(k,1)-100:genes.coordinates_concat(k,1)));
    %     promotorread(k)=sum(readPos_concat(genes.coordinates_concat(k,1)-100:genes.coordinates_concat(k,1)));
    else
        promotortn(k)=sum(tnPos_concat(genes.coordinates_concat(k,2):genes.coordinates_concat(k,2)+100));
    %     promotorread(k)=sum(readPos_concat(genes.coordinates_concat(k,2):genes.coordinates_concat(k,2)+100));
    end
end

% promotorreadpertn=promotorread./promotortn;



%% Determine which genes are essential

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Statistical Learning Matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
essentiality = zeros(geneCount,1);

for ii = 1:geneCount
    
    if sum(strcmp(essentialGeneList, T.GeneName{ii})) == 1
        essentiality(ii) = 1;
    end
end

%% form table

tnpergeneT = tnPerGene';
readpergeneT = readPerGene';
tndensityT = tnDensity;

% essentiality = char(essentiality);
essentialyesno=essentiality;

readnumb = 1 %TEMPORARY SO I CAN CONTINUE. VERY VERY WRONG! %READNUMB IS PROBABLY THE READ COUNT PER BP

readpergenepertn = readPerGene./(tnPerGene.*sum(readnumb));
readpergenepertnminten = readPerGeneMinTen./(tnPerGeneMinTen.*sum(readnumb));
tnFreeIntervalPerGeneLength = (double(tnFreeInterval)./geneLength);
intergenictn20kbup = intergenictn20kbpergene(1:6600,1); %FIX THIS
intergenictn20kbdown = intergenictn20kbpergene(1:6600,2); %FIX THIS
NI20kb = tndensityT./((intergenic_tndensity20kb(1:6600,1)+intergenic_tndensity20kb(1:6600,2))/2); %FIX THIS

% Select which data set you are working with

%EssentialityTable_WTI = table(tnpergeneminten', readpergeneminten', genelength, tndensityminten', tnfreeintervalpergenelength', promotortn', intergenictn20kbup, intergenictn20kbdown, NI20kb, essentialyesno);
%EssentialityTable_WTII = table(tnpergeneminten', readpergeneminten', genelength, tndensityminten', tnfreeintervalpergenelength', promotortn', intergenictn20kbup, intergenictn20kbdown, NI20kb, essentialyesno);

% TestTable_DplI = table(tnPerGeneMinTen', readPerGeneMinTen', geneLength, tndensityT, tnfreeintervalpergenelength', promotortn', intergenictn20kbup, intergenictn20kbdown, NI20kb);% 

%TestTable_DplIPsdII = table(tnpergeneminten', readpergeneminten', genelength, tndensityT, tnfreeintervalpergenelength', promotortn', intergenictn20kbup, intergenictn20kbdown, NI20kb);% 

dataTableOfFeatures = table(tnPerGeneMinTen, readPerGeneMinTen, geneLength, tndensityT, tnFreeIntervalPerGeneLength, intergenictn20kbup, intergenictn20kbdown, NI20kb, essentialyesno);% 
dataTableOfFeatures_yTW = dataTableOfFeatures;

%% Train classifier 
% [A, B, C] = trainClassifier(dataTableOfFeatures) USE BAGGED TREES !!!!
%%
% yfit = trainedModelWTStripped.predictFcn(dataTableOfFeatures);

classifier = trainedModel_bagged_211217;

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
[labels,score] = predict(kNNMdl,dataTableOfFeatures_yLIC);
Class_score_yLIC137_7 = score;

[labels,score] = predict(kNNMdl,dataTableOfFeatures_yTW);
Class_score_yTW001_4 = score;

%% check different essential genes in DplI and DplIPsdII with >95% accuracy

%sum(Predicted_DplIPsdII_WT1WT2training==Predicted_DplI_WT1WT2training)
%Predicted_ess_in_DplI_not_in_DplIPsdII = find((Predicted_DplIPsdII_WT1WT2training-Predicted_DplI_WT1WT2training)==-1);

Predicted_Essential_yLIC_95=(Class_score_yLIC137_7(:,2)>=0.95);
Predicted_Essential_yTW_95=(Class_score_yTW001_4(:,2)>=0.95);

Predicted_ess95_in_yLIC_not_in_yTW = find((Predicted_Essential_yTW_95-Predicted_Essential_yLIC_95)==-1);
Predicted_ess95_in_yTW_not_in_yLIC = find((Predicted_Essential_yTW_95-Predicted_Essential_yLIC_95)==1);

% Histogram of classification scores

% figure(1)
% hist(Class_score_yLIC137_7(:,1),100)
% xlabel('Classification score')
% ylabel('#genes')
% title('Distribution of Classification scores in yLIC137\_7')
% set (gca, 'Fontsize', 20)
% 
% figure(2)
% hist(Class_score_yTW001_4(:,1),100)
% xlabel('Classification score')
% ylabel('#genes')
% title('Distribution of Classification scores in yTW001\_4')
% set (gca, 'Fontsize', 20)

