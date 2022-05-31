%% This script, classifierApplication, uses the output of classifierWesselAdapted.m script. 
% It applies the supplied classifier to the supplied data and returns the
% scores. It also returns a matrix containing all essential genes in all
% different genetic background and a list of genes uniquely essential in
% one of the backgrounds (with high confidence)

%% load model
load('220224_WT_trainedModel.mat');
classifier = trainedModelWT_ab_220224;

%% load data table
dataset = 'allDataTables';
load(append(dataset,'.mat'));

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
predicted_WT = classifier.predictFcn(dataTableOfFeatures_wildType);
predicted_yLIC = classifier.predictFcn(dataTableOfFeatures_yLIC137_7);
predicted_yTW = classifier.predictFcn(dataTableOfFeatures_yTW001_4);

essdifference_yLIC_WT = ~(predicted_WT == predicted_yLIC);
essdifference_yLIC_yTW = ~(predicted_yTW == predicted_yLIC);


%% classification score
%USES fitcensemble ?
% temp = trainedModelWTStripped.ClassificationKNN;

%MAYBE MOVE THIS TO A LARGE MATRIX OR STRUCT?
kNNMdl = classifier.ClassificationEnsemble;
[labels_wildType, Class_score_wildType] = predict(kNNMdl,dataTableOfFeatures_wildType);
[labels_yLIC137_7, Class_score_yLIC137_7] = predict(kNNMdl,dataTableOfFeatures_yLIC137_7);
[labels_yLIC137_8, Class_score_yLIC137_8] = predict(kNNMdl,dataTableOfFeatures_yLIC137_8);
[labels_yTW001_4, Class_score_yTW001_4] = predict(kNNMdl,dataTableOfFeatures_yTW001_4);
[labels_yTW001_6, Class_score_yTW001_6] = predict(kNNMdl,dataTableOfFeatures_yTW001_6);
[labels_yWT03a_16, Class_score_yWT03a_16] = predict(kNNMdl,dataTableOfFeatures_yWT03a_16);
[labels_yWT03a_21, Class_score_yWT03a_21] = predict(kNNMdl,dataTableOfFeatures_yWT03a_21);
[labels_yWT04a_14, Class_score_yWT04a_14] = predict(kNNMdl,dataTableOfFeatures_yWT04a_14);
[labels_yWT04a_23, Class_score_yWT04a_23] = predict(kNNMdl,dataTableOfFeatures_yWT04a_23);

%% Find which genes are essential with >95% accuracy for all backgrounds

%sum(Predicted_DplIPsdII_WT1WT2training==Predicted_DplI_WT1WT2training)
%Predicted_ess_in_DplI_not_in_DplIPsdII = find((Predicted_DplIPsdII_WT1WT2training-Predicted_DplI_WT1WT2training)==-1);

%MAYBE MOVE THIS TO A LARGE MATRIX OR STRUCT?
predictedEssential_wildType=(Class_score_wildType(:,2)>=0.95);
predictedEssential_yLIC137_7_95=(Class_score_yLIC137_7(:,2)>=0.95);
predictedEssential_yLIC137_8_95=(Class_score_yLIC137_8(:,2)>=0.95);
predictedEssential_yTW001_4_95=(Class_score_yTW001_4(:,2)>=0.95);
predictedEssential_yTW001_6_95=(Class_score_yTW001_6(:,2)>=0.95);
predictedEssential_yWT03a_16_95=(Class_score_yWT03a_16(:,2)>=0.95);
predictedEssential_yWT03a_21_95=(Class_score_yWT03a_21(:,2)>=0.95);
predictedEssential_yWT04a_14_95=(Class_score_yWT04a_14(:,2)>=0.95);
predictedEssential_yWT04a_23_95=(Class_score_yWT04a_23(:,2)>=0.95);

%% Create a table (matrix) containing all the essential genes in different genetic backgrounds
load('geneList.mat')

essentialityMatrix = [dataTableOfFeatures_wildType.essentialGeneList, labels_wildType,labels_yLIC137_7,labels_yLIC137_8,labels_yTW001_4,labels_yTW001_6,labels_yWT03a_16,labels_yWT03a_21,labels_yWT04a_14,labels_yWT04a_23]; %The first column are the annotated essential genes in wildt-type. The second column is the predicted essential genes in wild-type
essentialAtLeastOnce = (sum(essentialityMatrix,2)>0);
essentialityMatrixReduced = [geneList(essentialAtLeastOnce), num2cell([dataTableOfFeatures_wildType.essentialGeneList(essentialAtLeastOnce), labels_yLIC137_7(essentialAtLeastOnce), ... 
    labels_yLIC137_8(essentialAtLeastOnce),labels_yTW001_4(essentialAtLeastOnce),labels_yTW001_6(essentialAtLeastOnce),labels_yWT03a_16(essentialAtLeastOnce), ...
    labels_yWT03a_21(essentialAtLeastOnce),labels_yWT04a_14(essentialAtLeastOnce),labels_yWT04a_23(essentialAtLeastOnce)])];


%%
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

%% Plot heatmap of classScores
A = 1000;

discreteTable(:,1) = classScoreTable(:,1);
discreteTable(:,2) = discretize(classScoreTable.classScores1,linspace(0,1,11));
discreteTable(:,3) = discretize(classScoreTable.classScores2,linspace(0,1,11));
discreteTable = discreteTable/10;
xvalues = {'yLIC137','yTW001'};
yvalues = geneList(1:A);

% h = heatmap(xvalues,yvalues,discreteTable);
h = heatmap(xvalues,yvalues,discreteTable(1:A,:));

h.Title = 'classification scores per gene';
h.XLabel = 'Genetic background';
h.YLabel = 'Genes';
h.GridVisible = 'off';

