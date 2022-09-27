function [PredLabels,PredLabels_ALL] = SemiSupervised_Graph_Classifier(Laplace_Mat, initial_Labels, TestSamLoc, alpha)
N=size(initial_Labels,1);
UU = (1-alpha) * speye(N) + alpha * Laplace_Mat;
F = UU \ initial_Labels; % classification function F = (I - \alpha S)^{-1}Y
[~, PredLabels] = max(F, [], 2); %simply checking which of elements is largest in each row
PredLabels_ALL=PredLabels;
PredLabels=PredLabels(TestSamLoc);