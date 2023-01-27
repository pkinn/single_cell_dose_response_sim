function [nInh] = calculatenInhThresh(params,thresh, doseVec)
%calculatenInhThresh Calculate the number of cells inhibited to a threshold
%at each dose
%   Function takes in:
%   params: [nCells x 4] matrix with dr params
%   thresh: scalar representing a percentage of inhibition required to be
%   called inhibited
%   doseVec: vector of doses to use
%   Function then finds the first dose at which each single cell is below
%   that cells individual threshold: thresh*(activity with no drug)
%
nCells = size(params,1);
%Calculate the activity with no drug
for ii = 1:nCells
    noDrugAct(ii) = returnDR(params(ii,:),0);
end
%Calculate the threshold for each cell to be inhibited
scThresh = thresh*noDrugAct;

%Calculate dose response
for ii = 1:nCells
    dr(ii,:) = returnDR(params(ii,:), doseVec);
    drThresh(ii,:) = dr(ii,:)<=scThresh(ii);
end

nInh = sum(drThresh,1);

end