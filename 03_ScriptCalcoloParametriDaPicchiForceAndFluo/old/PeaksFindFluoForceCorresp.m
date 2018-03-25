%
%check corresponding force-fluo
%the field CorrespondanceForceFluoPeaks is appended to dataGCamp
%CorrespondanceForceFluoPeaks is a cell where each element correspond to a
%fluo ROI (following the order specified in Info) -> inside each elementm, there are 2 columns:
% [index of force peaks to take      corresponding index of ROI fluo peaks to take]
%-> the index are referred to the corresponding lists in dataGCamp.PeaksPar
%The correspondance bw force peak and fluo peak is done based on the starting time point of the peaks
%
function dataGCamp = PeaksFindFluoForceCorresp(dataGCamp,ForceCol,StatusCheck)

%size
[sx sy] = size(dataGCamp.PeaksPar);
%num ROI fluo
nROI = sy-1;


%number of the trial of force
nTr  = dataGCamp.PeaksPar{1,ForceCol}(:,1);
%status of the trials of force
stTr = dataGCamp.PeaksPar{1,ForceCol}(:,2);
%time start of the peak of force
stPeak = dataGCamp.PeaksPar{1,ForceCol}(:,3);

%index status of the trials selected of force
stTr_sel_index = find( stTr == StatusCheck(1) | stTr == StatusCheck(2) ); 
lenstTr_sel = length(stTr_sel_index);

startROI = [];
for j=1:nROI   
    
    [xSR ySR]  = size(startROI);  
%     stTr_FluoROIsel = dataGCamp.PeaksPar{1,j+1}(:,2);
%     stTr_FluoROIsel_index = find( stTr_FluoROIsel == StatusCheck(1) | stTr_FluoROIsel == StatusCheck(2) );
%     arrayToAdd = dataGCamp.PeaksPar{1,j+1}(stTr_FluoROIsel_index,3);
    arrayToAdd = dataGCamp.PeaksPar{1,j+1}(:,3);
    xAA        = size(arrayToAdd,1);
    
    %to set different size of the ROI vectors
    diffx = xSR-xAA;
    if diffx>=0
        arrayToAdd = [arrayToAdd ; ones(abs(diffx),1)*NaN];
    else
        startROI   = [startROI   ; ones(abs(diffx),ySR)*NaN];
    end
    %time start of the peak of all ROI fluo
    startROI = [startROI arrayToAdd];
end
    

StoreMin       = [];
StoreMin_Index = [];
for i=1:lenstTr_sel
   
    %current start point
    curr_stPeak = stPeak(stTr_sel_index(i));
    
    %min diff start force peak vs fluo peak (for each ROI)
    [mindiff i_mindiff] = min(abs(curr_stPeak-startROI));
    
    StoreMin       = [StoreMin; mindiff];
    StoreMin_Index = [StoreMin_Index ; i_mindiff];
    
    
end

%check (remove the occurence of the same peak)
CorrespondanceForceFluoPeaks = [];
for i=1:nROI 
    %single occurences
    singleOccur = unique(StoreMin_Index(:,i));
    
    occ_ind_Force = [];
    occ_ind_ROI   = [];
    for j=1:length(singleOccur)
        %current occurence
        occ_ind = find(StoreMin_Index(:,i)==singleOccur(j));        
        %more than one occurence of the same peak
        if length(occ_ind)>1
             %occur to mantain
             [occurOK i_occurOK] = min(StoreMin(occ_ind,i)); 
              occ_ind = occ_ind(i_occurOK);
        end
        occ_ind_Force = [occ_ind_Force ;   occ_ind];
        occ_ind_ROI   = [occ_ind_ROI   ;   StoreMin_Index(occ_ind,i)];
    end
    %index of the force peaks and ROI Fluo peaks -> to take account for
    ForceToKeep   = stTr_sel_index(occ_ind_Force);
    FluoROIToKeep = occ_ind_ROI;
    CorrespondanceForceFluoPeaks{1,i} = [ForceToKeep FluoROIToKeep];
    
end

dataGCamp.CorrespondanceForceFluoPeaks = CorrespondanceForceFluoPeaks;
end
