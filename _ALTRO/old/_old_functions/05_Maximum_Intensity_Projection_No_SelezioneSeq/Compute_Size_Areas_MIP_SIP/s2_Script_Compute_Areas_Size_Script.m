%%% Scrpt to compute extentof the ROI selected with different threshold in the MIP and SIP

MIP_REHAB1   = [];
MIP_REHAB4   = [];
MIP_STROKE   = [];
MIP_CONTROL  = [];

aSIP_REHAB1  = [];
aSIP_REHAB4  = [];
aSIP_STROKE  = [];
aSIP_CONTROL = [];

pSIP_REHAB1  = [];
pSIP_REHAB4  = [];
pSIP_STROKE  = [];
pSIP_CONTROL = [];


iter = zeros(4,1);


%%%%%%%%%
CurrDir = cd;
cd .. 
PrevDir = cd;
Dir_Save_IM       = [CurrDir,'\Store_Figures'];
Dir_data          = [CurrDir,'\Store_MAT'];
MIP_tech_List     = {'max','sumA','sumP'};
LenMIP            = length(MIP_tech_List);
cd(CurrDir)
%%%%%%%%%

ListAnimal = dir(Dir_data);

for i=3:length(ListAnimal)
        
    NameAn = ListAnimal(i,1).name;
    load([Dir_data,'\',NameAn])
    
    sz = cellfun(@length,RegionBound_Index);
    
    sz(sz==0) = NaN;
    
    sz_MIP   = sz(1,1:5);
    sz_aSIP  = sz(2,1:5);
    sz_pSIP  = sz(3,1:5);
    
    if strfind(NameAn,'stroke_BoNT')
        
        %%%% Week_1
        if strfind(NameAn,'Week_1')
            
            %rehab 1
            MIP_REHAB1  = [MIP_REHAB1;  sz_MIP];
            aSIP_REHAB1 = [aSIP_REHAB1; sz_aSIP];
            pSIP_REHAB1 = [pSIP_REHAB1; sz_pSIP];
            
        elseif strfind(NameAn,'Week_4')
            
            %rehab 4
            MIP_REHAB4  = [MIP_REHAB4;  sz_MIP];
            aSIP_REHAB4 = [aSIP_REHAB4; sz_aSIP];
            pSIP_REHAB4 = [pSIP_REHAB4; sz_pSIP];
        end
        
        
        %%%% stroke
    elseif strfind(NameAn,'stroke')
        
        %stroke
        MIP_STROKE  = [MIP_STROKE;  sz_MIP];
        aSIP_STROKE = [aSIP_STROKE; sz_aSIP];
        pSIP_STROKE = [pSIP_STROKE; sz_pSIP];        
        
        
        %%%% control
    elseif strfind(NameAn,'control')
        
        %control
        MIP_CONTROL  = [MIP_CONTROL;  sz_MIP];
        aSIP_CONTROL = [aSIP_CONTROL; sz_aSIP];
        pSIP_CONTROL = [pSIP_CONTROL; sz_pSIP];      
        
    end
    
    clear sz_MIP sz_aSIP sz_pSIP
    
end



m_MIP_REHAB1   = nanmean(MIP_REHAB1);
s_MIP_REHAB1   = nanstd(MIP_REHAB1,[])/sqrt(size(MIP_REHAB1,2));
m_MIP_REHAB4   = nanmean(MIP_REHAB4);
s_MIP_REHAB4   = nanstd(MIP_REHAB4,[])/sqrt(size(MIP_REHAB4,2));
m_MIP_STROKE   = nanmean(MIP_STROKE);
s_MIP_STROKE   = nanstd(MIP_STROKE,[])/sqrt(size(MIP_STROKE,2));
m_MIP_CONTROL  = nanmean(MIP_CONTROL);
s_MIP_CONTROL  = nanstd(MIP_CONTROL,[])/sqrt(size(MIP_CONTROL,2));

m_aSIP_REHAB1   = nanmean(aSIP_REHAB1);
s_aSIP_REHAB1   = nanstd(aSIP_REHAB1,[])/sqrt(size(aSIP_REHAB1,2));
m_aSIP_REHAB4   = nanmean(aSIP_REHAB4);
s_aSIP_REHAB4   = nanstd(aSIP_REHAB4,[])/sqrt(size(aSIP_REHAB4,2));
m_aSIP_STROKE   = nanmean(aSIP_STROKE);
s_aSIP_STROKE   = nanstd(aSIP_STROKE,[])/sqrt(size(aSIP_STROKE,2));
m_aSIP_CONTROL  = nanmean(aSIP_CONTROL);
s_aSIP_CONTROL  = nanstd(aSIP_CONTROL,[])/sqrt(size(aSIP_CONTROL,2));

m_pSIP_REHAB1   = nanmean(pSIP_REHAB1);
s_pSIP_REHAB1   = nanstd(pSIP_REHAB1,[])/sqrt(size(pSIP_REHAB1,2));
m_pSIP_REHAB4   = nanmean(pSIP_REHAB4);
s_pSIP_REHAB4   = nanstd(pSIP_REHAB4,[])/sqrt(size(pSIP_REHAB4,2));
m_pSIP_STROKE   = nanmean(pSIP_STROKE);
s_pSIP_STROKE   = nanstd(pSIP_STROKE,[])/sqrt(size(pSIP_STROKE,2));
m_pSIP_CONTROL  = nanmean(pSIP_CONTROL);
s_pSIP_CONTROL  = nanstd(pSIP_CONTROL,[])/sqrt(size(pSIP_CONTROL,2));



figure('Name','Area Size')
subplot(311)
hold on
barwitherr([s_MIP_CONTROL' s_MIP_STROKE' s_MIP_REHAB1' s_MIP_REHAB4' ],[m_MIP_CONTROL' m_MIP_STROKE' m_MIP_REHAB1' m_MIP_REHAB4'])
xlim([0 6])
set(gca,'XTick',1:5)
xlabel('Threshold')
ylabel('Area [pixels]')
subplot(312)
hold on
barwitherr([s_aSIP_CONTROL' s_aSIP_STROKE' s_aSIP_REHAB1' s_aSIP_REHAB4' ],[m_aSIP_CONTROL' m_aSIP_STROKE' m_aSIP_REHAB1' m_aSIP_REHAB4'])
xlim([0 6])
set(gca,'XTick',1:5)
xlabel('Threshold')
ylabel('Area [pixels]')
subplot(313)
hold on
barwitherr([s_pSIP_CONTROL' s_pSIP_STROKE' s_pSIP_REHAB1' s_pSIP_REHAB4' ],[m_pSIP_CONTROL' m_pSIP_STROKE' m_pSIP_REHAB1' m_pSIP_REHAB4'])
xlim([0 6])
set(gca,'XTick',1:5)
xlabel('Threshold')
ylabel('Area [pixels]')






