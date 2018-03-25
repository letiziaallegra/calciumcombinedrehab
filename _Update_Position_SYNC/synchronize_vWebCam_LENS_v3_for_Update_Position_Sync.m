function [prove] = synchronize_vWebCam_LENS_v3_for_Update_Position_Sync(res,fileTXT, targetPos)
% sincronizza solo il primo evento assoluot di switch con il primo inizio assoluto della fase 4 (posizione raggiunta)
fs = 25;
Ts=1/fs;
overHome = 2; %mm

% units transformation

CurrDir = pwd;

%da Pos.mat (= filenamePOS_RAW) recupero le due colonne, res.X e res.LEDgrayscale
pos  = res.X(:,1);
len = length(pos);
t = 0:Ts:(len-1)*Ts;

%trovo il valore di posizione target
%il segnale del LED sarà una sequenza di 0 con alcuni 1 corrispondenti al raggiungimento della Pos Max e della Pos Target
%però mi interessa solo il primo 1 delle sequenze di 1
LED_value = res.LEDgrayscale;
LED_value_diff = diff(LED_value);
%starting and ending point of 1-sequence
LED_ON_StEn_index = [find(LED_value_diff ==1)+1  find(LED_value_diff ==-1)];

%
lenLED_ON = size(LED_ON_StEn_index,1);
list_pos_LED_valueON_index = [1:lenLED_ON]';

%le posizioni dispari nel vettore corrispondono al raggiungimento di PosMax
xodd = bitand(abs(list_pos_LED_valueON_index),1);
LED_valueON_index_PosMax = LED_ON_StEn_index(logical(xodd),:);

%le posizioni pari nel vettore corrispondono al raggiungimento di PosTarget
xeven=~xodd;
LED_valueON_index_PosTarget = LED_ON_StEn_index(logical(xeven),:);

ArrayPosMax = zeros(len,1);
for i=1:size(LED_valueON_index_PosMax,1)
    ArrayPosMax(LED_valueON_index_PosMax(i,1):LED_valueON_index_PosMax(i,2)) = 1;
end
ArrayPosTarget = zeros(len,1);
for i=1:size(LED_valueON_index_PosTarget,1)
    ArrayPosTarget(LED_valueON_index_PosTarget(i,1):LED_valueON_index_PosTarget(i,2)) = 1;
end
    

%posizione nei frames corrispondente a PosMax
pos_PosMax = mean(pos(logical(ArrayPosMax)));

%posizione nei frames corrispondente a PosTarget
pos_TargetPos = mean(pos(logical(ArrayPosTarget)));


%%trasformazione di scala
m = (18-targetPos)/(pos_PosMax-pos_TargetPos);
q = 18-m*pos_PosMax;

Pos = m*pos+q;

%filtraggio posizione
Pos = sgolayfilt(Pos,3,11);

% %da upperth elimino tutti i punti precedenti all'inizio del task che potrebbero
% %trovarsi oltre il valore di Home, a causa per esempio di un movimento
% %dell'animale prima dell'inizio
% %trovo il valore minimo di Pos e ne sommo un valore arbitrario basso overHome
% indexPos0 = find(Pos<min(Pos)+overHome);
% %poi da upperth (=indici in cui si è in pos target) rimuovo indici precedenti all'inizio del task
% ArrayPosTarget = ArrayPosTarget(ArrayPosTarget>indexPos0(end));



% sync
Tsound2 = zeros(len,1);
Tsound2(logical(ArrayPosTarget))=1; %mette a 1 tutti gli elementi upperth (quelli a 8mm)
                    % Tsound2 = T_LedON (nella tesi)

res.Ts=Ts;
res.t=t;
res.Tsound2=Tsound2;

    %file da CVI
prove = fileTXT;
len_f = length(prove);
t_f = 0:0.01:(len_f-1)*0.01;
Tsound2_f = zeros(len_f,1);
Tsound2_f(prove(:,7)==4)=1; %mette a 1 tutti gli elementi di status 4 
                            %(Tsound2 on) quindi LED on --> 8mm   ---->
                            % Tsound2_f = T_Sound2 (nella tesi)
    %
    
    
    %N.B. video 25 f/s --> 1 sample ogni 0.04 s
    %     CVI   100 Hz --> 1 sample ogni 0.01 s --> sottocampiono per
    %                                      riportarmi al sampling del video
    %                                      -->1 sample ogni 4 
    
    
% Tsound2_f_down = downsample(Tsound2_f,4);
% if(length(Tsound2_f_down)>length(Tsound2))
%     Tsound2_f_down_tr = Tsound2_f_down(1:length(Tsound2));
%     [cm icm]=max(xcorr(Tsound2,Tsound2_f_down_tr));
%     win = 2*len-1;
%     win_center = ceil(win/2);
%     delay = icm-win_center      %quanto un segnale deve essere anticipato o ritardato
%                                 % per renderlo uguale all'altro
% %     delay = abs(icm-win_center) NO
% else
%     Tsound2_f_down_tr = [Tsound2_f_down; zeros(length(Tsound2)-length(Tsound2_f_down),1)];
%     [cm icm]=max(xcorr(Tsound2,Tsound2_f_down_tr));
%     win = 2*len-1;
%     win_center = ceil(win/2);
%     delay = icm-win_center
% %     delay = abs(icm-win_center)NO
% end


Tsound2_f_down = downsample(Tsound2_f,4);


%differenza  - sincronizza solo il primo evento assoluto di switch con il primo inizio assoluto della fase 4 (posizione raggiunta)
Ones_Tsound2_f_down = find(Tsound2_f_down==1);
Ones_Tsound2 = find(Tsound2==1);
delay = Ones_Tsound2(1)-Ones_Tsound2_f_down(1);


if(delay>0)
    temp = Pos(delay+1:end);
    clear Pos t
    Pos = temp;
    len = length(Pos);
    t = Ts:Ts:(len-1)*Ts;
elseif(delay<0)
    Pos_New_first = [1:abs(delay)]' * 0; 
    
    Pos=cat(1, Pos_New_first, Pos);

end

Pos_filt = sgolayfilt(Pos,3,11);
Vel = derivative(Pos_filt,0.04);
Vel_filt = sgolayfilt(Vel,3,11);
Acc = derivative(Vel_filt,0.04);
Acc_filt = sgolayfilt(Acc,3,11);

Pos_res = resample(Pos_filt,4,1);
Vel_res = resample(Vel_filt,4,1);
Acc_res = resample(Acc_filt,4,1);
t_res = 0:0.01:(length(Pos_res)-1)*0.01;
    
res.t_res = t_res';
res.Pos=Pos_res;
res.Vel=Vel_res;
res.Acc=Acc_res;


if ( size(prove,1) > length(res.Pos))
    prove(length(res.Pos)+1:end,:)=[];
elseif ( size(prove,1) < length(res.Pos))
    res.Pos(size(prove,1)+1:end,:)=[];
    res.Vel(size(prove,1)+1:end,:)=[];
    res.Acc(size(prove,1)+1:end,:)=[];
end
prove(:,9:11) = [res.Pos res.Vel res.Acc];
% prove = prove(1:end-1,:);

% newfilename = strcat(filenameFORCE_RAW(1:end-4),'_sync.txt');
% if(fopen(newfilename)==-1)
%     save(newfilename, 'prove', '-ascii');
% else
%     fclose('all')
%     delete(newfilename);
%     save(newfilename, 'prove', '-ascii');
% end


return