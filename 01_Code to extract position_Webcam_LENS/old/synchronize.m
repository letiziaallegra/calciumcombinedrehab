function [res newfilename] = synchronize(filenamePOS_RAW,filenameFORCE_RAW, path,targetPos)
% function res = synchronize(res,filenameFORCE_RAW, path,targetPos)
fs = 25;
Ts=1/fs;

% units transformation

CurrDir = pwd;
cd(path)
load(filenamePOS_RAW)
cd(CurrDir)

% da Pos.mat (= filenamePOS_RAW) recupero le due colonne, res.X e res.LEDgrayscale
pos = sgolayfilt(res.X(:,1),3,11);
len = length(pos);
t = 0:Ts:(len-1)*Ts;
%trovo il valore di posizione target
mx = max(res.LEDgrayscale);
upperth = find(res.LEDgrayscale>mx-50); %<-- eventualmente cambiare (default mx-5)
pos_upperth = mean(pos(upperth));
%trovo il valore di posizione massima
mn = min(pos);
lowerth = find(pos<mn+3);
pos_lowerth = mean(pos(lowerth));

A = [pos_upperth 1; pos_lowerth 1];

%%targetPos = 8 or 10
par = A\[targetPos;18];
%par = A\[8;18];
%par = A\[10;18];
% par = A\[8;15];
Pos = par(1)*pos+par(2);




% sync
Tsound2 = zeros(len,1);
Tsound2(upperth)=1; %mette a 1 tutti gli elementi upperth (quelli a 8mm)
                    % Tsound2 = T_LedON (nella tesi)

res.Ts=Ts;
res.t=t;
res.Tsound2=Tsound2;

    %file da CVI
load(filenameFORCE_RAW)
prove = load(filenameFORCE_RAW, '-ascii' );
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
Tsound2_f_down = downsample(Tsound2_f,4);
if(length(Tsound2_f_down)>length(Tsound2))
    Tsound2_f_down_tr = Tsound2_f_down(1:length(Tsound2));
    [cm icm]=max(xcorr(Tsound2,Tsound2_f_down_tr));
    win = 2*len-1;
    win_center = floor(win/2);
    delay = icm-win_center      %quanto un segnale deve essere anticipato o ritardato
                                % per renderlo uguale all'altro
%     delay = abs(icm-win_center) NO
else
    Tsound2_f_down_tr = [Tsound2_f_down; zeros(length(Tsound2)-length(Tsound2_f_down),1)];
    [cm icm]=max(xcorr(Tsound2,Tsound2_f_down_tr));
    win = 2*len-1;
    win_center = floor(win/2);
    delay = icm-win_center
%     delay = abs(icm-win_center)NO
end


if(delay>0)
    temp = Pos(delay+1:end);
    clear Pos t
    Pos = temp;
    len = length(Pos);
%     t = 0:Ts:(len-1)*Ts;
    t = Ts:Ts:(len-1)*Ts;
elseif(delay<0)
    Pos_New_first = [1:abs(delay)]' * 0 %Pos(1,1);
%     Vel_New_first = [1:abs(delay)* res.Ts]' * 0 %Vel(1,1);
%     Acc_New_first = [1:abs(delay)* res.Ts]' * 0 %Acc_res;
    
    Pos=cat(1, Pos_New_first, Pos);
%     Vel=cat(1, Vel_New_first, Vel_res);
%     Acc=cat(1, Acc_New_first, Acc_res)

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

% prove(length(res.Pos)+1:end,:)=[];
if ( size(prove,1) > length(res.Pos))
    prove(length(res.Pos)+1:end,:)=[];
elseif ( size(prove,1) < length(res.Pos))
    res.Pos(size(prove,1)+1:end,:)=[];
    res.Vel(size(prove,1)+1:end,:)=[];
    res.Acc(size(prove,1)+1:end,:)=[];
end
prove(:,8:10) = [res.Pos res.Vel res.Acc];

newfilename = strcat(filenameFORCE_RAW(1:end-4),'_sync.txt');
if(fopen(newfilename)==-1)
    save(newfilename, 'prove', '-ascii');
else
    fclose('all')
    delete(newfilename);
    save(newfilename, 'prove', '-ascii');
end
% prova = [res.Pos res.Vel res.Acc];
% save('PROOOVA', 'prova', '-ascii');


return