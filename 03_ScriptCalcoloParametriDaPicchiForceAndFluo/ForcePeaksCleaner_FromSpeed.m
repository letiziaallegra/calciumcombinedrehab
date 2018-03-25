%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Funzione che mantiene i picchi di Forza, solo se in prossimità di
%   picchi di Velocità. Cancella gli altri da InfoPicchi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [StartStopPeak_Force_New PeakArray_New] = ForcePeaksCleaner_FromSpeed(StartStopPeak_Force, StartStopPeak_Speed,Fs,PeakArray, Force, Speed, pl_onoff)

num_SpeedPeaks = size(StartStopPeak_Speed,1);
num_ForcePeaks = size(StartStopPeak_Force,1);

%member or not
memb_St  = zeros(num_ForcePeaks,num_SpeedPeaks) ;
memb_En  = zeros(num_ForcePeaks,num_SpeedPeaks) ;
memb_Int = zeros(num_ForcePeaks,num_SpeedPeaks) ;

if strcmp(pl_onoff,'on')
    figure
    subplot(211)
    plot(Force)
    subplot(212)
    plot(Speed)
end

for ns=1:num_SpeedPeaks
    
    %Interval of the Speed Peak
    St_Sp = StartStopPeak_Speed(ns,1);
    En_Sp = StartStopPeak_Speed(ns,2);
    Interval_Sp = [St_Sp:En_Sp]';
    
    if strcmp(pl_onoff,'on')
        subplot(211)
        plot(St_Sp:En_Sp, Speed(St_Sp:En_Sp),'r')
    end
    
    for nf=1:num_ForcePeaks
        %Interval of the Force Peak
        St_F = StartStopPeak_Force(nf,1);
        En_F = StartStopPeak_Force(nf,2);
        Interval_F = [St_F:En_F]';
        
        %member or not
        memb_St(nf,ns) = sum(ismember(Interval_Sp,St_F));
        memb_En(nf,ns) = sum(ismember(Interval_Sp,En_F));
        
        %force peak broader than speed peak
        memb_Int(nf,ns) = sum(ismember(Interval_F,St_Sp));
        
        if ns==1
            if strcmp(pl_onoff,'on')
                subplot(212)
                hold on
                plot(St_F:En_F,Force(St_F:En_F),'r')
            end
        end
        
    end
    
end


index_ForcePeaks_OK_NO = sum([memb_St memb_En memb_Int],2);
%nuova matrice picchi forza
StartStopPeak_Force_New = StartStopPeak_Force(logical(index_ForcePeaks_OK_NO),:);


PeakArray_New = zeros(length(PeakArray),1);
for j=1:size(StartStopPeak_Force_New,1)
    PeakArray_New(StartStopPeak_Force_New(j,1):StartStopPeak_Force_New(j,2)) = 1;
end


end





%
% if(size(InfoPicchiForza,1) ~= 0    &&     size(InfoPicchiVelocita,1) ~= 0)
%
%    InfoPicchiForzeBuffer = zeros(size(InfoPicchiForza,1), size(InfoPicchiForza,2));
%    controlloV = 0;
%
%    for(i = 1:num_prove) %for num_prove
%
%         %trova righe in InfoPicchiVelocita corrispondenti alla prova in questione
%         RigaProofVelocita = find( InfoPicchiVelocita(:,2)== i);
%         %trova righe in InfoPicchiForza corrispondenti alla prova in questione
%         RigaProofForza = find( InfoPicchiForza(:,2)== i);
%
%         %ad ogni prova il indexPeakOK_PREVIOUS è messo a 0. Nello svolgersi della scrematura dei picchi di forza, bisogna fare attenzione
%         % che ad un picco di forza corrisponda un solo picco di velocità. A volte, è una cosa rara, a un picco di forza vanno a coincidere
%         % due picchi di velocità; in tal caso poichè bisogna mantenere la relazione uno a uno, il picco di velocità è cancellato se il
%         % precedente valore di forza (recuperabile dall'indexPeakOK_PREVIOUS) è uguale all'attuale indexPeakOK
%         indexPeakOK_PREVIOUS = 0;
%
%         indexSpeedOK_PREVIOUS = 0;
%
%         %scorre tutti i picchi di velocità della prova in questione
%         for(j = 1:length(RigaProofVelocita)) %for SpeedPeak
%
%             index = RigaProofVelocita(j);
%
%             %% Picco di Velocità %%
%
%             %ascissa inizio picco di velocità
%             ascissaStartV = InfoPicchiVelocita(index,3);
%                 ascissaStartV = round(ascissaStartV*Fs)/Fs;
%
%             %ascissa picco max di velocità
%             ascissaV = InfoPicchiVelocita(index,6);
%                 ascissaV = (round(ascissaV*Fs)/Fs);
%
%             %valore max picco di velocità
%             maxV = InfoPicchiVelocita(index,5);
%
%
%             %intervallo di durata picco di velocità
%             intervalloV = [InfoPicchiVelocita(index,3):...
%                            1/Fs:...
%                            InfoPicchiVelocita(index,3) + InfoPicchiVelocita(index,4)- 1/Fs];
%             intervalloV = (round(intervalloV*Fs)/Fs);
%             %%
%
%
%             %% Picchi di Forza %%
%
%             %ascisse picchi di Forza per la prova in questione
%             ascisseFORZE = InfoPicchiForza(RigaProofForza,6);
%
%             %ascisse INIZIO picchi di Forza per la prova in questione
%             ascisseStartFORZE = InfoPicchiForza(RigaProofForza,3);
%                 %arrotondamento per evitare problemi nel successivo
%                 %confronto
%                 ascisseStartFORZE = (round(ascisseStartFORZE*Fs)/Fs);
%
%             %ascisse FINE picchi di Forza per la prova in questione
%             ascisseEndFORZE = InfoPicchiForza(RigaProofForza,3)+...
%                               InfoPicchiForza(RigaProofForza,4)-...
%                               1/Fs;
%                 %arrotondamento per evitare problemi nel successivo
%                 %confronto
%                 ascisseEndFORZE = (round(ascisseEndFORZE*Fs)/Fs);
%             %%
%
%
%
%
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             %controllo che nella durata del primo picco di forza compaia
%             %un solo picco relativo e significativo di velocità
%                 ascissaStartVfirst = ascissaStartV + 0.0011;
%                 if(ascissaStartVfirst  ==  1/Fs + 0.0011)
%                     clear ascissaStartVfirst;
%                     y_V = cell2mat(  y_Vx(i)  );
%                     %prende i valori di V dell'intervallo
%                     valueControlloV = y_V( round(intervalloV*Fs));
%
%                     %trova picchi distanziati di almeno dist (numero di
%                     %elementi)
%               % N.B % per poter usare findpeaks_GUI devo utilizzare una soglia >= 0;
%                     % quindi sono costretto a considerare il valore delle V
%                     % nell'intervallo >0 (ecco perchè c'è un - davanti)
%                     if length(valueControlloV)>2
%                         [pks,locs] = findpeaks_GUI(-valueControlloV);
%                     else
%                         pks = 0;
%                     end
%
%                     if length(pks) > 1
%                         %faccio shiftare il vettore in su
%                         locsSHIFT = circshift(locs',-1)';
%                         pksSHIFT =  circshift(pks',-1)';
%
%                         distPeak = (locs - locsSHIFT)';
%                         dist = round(1/4* length(intervalloV));
%                         %tengo conto solo di picchi lontani tra loro di almeno dist
%                         %i valori 0 vogliono dire che il picco successivo è
%                         %vicino al picco attuale
%                         locsBINARY = (abs(distPeak) >= dist);
%                         index0_locsBINARY = find(locsBINARY == 0);
%
%                         diffAmpiezzaPeak = (pks - pksSHIFT)';
%                         %i valori 0 vogliono dire che il picco successivo è più
%                         %basso del picco attuale
%                         pksBINARY = abs(diffAmpiezzaPeak)>0;
%
%                         BINARY = ones(length(pks),1);
%                         if(pksBINARY(index0_locsBINARY) == 0)
%                            BINARY(index0_locsBINARY) = 0;
%                         elseif (index0_locsBINARY+1)<= length(BINARY)
%                            BINARY(index0_locsBINARY+1) = 0;
%                         end
%
%                         locs = locs .* BINARY' *(length(locsBINARY)>2) +...
%                                locs .* (length(locsBINARY)<=2);
%                         %cancello i picchi vicini più piccoli
%                         pks (find(locs==0))=[];
%                         locs(find(locs==0))=[];
%
%                         %ultimo controllo sui picchi rimasti; anche se
%                         %distanziati devono essere abbastanza grandi (rispetto alla media di tutti i picchi della prova in questione
% %                         BINARY_AMPLITUDE =     pks >= max(pks)/2;
%                         BINARY_AMPLITUDE =     pks >= -mean( InfoPicchiVelocita(RigaProofVelocita,5) )/2;
%                         if length(BINARY_AMPLITUDE)>0
%                             pks = pks.* BINARY_AMPLITUDE;
%                             locs(find(pks==0))=[];
%                             pks (find(pks==0))=[];
%                         end
%
%
%                        if length(locs)>1 &&...
%                           valueControlloV(1)<1/10*maxV &&...
%                           pks(1) == max(pks) &&...
%                           abs(pks(2)/pks(1))>= 20/100 %il secondo picco deve essere almeno il 20% del primo
%                                 %trova minimo (che dovrebbe separare i due o più semipicchi)
%                                 % che sia almeno a metà ampiezza media
%                            %N.B.% in questo caso lavoro con i valori negativi di velocità, quindi il mio massimo sarà il valore
%                                 % che più si avvicina a 0 (quello che posso chiamare minimo del picco)
%                                 [pks2,locs2] = findpeaks_GUI(valueControlloV);
%
%                                 if length(locs2)>0
% %                                     pks2BINARY = pks2>maxV/2 ;
%                                     pks2BINARY = pks2 >= mean( InfoPicchiVelocita(RigaProofVelocita,5) )/2;  
%                                     pks2 = pks2 .* pks2BINARY;
%                                     locs2(find(pks2==0))=[];
%                                     pks2 (find(pks2==0))=[];
%                                
%                                     if(length(locs2)>=1)
%                                         intervalloUtile = [locs(1):1:2*dist];
%                                         ascix = find(intervalloUtile == locs2(1));
%                                         InfoPicchiVelocita(index,5) = (length(ascix)>=1)* -pks(2) + (length(ascix)== 0)* -pks(1);
%                                         InfoPicchiVelocita(index,6) = (length(ascix)>=1)* locs(2)/Fs + (length(ascix)== 0)* locs(1)/Fs;
%                                         
%                                         %GLI EVENTUALI NUOVI VALORI SONO RIACQUISITI
%                                         %ascissa picco max di velocità
%                                         ascissaV = InfoPicchiVelocita(index,6);
%                                             ascissaV = (round(ascissaV*Fs)/Fs);
% 
%                                         %valore max picco di velocità
%                                         maxV = InfoPicchiVelocita(index,5);
%                                         
%                                         
%                                         
%                                     end
%                                 end
%                         end
%                     end
%                 end
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%                    
%             Control3 = 1;
%             CheckBef = 0;
%             %Controllo 1:
%             % picco di velocità deve iniziare dopo l'inizio del picco di forza
%             BeforeOK = find( (ascisseStartFORZE - ascissaStartV)<=0 );
%             
%            %Controllo 1
%             if ~isempty(BeforeOK) %if Controllo1
%                 
%                %trovo l'inizio del picco di forza più vicino all'inizio
%                %del picco di velocità
%                difFSP = (abs(ascisseStartFORZE(BeforeOK) - ascissaStartV));
%                indFok = find(difFSP == min(difFSP));
%                
%                %Controllo 2:
%                % inizio picco di velocità compare nella durata del picco di
%                % forza selezionato
%                
%                %intervallo di durata picco di forza (doppia)
%                percMore = 2;
%                intervalloF = [InfoPicchiForza(RigaProofForza(BeforeOK(indFok)),3):...
%                               1/Fs:...
%                               InfoPicchiForza(RigaProofForza(BeforeOK(indFok)),3) + InfoPicchiForza(RigaProofForza(BeforeOK(indFok)),4)...
%                               * percMore...
%                               - 1/Fs];
%                intervalloF = (round(intervalloF*Fs)/Fs); 
%                if  logical(length(   nonzeros( ismember(ascissaStartV, intervalloF) ) )) && (RigaProofForza(BeforeOK(indFok)) ~= indexPeakOK_PREVIOUS)
%                    
%                    %inserisco nel buffer i dati relativi a questo picco selezionato
%                    InfoPicchiForzeBuffer(RigaProofForza(BeforeOK(indFok)),:) = InfoPicchiForza(RigaProofForza(BeforeOK(indFok)),:);
%                    indexPeakOK_PREVIOUS = RigaProofForza(BeforeOK(indFok));
%                    indexSpeedOK_PREVIOUS = index;
%                    CheckBef = 1;
%                    Control3 = 1;
%                
%                else
%                    
%                    % passo alla condizione 3
%                     Control3 = 1;
% %                     InfoPicchiVelocita(index,:) = 0;
% %                     controlloV = 1;                    
%                end
%                
%             else %else if Controllo1
%                 Control3 = 1;
%                 
%             end %end if Controllo1         
%           
%            %non valido Controllo 1 o Controllo 2
%             % il max picco di velocità deve iniziare dopo l'inizio del picco di forza
%            if Control3 == 1 %if Control3
%                 
% %                 BeforeMaxOK = find( (ascisseStartFORZE - ascissaV)<=0 );
%                 BeforeMaxOK = find( (ascisseStartFORZE - ascissaV)<=0 );
%                 
%                 %Controllo 1bis
%                 if ~isempty(BeforeMaxOK) %if Controllo1bis
%                     
%                     %trovo l'inizio del picco di forza più vicino al max
%                     %del picco di velocità
%                     difFSP = (abs(ascisseStartFORZE(BeforeMaxOK) - ascissaV));
%                     indFok = find(difFSP == min(difFSP));
%                    
%                     %intervallo di durata picco di forza + una percentuale
%                     %variabile di aumento dell'intervallo di forza percMore
%                     % il doppio 
%                     %(7Jan2015)
%                     percMore = 2;
%                     intervalloF = [InfoPicchiForza(RigaProofForza(BeforeMaxOK(indFok)),3):...
%                                    1/Fs:...
%                                    InfoPicchiForza(RigaProofForza(BeforeMaxOK(indFok)),3) + InfoPicchiForza(RigaProofForza(BeforeMaxOK(indFok)),4)... 
%                                    * percMore...
%                                    - 1/Fs];
%                     intervalloF = (round(intervalloF*Fs)/Fs); 
%                     
%                     
%                     % il max del picco di velocità compare nella durata del picco di
%                     % forza selezionato
%                     if  logical(length(   nonzeros( ismember(ascissaV, intervalloF) ) ))
%                         
%                         %se indexPeakOK_PREVIOUS è diverso da 0
%                         if indexPeakOK_PREVIOUS ~=0 %if indPRE
%                             
%                             %se è stato trovato un picco precedente all'inizio del picco di velocità 
%                             if CheckBef == 1 %if CheckBef
%                                 
%                                 %%% se è lo stesso picco trovato nel controllo di prima, va bene 
%                                 if RigaProofForza(BeforeMaxOK(indFok)) == indexPeakOK_PREVIOUS
%                                     
%                                     %inserisco nel buffer i dati relativi a questo picco selezionato
%                                     InfoPicchiForzeBuffer(RigaProofForza(BeforeMaxOK(indFok)),:) = InfoPicchiForza(RigaProofForza(BeforeMaxOK(indFok)),:);
%                                     indexPeakOK_PREVIOUS = RigaProofForza(BeforeMaxOK(indFok));
%                                     indexSpeedOK_PREVIOUS = index;
%                                 
%                                 %%% altrimenti controllo che il nuovo picco trovato sia maggiore di quello trovato nel controllo precedente
%                                 elseif abs(InfoPicchiForzeBuffer(RigaProofForza(BeforeMaxOK(indFok)),5)) >  abs(InfoPicchiForzeBuffer(indexPeakOK_PREVIOUS,5))
%                                     
%                                     %inserisco nel buffer i dati relativi a questo picco selezionato
%                                     InfoPicchiForzeBuffer(indexPeakOK_PREVIOUS,:) = -1;
%                                     InfoPicchiForzeBuffer(RigaProofForza(BeforeMaxOK(indFok)),:) = InfoPicchiForza(RigaProofForza(BeforeMaxOK(indFok)),:);
%                                     indexPeakOK_PREVIOUS = RigaProofForza(BeforeMaxOK(indFok));
%                                     indexSpeedOK_PREVIOUS = index;
%                                     
%                                 else
%                                     %non fare nulla
%                                 end
%                                 
%                             else %else if CheckBef
%                                 
%                                 %%% se è lo stesso picco trovato prima
%                                 if RigaProofForza(BeforeMaxOK(indFok)) == indexPeakOK_PREVIOUS
%                                     
%                                     %%% controllo sul picco di velocità tra l'attuale e il precedente (uno verrà cancellato)
%                                      %maxV attuale
%                                     if abs(maxV) >= abs(InfoPicchiVelocita(indexSpeedOK_PREVIOUS,5))
%                                         
%                                         InfoPicchiVelocita(indexSpeedOK_PREVIOUS,:)=0;
%                                         controlloV = 1;
%                                     else
%                                         InfoPicchiVelocita(index,:)=0;
%                                         controlloV = 1;
%                                     end
%                                     
%                                 else
%                                     %inserisco nel buffer i dati relativi a questo picco selezionato
%                                     InfoPicchiForzeBuffer(RigaProofForza(BeforeMaxOK(indFok)),:) = InfoPicchiForza(RigaProofForza(BeforeMaxOK(indFok)),:);
%                                     indexPeakOK_PREVIOUS = RigaProofForza(BeforeMaxOK(indFok));
%                                     indexSpeedOK_PREVIOUS = index;                                    
%                                 end
%                               
%                                 
%                             end %end if CheckBef
%                                                              
%                         
%                         else %else if indPRE  
%                             
% %                             if (intervalloF(1)) == 1/Fs
% %                                 
% %                             else
%                                 %inserisco nel buffer i dati relativi a questo picco selezionato
%                                 InfoPicchiForzeBuffer(RigaProofForza(BeforeMaxOK(indFok)),:) = InfoPicchiForza(RigaProofForza(BeforeMaxOK(indFok)),:);
%                                 indexPeakOK_PREVIOUS = RigaProofForza(BeforeMaxOK(indFok));
%                                 indexSpeedOK_PREVIOUS = index;
% %                             end
%                         
%                         end %end if indPRE
%                      
%                         
%                     else
% 
%                         % cancella la riga relativa al picco di velocità non valido
%                         InfoPicchiVelocita(index,:) = 0;
%                         controlloV = 1;
% 
%                     end
%                 else %else if Controllo1bis
%                     
%                     % cancella la riga relativa al picco di velocità non valido
%                     InfoPicchiVelocita(index,:) = 0;
%                     controlloV = 1;
% 
%                 end %end if Controllo1bis
%                        
%             end %end %if Control3
%                   
%         end %for SpeedPeak
%    end %for num_prove
%    
% 
%    %%%%%%%%%%%% FORZE  controllo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    %%% Picchi che vanno totalmente eliminati perchè cadono in prossimità di
%    %%% un solo picco di velocità ma hanno un intensità troppo piccola
%    %%% rispetto a una soglia di attrito ThdMinF
%    if ~isempty(ThdMinF)
%        ClearTh = find( InfoPicchiForzeBuffer(:,5)> ThdMinF);
%        if ~isempty(ClearTh)
%             InfoPicchiForzeBuffer(ClearTh,:) = 0;
%        end
%    end
%    
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%    %%%%%%%%%%%%% FORZE primo controllo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    %%% Picchi che vanno totalmente eliminati perchè cadono in prossimità di
%    %%% un solo picco di velocità per cui è stato già rilevato un picco di
%    %%% forza
%    %cerca righe == -1
%    ClearZERO = find(InfoPicchiForzeBuffer(:,2)==-1);
%    delToT = ClearZERO;
%    %cancella righe == -1
%    InfoPicchiForzeBuffer(ClearZERO,:)=[];
%    
%    %numera nuovamente i picchi
%    InfoPicchiForzeBuffer(:,1)= [1:size(InfoPicchiForzeBuffer,1)]';
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%    
%    
%    %%%%%%%%%%%%% FORZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    %cerca righe == 0
%    ClearZERO = find(InfoPicchiForzeBuffer(:,2)==0);
%    %cancella righe == 0
%    InfoPicchiForzeBuffer(ClearZERO,:)=[];
%    
%    %numera nuovamente i picchi
%    InfoPicchiForzeBuffer(:,1)= [1:size(InfoPicchiForzeBuffer,1)]';
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  
% 
% 
%    
%    
%    
%    %%%%%%%%%%%%% VELOCITA'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    controlloV = 0;
%    if(controlloV == 1)
%        %cerca righe == 0
%        ClearZERO = find(InfoPicchiVelocita(:,2)==0);
%        %cancella righe == 0
%        InfoPicchiVelocita(ClearZERO,:)=[];
% 
%        %numera nuovamente i picchi
%        InfoPicchiVelocita(:,1)= [1:size(InfoPicchiVelocita,1)]';
%    end
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%    InfoPicchiVelocita = InfoPicchiVelocitaOLD;
% end
%    
%    
%    
%         
%         