
function PLOTsync()
[datafile datapath] = ...
        uigetfile({'*.txt','Data Files (*.txt)'; ...
        '*.*', 'All Files (*.*)'}, 'Select the Data File')

cd(datapath)

txt = load(datafile);
% datafile(length(datafile)-3:end)=[];
% txt = datafile

%%FORCES%
figure
% plot(txt(:,1),-(sgolayfilt(txt(:,2),3,21)-median(sgolayfilt(txt(:,2),3,21))))
hold on
% plot(txt(:,1),-(sgolayfilt(txt(:,3),3,21)-median(sgolayfilt(txt(:,3),3,21))),'r')
% hold on
% plot(txt(:,1),(sgolayfilt(txt(:,4),3,21)-median(sgolayfilt(txt(:,4),3,21))),'g')

%%KINEMATICS%
% figure
plot(txt(:,1),txt(:,7)/50,'y')
hold on
plot(txt(:,1),txt(:,8),'g')
hold on
% plot(txt(:,1),sgolayfilt(txt(:,9),3,21)/50,'b')
% plot(txt(:,1),txt(:,9)/25,'c')

%%STATUS AND POSITION%
% hold on
% plot(txt(:,1),txt(:,7)+10,'b')
% hold on
% plot(txt(:,1),txt(:,8),'r')

end
