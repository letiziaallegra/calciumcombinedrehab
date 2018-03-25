%centroid_control
load('Centroid_control');

%centroid_stroke
load('Centroid_stroke');

%centroid_rehab
load('Centroid_rehab');


days_cs    = [1:5];
days_rehab = [16:20];

for d=1:5
    
    %control - median pos
    x_cc_m   = median(Centroid_control(  Centroid_control(:,1) == days_cs(d),2));
    x_cc_std = std(Centroid_control(  Centroid_control(:,1) == days_cs(d),2),[]);
    y_cc_m   = median(Centroid_control(  Centroid_control(:,1) == days_cs(d),3));
    y_cc_std = std(Centroid_control(  Centroid_control(:,1) == days_cs(d),3),[]);
    %stroke - median pos
    x_cs_m   = median(Centroid_stroke(  Centroid_stroke(:,1) == days_cs(d),2));
    x_cs_std = std(Centroid_stroke(  Centroid_stroke(:,1) == days_cs(d),2),[]);
    y_cs_m   = median(Centroid_stroke(  Centroid_stroke(:,1) == days_cs(d),3));
    y_cs_std = std(Centroid_stroke(  Centroid_stroke(:,1) == days_cs(d),3),[]);
    %rehab - median pos
    x_cr_m   = median(Centroid_rehab(  Centroid_rehab(:,1) == days_rehab(d),2));
    x_cr_std = std(Centroid_rehab(  Centroid_rehab(:,1) == days_rehab(d),2),[]);
    y_cr_m   = median(Centroid_rehab(  Centroid_rehab(:,1) == days_rehab(d),3));
    y_cr_std = std(Centroid_rehab(  Centroid_rehab(:,1) == days_rehab(d),3),[]);
    
    %all days
    x_cc_m_tot(d)   =  x_cc_m;
    x_cc_std_tot(d) =  x_cc_std;
    y_cc_m_tot(d)   =  y_cc_m;
    y_cc_std_tot(d)  = y_cc_std;
    
    x_cs_m_tot(d)   = x_cs_m;
    x_cs_std_tot(d) = x_cs_std;
    y_cs_m_tot(d)   = y_cs_m;
    y_cs_std_tot(d) = y_cs_std;
    
    x_cr_m_tot(d)   = x_cr_m;
    x_cr_std_tot(d) = x_cr_std;
    y_cr_m_tot(d)   = y_cr_m;
    y_cr_std_tot(d) = y_cr_std;    
    
end


%%%%%
%%% Figure 1 %%%%%%%%
H_Main_Cent_Scatter = figure('Name',['Main Centroids']);
hold on
for itre=1:3
    
    if itre==1
        subplot(131)
        x  = x_cc_m_tot;
        y  = y_cc_m_tot;
        xe = x_cc_std_tot;
        ye = y_cc_std_tot;
    elseif itre==2
        subplot(132)
        x  = x_cs_m_tot;
        y  = y_cs_m_tot;
        xe = x_cs_std_tot;
        ye = y_cs_std_tot;
    elseif itre==3
        subplot(133)
        x  = x_cr_m_tot;
        y  = y_cr_m_tot;
        xe = x_cr_std_tot;
        ye = y_cr_std_tot;
    end


    %scatter with errors
    nD = length(x);
    %Make these defaults later:
    dotColor = [1 0.3 0.3]; % conservative pink
    yeColor = [0, 0.4, 0.8]; % bright navy blue
    xeColor = [0.35, 0.35, 0.35]; % not-too-dark grey
    dotSize = 23;
    % set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    set(gca, 'FontSize', 10);
    hold all;
    
  
    
    for i = 1:nD
        plot([(x(i) - xe(i)) (x(i) + xe(i))], [y(i) y(i)], 'Color', xeColor);
        plot([x(i) x(i)], [(y(i) - ye(i)) (y(i) + ye(i))], 'Color', yeColor);
    end
  scatter(x, y, dotSize, repmat(dotColor, nD, 1));
    set(gca,'Ydir','reverse')
    t = [1:5]'; b = num2str(t); Na = cellstr(b);
    % t = [1:9,11:20]'; b = num2str(t); Na = cellstr(b);
    dx = 0.2; dy = 0.2;
    text(x+dx, y+dy, Na);
    xlim([0 512])
    ylim([0 512])
end



















