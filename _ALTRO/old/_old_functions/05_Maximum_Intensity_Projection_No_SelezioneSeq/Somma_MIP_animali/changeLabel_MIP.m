%change Label
function changeLabel_MIP(ResX,ResY,bregmaRef)

    DimX = 4.4;%[mm];
    DimY = 4.4;%[mm];

    %resize X 1tick == 0.500 mm
    xTick_mm      = 0.5;%[mm] 
    xTick_px      = (ResX*xTick_mm)/DimX;
    xTick_mm_List = [xTick_mm: xTick_mm: DimX];
    xTick_mm_List_Label = num2cell(xTick_mm_List)   
    set(gca,'XTick',[xTick_px:xTick_px:ResX])
    set(gca,'XTickLabel',xTick_mm_List_Label)
    
    
    %resize X 1tick == 0.500 mm
    yTick_mm            = 0.5;%[mm] 
    yTick_px            = (ResY*yTick_mm)/DimY;
    yTick_mm_List       = [yTick_mm: yTick_mm: DimY];
    
    if nargin>2
        if bregmaRef
            yTick_mm_List = yTick_mm_List -0.25;
        end
    end
        
    yTick_mm_List_Label = num2cell(yTick_mm_List);
    
    set(gca,'YTick',[yTick_px:yTick_px:ResY])
    set(gca,'YTickLabel',yTick_mm_List_Label)
        

end