function ret = savefigs_v2(Type)
% This function allows you to quickly save all currently open figures with
% a custom filename for each in multiple formats.  To use the function
% simply call savefigs with no arguments, then follow the prompts
%
% Upon execution this function will one-by-one bring each currently open
% figure to the foreground.  Then it will supply a text prompt in the main
% console window asking you for a filename.  It will save that figure to
% that filename in the .fig, .emf, .png, and .eps formats.  
%
% The formats that it saves in can be changed by commenting out or adding
% lines below.
%
% Copyright 2010 Matthew Guidry 
% matt.guidry ATT gmail DOTT com  (Email reformatted for anti-spam)

hfigs = get(0, 'children')                          %Get list of figures

for m = 1:length(hfigs)
    
    figure(hfigs(m))                                %Bring Figure to foreground
    
    set(hfigs(m),'position',get(0,'screensize'))
    
    filename = input('Filename? (0 to skip)\n', 's')%Prompt user
    if strcmp(filename, '0')                        %Skip figure when user types 0
        continue
    else
        for i=1:length(Type)
            
            if strcmp(Type{i},'MAT')
                saveas(hfigs(m), [filename '.fig']) %Matlab .FIG file
            elseif strcmp(Type{i},'EMF')
                saveas(hfigs(m), [filename '.emf']) %Windows Enhanced Meta-File (best for powerpoints)
            elseif strcmp(Type{i},'PNG')
                saveas(hfigs(m), [filename '.png']) %Standard PNG graphics file (best for web)
            elseif strcmp(Type{i},'EPS')
                eval(['print -depsc2 ' filename])   %Enhanced Postscript (Level 2 color) (Best for LaTeX documents)
            end
        end
        
    end
end