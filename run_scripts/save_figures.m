function ret = save_figures(filename)
hfigs = get(0, 'children');                       %Get list of figures

for m = 1:length(hfigs);
    figure(hfigs(m));                          %Bring Figure to foreground
    if strcmp(filename, '0');                    %Skip figure when user types 0
        continue;
    else
        saveas(hfigs(m), strcat(filename, '/', num2str(m), '.fig')); %Matlab .FIG file
        saveas(hfigs(m), strcat(filename, '/', num2str(m), '.png')); %Standard PNG graphics file (best for web)
    end
end