function filt_image = Med_Filter(image,max_winsize)
%     if filter == "ramf"
        filt_image = RAMF(image,max_winsize);
%     elseif filter == "samf"
%         filt_image = SAMF(image,max_winsize);
%     else
%         throw("Filter not found");
%     end
end

function filt_image = RAMF(image,max_winsize)
    [r,c,l] = size(image);
    filt_image = image;
    for x = 1:l
        %Median filtering per slice
        for i = 1:r
            for j = 1:c
                win_s = min([max_winsize,i-1,j-1,r-i,c-j]);
                %first level
                for W = 0:win_s
                    S = image(i-W:i+W,j-W:j+W,x);
                    xmed = median(S,'all');
                    xmin = min(S,[],'all');
                    xmax = max(S,[],'all');
                    
                    Tminus = xmed - xmin;
                    Tplus = xmax - xmed;
                    if (Tminus > 0 && Tplus > 0)
                        break
                    end
                end
                %second level
                Uminus = image(i,j,x)-xmin;
                Uplus = xmax-image(i,j,x);
                if not (Uminus > 0 && Uplus > 0)
                    filt_image(i,j,x) = xmed;
                end
            end
        end
    end
end
% function filt_image = SAMF(image,max_winsize)
%     [r,c,l] = size(image);
%     filt_image = image;
%     
%     for x = i:l
%         for i = 1:r
%             for j = 1:c
%                 
% end