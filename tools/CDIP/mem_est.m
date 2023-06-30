
function [mem_out,chk] = mem_est(sp_data)

    % mem_out.freq_bands = sp_data.freq;
    % mem_out.dir_bands = MEM_dir_bands

    for i = 1 : 1 : length(sp_data.freq)
        mem_out.freq(i) = sp_data.freq(i);
    %     mem_out.band_width(i) = sp_data.band_width(i); % should be difined in sp_data

        if (sp_data.dir(i) ~= -1) 
            [d,chk] = mem_cdip(sp_data.a1(i), sp_data.a2(i), sp_data.b1(i), sp_data.b2(i));
            % merge into 5 deg directional bins, multiply by 0.2 to get units of m^2/Hz-deg
        end

        for j = 5 : 5 : 355
            if sp_data.dir(i) ~= -1
                mem_out.ds(j/5,i) = 0.2*sp_data.ener_dens(i)*(d(j-2)+d(j-1)+d(j)+d(j+1)+d(j+2));
            else
                mem_out.ds(j/5,i) = 0;
            end
        end

        if sp_data.dir(i) ~= -1
            mem_out.ds(72,i) = 0.2*sp_data.ener_dens(i)*(d(358)+d(359)+d(360)+d(1)+d(2));
        else
            mem_out.ds(72,i) = 0;
        end
    end

    for k = 5 : 5 : 360
        mem_out.dir(k/5) = k;
    end 
end

