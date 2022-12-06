clear
close all

%%%%%%%%%%%%
% Commands %
%%%%%%%%%%%%
cmd_write_ew    = 1;                % write EW to file; 
cmd_write_rv    = 0;                % write RV to file
cmd_write_ew_rv = 0;
cmd_write_asym  = 0;                % write measurement of asymmetry of lines to file
cmd_subplot     = 0;                % Plot all panels in one
cmd_scatter     = 0;                % scatter plot of EW - RV
cmd_scatter2    = 0;                % scatter plot of ASYM - RV
cmd_eps         = 1 - cmd_subplot;  % Individual plots

%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up working region %
%%%%%%%%%%%%%%%%%%%%%%%%%
c           = 299792458.0;
H_ALPHA     = 6562.808;
Telluric    = 6564.2585;
width       = 5;                % tunable. 5 is better for most stars and larger numbers for stars with wide H-alpha absorption lines
SPACING     = 0.01;             % spacing of interpolated data (raw data spaced by ~0.04)
WAV_MIN     = H_ALPHA - width;
WAV_MAX     = H_ALPHA + width;
w_low       = 0.5;              % width to determine H-alpha region; not tunable for best performance
w_low_num   = w_low / SPACING;

%%%%%%%%%%%%%%%%%%%%
% Read directories %
%%%%%%%%%%%%%%%%%%%%
dir_list    = regexp(fileread('dir.list'), '\n', 'split');
N_dir       = size(dir_list, 2) - 1;
% flag_list   = regexp(fileread('scatter_flag.list'), '\n', 'split'); 
% N_flag      = size(flag_list, 2) - 1;


%%%%%%%%%%%%%%%%
% Read RV data %
% %%%%%%%%%%%%%%%%
% allvels         = importdata('allvels.dat');
% % allvels         = importdata('RV_76920.dat');
% allvels_JD      = allvels(:,1);
% allvels_RV      = allvels(:,2);
% allvels_d_RV    = allvels(:,3);

% for w           = [0.5, 0.8, 1.0, 1.5];                % width to calculate EW of H-alpha; tunable
% for w = [1];

%     std_RV      = zeros(1, N_dir);
%     std_EW      = zeros(1, N_dir);
%     std_ASYM    = zeros(1, N_dir);
% 
%     %%%%%%%%%%%%%%%%%%
%     % Loop for stars %
%     %%%%%%%%%%%%%%%%%%
w = [1];
for n_dir = 1 : N_dir  % loops for each star
    
    telluric_min    = Telluric;
    telluric_max    = Telluric;

    %%%%%%%%%%%%
    % Commands %
    %%%%%%%%%%%%
    cmd_publ = 1;

    %%%%%%%%%%%%%%%%%%%
    % File management %
    %%%%%%%%%%%%%%%%%%%
    if cmd_write_ew
        ew_name     = [dir_list{n_dir}(1:end-1), '.ew'];
        fid         = fopen(ew_name, 'w');
    end

    if cmd_write_rv     % collect RVs
        rv_name     = [dir_list{n_dir}(1:end-1), '.vels'];
        sys_name    = [dir_list{n_dir}(1:end-1), '.sys'];
        fid2        = fopen(rv_name, 'w'); 
        fid3        = fopen(sys_name, 'w'); 
    end

    if cmd_write_ew_rv
        ew_rv_name  = [dir_list{n_dir}(1:end-1), '.ewrv'];
        fid4        = fopen(ew_rv_name, 'w');
    end

    if cmd_write_asym
        asym_name   = [dir_list{n_dir}(1:end-1), '.asym'];
        fid5        = fopen(asym_name, 'w');
    end

    % h = figure;
    % set(h, 'Visible', 'off');
    % hold on
    output      = strrep(dir_list{n_dir}(1:end-1), '_', '');
    fprintf([output, '\n']);
%     list        = [dir_list{n_dir}, 'ewha_in.list'];
%     file_list   = regexp( fileread(list), '\n', 'split');
    file_list   = dir(['./', dir_list{n_dir}]);
    file_list  = file_list(4:end);

    % Remove calculations on 20100703, 20100704 and 20100130 because the flux is systematically offset %
    N       = size(file_list, 1);
    N_file  = size(file_list, 1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find the highest S/N spectrum as a template (in order 19) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For now, the shift of the spectrum is done in wavelength space, under
    % the assumption that the logarithm of small wavelength region is
    % linear. 

    max_SN      = 0;
    for n_file = 1 : N_file
        file    = [dir_list{n_dir}, file_list(n_file).name];
        A       = regexp( fileread(file), '\n', 'split');
        SN      = str2num(A{7});  % S/N of order 19
        SN      = SN(2)^0.5;
        if SN > max_SN
            max_SN  = SN;   % stores the highest S/N
            max_A   = A;    % stores the highest S/N spectrum file
        end
    end
    [~, ~, ~, max_int_nor]      = rep_ord(max_A, WAV_MIN, WAV_MAX);
    max_int_nor(isnan(max_int_nor)) = 1;    % will become 0 in cross correlation

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % cross correlate all other spectra with above --> shift %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h_stack         = figure;
    hold on
    h_bar           = waitbar(0, 'Stacking...');         % Progress bar
    sum_SN          = 0;
    wav_arr         = (WAV_MIN : SPACING : WAV_MAX)';
    int_tmp         = zeros(2 * width / SPACING + 1, 1);    

    %%%%%%%%%%%%%%%%%%%%%
    % Loops for spectra %
    %%%%%%%%%%%%%%%%%%%%%
    len         = 2 * width / SPACING + 1;      % size of wav, int, err, int_nor arrays
    wav         = zeros(len, N_file);
    int         = zeros(len, N_file);
    err         = zeros(len, N_file);
    int_nor     = zeros(len, N_file);
    bary_JD     = zeros(1, N_file);
    EW          = zeros(1, N_file);
    delta_EW    = zeros(1, N_file);
    LPD         = zeros(1, N_file);
    delta_LPD   = zeros(1, N_file);

    for n_file = 1 : N_file

        telluric    = Telluric;
        waitbar((n_file - 1)/N_file);
        label3                              = 0;
        file                                = [dir_list{n_dir}, file_list(n_file).name];
        A                                   = regexp( fileread(file), '\n', 'split');
        [wav(:,n_file), int(:,n_file), err(:,n_file), int_nor(:, n_file)] = rep_ord(A, WAV_MIN, WAV_MAX);
        int_nor(isnan(int_nor(:,n_file)), n_file)   = 1;                                        % fake NaN
        int(isnan(int(:,n_file)), n_file)           = 0;                                        % fake NaN
        err(isnan(err(:,n_file)), n_file)           = mean(err(~isnan(err(:,n_file)), n_file)); % fake NaN

        % get barycentric RV 
        bc_RV       = str2double(A{3}(32:end));      
        % disp(n_file); disp(bc_RV)

        % get measured RV 
        bary_JD(n_file)     = str2double(A{1}(9:end));
%         idx_RV              = abs( allvels_JD - bary_JD(n_file)' ) < 0.001;
        JD                  = bary_JD(n_file);
        if_rv               = ~isempty(JD);

        if 0
            wav(:,n_file)   = wav(:,n_file) * (bc_RV / c + 1);
            telluric        = telluric * (bc_RV / c + 1);
            telluric_min    = min(telluric_min, telluric);
            telluric_max    = max(telluric_max, telluric); 
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cross Correlation to %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% determine the shift  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%
        if 1
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % 1st cross correlation %
            %%%%%%%%%%%%%%%%%%%%%%%%%
            [int_cor,lag]       = xcorr(1 - max_int_nor, 1 - int_nor(:,n_file));
            shift1              = lag(int_cor == max(int_cor));
            % Shift the flux grid %
            if shift1 < 0           % shift to left
                int_nor(:, n_file)  = [int_nor(abs(shift1)+1 : end, n_file)', 1 - zeros(1, abs(shift1))]';
                int(:, n_file)      = [int(abs(shift1)+1 : end, n_file)', zeros(1, abs(shift1))]';
                err(:, n_file)      = [err(abs(shift1)+1 : end, n_file)', mean(err(:, n_file)) * ones(1, abs(shift1))]';
            elseif shift1 > 0       % shift to right
                int_nor(:, n_file)  = [1 - zeros(1, abs(shift1)), int_nor(1 : end - abs(shift1), n_file)']';
                int(:, n_file)      = [zeros(1, abs(shift1)), int(1 : end - abs(shift1), n_file)']';
                err(:, n_file)      = [mean(err(:, n_file)) * ones(1, abs(shift1)), err(1 : end - abs(shift1), n_file)']';
            end
            telluric                = telluric + shift1 * SPACING;

            %%%%%%%%%%%%%%%%%%%%%%%%%
            % 2nd cross correlation %
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % fitting a Gaussian for precision shift performs poorly, as shift2 may > 0.01,
            % which violates the intended subpixel precision
            % reason: lines are not symmetric

            [int_cor,lag]   = xcorr(1 - max_int_nor, 1 - int_nor(:,n_file));
            idx             = (-5 < lag) & (lag < 5);
            lag             = lag(idx);
            int_cor         = int_cor(idx);
%                 ft = fittype('-a*(x-b)^2+c', 'independent','x', 'coefficients',{'a','b','c'});
%                 f = fit( lag', int_cor, ft, 'StartPoint', [SPACING, 0, max(int_cor)] );
            f = fit( lag', int_cor, '-a*(x-b)^2+c', 'StartPoint', [SPACING, 0, max(int_cor)] );
            % disp(f);
            % disp(f.b);

            if abs(f.b) > 1
                disp('@@@@@@ Achtung!!! @@@@@@')
            end

            % Shift the wavelength grid %
            wav(:,n_file)   = wav(:,n_file) + f.b * SPACING;      % Achtung! Plus, aber nicht Minus.
            telluric        = telluric + f.b * SPACING;
            telluric_min    = min(telluric_min, telluric);
            telluric_max    = max(telluric_max, telluric);
        end

        % Plotting %
        axis([WAV_MIN WAV_MAX 0 1.1])
        plot(wav(:,n_file), int_nor(:,n_file))
        title(output)

        %%%%%%%%%%%%%%%%%%%%%%
        % Build the template %
        %%%%%%%%%%%%%%%%%%%%%%
        SN      = str2num(A{7});  % S/N of order 19
        SN      = SN(2)^0.5;        
        sum_SN              = sum_SN + SN^2;
        int_nor(:,n_file)   = interp1q(wav(:,n_file), int_nor(:,n_file), wav_arr);
        int_nor(isnan(int_nor(:,n_file)), n_file) = 1;
        int(:,n_file)       = interp1q(wav(:,n_file), int(:,n_file), wav_arr);
        int(isnan(int(:,n_file)), n_file) = 0;            
        int_tmp             = int_tmp + int_nor(:,n_file) * SN^2;
    end % for n_file = 1 : N_file

    % disp([telluric_min, telluric_max])

    waitbar(1)
    int_tmp = int_tmp / sum_SN;
    close(h_bar)
    hold off
    close(h_stack)

    n_nc        = 0;
    idx_rv      = 0;

    %%%%%%%%%%%%%%
    % Plot setup %
    %%%%%%%%%%%%%%
    h_show = figure;
    if cmd_subplot && cmd_publ
        subplot(3,2,1)
    elseif cmd_subplot
        subplot(2,2,1)
    else 
        ha = figure;
    end
    title(output)
    hold on
    axis([WAV_MIN WAV_MAX 0 1.1])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine the center of H-alpha lines %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MIN     = 10;                   % arbitrary value that MEAN cannot reach above
    for i = 1 : (len - w_low_num)   % 2Å width for w = 2
        MEAN    = mean(int_tmp(i:(i + w_low_num)));
        if MIN > MEAN
            MIN     = MEAN;
            idx_min = i;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine the working region %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wav_cen         = wav_arr(idx_min) + w_low/2;
    wav_left        = wav_arr(idx_min) + w_low/2 - w/2;
    wav_right       = wav_arr(idx_min) + w_low/2 + w/2;
    w_cutoff        = 0.2;
    telluric_left   = telluric_min - w_cutoff;
    telluric_right  = telluric_max + w_cutoff;
    % Eligible counts must lie in 
    % H-alpha region [wav_left, wav_right] but NOT in the 
    % telluric region [telluric_left, telluric_right]
%     idx_fit         = (wav_arr < wav_right) & (wav_arr > wav_left) & ~((wav_arr > telluric_left) & (wav_arr < telluric_right));@jzhao
    idx_fit         = (wav_arr < wav_right) & (wav_arr > wav_left);

    cal_scatter2 = 1;

    for n_file = 1 : N_file
        file            = [dir_list{n_dir}, file_list(n_file).name];

        %%%%%%%%%%%%%%%%%%%%%
        % Identify outliers %
        %%%%%%%%%%%%%%%%%%%%%            
        diff_abs            = int_nor(:, n_file) - int_tmp;
        diff_relative   = diff_abs ./ int_tmp;
        if_cosmic       = (abs(diff_relative) > 21 ./ err(:, n_file)) & idx_fit;  % 10 times sigma            
        % 20 is tunable. To be further tested. Failed for HD71933. To
        % be adjusted when applied to specific spectrum. 

        if_cosmic_org   = if_cosmic;
        if_cosmic_temp  = if_cosmic;

        % if_cosmic will become 1 where outlier region extends more
        % than 14 interpolated pixels (~3.5 pixels in CCD)
        dd = 8;         % tunable, to be further tested
        for d = -dd:dd
            if_cosmic   = circshift(if_cosmic_temp, d) & if_cosmic;
        end

        % recover those if_cosmic to 1 --> do NOT mask these region as
        % they do not appear to be cosmics
        if_cosmic_temp  = if_cosmic;
        for d = -dd:dd
            if_cosmic   = circshift(if_cosmic_temp, d) | if_cosmic;
        end
        if_cosmic_org(if_cosmic == 1) = 0;

        % idx_mask        = if_cosmic_org & (wav_arr < wav_arr(idx_min) + 0.55 + w/2) & (wav_arr > wav_arr(idx_min) + 0.45 - w/2);
        idx_mask        = if_cosmic_org;
        idx_mask0       = idx_mask;

        % if an outlier is identified, remove also the neighbouring
        % pixels, as the filter applies only the top part of the outlier
        for d = -6:6    % tunable, to be further tested
            idx_mask    = circshift(idx_mask0, d) | idx_mask;
        end

        idx_del     = (idx_mask == 1);
        wav_del     = wav_arr(idx_del);
        int_del     = int_nor(idx_del, n_file);
        wav_keep    = wav_arr(~idx_del);
        int_keep    = int_nor(~idx_del, n_file); 
        int_arr     = csaps(wav_keep, int_keep, 0.999999, wav_arr);
        diff2       = int_arr - int_tmp;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get asymmetry of the line (scaled by area) % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        idx_links               = (wav_arr < wav_cen) & idx_fit;
        idx_rechts              = (wav_arr > wav_cen) & idx_fit;

        if 1
            LPD(n_file)         = sum(diff2(idx_fit) .* (wav_cen - wav_arr(idx_fit)));
            delta_LPD(n_file)   = sqrt(sum((diff2(idx_fit) .* (wav_cen - wav_arr(idx_fit)) ./ err(idx_fit)).^2));
        end

        if 0
            LPD(n_file)        = sum(diff2(idx_fit));
            STD_LPD(n_file)    = std(diff2(idx_fit));
        end

        if sum(isnan(delta_LPD)) ~= 0
            cal_scatter2 = 0;
        end

        % Colour plots of stacked profile % 
        cc = jet(101);
        BJD_max = max(bary_JD);
        BJD_min = min(bary_JD);
        k = floor((bary_JD(n_file) - BJD_min) / (BJD_max - BJD_min) * 100+1);
        plot(wav_arr, int_arr, 'color', cc(k,:))
        plot(wav_arr, abs(diff_abs), 'k')
        plot(wav_del, int_del, 'x', 'markers', 8, 'color', cc(k,:))
        grey = [0.4,0.4,0.4];
        % plot([WAV_MIN, WAV_MAX], [cosmic, cosmic], '--', 'Color', grey, 'LineWidth', 2)      

        if w > 1
            plot([wav_left, wav_left], [0.7, 1.05], '--', 'Color', grey)
            plot([wav_right, wav_right], [0.7, 1.05], '--', 'Color', grey)
            plot([telluric_left, telluric_left], [0, 1.05], 'r--')
            plot([telluric_right, telluric_right], [0, 1.05], 'r--')
        elseif (w <= 1) 
            area([telluric_left telluric_right], [0 0], 'LineStyle', 'none', 'FaceColor','r', 'FaceAlpha', 0.02);
            area([telluric_left telluric_right], [1.05 1.05], 'LineStyle', 'none', 'FaceColor','r', 'FaceAlpha', 0.02);
            plot([wav_left, wav_left], [0, 0.6], '--', 'Color', grey)
            plot([wav_right, wav_right], [0, 0.6], '--', 'Color', grey)            
            %plot([telluric_left, telluric_left], [0, 1.05], 'r--')
            %plot([telluric_right, telluric_right], [0, 1.05], 'r--')                    
        end                



        %%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate H-alpha EW %
        %%%%%%%%%%%%%%%%%%%%%%%%
        wav_ha              = wav_arr(idx_fit);
        int_ha              = int_arr(idx_fit);
        err_ha              = err(idx_fit, n_file);
        EW(n_file)          = sum(1 - int_ha) * SPACING;
        delta_EW(n_file)    = sum(int_ha) * SPACING / mean(err_ha);
        % disp([EW(n_file), delta_EW(n_file)]);

        %%%%%%%%%%%%%%%%%
        % Write to file %
        %%%%%%%%%%%%%%%%%
        % Format: BJD, RV, delta_RV, EW, delta_EW
        % idx_write   = abs( allvels_JD - bary_JD(n_file) ) <
        % 0.000005; (original)
%         idx_write   = abs( allvels_JD - bary_JD(n_file) ) < 0.001;
        JD                  = bary_JD(n_file);
        label1      = isempty(JD);
        disp(label1);
        label2      = 1 - label1;
        label3      = label3 | label2;

        if label1    % cannot find corresponding RV --> label1
            % generate a fake RV as an identifier. 

            fprintf(['____', file, ': \t no RV found in allvels.dat\n']);
            n_nc    = n_nc + label1;

        else        % corresponding RV found --> label2 = 1

            if cmd_write_rv     % collect RVs
                fprintf(fid2, '%15.6f\t%10.2f\t%7.2f\n', bary_JD(n_file), allvels_RV(idx_write), allvels_d_RV(idx_write));
            end         

            idx_rv              = idx_rv + label2;
            bary_JD_2(idx_rv)   = bary_JD(n_file);
            LPD_2(idx_rv)       = LPD(n_file);
            STD_LPD_2(idx_rv)   = delta_LPD(n_file);
%             RV_2(idx_rv)        = allvels_RV(idx_write);
%             delta_RV_2(idx_rv)  = allvels_d_RV(idx_write);
            EW_2(idx_rv)        = EW(n_file);
            delta_EW_2(idx_rv)  = delta_EW(n_file);            

            if cmd_write_ew_rv
                fprintf(fid4, '%10.2f\t%10.7f\t%10.7f\n', allvels_RV(idx_write), EW(n_file), delta_EW(n_file));
            end                
        end

    end % for n_file = 1 : N_file (identify outliers)

    %%%%%%%%%%%%%%%%%%%
    % File management %
    %%%%%%%%%%%%%%%%%%%
    if cmd_write_ew
        % MEAN_EW = mean(EW);
        for n_file = 1 : N_file
            fprintf(fid, '%15.6f\t%10.7f\t%10.7f\n', bary_JD(n_file), EW(n_file), delta_EW(n_file));
        end
    end        

    if cmd_write_asym
        for n_file = 1 : N_file
            fprintf(fid5, '%15.6f\t%10.7f\t%10.7f\n', bary_JD(n_file), LPD(n_file), delta_LPD(n_file));
        end
    end               

    if cmd_write_rv
        fclose(fid2);
        fprintf(fid3, 'Data  {\n	RV[] "');
        fprintf(fid3, dir_list{n_dir}(1:end-1));
        fprintf(fid3, '.vels"\n}\n"');
        fprintf(fid3, dir_list{n_dir}(1:end-1));
        fprintf(fid3, '" {\n  Mass   1.00\n}\n');
        fclose(fid3);
    end

    if cmd_write_ew_rv
        fclose(fid4);
    end

    if ~label3 && cmd_write_rv
        delete(sys_name)
        delete(rv_name)
    end

    if ~label3 && cmd_write_ew_rv
        delete(ew_rv_name)
    end        

    if cmd_write_ew
        fclose(fid);
    end

    if cmd_write_asym
        fclose(fid5);
    end

    xlabel('Wavelength (\AA)', 'Interpreter', 'latex')
    ylabel('Nomalized flux', 'Interpreter', 'latex')
    % plot(wav_arr, int_tmp, 'b.')
    hold off 
    if cmd_eps
        out_eps = [output, '_a.eps'];
        print(out_eps, '-depsc')
        close(ha)
    end            


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Time series of activity proxies %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bary_JD = bary_JD - 2450000;
    if cmd_subplot && cmd_publ
        subplot(3,2,2)
    elseif cmd_subplot
        subplot(2,2,2)
    else
        hb = figure;
    end
    hold on
    yyaxis left
    errorbar(bary_JD, EW, delta_EW, 's', 'markers', 8, 'MarkerFaceColor', 'b')
    ylabel('EW (\AA)', 'Interpreter', 'latex')
    % for n_file = 1 : N_file
    %   plot(bary_JD(n_file), EW(n_file), '^', 'markers', 8, 'color', cc(n_file,:))
    % end
    y_min = min(EW - delta_EW);
    y_max = max(EW + 2 * delta_EW);
    y_min = 2 * y_min - y_max;
    ylim([y_min y_max])

    yyaxis right
    errorbar(bary_JD, LPD, delta_LPD, 'r.', 'markers', 20)
    ylabel('LPD', 'Interpreter', 'latex')
    % for n_file = 1 : N_file
    %     plot(bary_JD(n_file), LPD(n_file), '.', 'markers', 20, 'color', cc(n_file,:))
    % end
    y_min = min(LPD - 5 * delta_LPD);
    y_max = max(LPD + delta_LPD);
    y_max = 2 * y_max - y_min;
    ylim([y_min y_max])

    hold off
    if strcmp(output, 'HD29399') 
        legend('EW', 'LPD', 'Location', 'east')
    else 
        legend('EW', 'LPD', 'Location', 'northwest')
    end
    % legend('boxoff')

    xlim([min(bary_JD)-100 max(bary_JD)+100])
    xlabel('BJD-2450000', 'Interpreter', 'latex')
    title('Time series of activity proxies')
    if cmd_eps
        out_eps = [output, '_b.eps'];
        print(out_eps, '-depsc')
        close(hb);
    end                  


%     if find(ismember(flag_list, dir_list(n_dir)))
%         flag = 0;   % flagged
%     else
%         flag = 1;   % not flagged
%     end
    flag = 1;

    if 0
        star_name = strrep(dir_list{n_dir}, '_', '');
        star_name = char(strrep(star_name, '/', ''));    

        %%%%%%%%%%%%
        % RV vs EW %
        %%%%%%%%%%%%
        if cmd_subplot && cmd_publ
            subplot(3,2,3)
        elseif cmd_subplot
            subplot(2,2,3)
        else 
            hc = figure;
        end

        if cmd_publ
            corr_EW_RV_res(star_name, 1, 0)
        else
            herrorbar(EW_2, RV_2, delta_EW_2, 'r.')
            errorbar(EW_2, RV_2, delta_RV_2, 'b^', 'markers', 8)
        end
        title('RV vs EW')
        xlabel('EW (\AA)', 'Interpreter','latex')
        ylabel('RV (m/s)', 'Interpreter','latex')    
        if cmd_eps
            out_eps = [output, '_c.eps'];
            print(out_eps, '-depsc')
            close(hc);
        end                     

        %%%%%%%%%%%%%
        % RV vs LPD %
        %%%%%%%%%%%%%
        if cmd_subplot && cmd_publ
            subplot(3,2,4)
        elseif cmd_subplot
            subplot(2,2,4)
        else 
            hd= figure;
        end            

        if cmd_publ
            corr_EW_RV_res(star_name, 2, 0)
        else
            herrorbar(LPD_2, RV_2, STD_LPD_2, 'r.')
            errorbar(LPD_2, RV_2, delta_RV_2, 'b.', 'markers', 20)                
        end

        title('RV vs LPD')
        xlabel('LPD', 'Interpreter','latex')
        ylabel('RV (m/s)', 'Interpreter','latex')

        if cmd_eps
            out_eps = [output, '_d.eps'];
            print(out_eps, '-depsc')
            close(hd);
        end                   

        %%%%%%%%%%%%%%%%%%%%%
        % RV residual vs EW %
        %%%%%%%%%%%%%%%%%%%%% 
        if cmd_subplot && cmd_publ
            subplot(3,2,5)
        elseif cmd_publ
            he = figure;
        end
        if 0
            hold on 
            errorbar(RV_2, EW_2, delta_EW_2, '.', 'markers', 20)
            dataa = [RV_2', EW_2', delta_EW_2'];
            result = weightedfit(dataa);
            bb = result.slope;
            aa = result.Intercept;
            xx = linspace(min(RV_2), max(RV_2));
            yy = bb * xx + aa;
            % plot(xx, yy, '-')
            xlabel('RV(m/s)','Interpreter','latex')
            ylabel('EW (\AA)','Interpreter','latex')
            % legend(sprintf('%0.2e', bb), 'Location', 'Best')
            % legend('boxoff')
            title('EW vs RV')
            hold off
        end
        if cmd_publ
            corr_EW_RV_res(star_name, 1, 1)
            title('RV residual vs EW')
            xlabel('EW (\AA)', 'Interpreter', 'latex')
            ylabel('RV residual (m/s)', 'Interpreter', 'latex')
            if cmd_eps && cmd_publ
                out_eps = [output, '_e.eps'];
                print(out_eps, '-depsc')
                close(he);
            end                    
        end

        %%%%%%%%%%%%%%%%%%%%%%
        % RV residual vs LPD %
        %%%%%%%%%%%%%%%%%%%%%%
        if cmd_subplot && cmd_publ
            subplot(3,2,6)
        elseif cmd_publ
            hf= figure;
        end                
        if 0
            % h_diff = figure;
            hold on
            errorbar(RV_2, LPD_2, STD_LPD_2, '.', 'markers', 20)
            hold off
            % pic_name    = ['FLX_DIFF', dir_list{n_dir}(1:end-1), '_w', num2str(w), '.eps'];
            % print(pic_name, '-depsc')  
            % close(h_diff);
            xlabel('RV (m/s)', 'Interpreter','latex')
            ylabel('LPD', 'Interpreter','latex')  
            title('LPD vs RV ')
            if 0
                x_min = min(LPD_2 - STD_LPD_2 * 1.5);
                x_max = max(LPD_2 + STD_LPD_2 * 1.5);
                y_min = min(RV_2 - delta_RV_2 * 1.5);
                y_max = max(RV_2 + delta_RV_2 * 1.5);
                axis([x_min x_max y_min y_max])
            end
        end
        if cmd_publ
            corr_EW_RV_res(star_name, 2, 1)
            title('RV residual vs LPD')
            xlabel('LPD', 'Interpreter','latex')
            ylabel('RV residual (m/s)', 'Interpreter','latex')
            if cmd_eps && cmd_publ
                out_eps = [output, '_f.eps'];
                print(out_eps, '-depsc')
                close(hf);
            end                    
        end

        if 0
            %%%%%%%%%%%%%%%%%%
            subplot(3,2,[5,6])
            %%%%%%%%%%%%%%%%%%  
            errorbar(bary_JD_2, RV_2, delta_RV_2, 'b.')
            xlabel('BJD', 'Interpreter','latex')
            ylabel('RV (m/s)', 'Interpreter','latex')
            title('RV')
        end

    end

    %%%%%%%%%%%%%%%%%
    % Save to image %
    %%%%%%%%%%%%%%%%%
    if cmd_subplot && cmd_publ
        pdf_name    = ['EW_RV', dir_list{n_dir}(1:end-1), '_w', num2str(w), '.pdf'];
        print('-fillpage', pdf_name, '-dpdf')
    elseif cmd_subplot
        eps_name    = ['EW_RV', dir_list{n_dir}(1:end-1), '_w', num2str(w), '.eps'];
        print(eps_name, '-depsc')
    end
    close(h_show);
    fprintf('##### %d/%d RV found #####\n\n', N_file - n_nc, N_file);

    %%%%%%%%%%%%%%%%%
    % Scatter Plots %
    %%%%%%%%%%%%%%%%%
    if cmd_scatter
        if (flag ~= 0) && (idx_rv ~= 0)   % not flagged
            std_RV(n_dir)   = std(RV_2, 1./delta_RV_2.^2);
            std_EW(n_dir)   = std(EW_2, 1./delta_EW_2.^2);
        end
    end

    if cmd_scatter2
        if (flag ~= 0) && (idx_rv ~= 0) && cal_scatter2  % not flagged
            std_RV(n_dir)   = std(RV_2, 1./delta_RV_2.^2);
            std_ASYM(n_dir) = std(LPD_2, 1./(STD_LPD_2).^2);
        end
    end        

    clear bary_JD_2 RV_2 delta_RV_2 EW_2 delta_EW_2 LPD_2 STD_LPD_2
end

%%%%%%%%%%%%%%%%%
% Scatter Plots %
%%%%%%%%%%%%%%%%%
if cmd_scatter
    h_scatter = figure; 
    plot(std_RV, std_EW, '.')
    xlabel('RV scatter (m/s)','Interpreter','latex')
    ylabel('EW scatter (\AA)','Interpreter','latex') 
    % saveas(h_scatter, pic_name);
    pic_name = ['scatter_', num2str(w), '.eps'];
    print(pic_name, '-depsc')
    % close(h_scatter)
end

if cmd_scatter2
    h_scatter2 = figure; 
    plot(std_RV, std_ASYM, '.', 'markers', 20)
    xlabel('RV scatter (m/s)', 'Interpreter','latex')
    ylabel('ASYM scatter', 'Interpreter','latex') 
    % saveas(h_scatter, pic_name);
    pic_name = ['scatter2_', num2str(w), '.eps'];
    print(pic_name, '-depsc')
    % close(h_scatter2)
end
    
% end