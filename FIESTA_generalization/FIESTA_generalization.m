% Test for different spot and plage configurations 

%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
N_FILE      = 100;
t           = 1:N_FILE;
grid_size   = 0.1;
Fs          = 1/grid_size;
v0          = (-20 : grid_size : 20)';          % km/s

% spot/plage size
N_sampling  = 20;
S1          = (0.5 :0.5: 20)' / 100;
for i = 1:N_sampling
    area        = sqrt(2*S1(i));
    folder_name = ['CCF_PROT=25.05_i=90.00_lon=(180.0,0.0,0.0,0.0)_lat=(0.0,0.0,0.0,0.0)_size=(' , sprintf('%.4f',area), ',0.0000,0.0000,0.0000)'];
    DIR         = ['/Volumes/DataSSD/SOAP_2/outputs/' , folder_name];

% spot/plage lat
% N_sampling  = 18;
% lat         = (0: 5: 85)';
% for i = 1:N_sampling
%     folder_name = ['CCF_PROT=25.05_i=90.00_lon=(180.0,0.0,0.0,0.0)_lat=(', sprintf('%.1f',lat(i)), ',0.0,0.0,0.0)_size=(0.1000,0.0000,0.0000,0.0000)'];
%     DIR         = ['/Volumes/DataSSD/SOAP_2/outputs/' , folder_name];

% spot/plage S/N
% N_sampling  = 10;
% SN = (1000: 1000: 10000)';
% for i = 1:N_sampling
% A           = importdata(filename);    mkdir([DIR, '/CCF_noise', num2str(SN(i))])
%     sn = SN(i);
%     folder_name = 'CCF_PROT=25.05_i=90.00_lon=(180.0,0.0,0.0,0.0)_lat=(0.0,0.0,0.0,0.0)_size=(0.1000,0.0000,0.0000,0.0000)';
%     DIR         = ['/Volumes/DataSSD/SOAP_2/outputs/' , folder_name];    

    jitter      = importdata([DIR, '/RV.dat']) / 1000;      % activity induced RV [km/s]
    if N_FILE == 100
        jitter      = jitter';
    elseif N_FILE == 200
        jitter      = [jitter', jitter'];
    elseif N_FILE == 400
        jitter      = [jitter', jitter', jitter', jitter'];               % comment this out if not tesitng "planet + jitter"
    end
    idx         = (v0 >= -10) & (v0 <= 10);
    v1          = v0(idx);

    %%%%%%%%%%%%%%%%%%%
    % Calculate Power %
    %%%%%%%%%%%%%%%%%%%

    % estimate the size of array FFT_power
    filename    = [DIR, '/fits/CCF', num2str(1), '.dat'];
    A           = 1 - importdata(filename);
    A           = A(idx);
    A1          = A;
    [~, bb, yy] = FUNCTION_FFT(A, Fs);
    size1       = length(bb);
    FFT_power   = zeros(size1, N_FILE);
    Y           = zeros(size1, N_FILE);
    RV_noise    = zeros(1,N_FILE);
    RV_gauss    = zeros(N_FILE,1);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Stacked cross correlation function %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h = figure;
    hold on
    for n = 1:N_FILE

        filename    = [DIR, '/fits/CCF', num2str(mod(n-1,100)), '.dat'];
        A           = importdata(filename);

        % add noise
%         A           = A + normrnd(0, A.^0.5/SN(i));
%         dlmwrite([DIR, '/CCF_noise', num2str(SN(i)), '/', num2str(n-1), '.dat'], A)
        
        A           = 1 - A;
        A           = A(idx);
        

        
        % Plot the one without noise
        if mod(n,10) == 1
    %         plot(v1, A_spline - A1, 'k')
            plot(v1, A, 'k')
        end    

        % obtain line centroid 
        idx_fit     = (v1 >= -9) & (v1 <= 9);
        v_fit       = v1(idx_fit);
        A_fit       = A(idx_fit);
        f           = fit(v_fit, A_fit, 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [0.5, 0, 4, 0] );
        RV_gauss(n) = f.b;

        [FFT_frequency, FFT_power(:, n), Y(:, n)] = FUNCTION_FFT(A, Fs);

    end     
    hold off
    % title('Stacked cross correlation function')
    % ylim([-1.2/1000 1.2/1000])
    set(gca,'fontsize',20)
    xlabel('Velocity [km/s]')
    ylabel('Normalized intensity')
    saveas(gcf, [DIR, '/1-Line_Profile'],'png')
    close(h)

    % Determine the midpoint the equally divides the power spectrum %
    % cutoff_power= max(max(FFT_power)) * 0.0001; % used for publication for demonstration purpose
    % cutoff_power= max(max(FFT_power)) * 0.001;
    cutoff_power= max(max(FFT_power)) * 0.0005;
    % cutoff_power= max(max(FFT_power)) * 0.000005;
    f_max       = max(FFT_frequency(FFT_power(:,1) > cutoff_power));
    n           = abs(FFT_frequency) <= f_max;
    power_sum   = sum(FFT_power(n,1));

    % half power %
    if 1
        cum = 0;
        for i = 1:fix(sum(n)/2)
            cum = cum + FFT_power(size(FFT_power,1)/2+1+i,1);
            if cum > power_sum/4
                break
            end
        end
        f_HL = FFT_frequency(size(FFT_power,1)/2+1+i);
    %     f_HL = 0.5 * f_max;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % FT power in all epochs %
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    if 1
        h = figure; 
        hold on
        for n = 1:N_FILE
            plot(FFT_frequency, FFT_power(:, n), 'k')
        end 
        grey = [0.4,0.4,0.4];
        plot([-f_max, -f_max], [0, 0.6], '--', 'Color', grey)
        plot([f_max, f_max], [0, 0.6], '--', 'Color', grey)
        plot([-f_HL, -f_HL], [0, 0.6], '--', 'Color', grey)
        plot([f_HL, f_HL], [0, 0.6], '--', 'Color', grey)
        hold off
        text(0, 0.3, 'L', 'FontSize', 20, 'HorizontalAlignment','center')
        text((f_max+f_HL)/2, 0.3, 'H', 'FontSize', 20, 'HorizontalAlignment','center')
        text(-(f_max+f_HL)/2, 0.3, 'H', 'FontSize', 20, 'HorizontalAlignment','center')
        text((f_max+0.2)/2, 0.3, '\oslash', 'FontSize', 20, 'HorizontalAlignment','center')
        text(-(f_max+0.2)/2, 0.3, '\oslash', 'FontSize', 20, 'HorizontalAlignment','center')
        xlabel('\xi [s/km]')
        ylabel('Power')   
        xlim([-0.199 0.199])
        set(gca,'fontsize',20)
    %     saveas(gcf,'LPD2-FT_power','png')
        saveas(gcf, [DIR, '/2-FT_power'],'png')
        close(h)
    end


    %%%%%%%%%%%%%%%
    % Phase angle %
    %%%%%%%%%%%%%%%
    h = figure;
    plot(FFT_frequency, unwrap(angle(Y(:, 51))), '.')
    title('Phase angle (Rotation phase = 0.51)')
    xlabel('FT frequency (1 / velocity in wavelength)')
    ylabel('Phase angle [radian]')
    xlim([-0.15 0.15])
    saveas(gcf, [DIR, '/3-Phase_angle'],'png')
    close(h)


    %%%%%%%%%%%%%%%%%%%%%
    % Phase angle -> RV %
    %%%%%%%%%%%%%%%%%%%%%
    n       = abs(FFT_frequency) <= f_max;
    slope   = zeros(1,N_FILE);
    RV_FT   = zeros(1,N_FILE);
    RV_FT_err  = zeros(1,N_FILE);
    % wegihted_velocity = zeros(1,N_FILE);
    h = figure; 
    hold on
    for i = 1:N_FILE
        xx  = FFT_frequency(n);
        yy  = angle(Y(n, i)) - angle(Y(n, 1));
        if mod(i,10) == 1
            plot(xx, yy, 'k-')
        end   
        % Phase angle -> RV
        weight          = FFT_power(n,i)';
        [fitresult, gof]= createFit(xx, yy', weight);
        slope(i)        = fitresult.p1;
        RV_FT(i)        = -slope(i) / (2*pi);
        ci              = confint(fitresult,0.95);
        RV_FT_err(i)    = abs(diff(ci(:,1))*1000 / (4*pi));
    end
    hold on 
    % ymax = 0.0159;
    % ymax = 0.009;
    ymax = 0.011;
    plot([-f_max, -f_max], [-ymax, ymax], '--', 'Color', grey)
    plot([f_max, f_max], [-ymax, ymax], '--', 'Color', grey)
    plot([-f_HL, -f_HL], [-ymax, ymax], '--', 'Color', grey)
    plot([f_HL, f_HL], [-ymax, ymax], '--', 'Color', grey)
    hold off
    text(0, 0.005, 'L', 'FontSize', 20, 'HorizontalAlignment','center')
    text((f_max+f_HL)/2, 0.005, 'H', 'FontSize', 20, 'HorizontalAlignment','center')
    text(-(f_max+f_HL)/2, -0.005, 'H', 'FontSize', 20, 'HorizontalAlignment','center')
    text((f_max+0.2)/2, 0, '\oslash', 'FontSize', 20, 'HorizontalAlignment','center')
    text(-(f_max+0.2)/2, 0, '\oslash', 'FontSize', 20, 'HorizontalAlignment','center')
    hold off
    set(gca,'fontsize',20)
    % ylim([-ymax ymax])
    xlim([-0.199 0.199])
    xlabel('\xi [s/km]')
    ylabel('\Delta \phi [radian]')
    saveas(gcf, [DIR, '/4-Relative_phase_angle'] ,'png')
    % saveas(gcf,'LPD4-Relative_phase_angle','png')
    % saveas(gcf,'SLPD4-Relative_phase_angle','png')
    close(h)

    % Low-pass %
    nl      =  (abs(FFT_frequency) <= f_HL);
    RV_FTL  = zeros(1,N_FILE);
    RV_FTL_err  = zeros(1,N_FILE);
    h       = figure; 
    hold on
    for i = 1:N_FILE
        xx  = FFT_frequency(nl);
        yy  = angle(Y(nl, i)) - angle(Y(nl, 1));
        if mod(i,10) == 1
            plot(xx, yy, 'k-')
        end
        % Phase angle -> RV
        weight          = FFT_power(nl,i)';
        [fitresult, gof]= createFit(xx, yy', weight);
        slope(i)        = fitresult.p1;
        RV_FTL(i)       = -slope(i) / (2*pi);
        ci              = confint(fitresult,0.95);
        RV_FTL_err(i)   = abs(diff(ci(:,1))*1000 / (4*pi));    
    end
    hold off
    set(gca,'fontsize',20)
    xlim([0 f_HL])
    xlabel('\xi [s/km]')
    ylabel('\Delta \phi [radian]')
    title('Low-pass')
    saveas(gcf, [DIR, '/4-Relative_phase_angle_L'],'png')
    % saveas(gcf,'LPD4-Relative_phase_angle_L','png')
    % saveas(gcf,'SLPD4-Relative_phase_angle_L','png')
    close(h)

    % high-pass % 
    n       = (FFT_frequency >= f_HL) & (FFT_frequency <= f_max);
    % n       = (FFT_frequency <= -f_HL) & (FFT_frequency >= -f_max);
    RV_FTH  = zeros(1,N_FILE);
    RV_FTH_err  = zeros(1,N_FILE);
    h       = figure; 
    hold on
    for i = 1:N_FILE
        xx  = FFT_frequency(n);
        yy  = angle(Y(n, i)) - angle(Y(n, 1));
    %     if mod(i,10) == 1
            plot(xx, yy, 'k-')
    %     end
        % Phase angle -> RV
        weight          = FFT_power(n,i)';
        [fitresult, gof]= createFit(xx, yy', weight);
        slope(i)        = fitresult.p1;
        RV_FTH(i)       = -slope(i) / (2*pi);    
        ci              = confint(fitresult,0.95);
        RV_FTH_err(i)   = abs(diff(ci(:,1))*1000 / (4*pi));        
    end
    hold off
    set(gca,'fontsize',20)
    xlim([f_HL f_max])
    % xlim([-f_max -f_HL])
    xlabel('\xi [s/km]')
    ylabel('\Delta \phi [radian]')
    title('High-pass')
    saveas(gcf, [DIR, '/4-Relative_phase_angle_H'],'png')
    % saveas(gcf,'LPD4-Relative_phase_angle_H','png')
    % saveas(gcf,'SLPD4-Relative_phase_angle_H','png')
    close(h)

    % test
    % figure; plot(1:100, wegihted_velocity, 1:100, RV_FT*1000)

    %%%%%%%%%%%%%%%%
    % Minimization % 
    %%%%%%%%%%%%%%%%
    GG = (RV_gauss-mean(RV_gauss))*1000;
    XX = (RV_FT-mean(RV_FT))'*1000;
    YY = (RV_FTL-mean(RV_FTL))'*1000;
    ZZ = (RV_FTH-mean(RV_FTH))'*1000;

    dlmwrite( [DIR, '/GG.txt'], GG)
    dlmwrite( [DIR, '/XX.txt'], XX)
    dlmwrite( [DIR, '/YY.txt'], GG-YY)
    dlmwrite( [DIR, '/ZZ.txt'], ZZ-GG)
    
%     dlmwrite( [DIR, '/CCF_noise', num2str(sn), '/GG.txt'], GG)
%     dlmwrite( [DIR, '/CCF_noise', num2str(sn), '/XX.txt'], XX)
%     dlmwrite( [DIR, '/CCF_noise', num2str(sn), '/YY.txt'], GG-YY)
%     dlmwrite( [DIR, '/CCF_noise', num2str(sn), '/ZZ.txt'], ZZ-GG)
%     
    % plot(t, GG, '*', t, (jitter-mean(jitter))*1000, 'o')


    %     [fitresult, gof]= createFit(ZZ-GG, GG-YY, XX*0+1);
    [fitresult, gof]= createFit(ZZ-XX, XX-YY, XX*0+1);
    alpha = fitresult.p1'


    
end     


%%%%%%%%%%%%%%%%%%%
% ONLY LINE SHIFT %
%%%%%%%%%%%%%%%%%%%
if 0
    % ---------
    % Version 1
    % ---------
    % Compare with simulated RV % 
    h = figure;
        ax1 =subplot(20,1,1:10);
        xx = (v_planet_array-v_planet_array(1))*1000;
        yy1 = (RV_FT-RV_FT(1)) *1000;
        yy2 = (RV_gauss-RV_gauss(1))'*1000;
        [fitresult, gof]= createFit(xx, yy1, 1./(0.08^2+RV_FTH_err.^2*0).^2);
        fitresult
        hold on 
        scatter(xx, yy1, 'rs', 'MarkerFaceColor', 'none', 'MarkerFaceAlpha', 0.5)
        scatter(xx, yy2, 'bo', 'MarkerFaceColor', 'none', 'MarkerFaceAlpha', 0.5)
        p0 = plot(xx, xx, 'g--', 'LineWidth', 3); p0.Color(4)=0.8;
%         errorbar(xx, yy1, RV_FT_err, 'r.', 'MarkerSize', 0.1)
        hold off
        ylim([-0.5 10.1])
%         title('RV Recovery')
        ylabel('Output RV [m/s]')
%         daspect(ax1,[1 1 1])
%         set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'fontsize',14)
        legend({'RV_{FT}', 'RV_{Gaussian}', 'Output RV = Input RV'}, 'Location', 'northwest')

        positions = ax1.Position;
        ax2 = subplot(20,1,11:15);
        rms_gauss   = rms(yy2-xx - mean(yy1-xx))
        rms_FT      = rms(yy1-xx - mean(yy1-xx))  
        hold on 
        scatter(xx, yy1-xx- mean(yy1-xx), 'rs', 'MarkerFaceColor', 'none', 'MarkerFaceAlpha', 0.5)
        scatter(xx, yy2-xx- mean(yy2-xx), 'bo', 'MarkerFaceColor', 'none', 'MarkerFaceAlpha', 0.5)
        p0 = plot([min(xx), max(xx)], [0,0], 'k--', 'LineWidth', 3); p0.Color(4)=0.3;
        hold off
        ylim([-0.24 0.24])
        ylabel('Residual [m/s]')
        set(gca,'xticklabel',[])
        set(gca,'fontsize',14)
        
        positions = ax2.Position;
        ax2 = subplot(20,1,16:20);
        hold on 
        scatter(xx, yy2-yy1-mean(yy2-yy1), 'ko', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
%         errorbar(xx, yy1-xx, RV_FT_err, 'r.', 'MarkerSize', 0.1)
        p0 = plot([min(xx), max(xx)], [0,0], 'k--', 'LineWidth', 3); p0.Color(4)=0.3;
        hold off
        ylim([-0.01 0.01])
        ylabel('\Delta RV [m/s]')
        set(gca,'fontsize',14)
        xlabel('Input RV [m/s]')    
        
        saveas(gcf,'5-LINE_SHIFT_ONLY','png')
    close(h)
    
    
    % ---------
    % Version 2
    % ---------
    h = figure;
        ax1 =subplot(20,1,1:15);
        xx = (v_planet_array-v_planet_array(1))*1000;
        yy1 = (RV_FT-RV_FT(1)) *1000;
        yy2 = (RV_gauss-RV_gauss(1))'*1000;
        [fitresult, gof]= createFit(xx, yy1, 1./(0.08^2+RV_FTH_err.^2*0).^2);
        fitresult
        hold on 
        plot(xx, yy1-(mean(yy1)-5), 'ks', 'MarkerSize', 10)
%         plot(xx, yy2, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'black')        
        scatter(xx, yy2-(mean(yy2)-5), 15, 'ko', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
        p0 = plot(xx, xx, 'k-', 'LineWidth', 3); p0.Color(4)=0.2;
        hold off
        ylim([-0.5 10.1])
        ylabel('Output RV [m/s]')
%         daspect(ax1,[1 1 1])
%         set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'fontsize',20)
        legend({'RV_{FT}', 'RV_{Gaussian}', 'Output RV = Input RV'}, 'Location', 'northwest')

        positions = ax1.Position;
        ax2     = subplot(20,1,16:19);
        rms_gauss   = rms(yy2-xx - mean(yy1-xx))
        rms_FT      = rms(yy1-xx - mean(yy1-xx))  
        hold on 
%         scatter(xx, yy1-xx- mean(yy1-xx), 'rs', 'MarkerFaceColor', 'none', 'MarkerFaceAlpha', 0.5)
%         scatter(xx, yy2-xx- mean(yy2-xx), 'bo', 'MarkerFaceColor', 'none', 'MarkerFaceAlpha', 0.5)
        plot(xx, yy1-xx- mean(yy1-xx), 'ks', 'MarkerSize', 10)
        scatter(xx, yy2-xx- mean(yy2-xx), 15, 'ko', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
        p0 = plot([min(xx), max(xx)], [0,0], 'k-', 'LineWidth', 3); p0.Color(4)=0.2;
        hold off
        ylim([-0.6 0.6])
        ylabel('Residual [m/s]')
        set(gca,'fontsize',20)
        
        xlabel('Input RV [m/s]')    
        
        saveas(gcf,'5-LINE_SHIFT_ONLY','png')
    close(h)
    
    

    % high-pass and low-pass %
    h = figure;
        ax1 =subplot(20,1,1:15);
        xx = (v_planet_array-v_planet_array(1))*1000;
        yy1 = RV_FTH *1000;
        yy1 = yy1 - mean(yy1 - xx);
        yy2 = RV_FTL *1000;
        yy2 = yy2 - mean(yy2 - xx);
        [fitresult, gof]= createFit(xx, yy1, 1./(0.01^2+RV_FTH_err.^2));
        fitresult
        [fitresult, gof]= createFit(xx, yy2, 1./(0.01^2+RV_FTL_err.^2));
        fitresult        
        hold on 
%         scatter(xx, yy1, 'k*', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
%         scatter(xx, yy2, 'kD', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
        scatter(xx, yy1, 'k+')
        scatter(xx, yy2, 'kD')
        p0 = plot(xx, xx, 'k-', 'LineWidth', 3); p0.Color(4)=0.2;
        hold off
        ylim([-0.5 10.1])
%         title('RV Recovery')
        ylabel('Output RV [m/s]')
%         daspect(ax1,[1 1 1])
%         set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'fontsize',20)
        legend({'RV_{FT,H}', 'RV_{FT,L}', 'Output RV = Input RV'}, 'Location', 'northwest')

        positions = ax1.Position;
        ax2     = subplot(20,1,16:19);
        rms1    = rms(yy1-xx - mean(yy1-xx))  
        rms2    = rms(yy2-xx - mean(yy2-xx))
        hold on 
%         scatter(xx, yy1-xx, 'k*', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
%         scatter(xx, yy2-xx, 'kD', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
        scatter(xx, yy1-xx, 'k+')
        scatter(xx, yy2-xx, 'kD')
        p0 = plot([min(xx), max(xx)], [0,0], 'k-', 'LineWidth', 3); p0.Color(4)=0.2;
        hold off
        ylim([-0.6 0.6])
        xlabel('Input RV [m/s]')
        ylabel('Residual [m/s]')
        set(gca,'fontsize',20)
        saveas(gcf,'5-LINE_SHIFT_ONLY-HL','png')
    close(h)    
    
    
end    

%%%%%%%%%%%%%%%
% ONLY JITTER %
%%%%%%%%%%%%%%%
if 0
    % Compare with intrinsic line deformation RV % 
    h = figure; 
        t_alpha = 0.4;
        yyL = RV_FTL*1000;
        yyH = RV_FTH*1000;
        xx = (RV_gauss - RV_gauss(1))'*1000;
        if 0
            xx = xx(1:100);
            yyH = yyH(1:100);
            yyL = yyL(1:100);
            RV_FTL_err = RV_FTL_err(1:100);
            RV_FTH_err = RV_FTH_err(1:100);
        end
        hold on
        p1 = scatter(xx(1:100), yyL(1:100), 'kD');
        p1 = scatter(xx(1:100), yyH(1:100), 'k+');
%         p1 = scatter(xx(1:100), yyL(1:100), 'kD', 'MarkerFaceColor', 'k','MarkerEdgeColor','k', 'MarkerFaceAlpha',t_alpha,'MarkerEdgeAlpha',t_alpha);
%         p1 = scatter(xx(1:100), yyH(1:100), 'k*', 'MarkerFaceColor', 'k','MarkerEdgeColor','k', 'MarkerFaceAlpha',t_alpha,'MarkerEdgeAlpha',t_alpha);
%         p1 = scatter(xx(101:200), yyL(101:200), 'cD', 'MarkerFaceColor', 'c','MarkerEdgeColor','c', 'MarkerFaceAlpha',t_alpha,'MarkerEdgeAlpha',t_alpha);
%         p1 = scatter(xx(101:200), yyH(101:200), 'c*', 'MarkerFaceColor', 'c','MarkerEdgeColor','c', 'MarkerFaceAlpha',t_alpha,'MarkerEdgeAlpha',t_alpha);
        [fitresult_L, gof]= createFit(xx, yyL, 1./(1+RV_FTL_err*0).^2);
        fitresult_L        
        L1 = fitresult_L.p1;
        L2 = fitresult_L.p2;
        p3 = plot([min(xx), max(xx)], [L1*min(xx)+L2, L1*max(xx)+L2], 'k-', 'LineWidth', 2); p3.Color(4)=0.3;
        [fitresult_H, gof]= createFit(xx, yyH, 1./(1+RV_FTL_err*0).^2);
        fitresult_H        
        H1 = fitresult_H.p1;
        H2 = fitresult_H.p2;
        p4 = plot([min(xx), max(xx)], [H1*min(xx)+H2, H1*max(xx)+H2], 'k-', 'LineWidth', 2); p4.Color(4)=0.3;
%         xlim([-0.3 4.05])
%         ylim([-1.3 6.1])
%         xlim([-35 35])
%         xlabel('"Jitter" [m/s]')
        xlabel('Jitter (RV_{Gaussian}) [m/s]')
        ylabel('{\Phi}ESTA RV [m/s]')           
        legend({'RV_{FT,L}', 'RV_{FT,H}'}, 'Location', 'northwest')
        hold off
        set(gca,'fontsize', 18)
        corr(xx',yyL')
        corr(xx',yyH')
%         title('Correlations in RM process')
%         saveas(gcf,'Jitter_HD189733','png')
        saveas(gcf,'5-JITTER_ONLY_1_1','png')
    close(h)
    % p_fit =     0.7974   -0.0025 FOR A THREE-SPOT CONFIGRATION (0.8104)
    % p_fit =     0.7736   -0.2375 FOR A TWO-SPOT CONFIGRATION 

        
    if 0
        h = figure; 
            plot(xx, xx - yyL, 'o')
            p_fit = polyfit(xx, xx - yyL,1)
            title('FT with jitter only (residual)')
            xlabel('RV_{Gaussian} [m/s]')
            ylabel('RV_{Gaussian} - RV_{FT} (m/s)')
            saveas(gcf,'5-JITTER_ONLY_2','png')
        close(h)
        h = figure; 
            plot(xx, yyH , 'o')
            p_fit = polyfit(xx, xx - yyL,1)
            title('FT with jitter only (residual)')
            xlabel('RV_{Gaussian} [m/s]')
            ylabel('RV_{FT} - RV_{Gaussian} (m/s)')
%             saveas(gcf,'5-JITTER_ONLY_2','png')
        close(h)        
    end
    
    % TIME SERIES
    h = figure;
    ax1 = subplot(20,1,1:15);
        xx = (1:N_FILE)/N_FILE;
        yy1 = (RV_FT - mean(RV_FT))'*1000;
        yy2 = (RV_gauss - mean(RV_gauss))*1000;
        hold on 
        plot(xx, yy1, 'ks', 'MarkerSize', 10)
        scatter(xx, yy2, 18, 'ko', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)        
%         scatter(xx, yy1, 'rs', 'MarkerFaceColor', 'none', 'MarkerFaceAlpha', 0.5);
%         scatter(xx, yy2, 'bo', 'MarkerFaceColor', 'none', 'MarkerFaceAlpha', 0.5);
        hold off
%         title('Apparent RV of deformed line profile')
        legend({'RV_{FT}', 'RV_{Gaussian}'}, 'Location', 'northwest')
        ylabel('RV [m/s]')
        ylim([-2.5 2.5])
        set(gca,'fontsize',18)
        set(gca,'xticklabel',[])
    ax2 = subplot(20,1,16:19);
        hold on
        scatter(xx, yy1 - yy2, 'k*')
        p0 = plot([min(xx), max(xx)], [0,0], 'k--', 'LineWidth', 3); p0.Color(4)=0.3;
        hold off
        xlabel('Stellar rotation phase')
        ylabel('\Delta RV [m/s]')
        ylim([-0.0149 0.0149])
        set(gca,'fontsize',18)
    saveas(gcf,'5-JITTER_ONLY_3','png')
    close(h)
    
    % TIME SERIES 2 
    h = figure;
    ax1 = subplot(20,1,1:15);
        xx = (1:N_FILE)/N_FILE;
        yy = (RV_gauss - RV_gauss(1))*1000;
        hold on 
        plot(xx, yy, '--', 'color', [0.9100    0.4100    0.1700], 'LineWidth', 3);
        scatter(xx, (yyL-L2)/L1, 'kD', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5);
        scatter(xx, (yyH-H2)/H1, 'k*', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5);
        hold off
        legend({'Jitter', 'RV_{FT,L} / k_{L}', 'RV_{FT,H} / k_{H}'}, 'Location', 'northwest')
        ylabel('RV [m/s]')
        ylim([-1.5 4.5])
        set(gca,'xticklabel',[])
%         title('Fitting apparent RV of deformed line profile')
        set(gca,'fontsize',15)
    ax2 = subplot(20,1,16:20);
        hold on
        scatter(xx, (yyL-L2)/L1 - yy', 'kD', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5);
        scatter(xx, (yyH-H2)/H1 - yy', 'k*', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5);
        hold off
        rms((yyL-L2)/L1 - yy')
        rms((yyH-H2)/H1 - yy')
        xlabel('Stellar rotation phase')
        ylabel('Residual [m/s]')
        ylim([-0.75 0.75])
        set(gca,'fontsize',15)
    saveas(gcf,'5-JITTER_ONLY_4','png')    
    close(h)    
    
    % TIME SERIES 3 
    h = figure;
    ax1 = subplot(20,1,1:15);
        yy2 = (RV_FT - RV_FT(1))*1000;
        yy1 = (RV_gauss - RV_gauss(1))'*1000;
        hold on 
        scatter(yy1, yy2, 18, 'ko', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)        
        [fitresult, gof]= createFit(yy1, yy2, 1./(1+RV_FTL_err*0).^2);
        hold off
        ylabel('RV_{FT} [m/s]')
        ylim([-0.2 4.05])
        xlim([-0.2 4.05])
        set(gca,'fontsize',18)
        set(gca,'xticklabel',[])
    ax2 = subplot(20,1,16:19);
        hold on
        scatter(yy1, yy1 - yy2, 18, 'ko', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
%         p0 = plot([min(xx), max(xx)], [0,0], 'k--', 'LineWidth', 3); p0.Color(4)=0.3;
        hold off
        xlabel('Jitter (RV_{Gaussian}) [m/s]')
        ylabel('Residual [m/s]')
        ylim([-0.0149 0.0149])
        xlim([-0.2 4.05])
        set(gca,'fontsize',18)
    saveas(gcf,'5-JITTER_ONLY_5','png')    
    close(h)        
    
    
    
    
    
    
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%
% LINE SHIFT AND JITTER %
%%%%%%%%%%%%%%%%%%%%%%%%%

% Compare the total input RV and with the recovered RV
if 0
    s_len = 2;
    % TIME SERIES
    h = figure; 
    ax1 = subplot(3,1,1);
        t    = (1:N_FILE)';
        yyL  = RV_FTL' * 1000;
        yyH  = RV_FTH' * 1000;
%         yy2  = (RV_gauss - RV_gauss(1)) * 1000;
        yy2  = (RV_FT - RV_FT(1))' * 1000;
        hold on
        scatter(t/100, yyL, 'kD');
        scatter(t/100, yyH, 'k+');
        scatter(t/100, yy2, 15, 'bo', 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5);
        hold off
        ylabel('RV [m/s]')    
        legend({'RV_{FT,L}', 'RV_{FT,H}', 'RV_{Gaussian}'}, 'Location', 'north')
        set(gca,'fontsize', 20)
%         set(gca,'xticklabel',[])
        ylim([-5 6])

    ax2 = subplot(3,1,2); 
        rv_L        = yy2 - yyL;
        rv_H        = yyH - yy2;
        rv_HL       = 0.5*rv_L + 0.5*rv_H *alpha;
        
        t_smooth    = linspace(1,N_FILE, 1000)';
        y_smooth1    = FUNCTION_GAUSSIAN_SMOOTHING(t, rv_L, t_smooth, s_len);
        y_smooth11   = FUNCTION_GAUSSIAN_SMOOTHING(t, rv_L, t, s_len);
        xx2 = (jitter- jitter(1))' * 1000;
        p_fit1 = polyfit(xx2, rv_L, 1)
        hold on
        jitter_model1 = (rv_L-p_fit1(2))/p_fit1(1);
        c = @cmu.colors;
        plot(t/100, xx2, '--', 'color', [0.9100    0.4100    0.1700], 'LineWidth', 3)
        scatter(t/100, jitter_model1, 'kD')

        y_smooth3    = FUNCTION_GAUSSIAN_SMOOTHING(t, rv_HL, t_smooth, s_len);
        y_smooth33    = FUNCTION_GAUSSIAN_SMOOTHING(t, rv_HL, t, s_len);
        p_fit3 = polyfit(xx2, rv_HL, 1)
        jitter_model3 = (rv_HL-p_fit3(2))/p_fit3(1);
        scatter(t/100, jitter_model3, 20, 'k^')
        
        y_smooth2    = FUNCTION_GAUSSIAN_SMOOTHING(t, rv_H, t_smooth, s_len);
        y_smooth22    = FUNCTION_GAUSSIAN_SMOOTHING(t, rv_H, t, s_len);
        p_fit2 = polyfit(xx2, rv_H, 1)
        jitter_model2 = (rv_H-p_fit2(2))/p_fit2(1);
        scatter(t/100, jitter_model2, 20, 'k+')
        
%         jitter_y_val    = importdata('jitter_y_val.txt');
%         p_fit4 = polyfit(xx2, jitter_y_val, 1)
%         jitter_model4 = (jitter_y_val-p_fit4(2))/p_fit4(1);
%         plot(t, jitter_model4)

        plot1 = plot(t_smooth/100, (y_smooth1-p_fit1(2))/p_fit1(1), 'k', 'LineWidth', 2);
        plot1.Color(4) = 0.2;        
        plot2 = plot(t_smooth/100, (y_smooth2-p_fit2(2))/p_fit2(1), 'k', 'LineWidth', 2);
        plot2.Color(4) = 0.2;
        plot3 = plot(t_smooth/100, (y_smooth3-p_fit3(2))/p_fit3(1), 'k', 'LineWidth', 2);
        plot3.Color(4) = 0.2;        
%             pbaspect(ax2,[5 1 1])
        hold off
        ylabel('Jitter [m/s]')
        ylim([-5 6])
        legend({'Input jitter', 'Model (w_1=1)', 'Model (w_1=0.5)', 'Model (w_1=0)'}, 'Location', 'north')
        set(gca,'fontsize', 20)
%         set(gca,'xticklabel',[])

    ax3 = subplot(3,1,3);
        hold on 
        scatter(t/100, jitter_model1 - xx2, 'kD')
        scatter(t/100, jitter_model2 - xx2, 'k+')
        scatter(t/100, jitter_model3 - xx2, 'k^')
        plot1 = plot(t/100, (y_smooth11-p_fit1(2))/p_fit1(1) - xx2, 'k', 'LineWidth', 2);
        plot1.Color(4) = 0.2; 
        plot2 = plot(t/100, (y_smooth22-p_fit2(2))/p_fit2(1) - xx2, 'k', 'LineWidth', 2);
        plot2.Color(4) = 0.2;         
        plot3 = plot(t/100, (y_smooth33-p_fit3(2))/p_fit3(1) - xx2, 'k', 'LineWidth', 2);
        plot3.Color(4) = 0.2;                 
        hold off
        ylim([-5 6])
        xlabel('Stellar rotation phase')
        ylabel('Residual [m/s]')
        set(gca,'fontsize', 20)
        saveas(gcf,'5-PLANET_AND_JITTER2','png')

        rms(xx2 - mean(xx2))    %1.2242
        
        rms(jitter_model1 - xx2) %0.6948
        rms((y_smooth11-p_fit1(2))/p_fit1(1) - xx2) %0.5358
        rms(jitter_model2 - xx2) %0.7755
        rms((y_smooth22-p_fit2(2))/p_fit2(1) - xx2) %0.6168
        rms(jitter_model3 - xx2) %0.7172
        rms((y_smooth33-p_fit3(2))/p_fit3(1) - xx2) %0.5701
        
        rms(jitter_model4 - xx2) 
        
        rms((jitter_model1 + jitter_model2)/2 - xx2) % 0.7191
        rms(((y_smooth22-p_fit2(2))/p_fit2(1) + (y_smooth11-p_fit1(2))/p_fit1(1))/2 - xx2) %0.5720
    close(h) 

    
    % Obtain kl, kh %
    yyL  = RV_FTL' * 1000;
    yyH  = RV_FTH' * 1000;
    yyG  = (RV_gauss - RV_gauss(1)) * 1000;
    plot(yyG-yyL, yyH-yyG, '.')
    p_fit = polyfit(yyG-yyL, yyH-yyG, 1)
    
    plot(yyG-yyL, yyH-yyL, '.')
    p_fit = polyfit(yyG-yyL, yyH-yyL, 1)    
    
    
    
    % REAL JITTER VS SCALED JITTER % 
    if 0
        h = figure; 
        hold on            
        plot(xx2, polyval(p_fit,xx2), '-')
        (sum((rv_d - polyval(p_fit,xx2)).^2)/N_FILE)^0.5
        (sum((rv_d - polyval(p_fit,xx2)).^2)/N_FILE)^0.5 / p_fit(1)
        hold off
        saveas(gcf,'7-Jitter_scaling','png')
        close(h) 
        % p_fit =    0.2019    0.0023 FOR A THREE-SPOT CONFIGRATION (0.1867)
        % p_fit =    0.2259    0.2377 FOR A TWO-SPOT CONFIGRATION
    end
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
    % NOT USED FOR NOW %
    figure; 
    xx = (jitter - v_planet_array(1:100)' - jitter(51))*1000;
    yy = slope/(-2*pi)*1000;
    plot(xx, yy' - xx, 'o')
    title('Recovered RV vs total RV')
    xlabel('input RV including jitter (m/s)')
    ylabel('residual (m/s)')    



    RVT_FT  = slope/(-2*pi)*1000;
    RVT     = (jitter - v_planet_array(1:100)' - jitter(51))*1000;

    figure; plot(t, RVT, t, RVT_FT')
    figure; plot(t, RVT - RVT_FT', t, (jitter-jitter(51))*1000)

    figure; plot(t, (RVT_FT' - RVT*0.7736)/(1-0.7736))




    x = (1:100)';
    y = (slope/(-2*pi)*1000)';
    rv = (jitter - v_planet_array(1:100)' - jitter(51))*1000;
    % dlmwrite('RV_tot.txt', rv)
    f = fit(x, y, 'a * sin(x/100*k+phi) * (1-m) + b + m * rv', 'StartPoint', [3 7 0.1 0.9 0 0]);






    ft = fittype('rv_ft(x, a, k, phi, m, b)');
    f = fit(t, RVT, ft);



    a0 = [3 7 0.1 0.9 0]; %?????a??????????
    options = statset('Jacobian', 'on');
    [a,r,J,~,msE] = nlinfit(x, y, @rv_ft, a0, options);%??
    [ypred,delta] = nlpredci(@rv_ft,x,a,r,'Jacobian',J,'predopt','observation','simopt','on');%??????????



    % test 

    x = (1:10)';
    y = x.^2+x/10;

    f = fit(x, y, 'a*x.^2+b', 'StartPoint', [0.01, 2]);
end



%%%%%%%%%%
% la fin %
%%%%%%%%%%

if 0


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct a template in FT space %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FFT_tpl = mean(FFT_power');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Correct each observation %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Y_c = zeros(size(Y));
    Y_i = zeros(size(Y));
    iFFT_power = zeros(size(FFT_power));

    % Scale the power, and keep the phase

    for n = 1:N_FILE
        for i = 1:length(FFT_tpl)
            Y_c(i,n) = Y(i,n) / abs(Y(i,n)) *  FFT_tpl(i);
        end
    end

    for n = 1:N_FILE
        [iFFT_power(:, n), Y_i(:,n)] = FUNCTION_iFFT(Y_c(:,n));
    end

    figure;
    hold on
    for n = 1:N_FILE
        plot(FFT_frequency, iFFT_power(:, n) - iFFT_power(:, 1), '.')
    end
    hold off 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Examine the corrected line centre %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i_planet = zeros(1,N_FILE);
    for n = 1:N_FILE
        % idx     = iFFT_power(:, n) > 0.01;
        % f_power = fit( FFT_frequency(idx)'-0.5, iFFT_power(idx, n), 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [1, 0.5, 1, 0] );
        idx = (FFT_frequency<1) & (FFT_frequency>0);
        f_power = fit( FFT_frequency(idx)', iFFT_power(idx, n), 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [0.5, 0.5, 0.18, 0] );
        i_planet(n) = f_power.b;goo
    end

    % test %
    figure; plot(f_power, FFT_frequency(:), iFFT_power(:, n))

    v_planet    = 2 * sin(mod(MJD, 100)/100*7*2*pi + 1) * 0.001;
    figure; plot(RV_noise, i_planet, '.', 'markers',12)


end