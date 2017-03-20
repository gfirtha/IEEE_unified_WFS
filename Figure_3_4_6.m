%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This work is supplementary material for the article                     %
%                                                                         %
% Gergely Firtha, Péter Fiala, Frank Schultz and Sascha Spors             %
%   Improved Referencing Schemes for 2.5D Wave Field Synthesis            %
%    Driving Functions                                                    %
% IEEE / ACM Transactions on Audio, Speech and Language Processing, 2017  %
%                                                                         %
%                                                                         %
% (c) 2017 by Gergely Firtha                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all

% simulation parameters
speed_of_sound = 343.1;
f = 1.5e3;
omega = 2*pi*f;
k = omega/speed_of_sound;

% create simulation mesh
dx = 0.025;
y = (dx:dx:3)';
x0 = (-70:dx:70)';
[X,Y] = meshgrid(x0,y);

% virtual field parameters
xs = [0,-1];                                                % virtual source position
alpha = 60*pi/180;                                          % virtual plane wave direction

% target sound fields
r_s = sqrt( ( X - xs(1) ) .^2 + ( Y - xs(2) ).^2 );
P_ref_ls =   - 1i/4 * besselh(0,2,k*r_s);                   % virtual line source
P_ref_pw =   exp( -1i*k*(cos(alpha)*X + sin(alpha)*Y) );    % virtual plane wave

% Calculate driving functions
dref = 2;                                                   % reference distance / line
ref_type = 'ref_line';
r0 = sqrt( ( x0 - xs(1) ).^2 + xs(2)^2 );

switch ref_type
    % Constant referencing function: Figure 3 is generated
    case 'constant'
        d_ls = dref;
        d_pw = dref;
    case 'ref_line'
        % Referencing on a parallel line: Figure 4 is generated
        d_ls = dref*r0./abs(xs(2));
        d_pw = dref/sin(alpha);
    case 'circle'
        % Referencing on a reference circle: Figure 6 is generated
        Rref = abs(xs(2)) + dref;
        d_ls = Rref-r0;
        d_pw = 0;
end
D_ls = sqrt(2*pi/(1i*k))*sqrt(d_ls)*0.5*1i*k*xs(2)./r0.*besselh(1,2,k*r0);
D_pw = sqrt(8*pi/(1i*k))*sqrt(d_pw)*1i*k*sin(alpha).*exp(-1i*k*cos(alpha)*x0);


% Calculate synthesized field based on Eq. (35) evaluted by numerical
% convolution
G0  =  1/(4*pi) * exp(-1i*k*sqrt( X.^2 + Y.^2 ))./sqrt( X.^2 + Y.^2 );
P_synth_ls = zeros(size(X));
P_synth_pw = zeros(size(X));
wb = waitbar(0,'Calculating synthesized field');
for n = 1:length(y)
    waitbar(n/length(y),wb);
    % Eq. (35)
    P_synth_ls(n,:)= conv(D_ls,G0(n,:),'same')*dx;
    P_synth_pw(n,:)= conv(D_pw,G0(n,:),'same')*dx;
end
close(wb);

% Calculate Positions of Correct Synthesis
% Eq. (30)
x_ref_ls = x0 + d_ls.*(x0-xs(1))./r0;
y_ref_ls = -d_ls.*xs(2)./r0;
% Eq. (28)
x_ref_pw = x0 + d_pw.*cos(alpha);
y_ref_pw = d_pw.*sin(alpha);

%% Draw results
x_lim = 6;
i = find(abs(x0)<= x_lim);
ftsize = 20;

switch ref_type
    case {'constant','ref_line'}
        f = figure('Units','points','Position',[200,50,600,420]);
        set(gcf,'Units','normalized');
        
        p1 = axes('Units','normalized','Position',[ 0.068 0.575 .87 .4 ]);
        pcolor(x0(i),y,20*log10(abs(P_ref_pw(:,i)-P_synth_pw(:,i))));shading interp;axis equal tight
        c = colorbar;
        ylabel(c,'[dB]' , 'Interpreter', 'LaTex' , 'FontSize', ftsize)
        caxis([-30,20])
        xlim([-x_lim,x_lim])
        hold on
        plot(x_ref_pw,0*x_ref_pw+y_ref_pw,'--w','LineWidth',1)
        ylim([y(1),y(end)])
        xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
        ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
        title('$20\mathrm{log}_{10}\left( P_{\mathrm{synth,pw}}(\mathbf{x},\omega)-P_{\mathrm{ref,pw}}(\mathbf{x},\omega) \right)$'...
            , 'Interpreter', 'LaTex' , 'FontSize', ftsize);
        set(gca,'FontName','Times New Roman');
        
        p2 = axes('Units','normalized','Position',[ 0.068 0.075 .87 .4 ]);
        pcolor(x0(i),y,20*log10(abs(P_ref_ls(:,i)-P_synth_ls(:,i))));shading interp;axis equal tight
        c = colorbar;
        ylabel(c,'[dB]' , 'Interpreter', 'LaTex' , 'FontSize', ftsize)
        caxis([-70,0])
        hold on
        plot(x_ref_ls,y_ref_ls,'--w','LineWidth',1)
        xlim([-x_lim,x_lim])
        ylim([y(1),y(end)])
        xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
        ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
        title('$20\mathrm{log}_{10}\left( P_{\mathrm{synth,ls}}(\mathbf{x},\omega)-P_{\mathrm{ref,ls}}(\mathbf{x},\omega) \right)$'...
            , 'Interpreter', 'LaTex' , 'FontSize', ftsize);
        set(gca,'FontName','Times New Roman');
        allAxesInFigure = findall(f,'type','axes');
        b = get(gca,'XTickLabel');
        set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-3);
        set(gcf,'PaperPositionMode','auto');
        
        %print -dpng fixed_referencing -r300
        
    case 'circle'
        
        f = figure('Units','points','Position',[200,50,600,210]);        
        set(gcf,'Units','normalized');
        p1 = axes('Units','normalized','Position',[ 0.068 0.075 .87 .9 ]);
        pcolor(x0(i),y,20*log10(abs(P_ref_ls(:,i)-P_synth_ls(:,i))));shading interp;axis equal tight
        
        c = colorbar;
        ylabel(c,'[dB]' , 'Interpreter', 'LaTex' , 'FontSize', ftsize)
        caxis([-70,0])
        hold on
        plot(x_ref_ls,y_ref_ls,'--w','LineWidth',1)
        xlim([-6,6])
        ylim([y(1),y(end)])
        xlabel( '$x \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
        ylabel( '$y \rightarrow [\mathrm{m}]$' , 'Interpreter', 'LaTex' , 'FontSize', ftsize );
        title('$20\mathrm{log}_{10}\left( P_{\mathrm{synth,ls}}(\mathbf{x},\omega)-P_{\mathrm{ref,ls}}(\mathbf{x},\omega) \right)$'...
            , 'Interpreter', 'LaTex' , 'FontSize', ftsize);
        set(gca,'FontName','Times New Roman');
        allAxesInFigure = findall(f,'type','axes');
        b = get(gca,'XTickLabel');
        set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-3);
        set(gcf,'PaperPositionMode','auto');
        %print -dpng circle_referencing -r300
end