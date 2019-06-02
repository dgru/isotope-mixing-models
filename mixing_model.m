
%% Script to Model Mixing of Two Endmembers
    % created by DSG and ARW on 10/08/2015
    % Used for modeling hydrate mixing or mixing microbial and themogenic methane

%% set some things up

% VPDB and VSMOW values for reference
R_13C_VPDB = 0.0112372;
R_2H_VSMOW = 155.76*10^-6;

f = 0:0.05:1; % fraction of mixing, 0 to 1

% d13C_1 = -50; % delta values of 13C methane
% R13C_1 = (d13C_1./1000 + 1).*0.0112372;  % 0.0112372 = PDB
% d13C_2 = -50; % delta values of 13C methane
% R13C_2 = (d13C_2./1000 + 1).*0.0112372;  % 0.0112372 = PDB

d13C_1 = -90; % delta values of 13C methane
R13C_1 = (d13C_1./1000 + 1).*0.0112372;  % 0.0112372 = PDB
d13C_2 = -30; % delta values of 13C methane
R13C_2 = (d13C_2./1000 + 1).*0.0112372;  % 0.0112372 = PDB

% dD_1 = -140;
% RD_1 = (dD_1./1000 + 1).*155.76*10^-6;
% dD_2 = -140;
% RD_2 = (dD_2./1000 + 1).*155.76*10^-6;

dD_1 = -190;
RD_1 = (dD_1./1000 + 1).*155.76*10^-6;
dD_2 = -20;
RD_2 = (dD_2./1000 + 1).*155.76*10^-6;

% Delta 13CH3D values
% Change these to the two endmember D13CH3D values
capD13CH3D_1 = 5.5; % 0.0 permil, replace with endmember 1 (biogenic)
capD13CH3D_2 = 2.0; % 0.0 permil, replace with endmember 2 (thermogenic)

% Convert D13CH3D to R_13CH3D
R_13CH3D_1 = R13C_1*RD_1*(capD13CH3D_1/1000+1);
R_13CH3D_2 = R13C_2*RD_2*(capD13CH3D_2/1000+1);

% Compute d13CH3D
d13CH3D_1 = 1000*(R_13CH3D_1/(R_13C_VPDB*R_2H_VSMOW)-1);
d13CH3D_2 = 1000*(R_13CH3D_2/(R_13C_VPDB*R_2H_VSMOW)-1);


% Calculate 13C Mixing of 2 End Members and Plot
    C_Mix = (1-f).*(R13C_1)+(f).*(R13C_2);  % mixing equation with R values
    C_Mix_Delta = (C_Mix/0.0112372 - 1)*1000;
    
    figure(1);
    hold on;
    scatter(f,C_Mix);
        [a,b] = polyfit(f,C_Mix,1);
        C_fit = a(1).*f + a(2);
    plot(f,C_fit,'k--','linewidth',1);
    
% Calculate D Mixing of 2 End Members and Plot
    D_Mix = (1-f).*(RD_1)+(f).*(RD_2);  % mixing equation
    D_Mix_Delta = (D_Mix/(155.76*10^-6) - 1)*1000;
    
    figure(2);
    hold on;
    scatter(f,D_Mix);
        [a,b] = polyfit(f,D_Mix,1);
        D_fit = a(1).*f + a(2);
    plot(f,D_fit,'k--','linewidth',1);

% Calculate 13CH3D Mixing of 2 End Members and Plot
    CD_Mix = (1-f).*(R_13CH3D_1)+(f).*(R_13CH3D_2);
    CD_Mix_Delta = 1000*(CD_Mix/(R_13C_VPDB*R_2H_VSMOW)-1); % little delta
    CD_Mix_capD = 1000*(CD_Mix./(C_Mix.*D_Mix)-1); % Capital delta 13CH3D
    
% Calculate 13CH3D and Plot

% stochastic distributions
%    s_12CH4 = 9.88*10^-1;
%    s_12CH3D = 6.16*10^-4;
%    s_13CH4 = 1.11*10^-2;
%    s_13CH3D = 6.93*10^-6;

%dCH3D = (s_12CH3D/s_12CH4./(4.*D_fit)-1)*1000;
%d13CH4 = (s_13CH4/s_12CH4./C_fit-1)*1000;
%d13CH3D = (s_13CH3D/s_12CH4./(4.*C_fit.*D_fit)-1)*1000; 

%D18 = log(d13CH3D) - log(d13CH4); 
D18 = CD_Mix_capD;
figure(3);
hold on;
scatter(f,D18);

%% Make plot

figure(4);
hold on;

subplot(1,3,1);
    h1 = scatter(f,C_Mix_Delta,'k','filled');
    xlabel('fraction')
    ylabel('\delta13C')
    axis square
    box on

subplot(1,3,2);
    h2 = scatter(f,D_Mix_Delta,'k','filled');
    xlabel('fraction')
    ylabel('\deltaD')
    axis square
    box on

subplot(1,3,3);
    x = D18;
    c = linspace(1,10,length(x));
    
    h3 = scatter(f,D18,[],c,'filled');
    xlabel('fraction')
    ylabel('\Delta13CH_3D')
    axis square
    box on

    
mymap = [   0.254901960784314,  0.713725490196078,  0.768627450980392
            0.4,                0.6,                0.4
            0.8,                0.5                 0.2
            1,                  0.498039215686275,  0];
        
    colormap jet;
    axis([0 1 1 7]);

suptitle('Mixing Model of 13C, D, and 13CH3D');
