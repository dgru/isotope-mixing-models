% Written by SO and DG for Gruen et al., 2018 GCA.

% Simulation and Parameters

        %13C/12C ratio of CO2
            R13 = 0.01; % for PDB
            R13i = 0.01*(1+10.9/1000); % 10.9 is estimated initial d13C of CO2

        % D/H ratio of H source
            % assume this is only water
            % assume that the water is -50 permil
            % MATH: (50 permil = 5%) 0.94*150e-6
            R2 = 150e-6; 
            R2w = 0.95*150e-6; % 141e-6;

        % clumped isotope effect for 13CH3D
            K = 1.0042; % 1.0035; % chosen to fit data
            K2 = 1.0032;

        N = 500; % number of steps for simulation

        % 13C/12C fractionation factor
        % assume = 0.98
        % assume this stays constant for duration
        a13 = linspace(0.9704,0.9704,N); % chosen to fit data

        % D/H fractionation factor 
            % assume this is one value for first 100 steps
            % changes for next 400 steps
            a2 = [linspace(0.69,0.57,100) linspace(0.57,0.57,400)];
            % a2 = [linspace(0.67,0.55,500)]; 

        % Valentine 13C alpha 1.023 - 1.064 (Valentine 2004 GCA)
        % Valentine 2D alpha 1.16 - 1.43
        % a2Vmax = 1/1.43
        % a2Vmin = 1/1.16
        % a2 = [linspace(a2Vmax,a2Vmax,500)];

        a = a13; 
        b = 4*a2*R2w; 
        c = 4*a13.*a2*R2w*K; 
        c2 = 4*a13.*a2*R2w*K2; 

        k = 0.01; % arbitarary rate constant
        % k = 0.03; % arbitarary rate constant

            CO2i = [1/(1+R13i) R13i/(1+R13i)]; % initial amount of 12CO2 and 13CO2
            CH4i = [0 0 0 0 0]; % 12CH4, 12CH3D, 13CH4, 13CH3D

            CO2 = CO2i;
            CH4 = CH4i;

        for n=1:500
            CO2(n+1,:) = CO2(n,:) + [-k*(1+b(n)) -k*(a(n)+c(n))].*CO2(n,:); 
            CH4(n+1,:) = CH4(n,:) + [k b(n)*k a(n)*k c(n)*k c2(n)*k].*[CO2(n,1) CO2(n,1) CO2(n,2) CO2(n,2) CO2(n,2)];
        end


            d13CO2 = (CO2(:,2)./CO2(:,1)/R13-1)*1000;
            dCH3D = (CH4(:,2)./CH4(:,1)/(4*R2)-1)*1000;
            d13CH4 = (CH4(:,3)./CH4(:,1)/R13-1)*1000;
            d13CH3D = (CH4(:,4)./CH4(:,1)/(4*R13*R2)-1)*1000; 
            d13CH3D2 = (CH4(:,5)./CH4(:,1)/(4*R13*R2)-1)*1000; 

            ln13CH4 = log(CH4(:,3)./CH4(:,1)/R13)*1000;
            lnCH3D = log(CH4(:,2)./CH4(:,1)/(4*R2))*1000;
            ln13CH3D = log(CH4(:,4)./CH4(:,2)/R13)*1000;
            ln13CH3D2 = log(CH4(:,5)./CH4(:,2)/R13)*1000;
            D18 = ln13CH3D - ln13CH4;
            D182 = ln13CH3D2 - ln13CH4;

            F = 1-CO2(:,1)./CO2(1,1);

    % assign data into a new matrix called data1
    % FYI 13C = column 4, error5, dD6, error7, 13CH3D8, error9, percent12 
       data1 = [Clumped_Data(23,:)
                Clumped_Data(26,:)
                Clumped_Data(29,:)
                Clumped_Data(31,:)];

    % Plot
        figure(2);

        % 13C
            subplot(1,3,1);
            hold on;
            box on;
            plot(F, d13CH4,'r--','linewidth',1);
            h1 = scatter(data1(:,12),data1(:,4),'MarkerEdgeColor','k','MarkerFaceColor','r','SizeData',90);
            errorbar(data1(:,12),data1(:,4),data1(:,5),'ro'); % y error bars
            % he1 = herrorbar(data1(:,1),data1(:,2),data1(:,8),'ro'); % x error bars
            ylabel('\delta^{13}C_{CH4} (permil)','fontsize',16)
            set(gca,'fontsize',16)
            % axis ([0 1 -30 5])
            hold off;
            axis equal square;

        % D
            subplot(1,3,2)
            hold on;
            box on;
            plot(F, dCH3D, 'b--','linewidth',1);
            h2 = scatter(data1(:,12),data1(:,6),'MarkerEdgeColor','k','MarkerFaceColor','b','SizeData',90);
            errorbar(data1(:,12),data1(:,6),data1(:,7),'bo');
            % herrorbar(data1(:,1),data1(:,6),data1(:,8),'bo');
            ylabel('\deltaD_{CH4} (permil)','fontsize',16)
            xlabel('Reacted Fraction','fontsize',16)
            set(gca,'fontsize',16)
            % axis ([0 1 -400 -300])
            hold off;
            axis equal square;

        % 13CH3D
            subplot(1,3,3);
            hold on;
            box on;
            plot(F, D18,'k--','linewidth',1);
            plot(F, D182,'k-.','linewidth',1);
            h3 = scatter(data1(:,12),data1(:,8),'MarkerEdgeColor','k','MarkerFaceColor','k','SizeData',90);
            errorbar(data1(:,12),data1(:,8),data1(:,9),'ko');
            %herrorbar(data1(:,1),data1(:,4),data1(:,8),'ko');
            ylabel('\Delta^{13}CH_3D (permil)','fontsize',16)
            set(gca,'fontsize',16)
            axis ([0 1 0 7])
            axis equal square;
