%%% Project: Sprat black sea

% Research question: how will +1.5 deg C and +5 deg C warming scenarios
% affect the MSY of sprat fisheries in the Black Sea?

%Warming function
T=@(t) 15+t/20;

%lower warming IPCC
R=@(t) 15+t*1.5/100;

%von Bertalanffy growth sprat
ag=linspace(0,7,10000);
lv=13.8*(1-exp(-0.42*(ag+1.09)));
plot(ag,lv,'m','linewidth', 2,'color',[0.4940, 0.1840, 0.5560]);
xlabel('Length (cm)');
ylabel('Age (years)');

% Model population - 1.5 degree warming
% Time increment: 1 year
% Model over 100 years

yhis = zeros(6,100);
yhis(:,1) = [3*10^12;3*10^9;3*10^8;3*10^8;3*10^7;3*10^7]; % column vector

for p = 1:99
        S2=yhis(3,p)+yhis(4,p)+yhis(5,p)+yhis(6,p);
        T2=@(t) 15+t*1.5/100;
        T2(p);
        Teff2=@(T2) max(0,1/5*heaviside(-T2+17.219).*(exp(0.2.*(T2-8))-1)+1/5*heaviside(T2-17.219).*4.*log(-(T2-21)));
        Teff2(T2(p));
        f22= (6000*exp(-0.0000000003*S2))*Teff2(T2(p));
        f32= (8000*exp(-0.0000000003*S2))*Teff2(T2(p));
        f42= (10000*exp(-0.0000000003*S2))*Teff2(T2(p));
        f52= (12000*exp(-0.0000000003*S2))*Teff2(T2(p));
        m232=0.6;
        m342=0.3;
        m452=0.4;
        m552=0.8;
        P232=exp(-m232);
        P342=exp(-m342);
        P452=exp(-m452);
        P552=exp(-m552);
        P012=0.001;
        P122=0.1;
        L2 = [0 0 f22 f32 f42 f52; P012 0 0 0 0 0; 0 P122 0 0 0 0; 0 0 P232 0 0 0; 0 0 0 P342 0 0; 0 0 0 0 P452 P552];
        yhis(:,p+1) = L2*yhis(:,p);
 end

% Plot
plot(log10(yhis(1,:)),'o', 'color', '#52b2bf','MarkerSize',4);
hold on
plot(log10(yhis(2,:)),'o','color','#4B9CD3','MarkerSize',6);
plot(log10(yhis(3,:)),'-.','color','#A1CAF1','linewidth',1.7);
plot(log10(yhis(4,:)),'--','Color','#6CA0DC','linewidth',1.5);
plot(log10(yhis(5,:)),'color','#0F52BA','linewidth',1.5);
plot(log10(yhis(6,:)),'color','#082567','linewidth',1.5);
xlabel('Years')
ylabel('# of individuals (log10 scale)')

yhis = zeros(6,100);
yhis(:,1) = [3*10^12;3*10^9;3*10^8;3*10^8;3*10^7;3*10^7]; % column vector

% Model population - 5 degree warming
% Time increment: 1 year
% Model over 100 years
 for p = 1:99
        S2=yhis(3,p)+yhis(4,p)+yhis(5,p)+yhis(6,p);
        T2=@(t) 15+t*5/100;
        T2(p);
        Teff2=@(T2) max(0,1/5*heaviside(-T2+17.219).*(exp(0.2.*(T2-8))-1)+1/5*heaviside(T2-17.219).*4.*log(-(T2-21)));
        Teff2(T2(p));
        f22= (6000*exp(-0.0000000003*S2))*Teff2(T2(p));
        f32= (8000*exp(-0.0000000003*S2))*Teff2(T2(p));
        f42= (10000*exp(-0.0000000003*S2))*Teff2(T2(p));
        f52= (12000*exp(-0.0000000003*S2))*Teff2(T2(p));
        m232=0.6;
        m342=0.3;
        m452=0.4;
        m552=0.8;
        P232=exp(-m232);
        P342=exp(-m342);
        P452=exp(-m452);
        P552=exp(-m552);
        P012=0.001;
        P122=0.1;
        L2 = [0 0 f22 f32 f42 f52; P012 0 0 0 0 0; 0 P122 0 0 0 0; 0 0 P232 0 0 0; 0 0 0 P342 0 0; 0 0 0 0 P452 P552];
        yhis(:,p+1) = L2*yhis(:,p);
  end


plot(log10(yhis(1,:)),'o', 'color', '#FA8072','MarkerSize',4) %[10:end]
hold on
plot(log10(yhis(2,:)),'o','color','#FF2800','MarkerSize',6)
plot(log10(yhis(3,:)),'-.','color','#D21F3C','linewidth',1.5)
plot(log10(yhis(4,:)),'--','Color','r','linewidth',1.5)
plot(log10(yhis(5,:)),'color','#B80F0A','linewidth',1.5)
plot(log10(yhis(6,:)),'color','#8B0000','linewidth',1.5)
xlabel('Years','FontSize',16)
ylabel('Number of individuals (log10 scale)','FontSize',16)
set(gca,'fontsize',14)

legend('fontsize',14)

%lgd = legend('Warming of 1.5°C','Warming of 5°C','fontsize',14)
%lgd.NumColumns = 8


% Fecundity plot: effect of temperature on fecundity
i=6:0.01:21;
Teff=max(0,1/5*heaviside(-i+17.219).*(exp(0.2.*(i-8))-1)+1/5*heaviside(i-17.219).*4.*log(-(i-21)));
plot(i,Teff,'linewidth',3,'Color','#7E2F8E')
set(gca,'fontsize',18)
xlabel('Annual mean temperature (°C)','fontsize',18)
ylabel('Effect on fecundity (unitless)','fontsize',18)


% Warming scenarios
x=0:0.1:100;
w1=15+x/20;
w2=15+x*1.5/100;
hold on
plot(x,w1,'red','linewidth',2)
plot(x,w2,'blue','linewidth',2)
set(gca,'fontsize',14)
ylabel('Annual mean temperature (°C)','fontsize',16)
xlabel('Time (years)','fontsize',16)
legend('Warming of 5°C','Warming of 1.5°C','fontsize',16)


%%% MSY plot for the two warming scenarios
clf

[f,yield1] = runfishing(1.5/100);

plot(f,yield1,'blue.','MarkerSize',12)
set(gca,'fontsize',14)
hold on

[f,yield2] = runfishing(1/20);

plot(f,yield2,'red.','MarkerSize',12)
set(gca,'fontsize',14)
hold off

xlabel('Fishing pressure (unitless)','fontsize',16)
ylabel('Catch at year 2100 (number of sprats)','fontsize',16)
legend('Warming of 1.5°C','Warming of 5°C','fontsize',14,'Box','off')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%       IF RUNFISHING IS MISSING, THE CODE BELOW CAN BE USED       %%%%
%%%%%       TO CREATE IT. CREATE NEW FUNCTION SCRIPT AND COPY IN.      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ftot2, yieldtot2] = runfishing(slope)

    % Empty vectors
    yieldtot2 = zeros(999, 1);
    ftot2 = zeros(999, 1);

    for z = 1:999
        yhis = zeros(6, 100);
        yhis(:, 1) = [3*10^12; 3*10^9; 3*10^8; 3*10^8; 3*10^7; 3*10^7]; % Initial conditions

        for p = 1:99
            S2 = yhis(3, p) + yhis(4, p) + yhis(5, p) + yhis(6, p);
            T2 = @(t22) 15 + t22 * slope;
            T2_value = T2(p);
            Teff2 = @(T2) max(0, 1/5 * heaviside(-T2 + 17.219) .* (exp(0.2 * (T2 - 8)) - 1) + ...
                                1/5 * heaviside(T2 - 17.219) .* 4 .* log(-(T2 - 21)));
            Teff2_value = Teff2(T2_value);

            f22 = (6000 * exp(-0.0000000003 * S2)) * Teff2_value;
            f32 = (8000 * exp(-0.0000000003 * S2)) * Teff2_value;
            f42 = (10000 * exp(-0.0000000003 * S2)) * Teff2_value;
            f52 = (12000 * exp(-0.0000000003 * S2)) * Teff2_value;

            m232 = 0.6;
            m342 = 0.3;
            m452 = 0.4;
            m552 = 0.8;

            F232 = 0.0003 * z;
            F342 = 0.000325 * z;
            F452 = 0.00035 * z;
            F552 = 0.000375 * z;

            P232 = exp(-m232 - F232);
            P342 = exp(-m342 - F342);
            P452 = exp(-m452 - F452);
            P552 = exp(-m552 - F552);

            P012 = 0.001;
            P122 = 0.1;

            L2 = [0 0 f22 f32 f42 f52;
                  P012 0 0 0 0 0;
                  0 P122 0 0 0 0;
                  0 0 P232 0 0 0;
                  0 0 0 P342 0 0;
                  0 0 0 0 P452 P552];

            yhis(:, p+1) = L2 * yhis(:, p); % Update state
        end

        % Calculate yields
        yield12 = yhis(3, end) * (1 - P232) * F232 / (F232 + m232);
        yield22 = yhis(4, end) * (1 - P342) * F342 / (F342 + m342);
        yield32 = yhis(5, end) * (1 - P452) * F452 / (F452 + m452);
        yield42 = yhis(6, end) * (1 - P552) * F552 / (F552 + m552);
        yieldtot2(z) = yield12 + yield22 + yield32 + yield42;
        
        % Total yield
        ftot2(z) = (F232 + F342 + F452 + F552) / 4;
    end

end

