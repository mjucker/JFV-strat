function TePV=Te_analytic_public(lat,pres,T_strat,epsS,epsN,limlatSH,limlatNH,A0_SH,A1_SH,A0_NH,A1_NH,days)
% Routine written by Martin Jucker, version February 2013
% Computes analytic profile as a function of latitude, pressure, and day of the year, approximating radiatively found profiles
% 
% Inputs are:
% lat: latitude for output grid [deg]
% pres: pressure for output grid [hPa]
% T_strat: Minimum temperature for Held-Suarez (1994) tropospheric Te [K]. Good value: 200
% epsS/epsN: Seasonal cycle in tropospheric Te for Southern/Northern hemisphere [K]. Standard: 10/10, Good: 10/40
% limlatSH/limlatNH: Latitudes beyond which Te is constant over the poles [deg]. Good value: -80/80
% A0_SH/A0_NH: Te amplitudes at tropopause p=pt [K]. Good value=30/20
% A1_SH/A1_NH: Te amplitude of seasonal cycle of the polar vortices [K] at p=p1. Good value: 60/45
% days: Days of the year for which Te should be computed []. Standard: Mid-month for each month

% Some constants. p0 is the reference pressure for log-p coordinates
p0=1e3;
Rd=287.04;
cp=1004;
kappa=Rd/cp;
grav=9.81;
A0s=15; %equator-to-pole summer amplitude [K]

p1=1; %reference pressure for SC amplitude [hPa]: this is where the amplitude equals A_SH or A_NH
pt=100;%p0; %reference pressure where equinox meridional profile vanishes

t_length=length(days);

Te=zeros(length(lat),length(pres));

%% Tropical profile: polynome based, one below 1hPa, another above
J=find(pres < 1);
x=log(pres/p0);
poly=-0.537*x.^4 -9.65*x.^3 -60.6*x.^2 -174*x +19.8;
poly2=0.668*x.^3 +22.3*x.^2 +248*x +1.16e3;
for l=1:length(lat)
    Te(l,:)   = max(T_strat,poly);
    Te(l,J)   = max(poly(J),poly2(J));
end

%% Stratospheric Seasonal cycle: Equinox symmetric polynome, amplitudes linear in lat and lnp

% equinox equator to pole shape
poly = -1.96e-9*lat.^4                 -1.15e-5*lat.^2              +1; %North-South symmetric
x = log(pres/pt)/log(p1/pt);
P=(poly - 1)*x' + 1;
J = find(pres < p1);
P(:,J) = P(:,J(end))*ones(size(J))'; %Don't let TOA get too cold: leave P constant above p1
JJ = find(pres >= pt);
P(:,JJ) = ones(size(P(:,JJ))); %no heating below cooling effect of poly

% SC amplitude
NH=find(lat >= 0);
SH=find(lat <= 0);

gammaSH = (A1_SH-A0_SH)/(log(p1/p0)-log(pt/p0)); % amplitude is A0 at p=pt, and A1 at p=p1
gammaNH = (A1_NH-A0_NH)/(log(p1/p0)-log(pt/p0)); % amplitude is A0 at p=pt, and A1 at p=p1
gammaS = (A0s-0)/(log(p1/p0)-log(pt/p0)); % amplitude is 0 at p=pt, and A0s at p=p1 in summer

%Amplitde is linear in latitude and log-pressure
A = zeros(size(Te));
A(NH,:) = -abs(lat(NH))/90*(gammaNH*log(pres'/pt) + A0_NH);
A(SH,:) = -abs(lat(SH))/90*(gammaSH*log(pres'/pt) + A0_SH);
A(:,J) = A(:,J(end)+1)*ones(size(J))'; % choice B: constant amplitude above p1
%Amplitude linear in latitude and log-pressure, but different strength
%in summer
AS = -abs(lat)/90*(gammaS*log(pres'/pt) + A0s);
AS(:,J) = AS(:,J(end)+1)*ones(size(J))'; % choice B: constant amplitude above p1


 for d=1:t_length  
    cosNH = cos(2*pi*(days(d)-349)/365); %355
    cosSH = cos(2*pi*(days(d)-166.5)/365); %172
    if(cosNH>0)
        Te0(NH,:,d) = Te(NH,:).*P(NH,:) + A(NH,:)*cosNH;
    else
        Te0(NH,:,d) = Te(NH,:).*P(NH,:) + (AS(NH,:) - Te(NH,:).*(1-P(NH,:)))*cosNH;
    end
    if(cosSH>0)
        Te0(SH,:,d) = Te(SH,:).*P(SH,:) + A(SH,:)*cosSH;
    else
        Te0(SH,:,d) = Te(SH,:).*P(SH,:) + (AS(SH,:) - Te(SH,:).*(1-P(SH,:)))*cosSH;
    end
    %Te0(:,:,d)  = min(Te(:,:),Te0(:,:,d)); % remove latitudinal structure in summer hemisphere
end
Slat = find(lat < limlatSH); % Te constant poleward of limlat
for l=1:length(Slat)
    Te0(Slat(l),:,:) = Te0(Slat(end),:,:);
end
Nlat = find(lat > limlatNH);
for l=1:length(Nlat)
    Te0(Nlat(l),:,:) = Te0(Nlat(1),:,:);
end

%% add Held-Suarez troposphere
JJ = find(pres >= 100);

TePV=zeros(size(Te0));
TePV(:,1:JJ(1),:)=Te0(:,1:JJ(1),:);

T0=315;
delT=60;
delv=10;
p0=1e3;
kap  = 287.04/1004;
sigma=pres/p0;

for k=1:length(JJ)
    for l=1:length(lat)
        for d=1:t_length
            if(lat(l)<=0)
                eps=epsS*cos(2*pi*days(d)/365);
            else
                eps=epsN*cos(2*pi*days(d)/365);
            end
            TePV(l,JJ(k),d)=max(T_strat,(T0 - delT*(sin(lat(l)*pi/180).^2) - eps*sin(lat(l)*pi/180)- delv*log(sigma(JJ(k))).*(cos(lat(l)*pi/180).^4)).*(sigma(JJ(k))^kap));
        end
    end
end


for d=1:t_length
    figure;
    contourf(lat,pres,squeeze(TePV(:,:,d))',[150:10:350]);
    set(gca,'yscale','log');
    set(gca,'ydir','rev');
    set(gca,'xtick',[-90:30:90]);
    set(gca,'clim',[180,320]);
    colorbar;
    title(['day ',num2str(days(d))]);
end

end

