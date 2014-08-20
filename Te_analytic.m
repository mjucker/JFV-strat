function TePV=Te_analytic(lat,pres,days,limlatSH,limlatNH,A0_SH,A0_NH,A1_SH,A1_NH,T_strat,epsS,epsN,p_hsin,p_bdin,filename)
% Computes analytic profile as a function of latitude, pressure, and day of the year, approximating radiatively found profiles
% 
% Inputs are:
% STRATOSPHERE:
% lat: latitude for output grid [deg]
% pres: pressure for output grid [hPa]
% days: Days of the year for which Te should be computed []. Standard: Mid-month for each month
% limlatSH/limlatNH: Latitudes beyond which Te is constant over the poles [deg]. Good value: -80/80
% A0_SH/A0_NH: Te amplitudes at tropopause p=pt [K]. Good value=30/20
% A1_SH/A1_NH: Te amplitude of seasonal cycle of the polar vortices [K] at p=p1. Good value: 60/45
% TROPOSPHERE:
% T_strat: Minimum temperature for Held-Suarez (1994) tropospheric Te [K]. Good value: 200
% epsS/epsN: Seasonal cycle in tropospheric Te for Southern/Northern hemisphere [K]. Standard: 10/10, Good: 10/40
% p_hsin: for pressures above p_hsin [hPa, size latin x days], Held-Suarez
%   is used. Standard: 100hPa
% p_bdin: for pressures below p_bdin [hPa, size latin x days], Tin is used.
%   standard: 100hPa
%   for p_bdin < p < p_hsin, linear interpolation between the two
% filename, if present, is the name of the file containing final tau. 
%   For usage with FMS, this should be [some directory/]temp.nc.
%   If empty, no file will be created.
%
% Routine written by Martin Jucker, please refer to M. Jucker, S.
% Fueglistaler, and G.K. Vallis (2014): [JGR REFERENCE HERE]

% Some constants. p0 is the reference pressure for log-p coordinates
p0=1e3;
Rd=287.04;
cp=1004;
kappa=Rd/cp;
grav=9.81;
A0s=15; %equator-to-pole summer amplitude [K]

p1=1; %reference pressure for SC amplitude [hPa]: this is where the amplitude equals A_SH or A_NH
pt=100;%p0; %reference pressure where equinox meridional profile vanishes

%be sure inputs are in the shape we need
if(pres(end)<pres(1))
    pres = pres(end:-1:1);
end
lat=lat(:);pres=pres(:);

t_length=length(days);

Te=zeros(length(lat),length(pres));

%% Tropical profile: polynome based, one below 1hPa, another above
II=find(pres < 1);
x=log(pres/p0);
poly=-0.537*x.^4 -9.65*x.^3 -60.6*x.^2 -174*x +19.8;
poly2=0.668*x.^3 +22.3*x.^2 +248*x +1.16e3;
for l=1:length(lat)
    Te(l,:)   = max(T_strat,poly);
    Te(l,II)   = max(poly(II),poly2(II));
end

%% Stratospheric Seasonal cycle: Equinox symmetric polynome, amplitudes linear in lat and lnp

% equinox equator to pole shape
poly = -1.96e-9*lat.^4                 -1.15e-5*lat.^2              +1; %North-South symmetric
x = log(pres/pt)/log(p1/pt);
P=(poly - 1)*x' + 1;
II = find(pres < p1);
P(:,II) = P(:,II(end))*ones(size(II))'; %Don't let TOA get too cold: leave P constant above p1
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
A(:,II) = A(:,II(end)+1)*ones(size(II))'; % choice B: constant amplitude above p1
%Amplitude linear in latitude and log-pressure, but different strength
%in summer
AS = -abs(lat)/90*(gammaS*log(pres'/pt) + A0s);
AS(:,II) = AS(:,II(end)+1)*ones(size(II))'; % choice B: constant amplitude above p1

Te0 = zeros(size(Te,1),size(Te,2),t_length); %this will be the stratospheric Te
 for d=1:t_length  
    cosNH = cos(2*pi*(days(d      ))/365); %days(d)-349 for realistic solstice, but needs change in Held-Suarez as well
    cosSH = cos(2*pi*(days(d-182.5))/365); %days(d)-166.5 for realistic solstice
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
TeHS=zeros(size(Te0));

T0=315;
delT=60;
delv=10;
p0=1e3;
kap  = 287.04/1004;
sigma=pres/p0;

for l=1:length(lat)
    for d=1:t_length
        for k=1:length(pres)
            if(lat(l)<=0)
                eps=epsS*cos(2*pi*days(d)/365);
            else
                eps=epsN*cos(2*pi*days(d)/365);
            end
            TeHS(l,k,d)=max(T_strat,(T0 - delT*(sin(lat(l)*pi/180).^2) - eps*sin(lat(l)*pi/180)- delv*log(sigma(k)).*(cos(lat(l)*pi/180).^4)).*(sigma(k)^kap));
        end
    end
end

%% construct complete Te
TePV = zeros(size(Te0));
for d=1:t_length
    for j=1:length(lat)
        I=find(pres > p_hsin(j,d));
        TePV(j,I,d) = TeHS(j,I,d);
        II=find(pres <= p_hsin(j,d) & pres > p_bdin(j,d));
        beta = (pres(II)-p_hsin(j,d))/(p_bdin(j,d)-p_hsin(j,d));
        TePV(j,II,d) = beta'.*Te0(j,II,d) + (1-beta').*TeHS(j,II,d);
        III=find(pres <= p_bdin(j,d));
        TePV(j,III,d) = Te0(j,III,d);
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


%% Now write Te into a netCDF file to be used as input to FMS
if ( exist('filename','var'))

    lonb = linspace(0,359,129)';
    difflon = diff(lonb);
    lon = cumsum(difflon)-difflon(1)/2;
    dlat = lat(2)-lat(1);
    latb = lat - dlat/2;
    latb(end+1) = lat(end)+dlat/2;
    phalf = zeros(length(pres)+1,1);
    phalf(2:end-1) = pres(2:end) - diff(pres)/2;
    phalf(end) = 1e3;
    
    
    Te = zeros(length(lon),length(lat),length(pres),length(days));
    for l=1:length(lon)
        Te(l,:,:,:)=TePV;
    end
    
    ncid = netcdf.create(filename,'64BIT_OFFSET');
    % define dimensions and variables
    lond_id = netcdf.defDim(ncid,'lon',length(lon));
    lonbd_id= netcdf.defDim(ncid,'lonb',length(lonb));
    latd_id = netcdf.defDim(ncid,'lat',length(lat));
    latbd_id= netcdf.defDim(ncid,'latb',length(latb));
    pd_id   = netcdf.defDim(ncid,'pfull',length(pres'));
    phd_id  = netcdf.defDim(ncid,'phalf',length(phalf));
    td_id  =  netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));
    
    lon_id  = netcdf.defVar(ncid,'lon','double',lond_id);
    netcdf.putAtt(ncid,lon_id,'long_name','longitude');
    netcdf.putAtt(ncid,lon_id,'units','degrees_E');
    netcdf.putAtt(ncid,lon_id,'cartesian_axis','X');
    netcdf.putAtt(ncid,lon_id,'edges','lonb');
    lonb_id = netcdf.defVar(ncid,'lonb','double',lonbd_id);
    netcdf.putAtt(ncid,lonb_id,'long_name','longitude_edges');
    netcdf.putAtt(ncid,lonb_id,'units','degrees_E');
    netcdf.putAtt(ncid,lonb_id,'cartesian_axis','X');
    lat_id  = netcdf.defVar(ncid,'lat','double',latd_id);
    netcdf.putAtt(ncid,lat_id,'long_name','latitude');
    netcdf.putAtt(ncid,lat_id,'units','degrees_N');
    netcdf.putAtt(ncid,lat_id,'cartesian_axis','Y');
    netcdf.putAtt(ncid,lat_id,'edges','latb');
    latb_id = netcdf.defVar(ncid,'latb','double',latbd_id);
    netcdf.putAtt(ncid,latb_id,'long_name','latitude_edges');
    netcdf.putAtt(ncid,latb_id,'units','degrees_N');
    netcdf.putAtt(ncid,latb_id,'cartesian_axis','Y');
    p_id    = netcdf.defVar(ncid,'pfull','double',pd_id);
    netcdf.putAtt(ncid,p_id,'long_name','approx full pressure level');
    netcdf.putAtt(ncid,p_id,'units','hPa');
    netcdf.putAtt(ncid,p_id,'cartesian_axis','Z');
    netcdf.putAtt(ncid,p_id,'positive','down');
    netcdf.putAtt(ncid,p_id,'edges','phalf');
    ph_id   = netcdf.defVar(ncid,'phalf','double',phd_id);
    netcdf.putAtt(ncid,ph_id,'long_name','approx half pressure level');
    netcdf.putAtt(ncid,ph_id,'units','hPa');
    netcdf.putAtt(ncid,ph_id,'cartesian_axis','Z');
    t_id    = netcdf.defVar(ncid,'time','double',td_id);
    netcdf.putAtt(ncid,t_id,'long_name','time');
    netcdf.putAtt(ncid,t_id,'cartesian_axis','T');
    netcdf.putAtt(ncid,t_id,'units','days since 0000-00-00 00:00:00');
    %netcdf.putAtt(ncid,t_id,'calendar_type','NO_CALENDAR'); %THIRTY_DAY_MONTHS %THIRTY_DAY_MONTHS
    %netcdf.putAtt(ncid,t_id,'calendar','NO_CALENDAR'); %writetempin_time.m
    
    T_id  = netcdf.defVar(ncid,'temp','double',[lond_id,latd_id,pd_id,td_id]);
    netcdf.putAtt(ncid,T_id,'long_name','radiative relaxation temperature');
    netcdf.putAtt(ncid,T_id,'units','deg_K');
    netcdf.putAtt(ncid,T_id,'valid_range',[0,500]);
    netcdf.putAtt(ncid,T_id,'missing_value',-1.e+10);
    
    
    netcdf.endDef(ncid);
    % write fields
    netcdf.putVar(ncid,lon_id,lon);
    netcdf.putVar(ncid,lonb_id,lonb);
    netcdf.putVar(ncid,lat_id,lat);
    netcdf.putVar(ncid,latb_id,latb);
    netcdf.putVar(ncid,p_id,pres');
    netcdf.putVar(ncid,ph_id,phalf);
    netcdf.putVar(ncid,t_id,0,length(days),days);
    netcdf.putVar(ncid,T_id,[0,0,0,0],[size(Te,1),size(Te,2),size(Te,3),size(Te,4)],Te);
    
    
    netcdf.close(ncid)
    
end
end

