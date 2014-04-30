function tau = tau_analytic(lat,pres,days,taut,taups,taupn,lbroad,tau_strat,p_hsin,p_bdin,filename)
% lat: latitude for output grid [deg]
% pres: pressure for output grid [hPa]
% days: Days of the year for which Te should be computed []. Standard: Mid-month for each month
% taut: value for tau at the equator and 100hPa; [d], length of days. Good:
%   40
% taups: value for tau at the south pole and 100hPa; [d], length of days.
%   Good: 20
% taupn: value for tau at the north pole and 100hPa; [d], length of days.
%   Good: 20
% lbroad: width of the tropical Gaussian at 100hPa, [degrees]. Good: 30
% tau_strat: value for tau in the Held-Suarez stratosphere [d]. Good & standard: 40
% p_hsin: for pressures above p_hsin [hPa, size latin x days], Held-Suarez
%   is used. Standard: 100hPa
% p_bdin: for pressures below p_bdin [hPa, size latin x days], Tin is used.
%   standard: 100hPa
%   for p_bdin < p < p_hsin, linear interpolation between the two
% filename, if present, is the name of the file containing final tau. When
%   used as input file, needs to be called tau.nc
%
% Routine written by Martin Jucker, please refer to M. Jucker, S.
% Fueglistaler, and G.K. Vallis (2014): [JGR REFERENCE HERE]

%some constants
p0  = 1e3; %do not change this
pt = 1e2;
taumin = 5;

pres=pres(:)';
lat=lat(:);
taut=taut(:);
taups=taups(:);
taupn=taupn(:);
SH=find(lat<0);
NH=find(lat>=0);

x = log(pres/p0);
poly = [0.045 1.38 15.9 81.6 162]; %polynomial as fitted from radiative tau
polynom = polyval(poly,x);
[m,I] = min(polynom); %correct for poly exploding at top
polynom(1:I) = polynom(I)*ones(size(polynom(1:I)));
polymin = min(polynom);
plast = find(pres < pt,1,'last');
polyt  = polynom(plast)-polymin; %grid point just above pt
polynorm = (polynom - polymin)/polyt; %normalized in [0,1]
polynorm(plast:end)=polynorm(plast)*ones(size(polynorm(plast:end))); %keep it one below pt

%Pn=(taut-taumin)./taut*polynorm+taumin./taut*ones(size(polynorm)); %adjust vertical profile to be between taut and taumin in the tropics

taupt = zeros(length(taut),length(lat));
lats=lat(SH)';
latn=lat(NH)';
expfun = zeros(size(lat));
expfun(SH) = exp(-(lats/lbroad).^2);
expfun(NH) = exp(-(latn/lbroad).^2);
taupt(:,SH) = taups*ones(size(lats)) + (taut - taups)*expfun(SH)'; %meridional tau profile at p=pt
taupt(:,NH) = taupn*ones(size(latn)) + (taut - taupn)*expfun(NH)'; %meridional tau profile at p=pt


taus=zeros(length(lat),length(pres),length(taut));
for t=1:length(taut)
    taus(:,:,t) = (taupt(t,:) - taumin)'*polynorm + taumin;
end

%% compute Held-Suarez tau
ka=1/tau_strat;
ks=1/4;
sigma=pres/p0;
sigmab=0.7;

kv = (ks-ka)*max(0,(sigma-sigmab)/(1-sigmab));

tauhs=zeros(length(lat),length(pres));
for k=1:length(pres)
    kT = ka + kv(k)*cos(lat*pi/180).^4;
    tauhs(:,k) = 1./kT;
end


%% construct complete tau, including Held-Suarez tau
tau = zeros(size(taus));
for d=1:length(days)
    for j=1:length(lat)
        I=find(pres > p_hsin(j,d));
        tau(j,I,d) = tauhs(j,I);
        II=find(pres <= p_hsin(j,d) & pres > p_bdin(j,d));
        beta = (pres(II)-p_hsin(j,d))/(p_bdin(j,d)-p_hsin(j,d));
        tau(j,II,d) = beta.*taus(j,II,d) + (1-beta).*tauhs(j,II);
        III=find(pres <= p_bdin(j,d));
        tau(j,III,d) = taus(j,III,d);
    end
    
end


for d=1:size(tau,3)
    figure;
    contourf(lat,pres,squeeze(tau(:,:,d))',[1:10,20:5:50]);
    set(gca,'yscale','log');
    set(gca,'ydir','rev');
    set(gca,'xtick',[-90:30:90]);
    set(gca,'clim',[0,50]);
    colorbar;
    title(['timestep ',num2str(d)]);
end

%% write new netcdf file with input tau
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
    
    
    tdamp = zeros(length(lon),length(lat),length(pres),length(days));
    for l=1:length(lon)
        tdamp(l,:,:,:)=tau;
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
    
    T_id  = netcdf.defVar(ncid,'tau','double',[lond_id,latd_id,pd_id,td_id]);
    netcdf.putAtt(ncid,T_id,'long_name','radiative relaxation time');
    netcdf.putAtt(ncid,T_id,'units','seconds');
    netcdf.putAtt(ncid,T_id,'valid_range',[0,100]*86400);
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
    netcdf.putVar(ncid,T_id,[0,0,0,0],[size(tdamp,1),size(tdamp,2),size(tdamp,3),size(tdamp,4)],tdamp);
    
    
    netcdf.close(ncid)
    
end
end