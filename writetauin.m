function writetauin(latin,pin,days,tauin,tau_strat,p_hsin,p_bdin,filename,basefile)
% writes given 3D (latitude x pressure x time) relaxation time profile into
% netCDF format to be used as input to the FMS Newtonian cooling model
%
% latin: latitudes [degrees] of the grid points
% pin:  pressure levels [hPa] of the grid points
% days: Days of the year for which Te should be computed []. Standard: Mid-month for each month
% tauin: Relaxation time profiles on latin x pin x days grid
% Held-Suarez:
% tau_strat: maximum tau for Held-Suarez (1994) tropospheric tau [days]. Good value: 40
% p_hsin: for pressures above p_hsin [hPa, size latin x days], Held-Suarez
%   is used. Standard: 100hPa
% p_bdin: for pressures below p_bdin [hPa, size latin x days], Tin is used.
%   standard: 100hPa
%   for p_bdin < p < p_hsin, linear interpolation between the two
% GENERAL:
% filename: Output file name. As code input, must be called [dir]/tau.nc
% basefile:  FMS-like file to read output grid. If empty, use input grid
%
% Routine written by Martin Jucker, please make reference to M. Jucker, S.
% Fueglistaler, and G.K. Vallis (2013): "Maintenance of the Stratospheric Structure in an
%   Idealized General Circulation Model", J. Atmos. Sci. 70, 3341. DOI: 10.1175/JAS-D-12-0305.1

%be sure inputs are in the shape we need
if(length(pin) ~= size(tauin,2))
    tauin = permute(tauin,[2,1,3]);
end
if(pin(end)<pin(1))
    pin = pin(end:-1:1);
    tauin = tauin(:,end:-1:1,:);
end
latin=latin(:);pin=pin(:);
if(length(pin) ~= size(tauin,2))
    tauin = permute(tauin,[2,1,3]);
end

grav = 9.81;
Rd   = 287.04;
secperday = 86400;


%% read fms dimensions and interpolate tau onto them
if( exist('basefile','var') )
    if( exist(basefile,'file') )

        ncid = netcdf.open(basefile,'NC_NOWRITE');
        varid = netcdf.inqVarID(ncid,'lon');
        lon=netcdf.getVar(ncid,varid);
        varid = netcdf.inqVarID(ncid,'lonb');
        lonb=netcdf.getVar(ncid,varid);
        varid = netcdf.inqVarID(ncid,'lat');
        lat=netcdf.getVar(ncid,varid);
        varid = netcdf.inqVarID(ncid,'latb');
        latb=netcdf.getVar(ncid,varid);
        varid = netcdf.inqVarID(ncid,'pfull');
        pfull=netcdf.getVar(ncid,varid);
        pfull=pfull(:)';
        varid = netcdf.inqVarID(ncid,'phalf');
        phalf=netcdf.getVar(ncid,varid);

        netcdf.close(ncid);
    else
        disp(['The file ',basefile,' does not exist! Have to stop here'])
        return
    end
else
    dummy=linspace(0,360,3);lon=dummy(1:end-1)';
    lonb=dummy-0.5*dummy(2)';
    latb=linspace(-90,90,length(latin)+1)';lat=0.5*(3*latb(1)-latb(2))+cumsum(diff(latb));
    pfull=pin';
    %phalf=[0.5*(3*pfull(2)-pfull(1)),0.5*(3*pfull(2)-pfull(1))+cumsum(diff(pfull)),0.5*(3*pfull(end)-pfull(end-1))];
    phalf(1)=0.;
    for l=1:length(pfull)-1
        phalf(1+l) = 0.5*(pfull(l)+pfull(l+1));
    end
    phalf(length(pfull)+1)=0.5*(3*pfull(end)-pfull(end-1));
end

t_length=size(tauin,3);
if(~exist('days','var'))
    dayint = 365/t_length;
    days   = dayint/2:dayint:365;
else
    if(length(days)~=t_length)
        'length of days must be equal length of input tau'
        return
    end
end

tdamp = zeros(length(lon),length(lat),length(pfull),t_length);

%% compute Held-Suarez tau

tauhs = heldsuarez_tau('',tau_strat,lat,pfull);

%% interpolate inputs onto code grid
tau  = zeros(length(lat),length(pfull),t_length);  
p_hs = zeros(length(lat),t_length);
p_bd = zeros(length(lat),t_length);  
for d=1:t_length
    if(length(latin)>1)
        tau_tmp = zeros(length(lat),length(pin));
        for k=1:length(pin);
            tau_tmp(:,k) = interp1(latin,tauin(:,k,d),lat,'linear','extrap');
        end
%         K = find(pfull >= min(pin) & pfull <= max(pin));
        for j=1:length(lat)
%             tau(j,K,d) = interp1(pin,tau_tmp(j,:),pfull(K));
            tau(j,:,d) = interp1(pin,tau_tmp(j,:),pfull,'linear','extrap');
        end
        if(size(p_hsin,2)>1)
            p_hs(:,d) = interp1(latin,p_hsin(:,d),lat,'linear','extrap');
            p_bd(:,d) = interp1(latin,p_bdin(:,d),lat,'linear','extrap');
        else
            p_hs(:,d) = interp1(latin,p_hsin,lat,'linear','extrap');
            p_bd(:,d) = interp1(latin,p_bdin,lat,'linear','extrap');
        end
        % what to do with layers outside pin?
%         if(min(K)>1)
%             for kk=1:min(K)
%                 %tau(:,kk,d)=tau(:,min(K),d) + 10*(pfull(kk)-pfull(min(K)))/(pfull(1)-pfull(min(K)))*ones(length(lat),1);
%                 tau(:,kk,d)=tau(:,min(K),d);
%             end
%         end
%         if(max(K)<length(pfull))
%             for kk=max(K):length(pfull)
%                 tau(:,kk,d)=tau(:,max(K),d);
%             end
%         end
        clear tau_tmp
    else
        for j=1:length(lat)
            tau(j,:,d) = interp1(pin,tauin(:,d),pfull,'linear','extrap');
        end
    end
end
%% concatenate tau in troposphere (tauhs) stratosphere (tauin,tau)


for d=1:t_length
    for i=1:length(lon)
        for j=1:length(lat)
            tdamp(i,j,:,d) = tau(j,:,d);
        end
        for j=1:length(lat)
            I=find(pfull > p_hs(j,d));
            tdamp(i,j,I,d) = tauhs(j,I);
            J=find(pfull <= p_hs(j,d) & pfull > p_bd(j,d));
            beta = (pfull(J)-p_hs(j,d))/(p_bd(j,d)-p_hs(j,d));
            tdamp(i,j,J,d) = beta.*tau(j,J,d) + (1-beta).*tauhs(j,J);
        end
    end
end

%% visualize resulting damping times

figure;
vv = linspace(min(tdamp(:)),max(tdamp(:)),20);
contourf(lat,pfull,double(squeeze(mean(squeeze(mean(tdamp,4)),1)))',vv)
colorbar;
set(gca,'ydir','rev');
set(gca,'yscale','log');
title('tdamp [days]');
if(t_length > 1 && t_length <= 12)
    for t=1:t_length
        figure;
        ttmp = double(squeeze(mean(tdamp(:,:,:,t),1)));
        [h,c]=contourf(lat,pfull,ttmp',vv);
        clabel(h,c);
        colorbar;
        set(gca,'ydir','rev');
        set(gca,'yscale','log');
        title(['tau, t = ',num2str(days(t))]);
    end
end

disp(['tau_min = ',num2str(min(tdamp(:))),', tau_max = ',num2str(max(tdamp(:)))]);
%% convert to seconds

tdamp = tdamp*secperday;

%% write new netcdf file with input tau
ncid = netcdf.create(filename,'64BIT_OFFSET');
% define dimensions and variables
lond_id = netcdf.defDim(ncid,'lon',length(lon));
lonbd_id= netcdf.defDim(ncid,'lonb',length(lonb));
latd_id = netcdf.defDim(ncid,'lat',length(lat));
latbd_id= netcdf.defDim(ncid,'latb',length(latb));
pd_id   = netcdf.defDim(ncid,'pfull',length(pfull));
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
netcdf.putAtt(ncid,T_id,'valid_range',[0,100]*secperday);
netcdf.putAtt(ncid,T_id,'missing_value',-1.e+10);


netcdf.endDef(ncid);
% write fields
netcdf.putVar(ncid,lon_id,lon);
netcdf.putVar(ncid,lonb_id,lonb);
netcdf.putVar(ncid,lat_id,lat);
netcdf.putVar(ncid,latb_id,latb);
netcdf.putVar(ncid,p_id,pfull);
netcdf.putVar(ncid,ph_id,phalf);
netcdf.putVar(ncid,t_id,0,length(days),days);
netcdf.putVar(ncid,T_id,[0,0,0,0],[size(tdamp,1),size(tdamp,2),size(tdamp,3),size(tdamp,4)],tdamp);


netcdf.close(ncid)


end
