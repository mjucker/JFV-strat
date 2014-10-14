function writetempin(latin,pin,days,Tin,T_strat,epsS,epsN,p_hsin,p_bdin,filename,basefile)
% writes given 3D (latitude x pressure x time) temperature profile into
% netCDF format to be used as input to the FMS Newtonian cooling model
%
% latin: latitudes [degrees] of the grid points
% pin:  pressure levels [hPa] of the grid points
% days: Days of the year for which Te should be computed []. Standard: Mid-month for each month
% Tin: Relaxation temperature profiles on latin x pin x days grid
% Held-Suarez:
% T_strat: Minimum temperature for Held-Suarez (1994) tropospheric Te [K]. Good value: 200
% epsS/epsN: Seasonal cycle in tropospheric Te for Southern/Northern hemisphere [K]. Standard: 10/10, Good: 10/40
% p_hsin: for pressures above p_hsin [hPa, scalar or size latin x days], Held-Suarez
%   is used. Standard: 100hPa
% p_bdin: for pressures below p_bdin [hPa, scalar or size latin x days], Tin is used.
%   standard: 100hPa
%   for p_bdin < p < p_hsin, linear interpolation between the two
% GENERAL:
% filename: Output file name. As code input, must be called [dir]/temp.nc
% basefile:  FMS-like file to read output grid. If empty, use input grid
%
% Routine written by Martin Jucker, please make reference to M. Jucker, S.
% Fueglistaler, and G.K. Vallis (2013): "Maintenance of the Stratospheric Structure in an
%   Idealized General Circulation Model", J. Atmos. Sci. 70, 3341. DOI: 10.1175/JAS-D-12-0305.1


%be sure inputs are in the shape we need
if(length(pin) ~= size(Tin,2))
    Tin = permute(Tin,[2,1,3]);
end
if(pin(end)<pin(1))
    pin = pin(end:-1:1);
    Tin = Tin(:,end:-1:1,:);
end
latin=latin(:);pin=pin(:);
if(length(p_hsin) == 1)
    p_hsin = p_hsin*ones(length(latin),length(days));
end
if(length(p_bdin) == 1)
    p_bdin = p_bdin*ones(length(latin),length(days));
end

grav = 9.81;
Rd   = 287.04;
T_range = [100,400]; 


%% read fms dimensions and interpolate Teq onto them (if file exists)
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
    phalf(1)=0.;
    for l=1:length(pfull)-1
        phalf(1+l) = 0.5*(pfull(l)+pfull(l+1));
    end
    phalf(length(pfull)+1)=0.5*(3*pfull(end)-pfull(end-1));
end

t_length=size(Tin,3);
if(~exist('days','var'))
    dayint = 365/t_length;
    days   = dayint/2:dayint:365;
else
    if(length(days)~=t_length)
        'length of days must be equal length of input Temperature'
        return
    end
end

%% add Held-Suarez troposphere on output grid
TeHS=zeros(length(lat),length(pfull),t_length);

T0=315;
delT=60;
delv=10;
p0=1e3;
kap  = 287.04/1004;
sigma=pfull/p0;

for l=1:length(lat)
    for d=1:t_length
        for k=1:length(pfull)
            if(lat(l)<=0)
                eps=epsS*cos(2*pi*days(d)/365);
            else
                eps=epsN*cos(2*pi*days(d)/365);
            end
            TeHS(l,k,d)=max(T_strat,(T0 - delT*(sin(lat(l)*pi/180).^2) - eps*sin(lat(l)*pi/180)- delv*log(sigma(k)).*(cos(lat(l)*pi/180).^4)).*(sigma(k)^kap));
        end
    end
end

%% interpolate input T onto output grid

TePV = zeros(length(lon),length(lat),length(pfull),t_length);
JJ = find(lat >= min(latin) & lat <= max(latin)); %avoid extrapolating
% K = find(pfull >= min(pin) & pfull <= max(pin)); %avoid extrapolating

for i=1:length(lon)
for d=1:t_length
    if(length(latin)>1)
        T_s_tmp = zeros(length(lat),length(pin));
        for k=1:length(pin);
            T_s_tmp(JJ,k) = interp1(latin,Tin(:,k,d),lat(JJ),'linear','extrap');
            T_s_tmp(1:JJ(1),k) = T_s_tmp(JJ(1),k);
            T_s_tmp(JJ(end):end,k) = T_s_tmp(JJ(end),k);
        end
        for j=1:length(lat)
%             TePV(i,j,K,d) = interp1(pin,T_s_tmp(j,:),pfull(K));
            TePV(i,j,:,d) = interp1(pin,T_s_tmp(j,:),pfull,'linear','extrap');
        end
%         TePV(i,j,K(end):end,d) = TePV(i,j,K(end),d);
        clear T_s_tmp
        %treat outside of input domain
        %Ttopmean = mean(TePV(:,min(K),d),1);
%         if(T_strat>0 && min(K)>1)
%             for kk=1:min(K)
%                 TePV(i,:,kk,d)=squeeze(TePV(i,:,min(K),d))' - 10*(pfull(kk)-pfull(min(K)))/(pfull(1)-pfull(min(K)))*ones(length(lat),1); %decrease by 10K
%                 %T_s(:,kk,d)=T_s(:,min(K),d) + (pfull(kk)-pfull(min(K)))/(pfull(1)-pfull(min(K)))*(Ttopmean - T_s(:,min(K),d)); %towards latitudinal mean on top layer
%                 %T_s(:,kk,d) = T_s(:,min(K),d); %no vertical gradient
%             end
%         end
    else
        for j=1:length(lat)
            TePV(i,j,:,d) = interp1(pin,Tin(:,d),pfull);
        end
    end
end
end
p_hs = zeros(length(lat),t_length);
p_bd = zeros(length(lat),t_length);
for d=1:t_length
    p_hs(JJ,d) = interp1(latin,p_hsin(:,d),lat(JJ),'linear','extrap');
    p_bd(JJ,d) = interp1(latin,p_bdin(:,d),lat(JJ),'linear','extrap');
    p_hs(1:JJ(1),d) = p_hs(JJ(1),d);
    p_bd(1:JJ(1),d) = p_bd(JJ(1),d);
    p_hs(JJ(end):end,d) = p_hs(JJ(end),d);
    p_bd(JJ(end):end,d) = p_bd(JJ(end),d);
end

%% construct complete Te
Te = zeros(size(TePV));
TeHS = permute(repmat(TeHS,[1,1,1,length(lon)]),[4,1,2,3]);
if(T_strat == 0)
    TeHS=zeros(size(TeHS));
    T_range=[-1e10,1e10];
end
for d=1:t_length
    for j=1:length(lat)
        I=find(pfull > p_hs(j,d));
        Te(:,j,I,d) = TeHS(:,j,I,d);
        II=find(pfull <= p_hs(j,d) & pfull > p_bd(j,d));
        beta = (pfull(II)-p_hs(j,d))/(p_bd(j,d)-p_hs(j,d));
        beta = reshape(permute(repmat(beta,[1,length(lon)]),[2,1]),[length(lon),1,length(II),1]);
        Te(:,j,II,d) = beta.*TePV(:,j,II,d) + (1-beta).*TeHS(:,j,II,d);
        III=find(pfull <= p_bd(j,d));
        Te(:,j,III,d) = TePV(:,j,III,d);
    end
end



%keep Te in range
I = find(Te(:) < T_range(1));
Te(I) = T_range(1);
I = find(Te(:) > T_range(2));
Te(I) = T_range(2);

figure;
vv = 100:5:320;
[h,c]=contourf(lat,pfull,double(squeeze(mean(squeeze(mean(Te,4)),1)))',vv);
clabel(h,c);
colorbar;
set(gca,'ydir','rev');
set(gca,'yscale','log');
title('Temporal mean Te');
if(t_length > 1 && t_length <= 12)
    for t=1:t_length
        figure;
        Ttmp = double(squeeze(mean(TePV(:,:,:,t),1)));
        pcolor(lat,pfull,Ttmp')
        shading interp
        hold on
        [h,c]=contour(lat,pfull,Ttmp',vv,'-k');
        clabel(h,c);
        colorbar;
        set(gca,'ydir','rev');
        set(gca,'yscale','log');
        title(['Te, t = ',num2str(t)]);
    end
end

%% write new netcdf file with input teq
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

T_id  = netcdf.defVar(ncid,'temp','double',[lond_id,latd_id,pd_id,td_id]);
netcdf.putAtt(ncid,T_id,'long_name','equilibrium temperature');
netcdf.putAtt(ncid,T_id,'units','deg_K');
netcdf.putAtt(ncid,T_id,'valid_range',[50,400]);
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
netcdf.putVar(ncid,T_id,[0,0,0,0],[size(Te,1),size(Te,2),size(Te,3),size(Te,4)],Te);


netcdf.close(ncid)


end