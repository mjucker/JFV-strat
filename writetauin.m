function writetauin(latin,pin,days,tauin,tau_strat,p_hsin,p_bdin,filename,directory)
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
% directory: directory containing file 'atmos_average.nc' to read FMS
%   dimensions from. If empty, use default T42 grid
%
% Routine written by Martin Jucker, please make reference to M. Jucker, S.
% Fueglistaler, and G.K. Vallis (2013): "Maintenance of the Stratospheric Structure in an
%   Idealized General Circulation Model", J. Atmos. Sci. 70, 3341. DOI: 10.1175/JAS-D-12-0305.1

grav = 9.81;
Rd   = 287.04;
secperday = 86400;


%% read fms dimensions and interpolate tau onto them
if(exist('directory','var') && exist([directory,'/atmos_average.nc'],'file'))
    file=[directory,'/atmos_average_pfull.nc'];

    ncid = netcdf.open(file,'NC_NOWRITE');
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
    lon=[0    2.8125    5.6250    8.4375   11.2500   14.0625   16.8750   19.6875   22.5000   25.3125   28.1250   30.9375   33.7500 ...
   36.5625   39.3750   42.1875   45.0000   47.8125   50.6250   53.4375   56.2500   59.0625   61.8750   64.6875   67.5000   70.3125 ...
   73.1250   75.9375   78.7500   81.5625   84.3750   87.1875   90.0000   92.8125   95.6250   98.4375  101.2500  104.0625  106.8750 ...
  109.6875  112.5000  115.3125  118.1250  120.9375  123.7500  126.5625  129.3750  132.1875  135.0000  137.8125  140.6250  143.4375 ...
  146.2500  149.0625  151.8750  154.6875  157.5000  160.3125  163.1250  165.9375  168.7500  171.5625  174.3750  177.1875  180.0000 ...
  182.8125  185.6250  188.4375  191.2500  194.0625  196.8750  199.6875  202.5000  205.3125  208.1250  210.9375  213.7500  216.5625 ...
  219.3750  222.1875  225.0000  227.8125  230.6250  233.4375  236.2500  239.0625  241.8750  244.6875  247.5000  250.3125  253.1250 ...
  255.9375  258.7500  261.5625  264.3750  267.1875  270.0000  272.8125  275.6250  278.4375  281.2500  284.0625  286.8750  289.6875 ...
  292.5000  295.3125  298.1250  300.9375  303.7500  306.5625  309.3750  312.1875  315.0000  317.8125  320.6250  323.4375  326.2500 ...
  329.0625  331.8750  334.6875  337.5000  340.3125  343.1250  345.9375  348.7500  351.5625  354.3750  357.1875  ]';
    lonb=[-1.4062    1.4062    4.2188    7.0312    9.8438   12.6562   15.4687   18.2812   21.0938   23.9062   26.7188   29.5312   32.3438 ...
   35.1562   37.9688   40.7812   43.5938   46.4062   49.2188   52.0312   54.8438   57.6562   60.4688   63.2812   66.0938   68.9062 ...
   71.7188   74.5312   77.3438   80.1562   82.9688   85.7812   88.5938   91.4062   94.2188   97.0312   99.8438  102.6562  105.4688 ...
  108.2812  111.0938  113.9062  116.7188  119.5312  122.3437  125.1563  127.9688  130.7812  133.5938  136.4062  139.2188  142.0312 ...
  144.8438  147.6562  150.4688  153.2812  156.0938  158.9062  161.7188  164.5312  167.3438  170.1562  172.9688  175.7812  178.5938 ...
  181.4062  184.2188  187.0312  189.8438  192.6562  195.4688  198.2812  201.0938  203.9062  206.7188  209.5312  212.3438  215.1562 ...
  217.9688  220.7812  223.5938  226.4062  229.2187  232.0312  234.8438  237.6562  240.4688  243.2812  246.0937  248.9063  251.7188 ...
  254.5312  257.3438  260.1562  262.9688  265.7812  268.5938  271.4062  274.2188  277.0312  279.8438  282.6562  285.4688  288.2812 ...
  291.0938  293.9062  296.7188  299.5312  302.3438  305.1562  307.9688  310.7812  313.5938  316.4062  319.2188  322.0312  324.8438 ...
  327.6562  330.4688  333.2812  336.0938  338.9062  341.7188  344.5312  347.3438  350.1562  352.9688  355.7812  358.5938  ]';
    lat=[-87.8638  -85.0965  -82.3129  -79.5256  -76.7369  -73.9475  -71.1578  -68.3678  -65.5776  -62.7874  -59.9970  -57.2066  -54.4162 ...
  -51.6257  -48.8352  -46.0447  -43.2542  -40.4636  -37.6731  -34.8825  -32.0919  -29.3014  -26.5108  -23.7202  -20.9296  -18.1390 ...
  -15.3484  -12.5578   -9.7671   -6.9765   -4.1859   -1.3953    1.3953    4.1859    6.9765    9.7671   12.5578   15.3484   18.1390 ...
   20.9296   23.7202   26.5108   29.3014   32.0919   34.8825   37.6731   40.4636   43.2542   46.0447   48.8352   51.6257   54.4162 ...
   57.2066   59.9970   62.7874   65.5776   68.3678   71.1578   73.9475   76.7369   79.5256   82.3129   85.0965   87.8638   ]';
    latb=[-90.0000  -86.5777  -83.7570  -80.9550  -78.1583  -75.3639  -72.5707  -69.7781  -66.9860  -64.1942  -61.4026  -58.6111  -55.8198 ...
  -53.0285  -50.2373  -47.4462  -44.6551  -41.8640  -39.0730  -36.2820  -33.4910  -30.7000  -27.9091  -25.1181  -22.3272  -19.5363 ...
  -16.7454  -13.9545  -11.1636   -8.3727   -5.5818   -2.7909   -0.0000    2.7909    5.5818    8.3727   11.1636   13.9545   16.7454 ...
   19.5363   22.3272   25.1181   27.9091   30.7000   33.4910   36.2820   39.0730   41.8640   44.6551   47.4462   50.2373   53.0285 ...
   55.8198   58.6111   61.4026   64.1942   66.9860   69.7781   72.5707   75.3639   78.1583   80.9550   83.7570   86.5777   90.0000 ]';
    pfull=[0.0071    0.0324    0.0738    0.1496    0.2782    0.4832    0.7946    1.2495    1.8924    2.7768    3.9656    5.5303    7.5544 ...
   10.1337   13.3701   17.3855   22.3174   28.3123   35.5326   44.1505   54.3442   66.3102   80.2926   96.5594  115.3391  136.8651 ...
  161.5270  189.6748  221.5448  257.5338  298.0871  343.5376  394.3759  451.0932  514.0716  583.8445  661.2804  747.0481  841.4656 ...
  945.1800  ]';
    phalf=[0    0.0192    0.0478    0.1032    0.2013    0.3628    0.6144    0.9895    1.5290    2.2808    3.3042    4.6660    6.4422 ...
    8.7241   11.6120   15.2090   19.6570   25.0880   31.6640   39.5470   48.9200   59.9553   72.8750   87.9460  105.4370  125.5330 ...
  148.5190  174.8940  204.8500  238.6700  276.8700  319.8200  367.8140  421.5480  481.2980  547.5570  620.9000  702.5000  792.5000 ...
  891.4000 1000.0000  ]';
end

t_length=size(tauin,3);
if(~exist('days','var'))
    dayint = 360/t_length;
    days   = dayint/2:dayint:360;
else
    if(length(days)~=t_length)
        'length of days must be equal length of input tau'
        return
    end
end

tdamp = zeros(length(lon),length(lat),length(pfull),t_length);

%% compute Held-Suarez tau

ka=1/tau_strat;
ks=1/4;
p0=1e3;
sigma=pfull/p0;
sigmab=0.7;

kv = (ks-ka)*max(0,(sigma-sigmab)/(1-sigmab));

tauhs=zeros(length(lat),length(pfull));
for k=1:length(pfull)
    kT = ka + kv(k)*cos(lat*pi/180).^4;
    tauhs(:,k) = 1./kT;
end

%% interpolate inputs onto code grid
tau  = zeros(length(lat),length(pfull),t_length);  
p_hs = zeros(length(lat),t_length);
p_bd = zeros(length(lat),t_length);  
for d=1:t_length
    if(length(latin)>1)
        tau_tmp = zeros(length(lat),length(pin));
        for k=1:length(pin);
            tau_tmp(:,k) = interp1(latin,tauin(:,k,d),lat);
        end
        K = find(pfull >= min(pin));
        for j=1:length(lat)
            tau(j,K,d) = interp1(pin,tau_tmp(j,:),pfull(K));
        end
        if(size(p_hsin,2)>1)
            p_hs(:,d) = interp1(latin,p_hsin(:,d),lat);
            p_bd(:,d) = interp1(latin,p_bdin(:,d),lat);
        else
            p_hs(:,d) = interp1(latin,p_hsin,lat);
            p_bd(:,d) = interp1(latin,p_bdin,lat);
        end
        % what to do with layers above ECMWF?
        if(min(K)>1)
            for kk=1:min(K)
                %tau(:,kk,d)=tau(:,min(K),d) + 10*(pfull(kk)-pfull(min(K)))/(pfull(1)-pfull(min(K)))*ones(length(lat),1);
                tau(:,kk,d)=tau(:,min(K),d);
            end
        end
        clear tau_tmp
    else
        for j=1:length(lat)
            tau(j,:,d) = interp1(pin,tauin(:,d),pfull);
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
            tdamp(i,j,J,d) = beta'.*tau(j,J,d) + (1-beta').*tauhs(j,J);
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
