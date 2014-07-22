function [Teq,lat,pfull]=heldsuarez(directory,T_strat,epsS,epsN,lat,pfull)


    
%% read fms dimensions and interpolate Teq onto them
if(length(directory)>0)
    filename=[directory,'/atmos_average_pfull.nc'];
    ncid = netcdf.open(filename,'NC_NOWRITE');
    varid = netcdf.inqVarID(ncid,'lat');
    lat=netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'pfull');
    pfull=netcdf.getVar(ncid,varid);
end
if(~exist('lat','var'))
    latb=linspace(-90,90,length(latin)+1)';lat=0.5*(3*latb(1)-latb(2))+cumsum(diff(latb));
%     lat=[-87.8638  -85.0965  -82.3129  -79.5256  -76.7369  -73.9475  -71.1578  -68.3678  -65.5776  -62.7874  -59.9970  -57.2066  -54.4162 ...
%   -51.6257  -48.8352  -46.0447  -43.2542  -40.4636  -37.6731  -34.8825  -32.0919  -29.3014  -26.5108  -23.7202  -20.9296  -18.1390 ...
%   -15.3484  -12.5578   -9.7671   -6.9765   -4.1859   -1.3953    1.3953    4.1859    6.9765    9.7671   12.5578   15.3484   18.1390 ...
%    20.9296   23.7202   26.5108   29.3014   32.0919   34.8825   37.6731   40.4636   43.2542   46.0447   48.8352   51.6257   54.4162 ...
%    57.2066   59.9970   62.7874   65.5776   68.3678   71.1578   73.9475   76.7369   79.5256   82.3129   85.0965   87.8638   ]';
end
if(~exist('pfull','var'))
    pfull=logspace(-3,3,40)';
%     pfull=[0.0071    0.0324    0.0738    0.1496    0.2782    0.4832    0.7946    1.2495    1.8924    2.7768    3.9656    5.5303    7.5544 ...
%    10.1337   13.3701   17.3855   22.3174   28.3123   35.5326   44.1505   54.3442   66.3102   80.2926   96.5594  115.3391  136.8651 ...
%   161.5270  189.6748  221.5448  257.5338  298.0871  343.5376  394.3759  451.0932  514.0716  583.8445  661.2804  747.0481  841.4656 ...
%   945.1800  ]';
end

Teq = zeros(length(lat),length(pfull));

%% compute Held-Suarez Teq

T0=315;
delT=60;
delv=10;
p0=1e3;
kap  = 2./7.;
sigma=pfull/p0;
Js = find(lat<0);
Jn = find(lat>=0);

te=zeros(length(lat),length(pfull));
for k=1:length(pfull)
    Teq(Js,k)=max(T_strat,(T0 - delT*(sin(lat(Js)*pi/180).^2) - epsS*sin(lat(Js)*pi/180)- delv*log(sigma(k)).*(cos(lat(Js)*pi/180).^2)).*(sigma(k)^kap));
    Teq(Jn,k)=max(T_strat,(T0 - delT*(sin(lat(Jn)*pi/180).^2) - epsN*sin(lat(Jn)*pi/180)- delv*log(sigma(k)).*(cos(lat(Jn)*pi/180).^2)).*(sigma(k)^kap));
end


end