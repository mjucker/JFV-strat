function TePK = polvanikushner(lat,pfull,T_strat,epsS,epsN,gamma)
    
    pv_phi0=50.;
    pv_dphi=-10;
    pv_gamma=gamma*1e-3;
    pref=1e3;
    p_tropopause=0.1;
    
    rdgas = 287.04;
    grav  = 9.80;
    
    

    Te0=heldsuarez([],T_strat,epsS,epsN,lat,pfull);
    
    sat_p = [1.0000000000000,0.2233611050922,0.0540329501078,0.0085666783593,0.0010945601338, 0.0006606353133,0.0000390468337,0.0000000000001];
    sat_t = [288.150, 216.650, 216.650, 228.650, 270.650, 270.650, 214.650,186.946];
    sat_g = [-6.5e-3,0.0e-3, 1.0e-3, 2.8e-3,0.0e-3,-2.8e-3, -2.0e-3, 0.0e-3];
    t_tropopause=T_strat;
    pif = 1;
    
    lat_wgt =  0.5 * ( 1. - tanh((lat(:)/pif-pv_phi0)/pv_dphi) ) * ones(size(pfull(:)'));
    
    p_norm = repmat(pfull(:)',[length(lat),1])/pref;
    
    t_sat=zeros(size(p_norm));
    for l=1:length(sat_p)-1
        stratlev = p_norm < p_tropopause & p_norm > sat_p(l+1) & p_norm <= sat_p(l);
        t_sat( stratlev ) = sat_t(l)*(p_norm(stratlev)/sat_p(l)).^(-rdgas*sat_g(l)/grav);
    end
    t_pv = zeros(size(p_norm));
    stratlev = p_norm < p_tropopause;
    t_pv(stratlev) = t_tropopause*(p_norm(stratlev)/p_tropopause).^(rdgas*pv_gamma/grav);
    
    TePK=Te0;
    TePK(stratlev) = (1-lat_wgt(stratlev)).*t_sat(stratlev) + lat_wgt(stratlev).*t_pv(stratlev);
    
end
