# this compares the expressions that go into the Hamiltonian term by term.

assume(r>0);

tortoise := 1;
copysign := proc(a,b)
  if(b >= 0) then
    return +abs(a);
   else
    return -abs(a);
   end if
end proc;

log2 := x -> log[2](x);
log1p := x -> log(1+x);
invlog_2e := 1/log2(exp(1));
u := 1/r;
u2 := 1/r^2;
u3 := 1/r^3;
u4 := 1/r^4;
u5 := 1/r^5;
logu := log(u);

# some things I removed and did not just change
m1PlusetaKK := 1/invm1PlusetaKK;
nu := 1/2 * simplify(log(deltaT*rho2/Lambda), symbolic);
MU := 1/2 * log(rho2);  

simp := proc(x)
  return simplify(subs(B=sqrt(deltaT), x), symbolic)
end proc;

# original expressions
  orig_r2 := x_data[0]*x_data[0] + x_data[1]*x_data[1] + x_data[2]*x_data[2];
  orig_r := sqrt(r2);
  orig_nx := x_data[0] / r;
  orig_ny := x_data[1] / r;
  orig_nz := x_data[2] / r;   
     
  orig_sKerr_x := sigmaKerr_data[0];
  orig_sKerr_y := sigmaKerr_data[1];
  orig_sKerr_z := sigmaKerr_data[2];

  orig_sStar_x := sigmaStar_data[0];
  orig_sStar_y := sigmaStar_data[1];
  orig_sStar_z := sigmaStar_data[2];
     
  orig_a2 := sKerr_x*sKerr_x + sKerr_y*sKerr_y + sKerr_z*sKerr_z;
  orig_a := sqrt( a2 );

    orig_e3_x := sKerr_x / a;
    orig_e3_y := sKerr_y / a;
    orig_e3_z := sKerr_z / a;
    
  orig_costheta := e3_x*nx + e3_y*ny + e3_z*nz; 
    
  orig_xi2 :=1. - costheta*costheta;

  orig_xi_x := -e3_z*ny + e3_y*nz;
  orig_xi_y :=  e3_z*nx - e3_x*nz;
  orig_xi_z := -e3_y*nx + e3_x*ny;

  orig_vx := -nz*xi_y + ny*xi_z;
  orig_vy :=  nz*xi_x - nx*xi_z;
  orig_vz := -ny*xi_x + nx*xi_y;

  orig_w2 := r2 + a2;
  orig_rho2 := r2 + a2*costheta*costheta;

  orig_u := 1./r;
  orig_u2 := u*u;
  orig_u3 := u2*u;
  orig_u4 := u2*u2;
  orig_u5 := u4*u;

  orig_m1PlusetaKK := -1. + eta * coeffs_KK;
  # Eq. 5.75 of BB1 */
  orig_bulk := 1./(m1PlusetaKK*m1PlusetaKK) + (2.*u)/m1PlusetaKK + a2*u2;
  # Eq. 5.73 of BB1 */
  orig_logTerms := 1. + eta*coeffs_k0 + eta*log((1. + coeffs_k1*u + coeffs_k2*u2 + coeffs_k3*u3 + coeffs_k4*u4
                                              + coeffs_k5*u5 + coeffs_k5l*u5*log(u)));
  # Eq. 5.73 of BB1 */
  orig_deltaU := bulk*logTerms;
  # Eq. 5.71 of BB1 */
  orig_deltaT := r2*deltaU;
  # ddeltaU/du */
  orig_deltaU_u := 2.*(1./m1PlusetaKK + a2*u)*logTerms + 
	  bulk * (eta*(coeffs_k1 + u*(2.*coeffs_k2 + u*(3.*coeffs_k3 + u*(4.*coeffs_k4 + 5.*(coeffs_k5+coeffs_k5l*log(u))*u)))))
          / (1. + coeffs_k1*u + coeffs_k2*u2 + coeffs_k3*u3 + coeffs_k4*u4 + (coeffs_k5+coeffs_k5l*log(u))*u5);
  # ddeltaT/dr */
  orig_deltaT_r := 2.*r*deltaU - deltaU_u;
  # Eq. 5.39 of BB1 */
  orig_Lambda := (w2*w2 - a2*deltaT*xi2);
  # Eq. 5.83 of BB1, inverse */
  orig_D := 1. + log(1. + 6.*eta*u2 + 2.*(26. - 3.*eta)*eta*u3);
  # Eq. 5.38 of BB1 */
  orig_deltaR := deltaT*D;
  # See Hns below, Eq. 4.34 of Damour et al. PRD 62, 084011 (2000) */
  orig_qq := 2.*eta*(4. - 3.*eta);
  # See Hns below. In Sec. II D of BB2 b3 and bb3 coeffs are chosen to be zero. */
  orig_ww :=2.*a*r + coeffs_b3*eta*a2*a*u + coeffs_bb3*eta*a*u;

  # We need to transform the momentum to get the tortoise co-ord */
    orig_csi := sqrt( (deltaT * deltaR) )/ w2; # Eq. 28 of Pan et al. PRD 81, 084041 (2010) */


	  orig_prT := p_data[0]*nx + p_data[1]*ny + p_data[2]*nz;
      # p_data is BL momentum vector; tmpP is tortoise momentum vector */ 
      orig_tmpP[0] := p_data[0] - nx * prT * (csi - 1.)/csi;
      orig_tmpP[1] := p_data[1] - ny * prT * (csi - 1.)/csi;
      orig_tmpP[2] := p_data[2] - nz * prT * (csi - 1.)/csi;
  
  orig_pxir := (tmpP[0]*xi_x + tmpP[1]*xi_y + tmpP[2]*xi_z) * r;
  orig_pvr := (tmpP[0]*vx + tmpP[1]*vy + tmpP[2]*vz) * r;
  orig_pn := tmpP[0]*nx + tmpP[1]*ny + tmpP[2]*nz;
          
  orig_pr := pn;
  orig_pf := pxir;
  orig_ptheta2 := pvr * pvr / xi2;

  # Eqs. 5.36 - 5.46 of BB1 */
  # Note that the tortoise prT appears only in the quartic term, explained in Eqs. 14 and 15 of Tarrachini et al. */
  orig_Hns := sqrt(1. + prT*prT*prT*prT*qq*u2 + ptheta2/rho2 + pf*pf*rho2/(Lambda*xi2) + pr*pr*deltaR/rho2)
      / sqrt(Lambda/(rho2*deltaT)) + pf*ww/Lambda;
  
  # Eqs. 5.30 - 5.33 of BB1 */
  orig_B := sqrt(deltaT);
  orig_w := ww/Lambda;
  orig_nu := 0.5 * log(deltaT*rho2/Lambda);
  orig_MU := 0.5 * log(rho2);  
  # dLambda/dr */
  orig_Lambda_r := 4.*r*w2 - a2*deltaT_r*xi2;
     
  orig_ww_r :=2.*a - (a2*a*coeffs_b3*eta)*u2 - coeffs_bb3*eta*a*u2;
  # Eqs. 5.47a - 5.47d of BB1 */
  orig_BR := (-2.*deltaT + sqrt(deltaR)*deltaT_r)/(2.*sqrt(deltaR*deltaT));
  orig_wr := (-Lambda_r*ww + Lambda*ww_r)/(Lambda*Lambda);
  orig_nur := (r/rho2 + (w2 * (-4.*r*deltaT + w2*deltaT_r) ) / (2.*deltaT*Lambda) );
  orig_mur := (r/rho2 - 1./sqrt(deltaR));
  # Eqs. 5.47f - 5.47h of BB1 */
  orig_wcos := -2.*a2*costheta*deltaT*ww/(Lambda*Lambda);  
  orig_nucos := a2*costheta*w2*(w2-deltaT)/(rho2*Lambda);  
  orig_mucos := a2*costheta/rho2;
  # Eq. 5.52 of BB1, (YP) simplified */
  orig_Q := 1. + pvr*pvr/(rho2*xi2) + deltaT*rho2/Lambda*pxir*pxir/(B*B*xi2) + pn*pn*deltaR/rho2;
  orig_pn2 := pr * pr * deltaR / rho2;
  orig_pp := Q - 1.;

  # Eq. 5.68 of BB1, (YP) simplified for orig_aa :=bb=0. */
  orig_deltaSigmaStar_x :=eta*(-8.*sKerr_x - 36.*pn2*r*sKerr_x + 3.*pp*r*sKerr_x + 14.*sStar_x - 30.*pn2*r*sStar_x + 4.*pp*r*sStar_x)/(12.*r);

  orig_deltaSigmaStar_y :=eta*(-8.*sKerr_y - 36.*pn2*r*sKerr_y + 3.*pp*r*sKerr_y + 14.*sStar_y - 30.*pn2*r*sStar_y + 4.*pp*r*sStar_y)/(12.*r);

  orig_deltaSigmaStar_z :=eta*(-8.*sKerr_z - 36.*pn2*r*sKerr_z + 3.*pp*r*sKerr_z + 14.*sStar_z - 30.*pn2*r*sStar_z + 4.*pp*r*sStar_z)/(12.*r);


  # Now compute the additional 3.5PN terms. */
  # The following gauge parameters correspond to those given by 
  # * Eqs. (69) and (70) of BB2 (aaa _ a0, bbb _ b0).
  # * In SEOBNRv1 model, we chose to set all of them to zero,
  # * described between Eqs. (3) and (4).
  # */
  # Eq. 52 of BB2, (YP) simplified for zero gauge parameters */    
  orig_sMultiplier1 := -(2.*eta*(-353. + 27.*eta) + 2.*(103.*eta - 60.*eta*eta)*pp*r 
               + 120.*(-3.*eta*eta)*pn2*pn2*r*r + (eta*(23. + 3.*eta))*pp*pp*r*r 
               + 6.*pn2*r*(- 47.*eta + 54.*eta*eta + (- 16.*eta + 21.*eta*eta)*pp*r))
               / (72.*r*r);                        
  # Eq. 52 of BB2, (YP) simplified for zero gauge parameters */       
  orig_sMultiplier2 := (-16.*(7.*eta*(8. + 3.*eta)) + 4.*(- 109.*eta + 51.*eta*eta)*pp*r 
               + 810.*eta*eta*pn2*pn2*r*r - 45.*eta*pp*pp*r*r 
               - 6.*pn2*r*(16.*eta + 147.*eta*eta + (- 6.*eta + 39.*eta*eta)*pp*r))
               / (144.*r*r);
  # Eq. 52 of BB2 */                     
  orig_deltaSigmaStar_x := deltaSigmaStar_x +sMultiplier1*sigmaStar_data[0] + sMultiplier2*sigmaKerr_data[0];
  orig_deltaSigmaStar_y := deltaSigmaStar_y +sMultiplier1*sigmaStar_data[1] + sMultiplier2*sigmaKerr_data[1];
  orig_deltaSigmaStar_z := deltaSigmaStar_z +sMultiplier1*sigmaStar_data[2] + sMultiplier2*sigmaKerr_data[2];

  # And now the (calibrated) 4.5PN term */
  orig_deltaSigmaStar_x := deltaSigmaStar_x +coeffs_d1 * eta * sigmaStar_data[0] / (r*r*r);
  orig_deltaSigmaStar_y := deltaSigmaStar_y +coeffs_d1 * eta * sigmaStar_data[1] / (r*r*r);
  orig_deltaSigmaStar_z := deltaSigmaStar_z +coeffs_d1 * eta * sigmaStar_data[2] / (r*r*r);
  orig_deltaSigmaStar_x := deltaSigmaStar_x +coeffs_d1v2 * eta * sigmaKerr_data[0] / (r*r*r);
  orig_deltaSigmaStar_y := deltaSigmaStar_y +coeffs_d1v2 * eta * sigmaKerr_data[1] / (r*r*r);
  orig_deltaSigmaStar_z := deltaSigmaStar_z +coeffs_d1v2 * eta * sigmaKerr_data[2] / (r*r*r);


  orig_sx := sStar_x + deltaSigmaStar_x;
  orig_sy := sStar_y + deltaSigmaStar_y;
  orig_sz := sStar_z + deltaSigmaStar_z;     
     
     
  orig_sxi := sx*xi_x + sy*xi_y + sz*xi_z;
  orig_sv := sx*vx + sy*vy + sz*vz;
  orig_sn := sx*nx + sy*ny + sz*nz; 
     
  orig_s3 := sx*e3_x + sy*e3_y + sz*e3_z;  
  # Eq. 3.45 of BB1, second term */        
  orig_Hwr := (exp(-3.*MU - nu)*sqrt(deltaR)*(exp(2.*(MU + nu))*pxir*pxir*sv - B*exp(MU + nu)*pvr*pxir*sxi + 
        B*B*xi2*(exp(2.*MU)*(sqrt(Q) + Q)*sv + pn*pvr*sn*sqrt(deltaR) - pn*pn*sv*deltaR)))/(2.*B*(1. + sqrt(Q))*sqrt(Q)*xi2);
  # Eq. 3.45 of BB1, third term */     
  orig_Hwcos := (exp(-3.*MU - nu)*(sn*(-(exp(2.*(MU + nu))*pxir*pxir) + B*B*(pvr*pvr - exp(2.*MU)*(sqrt(Q) + Q)*xi2)) - 
        B*pn*(B*pvr*sv - exp(MU + nu)*pxir*sxi)*sqrt(deltaR)))/(2.*B*(1. + sqrt(Q))*sqrt(Q));
  # Eq. 3.44 of BB1, leading term */     
  orig_HSOL := (exp(-MU + 2.*nu)*(-B + exp(MU + nu))*pxir*s3)/(B*B*sqrt(Q)*xi2);
  # Eq. 3.44 of BB1, next-to-leading term */
  orig_HSONL := (exp(-2.*MU + nu)*(-(B*exp(MU + nu)*nucos*pxir*(1. + 2.*sqrt(Q))*sn*xi2) + 
        (-(BR*exp(MU + nu)*pxir*(1. + sqrt(Q))*sv) + B*(exp(MU + nu)*nur*pxir*(1. + 2.*sqrt(Q))*sv + B*mur*pvr*sxi + 
        B*sxi*(-(mucos*pn*xi2) + sqrt(Q)*(mur*pvr - nur*pvr + (-mucos + nucos)*pn*xi2))))*sqrt(deltaR)))/(B*B*(sqrt(Q) + Q)*xi2);   
  # Eq. 3.43 and 3.45 of BB1 */
  orig_Hs := w*s3 + Hwr*wr + Hwcos*wcos + HSOL + HSONL;
  # Eq. 5.70 of BB1, last term */   
  orig_Hss := -0.5*u3 * (sx*sx + sy*sy + sz*sz - 3.*sn*sn);
  # Eq. 5.70 of BB1 */
  orig_H := Hns + Hs + Hss;

  # Add the additional calibrated term */
  orig_H := H + coeffs_dheffSS * eta * (sKerr_x*sStar_x + sKerr_y*sStar_y + sKerr_z*sStar_z) / (r*r*r*r);
  # One more calibrated term proportional to S1^2+S2^2. Note that we use symmetric expressions of m1,m2 and S1,S2 */
  orig_H := H + coeffs_dheffSSv2 * eta / (r*r*r*r)
                         * (s1Vec_data[0]*s1Vec_data[0] + s1Vec_data[1]*s1Vec_data[1] + s1Vec_data[2]*s1Vec_data[2]
                           +s2Vec_data[0]*s2Vec_data[0] + s2Vec_data[1]*s2Vec_data[1] + s2Vec_data[2]*s2Vec_data[2]);
  # Real Hamiltonian given by Eq. 2, ignoring the constant -1. */
  orig_Hreal := sqrt(1. + 2.*eta *((H) - 1.));

# current expression
  current_r2 := x_data[0]*x_data[0] + x_data[1]*x_data[1] + x_data[2]*x_data[2];
  current_r := sqrt(r2);
  current_u  := 1./r;
  current_u2 := u*u;
  current_u3 := u2*u;
  current_u4 := u2*u2;
  current_u5 := u4*u;

  current_nx := x_data[0] *u;
  current_ny := x_data[1] *u;
  current_nz := x_data[2] *u;   
     
  current_sKerr_x := sigmaKerr_data[0];
  current_sKerr_y := sigmaKerr_data[1];
  current_sKerr_z := sigmaKerr_data[2];

  current_sStar_x := sigmaStar_data[0];
  current_sStar_y := sigmaStar_data[1];
  current_sStar_z := sigmaStar_data[2];
     
  current_a2 := sKerr_x*sKerr_x + sKerr_y*sKerr_y + sKerr_z*sKerr_z;
  current_a := sqrt( a2 );

  # RH: the simplest way to get rid of the "if" is to set
  # sKerr_y += copysign(sKerr_y, 1e-15)
  # a += 1e-15
  # such that for current_sKerr_x :==sKerr_y==sKerr_z==0. the devision takes care of
  # getting the correct value and 1e-15 is small enough that we don't are
  # about it otherwise
    inva := 1./a;
    current_e3_x := sKerr_x * inva;
    current_e3_y := sKerr_y * inva;
    current_e3_z := sKerr_z * inva;

  current_costheta := e3_x*nx + e3_y*ny + e3_z*nz; 
    
  current_xi2 :=1. - costheta*costheta;

  current_xi_x := -e3_z*ny + e3_y*nz;
  current_xi_y :=  e3_z*nx - e3_x*nz;
  current_xi_z := -e3_y*nx + e3_x*ny;

  current_vx := -nz*xi_y + ny*xi_z;
  current_vy :=  nz*xi_x - nx*xi_z;
  current_vz := -ny*xi_x + nx*xi_y;

  current_w2 := r2 + a2;
  current_rho2 := r2 + a2*costheta*costheta;

  invm1PlusetaKK := 1./(-1. + eta * coeffs_KK);
  # Eq. 5.75 of BB1 */
  current_bulk := invm1PlusetaKK*(invm1PlusetaKK + (2.*u)) + a2*u2;
  # Eq. 5.73 of BB1 */
  # use ln(u) = log_2(u)/log_2(e) and the fact that log2 is faster than ln
  # this relies on the compiler evaluating the expression at compile time.
  # which apparently not all do so in stead of 1./log2(exp(1.)) I use the
  # result returned by Maple.
  current_logu := log2(u)*invlog_2e;
  current_logTerms := 1. + eta*coeffs_k0 + eta*log1p(coeffs_k1*u + coeffs_k2*u2 + coeffs_k3*u3 + coeffs_k4*u4
                                             + coeffs_k5*u5 + coeffs_k5l*u5*logu);
  # Eq. 5.73 of BB1 */
  current_deltaU := bulk*logTerms;
  # Eq. 5.71 of BB1 */
  current_deltaT := r2*deltaU;
  # ddeltaU/du */
  current_deltaU_u := 2.*(invm1PlusetaKK + a2*u)*logTerms +
	  bulk * (eta*(coeffs_k1 + u*(2.*coeffs_k2 + u*(3.*coeffs_k3 + u*(4.*coeffs_k4 + 5.*(coeffs_k5+coeffs_k5l*logu)*u)))))
          / (1. + coeffs_k1*u + coeffs_k2*u2 + coeffs_k3*u3 + coeffs_k4*u4 + (coeffs_k5+coeffs_k5l*logu)*u5);
  # ddeltaT/dr */
  current_deltaT_r := 2.*r*deltaU - deltaU_u;
  # Eq. 5.39 of BB1 */
  current_Lambda := w2*w2 - a2*deltaT*xi2;
  # RH: this is horrible, but faster than 3 divisions
  invrho2xi2Lambda := 1./(rho2*xi2*Lambda);
  invrho2 := xi2 * (Lambda*invrho2xi2Lambda);
  invxi2 := rho2 * (Lambda*invrho2xi2Lambda);
  invLambda := xi2*rho2*invrho2xi2Lambda;
  # Eq. 5.83 of BB1, inverse */
  current_D := 1. + log2(1. + 6.*eta*u2 + 2.*(26. - 3.*eta)*eta*u3)*invlog_2e;
  # Eq. 5.38 of BB1 */
  current_deltaR := deltaT*D;
  # See Hns below, Eq. 4.34 of Damour et al. PRD 62, 084011 (2000) */
  current_qq := 2.*eta*(4. - 3.*eta);
  # See Hns below. In Sec. II D of BB2 b3 and bb3 coeffs are chosen to be zero. */
  current_ww :=2.*a*r + coeffs_b3*eta*a2*a*u + coeffs_bb3*eta*a*u;

  # We need to transform the momentum to get the tortoise co-ord */
  # Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
  # RH: this assumes that tortoise can be 0 or 1 or 2.
  current_csi := sqrt( deltaT * deltaR )/ w2;
  # non-unity only for current_tortoise :==1
  csi1 := 1.0 + (1.-(1.-tortoise)) * (csi - 1.0); 
  # non-unity only for current_tortoise :==2
  csi2 := 1.0 + (0.5-copysign(0.5, 1.5-tortoise)) * (csi - 1.0); 

  current_prT := (p_data[0]*nx + p_data[1]*ny + p_data[2]*nz)*csi2;
  # p_data is BL momentum vector; tmpP is tortoise momentum vector */
  current_tmpP[0] := p_data[0] - nx * prT * ((csi1 - 1.)/csi1);
  current_tmpP[1] := p_data[1] - ny * prT * ((csi1 - 1.)/csi1);
  current_tmpP[2] := p_data[2] - nz * prT * ((csi1 - 1.)/csi1);
  
  current_pxir := (tmpP[0]*xi_x + tmpP[1]*xi_y + tmpP[2]*xi_z) * r;
  current_pvr := (tmpP[0]*vx + tmpP[1]*vy + tmpP[2]*vz) * r;
  current_pn := tmpP[0]*nx + tmpP[1]*ny + tmpP[2]*nz;
          
  current_pr := pn;
  current_pf := pxir;
  current_ptheta2 := pvr * pvr *invxi2;

  # Eqs. 5.36 - 5.46 of BB1 */
  # Note that the tortoise prT appears only in the quartic term, explained in Eqs. 14 and 15 of Tarrachini et al. */
  current_Hns := sqrt((1. + prT*prT*prT*prT*qq*u2 + ptheta2*invrho2 + pf*pf*rho2*invLambda*invxi2 + pr*pr*deltaR*invrho2)
             * (rho2*deltaT) * invLambda) + pf*ww*invLambda;
  
  # Eqs. 5.30 - 5.33 of BB1 */
  current_B := sqrt(deltaT);
  # RH: this is horrible but faster than 3 divisions
  sqrtdeltaT := sqrt(deltaT);
  sqrtdeltaR := sqrt(deltaR);
  invdeltaTsqrtdeltaTsqrtdeltaR := 1./(sqrtdeltaT*deltaT*sqrtdeltaR);
  invdeltaT := sqrtdeltaT*(sqrtdeltaR*invdeltaTsqrtdeltaTsqrtdeltaR);
  invsqrtdeltaT := deltaT*(sqrtdeltaR*invdeltaTsqrtdeltaTsqrtdeltaR);
  invsqrtdeltaR := deltaT*sqrtdeltaT*invdeltaTsqrtdeltaTsqrtdeltaR;
  current_w := ww*invLambda;
  expnu := simplify(sqrt(deltaT*rho2*invLambda), symbolic);
  expMU := simplify(sqrt(rho2), symbolic);
  # RH: this is horrible but faster than 2 divisions
  invexpnuexpMU := 1./(expnu*expMU);
  invexpnu := expMU*invexpnuexpMU;
  invexpMU := expnu*invexpnuexpMU;
  # dLambda/dr */
  current_Lambda_r := 4.*r*w2 - a2*deltaT_r*xi2;
     
  current_ww_r :=2.*a - (a2*a*coeffs_b3*eta)*u2 - coeffs_bb3*eta*a*u2;
  # Eqs. 5.47a - 5.47d of BB1 */
  current_BR := (-deltaT*invsqrtdeltaR + deltaT_r*0.5)*invsqrtdeltaT;
  current_wr := (-Lambda_r*ww + Lambda*ww_r)*(invLambda*invLambda);
  current_nur := (r*invrho2 + (w2 * (-4.*r*deltaT + w2*deltaT_r) ) * 0.5*invdeltaT*invLambda );
  current_mur := (r*invrho2 - invsqrtdeltaR);
  # Eqs. 5.47f - 5.47h of BB1 */
  current_wcos := -2.*(a2*costheta)*deltaT*ww*(invLambda*invLambda);
  current_nucos := (a2*costheta)*w2*(w2-deltaT)*(invrho2*invLambda);
  current_mucos := (a2*costheta)*invrho2;
  # Eq. 5.52 of BB1, (YP) simplified */
  current_Q := 1. + pvr*pvr*invrho2*invxi2 + deltaT*rho2*invLambda*pxir*pxir*invdeltaT*invxi2 + pn*pn*deltaR*invrho2;
  current_pn2 := pr * pr * deltaR * invrho2;
  current_pp := Q - 1.;

  # Eq. 5.68 of BB1, (YP) simplified for current_aa :=bb=0. */
  current_deltaSigmaStar_x :=eta*((-8.*sKerr_x - 36.*pn2*r*sKerr_x) + (3.*pp*r*sKerr_x + 14.*sStar_x) + (- 30.*pn2*r*sStar_x + 4.*pp*r*sStar_x))*(1./12.)*u;

  current_deltaSigmaStar_y :=eta*((-8.*sKerr_y - 36.*pn2*r*sKerr_y) + (3.*pp*r*sKerr_y + 14.*sStar_y) + (- 30.*pn2*r*sStar_y + 4.*pp*r*sStar_y))*(1./12.)*u;

  current_deltaSigmaStar_z :=eta*((-8.*sKerr_z - 36.*pn2*r*sKerr_z) + (3.*pp*r*sKerr_z + 14.*sStar_z) + (- 30.*pn2*r*sStar_z + 4.*pp*r*sStar_z))*(1./12.)*u;


  # Now compute the additional 3.5PN terms. */
  # The following gauge parameters correspond to those given by 
  # * Eqs. (69) and (70) of BB2 (aaa _ a0, bbb _ b0).
  # * In SEOBNRv1 model, we chose to set all of them to zero,
  # * described between Eqs. (3) and (4).
  # */
  # Eq. 52 of BB2, (YP) simplified for zero gauge parameters */    
  current_sMultiplier1 := -(2.*eta*(-353. + 27.*eta) + 2.*(103.*eta - 60.*eta*eta)*pp*r 
               + (120.*(-3.))*(eta*eta)*(pn2*pn2)*(r*r) + (eta*(23. + 3.*eta))*(pp*pp)*(r*r )
               + 6.*pn2*r*(- 47.*eta + 54.*(eta*eta) + (- 16.*eta + 21.*(eta*eta))*pp*r))
               * (1./72.) * u2;
  # Eq. 52 of BB2, (YP) simplified for zero gauge parameters */       
  current_sMultiplier2 := (-16.*(7.*eta*(8. + 3.*eta)) + 4.*(- 109.*eta + 51.*eta*eta)*pp*r 
               + 810.*(eta*eta)*(pn2*pn2)*(r*r) - 45.*eta*(pp*pp)*(r*r)
               - 6.*pn2*r*(16.*eta + 147.*eta*eta + (- 6.*eta + 39.*(eta*eta))*pp*r))
               * (1./144.) * u2;
  # Eq. 52 of BB2 */                     
  current_deltaSigmaStar_x := deltaSigmaStar_x +sMultiplier1*sigmaStar_data[0] + sMultiplier2*sigmaKerr_data[0];
  current_deltaSigmaStar_y := deltaSigmaStar_y +sMultiplier1*sigmaStar_data[1] + sMultiplier2*sigmaKerr_data[1];
  current_deltaSigmaStar_z := deltaSigmaStar_z +sMultiplier1*sigmaStar_data[2] + sMultiplier2*sigmaKerr_data[2];

  # And now the (calibrated) 4.5PN term */
  current_deltaSigmaStar_x := deltaSigmaStar_x +coeffs_d1 * eta * sigmaStar_data[0] * u3;
  current_deltaSigmaStar_y := deltaSigmaStar_y +coeffs_d1 * eta * sigmaStar_data[1] * u3;
  current_deltaSigmaStar_z := deltaSigmaStar_z +coeffs_d1 * eta * sigmaStar_data[2] * u3;
  current_deltaSigmaStar_x := deltaSigmaStar_x +coeffs_d1v2 * eta * sigmaKerr_data[0] * u3;
  current_deltaSigmaStar_y := deltaSigmaStar_y +coeffs_d1v2 * eta * sigmaKerr_data[1] * u3;
  current_deltaSigmaStar_z := deltaSigmaStar_z +coeffs_d1v2 * eta * sigmaKerr_data[2] * u3;


  current_sx := sStar_x + deltaSigmaStar_x;
  current_sy := sStar_y + deltaSigmaStar_y;
  current_sz := sStar_z + deltaSigmaStar_z;     
     
     
  current_sxi := sx*xi_x + sy*xi_y + sz*xi_z;
  current_sv := sx*vx + sy*vy + sz*vz;
  current_sn := sx*nx + sy*ny + sz*nz; 
     
  current_s3 := sx*e3_x + sy*e3_y + sz*e3_z;  
  # Eq. 3.45 of BB1, second term */        
  sqrtQ := sqrt(Q);
  inv2B1psqrtQsqrtQ := 1./(2.*B*(1. + sqrtQ)*sqrtQ);
  current_Hwr := ((invexpMU*invexpMU*invexpMU*invexpnu)*sqrtdeltaR*((expMU*expMU*expnu*expnu)*pxir*pxir*sv - B*(expMU*expnu)*pvr*pxir*sxi +
                                                       B*B*xi2*((expMU*expMU)*(sqrtQ + Q)*sv + pn*pvr*sn*sqrtdeltaR - pn*pn*sv*deltaR)))*inv2B1psqrtQsqrtQ*invxi2;
  # Eq. 3.45 of BB1, third term */     
  current_Hwcos := ((invexpMU*invexpMU*invexpMU*invexpnu)*(sn*(-(expMU*expMU*expnu*expnu*pxir*pxir) + B*B*(pvr*pvr - (expMU*expMU)*(sqrtQ + Q)*xi2)) -
                                            B*pn*(B*pvr*sv - (expMU*expnu)*pxir*sxi)*sqrtdeltaR))*inv2B1psqrtQsqrtQ;
  # Eq. 3.44 of BB1, leading term */     
  current_HSOL := ((expnu*expnu*invexpMU)*(-B + (expMU*expnu))*pxir*s3)/(deltaT*sqrtQ)*invxi2;
  # Eq. 3.44 of BB1, next-to-leading term */
  current_HSONL := ((expnu*(invexpMU*invexpMU))*(-(B*expMU*expnu*nucos*pxir*(1. + 2.*sqrtQ)*sn*xi2) +
        (-(BR*(expMU*expnu)*pxir*(1. + sqrtQ)*sv) + B*((expMU*expnu)*nur*pxir*(1. + 2.*sqrtQ)*sv + B*mur*pvr*sxi + 
        B*sxi*(-(mucos*pn*xi2) + sqrtQ*(mur*pvr - nur*pvr + (-mucos + nucos)*pn*xi2))))*sqrtdeltaR))*invxi2/(deltaT*(sqrtQ + Q));
  # Eq. 3.43 and 3.45 of BB1 */
  current_Hs := w*s3 + Hwr*wr + Hwcos*wcos + HSOL + HSONL;
  # Eq. 5.70 of BB1, last term */   
  current_Hss := -0.5*u3 * (sx*sx + sy*sy + sz*sz - 3.*sn*sn);
  # Eq. 5.70 of BB1 */
  current_H := Hns + Hs + Hss;

  # Add the additional calibrated term */
  current_H := H+coeffs_dheffSS * eta * (sKerr_x*sStar_x + sKerr_y*sStar_y + sKerr_z*sStar_z) *u4;
  # One more calibrated term proportional to S1^2+S2^2. Note that we use symmetric exp2ressions of m1,m2 and S1,S2 */
  current_H := H+coeffs_dheffSSv2 * eta * u4
                         * (s1Vec_data[0]*s1Vec_data[0] + s1Vec_data[1]*s1Vec_data[1] + s1Vec_data[2]*s1Vec_data[2]
                           +s2Vec_data[0]*s2Vec_data[0] + s2Vec_data[1]*s2Vec_data[1] + s2Vec_data[2]*s2Vec_data[2]);
  # Real Hamiltonian given by Eq. 2, ignoring the constant -1. */
  current_Hreal := sqrt(1. + 2.*eta *(H - 1.));

# now compare the terms
# gawk '$1~/current_/{print gensub("current","diff","",$1)" := simp("$1" - "gensub("current","orig","",$1)");"}' | sort -u 
# shows all defined terms
diff_B := simp(current_B - orig_B);
diff_BR := simp(current_BR - orig_BR);
diff_D := simp(current_D - orig_D);
diff_H := simp(current_H - orig_H);
diff_HSOL := simp(current_HSOL - orig_HSOL);
diff_HSONL := simp(current_HSONL - orig_HSONL);
diff_Hns := simp(current_Hns - orig_Hns);
diff_Hreal := simp(current_Hreal - orig_Hreal);
diff_Hs := simp(current_Hs - orig_Hs);
diff_Hss := simp(current_Hss - orig_Hss);
diff_Hwcos := simp(current_Hwcos - orig_Hwcos);
diff_Hwr := simp(current_Hwr - orig_Hwr);
diff_Lambda := simp(current_Lambda - orig_Lambda);
diff_Lambda_r := simp(current_Lambda_r - orig_Lambda_r);
diff_Q := simp(current_Q - orig_Q);
diff_a := simp(current_a - orig_a);
diff_a2 := simp(current_a2 - orig_a2);
diff_bulk := simp(current_bulk - orig_bulk);
diff_costheta := simp(current_costheta - orig_costheta);
diff_csi := simp(current_csi - orig_csi);
diff_deltaR := simp(current_deltaR - orig_deltaR);
diff_deltaSigmaStar_x := simp(current_deltaSigmaStar_x - orig_deltaSigmaStar_x);
diff_deltaSigmaStar_y := simp(current_deltaSigmaStar_y - orig_deltaSigmaStar_y);
diff_deltaSigmaStar_z := simp(current_deltaSigmaStar_z - orig_deltaSigmaStar_z);
diff_deltaT := simp(current_deltaT - orig_deltaT);
diff_deltaT_r := simp(current_deltaT_r - orig_deltaT_r);
diff_deltaU := simp(current_deltaU - orig_deltaU);
diff_deltaU_u := simp(current_deltaU_u - orig_deltaU_u);
diff_e3_x := simp(current_e3_x - orig_e3_x);
diff_e3_y := simp(current_e3_y - orig_e3_y);
diff_e3_z := simp(current_e3_z - orig_e3_z);
diff_invm1PlusetaKK := simp(current_invm1PlusetaKK - orig_invm1PlusetaKK);
diff_logTerms := simp(current_logTerms - orig_logTerms);
diff_logu := simp(current_logu - orig_logu);
diff_mucos := simp(current_mucos - orig_mucos);
diff_mur := simp(current_mur - orig_mur);
diff_nucos := simp(current_nucos - orig_nucos);
diff_nur := simp(current_nur - orig_nur);
diff_nx := simp(current_nx - orig_nx);
diff_ny := simp(current_ny - orig_ny);
diff_nz := simp(current_nz - orig_nz);
diff_pf := simp(current_pf - orig_pf);
diff_pn := simp(current_pn - orig_pn);
diff_pn2 := simp(current_pn2 - orig_pn2);
diff_pp := simp(current_pp - orig_pp);
diff_pr := simp(current_pr - orig_pr);
diff_prT := simp(current_prT - orig_prT);
diff_ptheta2 := simp(current_ptheta2 - orig_ptheta2);
diff_pvr := simp(current_pvr - orig_pvr);
diff_pxir := simp(current_pxir - orig_pxir);
diff_qq := simp(current_qq - orig_qq);
diff_r := simp(current_r - orig_r);
diff_r2 := simp(current_r2 - orig_r2);
diff_rho2 := simp(current_rho2 - orig_rho2);
diff_s3 := simp(current_s3 - orig_s3);
diff_sKerr_x := simp(current_sKerr_x - orig_sKerr_x);
diff_sKerr_y := simp(current_sKerr_y - orig_sKerr_y);
diff_sKerr_z := simp(current_sKerr_z - orig_sKerr_z);
diff_sMultiplier1 := simp(current_sMultiplier1 - orig_sMultiplier1);
diff_sMultiplier2 := simp(current_sMultiplier2 - orig_sMultiplier2);
diff_sStar_x := simp(current_sStar_x - orig_sStar_x);
diff_sStar_y := simp(current_sStar_y - orig_sStar_y);
diff_sStar_z := simp(current_sStar_z - orig_sStar_z);
diff_sn := simp(current_sn - orig_sn);
diff_sv := simp(current_sv - orig_sv);
diff_sx := simp(current_sx - orig_sx);
diff_sxi := simp(current_sxi - orig_sxi);
diff_sy := simp(current_sy - orig_sy);
diff_sz := simp(current_sz - orig_sz);
diff_tmpP[0] := simp(current_tmpP[0] - orig_tmpP[0]);
diff_tmpP[1] := simp(current_tmpP[1] - orig_tmpP[1]);
diff_tmpP[2] := simp(current_tmpP[2] - orig_tmpP[2]);
diff_vx := simp(current_vx - orig_vx);
diff_vy := simp(current_vy - orig_vy);
diff_vz := simp(current_vz - orig_vz);
diff_w := simp(current_w - orig_w);
diff_w2 := simp(current_w2 - orig_w2);
diff_wcos := simp(current_wcos - orig_wcos);
diff_wr := simp(current_wr - orig_wr);
diff_ww := simp(current_ww - orig_ww);
diff_ww_r := simp(current_ww_r - orig_ww_r);
diff_xi2 := simp(current_xi2 - orig_xi2);
diff_xi_x := simp(current_xi_x - orig_xi_x);
diff_xi_y := simp(current_xi_y - orig_xi_y);
diff_xi_z := simp(current_xi_z - orig_xi_z);
