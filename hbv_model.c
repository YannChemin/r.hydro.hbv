#include<math.h>
#include<sys/param.h> /*MIN() && MAX()*/

int hbv_model(int n,int t,int b,float **p,float **po,float pcalt,float *dep,float **ta,float **tm,float tcalt,float *det,float tt,float tti,float **s,float *ffo,float foscf,float *ffi,float sfcf,float **r,float rfcf,float **sm,float cfmax,float dttm,float **ssp,float cfr,float focfmax,float **smw,float **inp,float **sr,float whc,float **qd,float **ssm,float **fc,float **qin,float **beta,float **etp,float **etpo,float ecevpfo,float ecalt,float *dee,float **lp,float **qc,float **cflux,float **qf,float **kf,float **ssw,float **alpha,float **qs,float **ks,float **sgw,float **qt,float *area,float tcon,float **sgwx,float sgwmax,float **perc,float **eta){
	// snow pack balance components (mm/d)
	p[b][t]		= (1+0)*po[b][t]*(1+pcalt*dep[b]);		// Error
	ta[b][t]	= tm[b][t]+0-tcalt*det[b] ;			// Error
	if( ta[b][t] < (tt-tti/2)){
		s[b][t]	= p[b][t]*(ffo[b]*foscf+ffi[b])*sfcf;
	} else if (ta[b][t] >= (tt-tti/2) && ta[b][t] < (tt+tti/2)){
		s[b][t]	= p[b][t] * ((tt+tti/2)-ta[b][t])/tti*(ffo[b]*foscf+ffi[b])*sfcf;
		r[b][t]	= p[b][t] * ( ta[b][t] - (tt-tti/2 ) ) / tti * rfcf;
	}else{
		r[b][t]	= p[b][t] * rfcf;
	}
	sm[b][t]	= MIN(MAX(cfmax*(ta[b][t]-(tt+dttm)),0.),ssp[b][t]);
	sr[b][t]	= MIN(MAX(cfr*(ffo[b]*focfmax+ffi[b])*cfmax*((tt+dttm)-ta[b][t]),0.),smw[b][t]);
	inp[b][t]	= MAX(smw[b][t]+sm[b][t]+r[b][t]-sr[b][t]-whc*ssp[b][t],0.);
	// soil moisture balance components (mm/d)
	qd[b][t]	= MAX((inp[b][t]+ssm[b][t]-fc[b][n]),0.);
	// fc[b][n] (field capacity) is a Monte-Carlo-sampled parameter, not
	// a constant -- guard every division by it (and by lp[b][n]*fc[b][n]
	// below) against a sampled value of exactly 0, which would
	// otherwise produce inf/nan that then propagates forward forever
	// through ssm[b][t+1]'s recursion. fc<=0 has no valid soil-moisture
	// interpretation, so qin/qc are defined as 0 in that case (no
	// interflow/capillary flux out of a field capacity that doesn't
	// exist); confirmed needed in practice for lp (dataset=original's
	// own bounds allow lp=0 for several sub-basins -- see
	// docs/raster_options.md).
	qin[b][t]	= (fc[b][n] > 0.)
			? pow(ssm[b][t]/fc[b][n],beta[b][n]) * (inp[b][t]-qd[b][t])
			: 0.;
	etp[b][t]	= etpo[b][t] * (ffo[b] * ecevpfo+ffi[b]) * (1-ecalt*dee[b]);
	{
		float lpfc = lp[b][n] * fc[b][n];
		// as lp*fc -> 0+, etp*ssm/(lp*fc) -> +inf for any ssm > 0, so
		// the mathematical limit (and this guard's result) is
		// eta = etp; only a simultaneously zero ssm (no moisture to
		// evaporate at all) forces eta = 0 regardless of lp/fc.
		eta[b][t]	= (lpfc > 0.)
				? MIN( etp[b][t], ( etp[b][t]*ssm[b][t] / lpfc ) )
				: (ssm[b][t] > 0. ? etp[b][t] : 0.);
	}
	qc[b][t]	= (fc[b][n] > 0.)
			? cflux[b][n] * ( fc[b][n]-ssm[b][t] ) / fc[b][n]
			: 0.;
	// surface water balance components (mm/d)
	qf[b][t]	= kf[b][n] * pow( ssw[b][t], 1+alpha[b][n] );
	// ground water balance components (mm/d)
	qs[b][t]	= ks[b][n] * sgw[b][t];
	// total discharge (m3/s)
	qt[b][t]	= ( qs[b][t] + qf[b][t] ) * area[b] / tcon;
	// snow pack balance (mm)
	ssp[b][t+1]	= ssp[b][t] + s[b][t] + sr[b][t] - sm[b][t];
	smw[b][t+1]	= smw[b][t] + r[b][t] - sr[b][t] + sm[b][t] - inp[b][t];
	// surface water balance (mm)
	if(sgw[b][t] >= sgwmax){
		ssw[b][t+1] = MAX(ssw[b][t]+MAX((qd[b][t]+qin[b][t]),0.)-qf[b][t]-MIN(ssw[b][t],qc[b][t]),0.);
	}else{
		ssw[b][t+1] = MAX(ssw[b][t]+MAX((qd[b][t]+qin[b][t]-perc[b][n]),0.)-qf[b][t]-MIN(ssw[b][t],qc[b][t]),0.);
	}
	//soil moisture balance (mm)
	if(ssw[b][t+1]==0.){
		qc[b][t] = ssw[b][t]+MAX( (qd[b][t]+qin[b][t]-perc[b][n]),0.)-qf[b][t];
	}else{
		qc[b][t] = MIN( ssw[b][t],qc[b][t] );
	}
	// unlike ssp/ssw just above, ssm was never clamped to >= 0 here --
	// a negative ssm[b][t] then feeds pow(ssm[b][t]/fc[b][n],
	// beta[b][n]) above (qin's formula) with a negative base and a
	// non-integer (sampled) exponent, which is mathematically
	// undefined and returns nan in C, permanently corrupting every
	// subsequent timestep's state through this same recursion.
	// Confirmed as the actual cause of a real nan run (343 of 916 days
	// nan in one dataset=original ETout report) that guarding the
	// lp[b][n]*fc[b][n] division in eta's formula, above, did not by
	// itself fix.
	ssm[b][t+1] = MAX(ssm[b][t]+inp[b][t]-qd[b][t]-qin[b][t]+qc[b][t]-eta[b][t], 0.);
	//ground water balance (mm)
	if(sgw[b][t] >= sgwmax){
		sgw[b][t+1] = (1-ks[b][n])*sgw[b][t];
	}else{
		sgwx[b][t+1] = (1-ks[b][n])*sgw[b][t]+MIN( qd[b][t]+qin[b][t] , perc[b][n] );
	}
	sgw[b][t+1] = MIN(sgwx[b][t+1] , sgwmax);
	return (1);
}
