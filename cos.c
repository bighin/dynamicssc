#include <math.h>
#include <complex.h>
#include "cubature/cubature.h"

#include "auxx.h"
#include "config.h"
#include "dsc.h"
#include "bigpsi.h"

struct across_container_t
{
	struct interpolation_t *intre1;
	struct interpolation_t *intim1;

	struct interpolation_t *intre2;
	struct interpolation_t *intim2;
	
	int L,Lprime;
	
	struct configuration_t *config;
};

int fAcross(unsigned ndim,const double *x,void *fdata,unsigned fdim,double *fval)
{
        /* The context from which we read the global variables */

        struct across_container_t *container=(struct across_container_t *)(fdata);
	struct configuration_t *config=container->config;

        /* The integration variables and other auxiliary variables */

        double k;
	double complex f;
	
	k=x[0];

	f=conj(get_point(container->intre1,k)+I*get_point(container->intim1,k));
	f*=get_point(container->intre2,k)+I*get_point(container->intim2,k);
	f/=fscale(k,container->L,config)*fscale(k,container->Lprime,config);

	fval[0]=creal(f);
	fval[1]=cimag(f);

	return 0;
}

double complex Across(struct bigpsi_t *psi,int L,int Lprime,int n,struct configuration_t *config)
{
	struct across_container_t container;
	double *x,*y1re,*y1im,*y2re,*y2im;

	double xmin,xmax,res[2],err[2];
	int c;

	extern double relError;
	extern size_t maxEval;

	int offsetL=L*(2+10*config->gridpoints);
	int offsetLprime=Lprime*(2+10*config->gridpoints);

	x=malloc(sizeof(double)*config->gridpoints);
	y1re=malloc(sizeof(double)*config->gridpoints);
	y1im=malloc(sizeof(double)*config->gridpoints);
	y2re=malloc(sizeof(double)*config->gridpoints);
	y2im=malloc(sizeof(double)*config->gridpoints);

        for(c=0;c<config->gridpoints;c++)
	{
                double gridstep=config->cutoff/config->gridpoints;
                double k=c*gridstep;
		int noffset;

		double complex phase1,phase2;

                x[c]=k;

		switch(n)
		{
			case -2:
			phase1=timephase(-(L*(L+1.0f)+omegak(k)-2.0f),psi->t,config);
			phase2=timephase(-(Lprime*(Lprime+1.0f)+omegak(k)-2.0f),psi->t,config);
			noffset=0;
			break;

			case -1:
			phase1=timephase(-(L*(L+1.0f)+omegak(k)+4.0f),psi->t,config);
			phase2=timephase(-(Lprime*(Lprime+1.0f)+omegak(k)+4.0f),psi->t,config);
			noffset=2;
			break;

			case 0:
			phase1=timephase(-(L*(L+1.0f)+omegak(k)+6.0f),psi->t,config);
			phase2=timephase(-(Lprime*(Lprime+1.0f)+omegak(k)+6.0f),psi->t,config);
			noffset=4;
			break;

			case 1:
			phase1=timephase(-(L*(L+1.0f)+omegak(k)+4.0f),psi->t,config);
			phase2=timephase(-(Lprime*(Lprime+1.0f)+omegak(k)+4.0f),psi->t,config);
			noffset=6;
			break;

			case 2:
			phase1=timephase(-(L*(L+1.0f)+omegak(k)-2.0f),psi->t,config);
			phase2=timephase(-(Lprime*(Lprime+1.0f)+omegak(k)-2.0f),psi->t,config);
			noffset=8;
			break;
			
			default:
			phase1=phase2=0.0f;
			noffset=0;
			fprintf(stderr,"Bug in Across()\n");
			exit(0);
		}

                y1re[c]=creal(phase1*(psi->y[offsetL+2+10*c+noffset]+I*psi->y[offsetL+2+10*c+noffset+1]));
                y1im[c]=cimag(phase1*(psi->y[offsetL+2+10*c+noffset]+I*psi->y[offsetL+2+10*c+noffset+1]));

                y2re[c]=creal(phase2*(psi->y[offsetLprime+2+10*c+noffset]+I*psi->y[offsetLprime+2+10*c+noffset+1]));
                y2im[c]=cimag(phase2*(psi->y[offsetLprime+2+10*c+noffset]+I*psi->y[offsetLprime+2+10*c+noffset+1]));
	}

        container.intre1=init_interpolation(x,y1re,config->gridpoints);
        container.intim1=init_interpolation(x,y1im,config->gridpoints);
        container.intre2=init_interpolation(x,y2re,config->gridpoints);
        container.intim2=init_interpolation(x,y2im,config->gridpoints);
	
	container.L=L;
	container.Lprime=Lprime;
	container.config=config;

	xmin=0.0f;
	xmax=config->cutoff;

	hcubature(2,fAcross,&container,1,&xmin,&xmax,maxEval,0,10*relError,ERROR_INDIVIDUAL,res,err);

	fini_interpolation(container.intre1);
	fini_interpolation(container.intim1);
	fini_interpolation(container.intre2);
	fini_interpolation(container.intim2);

	if(x)	free(x);
	if(y1re) free(y1re);
	if(y1im) free(y1im);
	if(y2re) free(y2re);
	if(y2im) free(y2im);

	return res[0]+I*res[1];
}

/*
	Please refer to my notes, here's a precalculated version of f_{\lambda 0},
	i.e. the coefficients in the spherical harmonics expansion of the alignment
	cosine operator.
*/

double ftable[]=
{
   1.7724538509055159,0.,0.9908318244015027,0.,-0.22155673136318949,0.,
   0.09985426941927254,0.,-0.057093862842314214,0.,0.0370161637780817,0.,
   -0.025963679456623767,0.,0.019225070270110573,0.,-0.01481142194011715,0.,
   0.011762561565771736,0.,-0.009567957034891546,0.,0.00793552938587553,0.,
   -0.006688267492838287,0.,0.0057137776188041075,0.,-0.00493789192722724,0.,
   0.004310055435583415,0.,-0.0037948419927762893,0.,0.0033668266630515983,0.,
   -0.0030073777743079306,0.,0.0027025887919370504,0.,-0.0024419094292527006,0.,
   0.0022172182050591944,0.,-0.002022180887315326,0.,0.0018517983401647546,0.,
   -0.0017020824461301998,0.,0.0015698202494509703,0.,-0.001452399897833141,0.,
   0.0013476805426497135,0.,-0.001253893950288784,0.,0.001169569287323113,0.,
   -0.0010934750439840373,0.,0.001024573773244495,0.,-0.000961986512008685,0.,
   0.0009049645873477925,0.,-0.0008528671062057122,0.,0.0008051428557483715,0.,
   -0.0007613156535040531,0.,0.0007209724157007564,0.,-0.0006837533822523379,0.,
   0.0006493440640782344,0.,-0.0006174685744283888,0.,0.000587884078864593,0.,
   -0.0005603761544462608,0.,0.0005347548917807982,0.,-0.0005108516070709821,0.,
   0.0004885160574433798,0.,-0.000467614073396409,0.,0.0004480255384542761,0.,
   -0.0004296426590261236,0.,0.00041236847778593455,0.,-0.000396115592171254,0.,
   0.0003808050462800969,0.,-0.00036636536985956513,0.,0.00035273174248615536,0.,
   -0.0003398452646390367,0.,0.00032765232032265626,0.,-0.0003161040183291965,0.,
   0.0003051557012439238,0.,-0.0002947665129663062,0.,0.00028489901691005774,0.,
   -0.0002755188582066143,0.,0.0002665944642097673,0.,-0.0002580967784172574,0.,
   0.00024999902361484786,0.,-0.000242276490631557,0.,0.00023490634958918229,0.,
   -0.00022786748094959472,0.,0.00022114032402157868,0.,-0.00021470674089515287,0.,
   0.00020854989403356957,0.,-0.0002026541359783689,0.,0.00019700490981665601,0.,
   -0.00019158865922691567,0.,0.00018639274706415702,0.,-0.00018140538157032489,0.,
   0.00017661554940453868,0.,-0.00017201295478218258,0.,0.00016758796409418602,0.,
   -0.00016333155544969368,0.,0.00015923527264816772,0.,-0.00015529118314203184,0.,
   0.00015149183959928977,0.,-0.0001478302447180367,0.,0.0001442998189821943,0.,
   -0.00014089437108079508,0.,0.00013760807074228988,0.,-0.0001344354237611405,0.,
   0.0001313712490168061,0.,-0.0001284106573055064,0.,0.0001255490318231585,0.,
   -0.00012278201015391234,0.,0.0001201054676329925,0.,-0.00011751550196529451,0.,
   0.00011500841899256702,0.,-0.00011258071951219113,0.,0.00011022908705968641,0.,
   -0.00010795037657524743,0.,0.00010574160388195363,0.,-0.00010359993590989046,0.,
   0.00010152268160635366,0.,-0.00009950728347765354,0.,0.00009755130971285613,0.,
   -0.00009565244684414804,0.,0.0000938084929024433,0.,-0.00009201735103040629,0.,
   0.00009027702351828418,0.,-0.00008858560623086051,0.,0.0000869412833964882,0.,
   -0.00008534232273156457,0.,0.00008378707087599468,0.,-0.00008227394911717742,0.,
   0.00008080144938185798,0.,-0.00007936813047683893,0.,0.00007797261456104714,0.,
   -0.0000766135838328252,0.,0.00007528977741757211,0.,-0.00007399998844200294,0.,
   0.0000727430612823485,0.,-0.00007151788897477553,0.,0.0000703234107771898,0.,
   -0.00006915860987239174,0.,0.00006802251120329631,0.,-0.00006691417943160906,0.,
   0.00006583271701197745,0.,-0.00006477726237421116,0.,0.00006374698820669651,0.,
   -0.00006274109983461684,0.,0.0000617588336870424,0.,-0.000060799455847366675,0.,
   0.00005986226068195042,0.,-0.00005894656954218794,0.,0.00005805172953553696,0.,
   -0.000057177112361355875,0.,0.00005632211320767206,0.,-0.00005548614970526337,0.,
   0.0000546686609356752,0.,-0.00005386910649001778,0.,0.00005308696557559418,0.,
   -0.00005232173616760136,0.,0.000051572934203323756,0.,-0.00005084009281640412,0.,
   0.00005012276160892969,0.,-0.00004942050595921443,0.,0.00004873290636329036,0.,
   -0.00004805955780824537,0.,0.00004740006917565892,0.,-0.0000467540626734952,0.,
   0.000046121173294912635,0.,-0.00004550104830254204,0.,0.00004489334673687293,0.,
   -0.000044297738947468434,0.,0.000043713906145805336,0.,-0.00004314153997860662,0.,
   0.00004258034212060013,0.,-0.00004203002388569925,0.,0.00004149030585565898,0.,
   -0.00004096091752531593,0.,0.00004044159696357079,0.,-0.0000399320904893201,0.,
   0.00003943215236158877,0.,-0.0000389415444831563,0.,0.00003846003611700959,0.,
   -0.000037987403614991656,0.,0.000037523430158050435,0.,-0.00003706790550752439,0.,
   0.00003662062576693226,0.,-0.000036181393153762805,0.,0.00003575001578078794,0.,
   -0.00003532630744644735,0.,0.00003491008743387752,0.,-0.00003450118031817971,0.,
   0.0000340994157815435,0.,-0.00003370462843586173,0.,0.00003331665765249222,0.,
   -0.000032935347398838905,0.,0.00003256054608144192,0.,-0.000032192106395282226,0.,
   0.00003182988517902086,0.,-0.00003147374327590731,0.,0.00003112354540010468,0.,
   -0.00003077916000819182,0.,0.0000304404591756147,0.,-0.00003010731847787014,0.,
   0.000029779616876216118,0.,-0.000029457236607712406,0.,0.000029140063079405156,0.,
   -0.000028827984766477927,0.,0.000028520893114200012,0.,-0.000028218682443511312,0.,
   0.00002792124986009021,0.,-0.00002762849516675855,0.,0.000027340320779084406,0.,
   -0.000027056631644049903,0.,0.000026777335161657594,0.,-0.000026502341109354594,0.,
   0.00002623156156915938,0.,-0.000025964910857381352,0.,0.000025702305456828216,0.,
   -0.00002544366395140111,0.,0.00002518890696298178,0.,-0.000024937957090520544,0.,
   0.000024690738851237678,0.,-0.000024447178623854906,0.,0.000024207204593777233,0.,
   -0.00002397074670014899,0.,0.000023737736584711112,0.,-0.000023508107542390035,0.,
   0.000023281794473551464,0.,-0.000023058733837855308,0.,0.000022838863609650566,0.,
   -0.000022622123234851863,0.,0.000022408453589241484,0.,-0.000022197796938143516,0.,
   0.000021990096897418523,0.,-0.00002178529839572974,0.,0.000021583347638033492,0.,
   -0.000021384192070248775,0.,0.000021187780345062584,0.,-0.000020994062288829533,0.,
   0.00002080298886952582,0.,-0.000020614512165719455,0.,0.000020428585336519988,0.,
   -0.00002024516259247259,0.,0.000020064199167362757,0.,-0.000019885651290899155,0.,
   0.00001970947616224353,0.,-0.000019535631924357824,0.,0.000019364077639139695,0.,
   -0.00001919477326331898,0.,0.00001902767962508853,0.,-0.000018862758401443988,0.,
   0.000018699972096208025,0.,-0.000018539284018715508,0.,0.000018380658263136982,0.,
   -0.000018224059688418725,0.,0.000018069453898818433,0.,-0.000017916807225016403,0.,
   0.000017766086705782886,0.,-0.00001761726007018294,0.,0.000017470295720300834,0.,
   -0.000017325162714466803,0.,0.000017181830750969468,0.,-0.000017040270152237887,0.,
   0.00001690045184947799,0.,-0.000016762347367748288,0.,0.000016625928811460802,0.,
   -0.000016491168850293253,0.,0.00001635804070549935,0.,-0.000016226518136604256,0.,
   0.000016096575428473032,0.,-0.000015968187378740023
};

double flambda(int lambda)
{
	return ftable[lambda];
}

double complex old_fcostheta2d(struct bigpsi_t *psi,int L,int Lprime,int M,int lambda,struct configuration_t *config)
{
	double f=flambda(lambda);
	double complex gL,gLprime,res;
	int n;

	int offsetL=L*(2+10*config->gridpoints);
	int offsetLprime=Lprime*(2+10*config->gridpoints);

	gL=(psi->y[offsetL+0]+I*psi->y[offsetL+1])*timephase(-L*(L+1),psi->t,config);
	gLprime=(psi->y[offsetLprime+0]+I*psi->y[offsetLprime+1])*timephase(-Lprime*(Lprime+1),psi->t,config);

	res=conj(gL)*gLprime*sqrtf((2.0f*L+1.0f)/(2.0f*Lprime+1.0f))*sqrtf((2.0f*lambda+1.0f)/(4*M_PI));
	res*=cg(L,M,lambda,0,Lprime,M)*cg(L,0,lambda,0,Lprime,0);

	if((config->freeevolution)||(config->altcos))
		return f*res;

	for(n=-2;n<=2;n++)
	{
		double cgs;
		
		cgs=cg(L,M,lambda,0,Lprime,M)*cg(L,n,lambda,0,Lprime,n);
		
		if(fabs(cgs)>1e-10)
		{
			double complex tmp;

			tmp=Across(psi,L,Lprime,n,config)*sqrtf((2.0f*L+1.0f)/(2.0f*Lprime+1.0f))*sqrtf((2.0f*lambda+1.0f)/(4*M_PI));
			tmp*=cgs;
		
			res+=tmp;
		}
	}

	return f*res;
}

double complex old_costheta2d(struct bigpsi_t *psi,struct configuration_t *config)
{
	int L,Lprime,M,lambda;
	double complex total=0.0f;

	M=psi->params[0].M;

	for(L=0;L<config->maxl;L++)
	{
		for(Lprime=0;Lprime<config->maxl;Lprime++)
		{
			for(lambda=0;lambda<config->maxl;lambda+=2)
			{
				total+=old_fcostheta2d(psi,L,Lprime,M,lambda,config);
			}
		}
	}

	return total;
}

double complex fcostheta2d(struct bigpsi_t *psi,int L,int Lprime,int M,int lambda,struct configuration_t *config)
{
	double f=flambda(lambda);
	double complex gL,gLprime,tmp,res=0.0f;
	int n;

	int offsetL=L*(2+10*config->gridpoints);
	int offsetLprime=Lprime*(2+10*config->gridpoints);

	gL=(psi->y[offsetL+0]+I*psi->y[offsetL+1])*timephase(-L*(L+1),psi->t,config);
	gLprime=(psi->y[offsetLprime+0]+I*psi->y[offsetLprime+1])*timephase(-Lprime*(Lprime+1),psi->t,config);

	tmp=conj(gL)*gLprime*sqrtf((2.0f*L+1.0f)/(2.0f*Lprime+1.0f))*sqrtf((2.0f*lambda+1.0f)/(4*M_PI));
	tmp*=cg(L,M,lambda,0,Lprime,M)*cg(L,0,lambda,0,Lprime,0);
	res+=tmp;

	gL=(psi->y[offsetLprime+0]+I*psi->y[offsetLprime+1])*timephase(-Lprime*(Lprime+1),psi->t,config);
	gLprime=(psi->y[offsetL+0]+I*psi->y[offsetL+1])*timephase(-L*(L+1),psi->t,config);

	tmp=conj(gL)*gLprime*sqrtf((2.0f*Lprime+1.0f)/(2.0f*L+1.0f))*sqrtf((2.0f*lambda+1.0f)/(4*M_PI));
	tmp*=cg(Lprime,M,lambda,0,L,M)*cg(Lprime,0,lambda,0,L,0);
	res+=tmp;

	if((config->freeevolution)||(config->altcos))
		return f*res;

	for(n=-2;n<=2;n++)
	{
		double cgs;
		double complex across;
		int across_has_been_calculated=0;
		
		cgs=cg(L,M,lambda,0,Lprime,M)*cg(L,n,lambda,0,Lprime,n);
		
		if(fabs(cgs)>1e-10)
		{
			double complex tmp;

			across=Across(psi,L,Lprime,n,config);
			across_has_been_calculated=1;

			tmp=across*sqrtf((2.0f*L+1.0f)/(2.0f*Lprime+1.0f))*sqrtf((2.0f*lambda+1.0f)/(4*M_PI));
			tmp*=cgs;

			res+=tmp;
		}

		cgs=cg(Lprime,M,lambda,0,L,M)*cg(Lprime,n,lambda,0,L,n);

		if(fabs(cgs)>1e-10)
		{
			if(across_has_been_calculated==0)
				across=Across(psi,L,Lprime,n,config);

			across=conj(across);

			tmp=across*sqrtf((2.0f*Lprime+1.0f)/(2.0f*L+1.0f))*sqrtf((2.0f*lambda+1.0f)/(4*M_PI));
			tmp*=cgs;

			res+=tmp;
		}
	}

	return f*res;
}

double complex costheta2d(struct bigpsi_t *psi,struct configuration_t *config)
{
	int L,Lprime,M,lambda;
	double complex total=0.0f;

	M=psi->params[0].M;

	for(L=0;L<config->maxl;L++)
	{
		for(Lprime=0;Lprime<L;Lprime++)
		{
			for(lambda=0;lambda<config->maxl;lambda+=2)
			{
				total+=fcostheta2d(psi,L,Lprime,M,lambda,config);
			}
		}

		for(lambda=0;lambda<config->maxl;lambda++)
		{
			total+=old_fcostheta2d(psi,L,L,M,lambda,config);
		}
	}

	return total;
}

double complex fcosthetasquared(struct bigpsi_t *psi,int L,int Lprime,int M,int lambda,struct configuration_t *config)
{
	double f;
	double complex gL,gLprime,res;
	int n;

	int offsetL=L*(2+10*config->gridpoints);
	int offsetLprime=Lprime*(2+10*config->gridpoints);

	if(lambda==0)
		f=(1.0f/3.0f)*sqrt(4.0f*M_PI);

	if(lambda==2)
		f=(4.0/3.0)*sqrt(M_PI/5.0f);

	gL=(psi->y[offsetL+0]+I*psi->y[offsetL+1])*timephase(-L*(L+1),psi->t,config);
	gLprime=(psi->y[offsetLprime+0]+I*psi->y[offsetLprime+1])*timephase(-Lprime*(Lprime+1),psi->t,config);

	res=conj(gL)*gLprime*sqrtf((2.0f*L+1.0f)/(2.0f*Lprime+1.0f))*sqrtf((2.0f*lambda+1.0f)/(4*M_PI));
	res*=cg(L,M,lambda,0,Lprime,M)*cg(L,0,lambda,0,Lprime,0);

	if(config->freeevolution==true)
		return f*res;

	for(n=-2;n<=2;n++)
	{
		double cgs;
		
		cgs=cg(L,M,lambda,0,Lprime,M)*cg(L,n,lambda,0,Lprime,n);
		
		if(fabs(cgs)>1e-10)
		{
			double complex tmp;

			tmp=Across(psi,L,Lprime,n,config)*sqrtf((2.0f*L+1.0f)/(2.0f*Lprime+1.0f))*sqrtf((2.0f*lambda+1.0f)/(4*M_PI));
			tmp*=cgs;
		
			res+=tmp;
		}
	}

	return f*res;
}

double complex costhetasquared(struct bigpsi_t *psi,struct configuration_t *config)
{
	int L,Lprime,M,lambda;
	double complex total=0.0f;

	M=psi->params[0].M;

	for(L=0;L<config->maxl;L++)
	{
		for(Lprime=0;Lprime<config->maxl;Lprime++)
		{
			for(lambda=0;lambda<=2;lambda+=2)
			{
				total+=fcosthetasquared(psi,L,Lprime,M,lambda,config);
			}
		}
	}

	return total;
}
