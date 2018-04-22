#ifndef _conWENO5_h_
#define _conWENO5_h_

double WENONI31(double v1,double v2,double v3,double v4,double v5);
double WENONI32(double v1,double v2,double v3,double v4,double v5);
double WENONI33(double v5,double v4,double v3,double v2,double v1);

double WENONI31(double v1,double v2,double v3,double v4,double v5)
{
     double s1,s2,s3;
     double sum;
     double tmp1,tmp2;
     double w1,w2,w3;
     double value;
     double q15=sqrt(15.0);
     double eps=1e-14;
     
     tmp1=v1-2*v2+v3;
     tmp2=v1-4*v2+3*v3;
     s1=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2+eps;
     s1=(22*q15+9)/((3*q15-2)*40)/(s1*s1);

     tmp1=v2-2*v3+v4;
     tmp2=v2-v4;
     s2=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2+eps;
     s2=403./(655*s2*s2);

     tmp1=v3-2*v4+v5;
     tmp2=3*v3-4*v4+v5;
     s3=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2+eps;
     s3=(22*q15-9)/((3*q15+2)*40)/(s3*s3);

     sum=s1+s2+s3;
     w1=s1/sum;
     w2=s2/sum;
     w3=s3/sum;

     value=w1*((1./30-0.05*q15)*v1+(0.2*q15-1./15)*v2+(31./30-0.15*q15)*v3)
		 +w2*((1./30+0.05*q15)*v2+14./15*v3+(1./30-0.05*q15)*v4)
		 +w3*((31./30+0.15*q15)*v3-(1./15+0.2*q15)*v4+(1./30+0.05*q15)*v5);
     
     return value;
}

double WENONI32(double v1,double v2,double v3,double v4,double v5)
{
     double s1,s2,s3;
     double sum,eps=1e-20;
     double tmp1,tmp2;
     double w1,w2,w3;
     double value,p1,p2,p3;

     tmp1=v1-2*v2+v3;
     tmp2=v1-4*v2+3*v3;
     s1=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2+eps;


     tmp1=v2-2*v3+v4;
     tmp2=v2-v4;
     s2=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2+eps;
   

     tmp1=v3-2*v4+v5;
     tmp2=3*v3-4*v4+v5;
     s3=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2+eps;

     p1=(-v1+2*v2+23*v3)/24.0;
     p2=(-v2+26*v3-v4)/24.0;
     p3=(23*v3+2*v4-v5)/24.0;

    
     w1=9.0/(214*s1*s1);
     w2=98.0/(107*s2*s2);
     w3=9.0/(214*s3*s3);
 

     sum=w1+w2+w3;
     w1=w1/sum;
     w2=w2/sum;
     w3=w3/sum;

     value=107.0/40*(w1*p1+w2*p2+w3*p3);

     w1=9.0/(67*s1*s1);
     w2=49.0/(67*s2*s2);
     w3=9.0/(67*s3*s3);

     sum=w1+w2+w3;
     w1=w1/sum;
     w2=w2/sum;
     w3=w3/sum;

     value-=67.0/40*(w1*p1+w2*p2+w3*p3);

     return value;
}

double WENONI33(double v5,double v4,double v3,double v2,double v1)
{
     double s1,s2,s3;
     double sum;
     double tmp1,tmp2;
     double w1,w2,w3;
     double value;
     double q15=sqrt(15.0);
     double eps=1e-14;
     
     tmp1=v1-2*v2+v3;
     tmp2=v1-4*v2+3*v3;
     s1=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2+eps;
     s1=(22*q15+9)/((3*q15-2)*40)/(s1*s1);

     tmp1=v2-2*v3+v4;
     tmp2=v2-v4;
     s2=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2+eps;
     s2=403./(655*s2*s2);

     tmp1=v3-2*v4+v5;
     tmp2=3*v3-4*v4+v5;
     s3=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2+eps;
     s3=(22*q15-9)/((3*q15+2)*40)/(s3*s3);

     sum=s1+s2+s3;
     w1=s1/sum;
     w2=s2/sum;
     w3=s3/sum;

     value=w1*((1./30-0.05*q15)*v1+(0.2*q15-1./15)*v2+(31./30-0.15*q15)*v3)+w2*((1./30+0.05*q15)*v2+14./15*v3+(1./30-0.05*q15)*v4)+w3*((31./30+0.15*q15)*v3-(1./15+0.2*q15)*v4+(1./30+0.05*q15)*v5);
     return value;
}


#endif
