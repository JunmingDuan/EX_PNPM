double WENONI21(double v1,double v2,double v3);
double WENONI22(double v3,double v2,double v1);

/** 
 * 
 * 
 * @param v1 cell average u_{i-1}  
 * @param v2  
 * @param v3 cell average u_{i+1} 
 *  
 * @return  u point value at i-sqrt(3)/6 
 */
 
double WENONI21(double v1,double v2,double v3)
{
 
     double s1,s2;
     double sum;
     double w1,w2;
     double value;
     double q3=sqrt(3.0);
     double eps=1e-14;
 
     s1=(v1-v2)*(v1-v2)+eps;
 
     s2=(v2-v3)*(v2-v3)+eps;
 
     s1=0.5/(s1*s1);
 
     s2=0.5/(s2*s2);
 

     sum=s1+s2;
     w1=s1/sum; 
     w2=s2/sum;
 
     value=w1*(q3/6*v1+(1-q3/6)*v2)+w2*((1+q3/6)*v2-q3/6*v3);
 
     return value;
 
}
 
/** 
 * 
 * 
 * @param v3 cell average u_{i-1}  
 * @param v2 
 * @param v1 cell average u_{i+1}
 * 
 * @return  u point value at i+sqrt(3)/6 
 */
 
double WENONI22(double v3,double v2,double v1)
{
 
     double s1,s2;
     double sum;
     double w1,w2;
     double value;
     double q3=sqrt(3.0);
     double eps=1e-14;
 
     s1=(v1-v2)*(v1-v2)+eps;
 
     s2=(v2-v3)*(v2-v3)+eps;
 
     s1=0.5/(s1*s1);
     s2=0.5/(s2*s2);
 
     sum=s1+s2; 
     w1=s1/sum;
     w2=s2/sum;

     value=w1*(q3/6*v1+(1-q3/6)*v2)+w2*((1+q3/6)*v2-q3/6*v3);
 
     return value;
 
} 

