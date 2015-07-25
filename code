#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include <string.h>
#include <fstream.h>
#include <stdio.h>

double f1(double,double);
double f2(double,double);
double G(double,double);
double pi=3.141592653589793;
double t=120;
int main()
 {
	double lp,l1,l2,l3,lG1,lG2,kG1,kG2,D,LA,LB,LAA,LBB,l,lc,f,f1,F1,F2,fmn,fm1n1,Gmn,Gm1n1,gmn,gm1n1,tao,tao1,x1,x2,x3,x4,index1,index2;//l是正畴的长度，
	int m,n,m1,n1,p,q,p1,q1,x;
	FILE *fp;
	
	//带隙出现的位置是1310nm和1550nm两个波段
	
	
	l1=1.310;
	kG1=2*pi/(l1/2);
	lG1=2*pi/kG1;
	l2=1.550;
    kG2=2*pi/(l2/2);
    lG2=2*pi/kG2;
    

	//index1=1.435;       //两种金属折射率的实部
	//index2=1.669;

	f1=0;
	F1=0;
	F2=0;

	for(m=0;m<=3;m++)
	{
		for(n=0;n<=3;n++)
		{
           	if(m==0&&n==0)
				continue;
			for(m1=0;m1<=3;m1++)
			{
	            for(n1=0;n1<=3;n1++)
				{             				
					
					if(m1==0||n1==0)
						continue;
					if(m==0||n==0)
						continue;
					x=n*n1;                                    //x是一个判断量，使得n和n1不为0，以免方程无解
					if(x==0)
						continue;
				    
					tao=(m1*kG1-m*kG2)/(n*kG2-n1*kG1);         //计算tao 和 平均结构参数
					if(tao<0)
						continue;
					D=2*pi*(m+n*tao)/kG1;
					

					for(LB=0.01;LB<D-0.01;LB=LB+0.01)                //LB从2um开始取
					{
						LA=(D-LB)/tao;
					    if(LA>LB)
						{
							
						  for(l=0.01;l<=LB-0.01;l=l+0.01)
						  {
						                                   
							 Gmn=2*pi*(m+n*tao)/D;
							 x1=Gmn*l/2;
                             x2=pi*(1+tao)*(m*LA-n*LB)/D;
                             fmn=2*(1+tao)*(l/D)*sin(x1)/x1*sin(x2)/x2;
                            //  if(fmn<0.45)
							//	 continue;
							 
							 Gm1n1=2*pi*(m1+n1*tao)/D;
							 x3=Gm1n1*l/2;
							 x4=pi*(1+tao)*(m1*LA-n1*LB)/D;
					 		 fm1n1=2*(1+tao)*(l/D)*sin(x3)/x3*sin(x4)/x4;
                             f=fmn*fm1n1;						
							 if(f>f1)
							 {
								 f1=f;
								 F1=fmn;                                //F1,F2用来存储的得到的最大的傅立叶系数
                                 F2=fm1n1;
								 p=m;q=n;p1=m1;q1=n1;
								 LBB=LB;
								 LAA=LA;
								 tao1=tao;
                                 lc=l;
							     gmn=Gmn;
								 gm1n1=Gm1n1;
							 }
						  } 
                          
						}  
						  else
						  {
                            for(l=0.01;l<=LA-0.01;l=l+0.01)
							{
							                            
							   Gmn=2*pi*(m+n*tao)/D;							  
							   x1=Gmn*l/2;
                               x2=pi*(1+tao)*(m*LA-n*LB)/D;
                               fmn=2*(1+tao)*(l/D)*sin(x1)/x1*sin(x2)/x2;
                             //  if(fmn<0.45)
							 //	   continue;
							   Gm1n1=2*pi*(m1+n1*tao)/D;
						 	   x3=Gm1n1*l/2;
							   x4=pi*(1+tao)*(m1*LA-n1*LB)/D;
					 		   fm1n1=2*(1+tao)*(l/D)*sin(x3)/x3*sin(x4)/x4;
                               f=fmn*fm1n1;
							
						
							   if(f>f1)
							   {
							  	 f1=f;
								 F1=fmn;                                //F1,F2用来存储的得到的最大的傅立叶系数
                                 F2=fm1n1;
								 p=m;q=n;p1=m1;q1=n1;
								 LBB=LB;
								 LAA=LA;
							     tao1=tao;
							     lc=l;
							     gmn=Gmn;
								 gm1n1=Gm1n1;
							   }
							}					
						  }					
					}
				}
			}	
		}
	}
	fp=fopen("e:\\quasiperiod_structure_double_PBG.txt","w");
	fprintf(fp,"lG1=%f, lG2=%f\n",lG1,lG2);
	fprintf(fp,"kG1=%f, kG2=%f\n",gmn,gm1n1);
	fprintf(fp,"LA=%f, LB=%f\n",LAA,LBB);
	fprintf(fp,"m=%d, n=%d, fmn=%f\n",p,q,F1);
    fprintf(fp,"m1=%d, n1=%d, fm1n1=%f\n",p1,q1,F2);
	fprintf(fp,"tao=%f\n",tao1);
	fprintf(fp,"l=%f",lc);
	fclose(fp);
	
 }



double G(double lp,double l1)
{
	double np,n1,n2,l2,kp,k1,k2,kG;
	np=f1(lp,t);n1=f1(l1,t);
	l2=lp*l1/(l1-lp);
	n2=f1(l2,t);
    kp=2*pi*np/lp;k1=2*pi*n1/l1;k2=2*pi*n2/l2;
	kG=kp-k1-k2;
	return kG;
}




double f1(double w,double t)
{
   double b, c, z, n;
    b=0.000000026794*(t+273.15)*(t+273.15);
    c=0.000000016234*(t+273.15)*(t+273.15);
    z=w*w;
    n=sqrt(4.5284+(0.0072449+b)/(z-(0.2453+c)*(0.2453+c))+0.07769/(z-0.1838*0.1838)-0.02367* z);
    return n;
}


double f2(double w,double t1)
{
   double a1,a2,a3,a4,b1,b2,b3,c1,c,c2,n;
   a1=40582;
   a2=0.09921;
   a3=0.2109;
   a4=-0.0219;
   b1=0.00000022971;
   b2=0.000000052716;
   b3=-0.000000049143;
   c=88506.25;
   t1=t+273.15;
   c1=w*w-(a3+b3*(t1*t1-c)*(t1*t1-c));
   c2=a1+b1*(t1*t1-c)+(a2+b2*(t1*t1-c))/c1+a4*w*w;
   n=sqrt(c2);
   return n;
}
