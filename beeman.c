#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void beeman_predict(double *x, double *vx, double *ax ,double *a0x ,double *a1x, double *dt);
void beeman_correct(double *vx, double *ax ,double *a0x ,double *a1x, double *dt);
void beeman_update(double *a0x ,double *a1x);

/*
void beeman_predict(double *x, double *vx, double *ax ,double *a0x ,double *a1x, double *dt){
  (*x)+=((*vx) + (4.0*(*ax)-(*a1x))*(*dt)/6.0)*(*dt);
  (*a0x)=(*ax);
}

void beeman_correct(double *vx, double *ax ,double *a0x ,double *a1x, double *dt){
    (*vx)+=(2.0*(*ax)+5.0*(*a0x)-(*a1x))*(*dt)/6.0;
}

void beeman_update(double *a0x ,double *a1x){
    (*a1x)=(*a0x);
}
*/

void beeman_predict(double *x, double *vx, double *ax ,double *a0x ,double *a1x, double *dt){
  (*x)+=((*vx) + (4.0*(*ax)-(*a1x))*(*dt)/6.0)*(*dt);
  (*a0x)=(*ax);
}

void beeman_correct(double *vx, double *ax ,double *a0x ,double *a1x, double *dt){
  (*vx)+=(2.0*(*ax)+5.0*(*a0x)-(*a1x))*(*dt)/6.0;
}

void beeman_update(double *a0x ,double *a1x){
  (*a1x)=(*a0x);
}

double beeman_SHOtest(int Ns){

  const double Ti=0.0;
  const double Tf=10.0;

  int i;
  double p_x,p_vx,p_ax,p_a0x,p_a1x,t,dt;
  double err;
  //initial;
  t=Ti;
  p_x=0.0;
  p_vx=10.0;
  err=0.0;
  
  dt=(Tf-Ti)/((double)Ns);
  
  for(i=0;i<Ns;i++){
    beeman_predict(&p_x,&p_vx,&p_ax,&p_a0x,&p_a1x,&dt);

    //force calculus
    p_ax=-1.0*p_x;

    beeman_correct(&p_vx,&p_ax,&p_a0x,&p_a1x,&dt);

    beeman_update(&p_a0x,&p_a1x);
    
    t=t+dt;

/*    printf("%lf %lf %lf %lf\n",t,p_x,p_vx,10.0*sin(t));*/

  };

  err=fabs(p_x-10.0*sin(t))/10.0;
  return err;
};

double beeman_DHOtest(int Ns){

  const double Ti=0.0;
  const double Tf=10.0;

  int i;
  double p_x,p_vx,p_ax,p_a0x,p_a1x,t,dt,g;
  double err;
  //initial;
  t=Ti;
  p_x=0.0;
  p_vx=10.0;
  p_ax=0.0;
  p_a0x=0.0;
  p_a1x=0.0;
  err=0.0;
  
  g=0.5;
  
  dt=(Tf-Ti)/((double)Ns);
  
  for(i=0;i<Ns;i++){
    beeman_predict(&p_x,&p_vx,&p_ax,&p_a0x,&p_a1x,&dt);

    //force calculus
    p_ax=-1.0*p_x-2.0*g*p_vx;

    beeman_correct(&p_vx,&p_ax,&p_a0x,&p_a1x,&dt);

    beeman_update(&p_a0x,&p_a1x);
    
    t=t+dt;

/*    printf("%lf %lf %lf %lf\n",t,p_x,p_vx,(10.0/sqrt(1.0-g*g))*exp(-g*t)*sin(t*sqrt(1.0-g*g)));*/

  };

  err=fabs(p_x-10.0*exp(-t/2.0)*sin(t))/10.0;
  return err;
};

int main(){
  printf("Erro: %lf %lf %lf\n",beeman_oscillator_test(10),beeman_oscillator_test(100),beeman_oscillator_test(1000));


  beeman_dampedoscillator_test(1000);

  return 0;
}
