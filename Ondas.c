#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void copiar(double arr1[], double arr2[]);
int i,j,k,l,ii,iii,jj,kk,ll,mm,m,nt,nc=129, n1, n2, n3;
double T,t0,t1,t2,t3,L=0.64, c=250.0, tf=1.0; // m, m/s, s
double ci2(double tiempo);

int main(){
    
    
    double * x= malloc(nc*sizeof(double));
    double * yOld= malloc(nc*sizeof(double));
    double * yNew= malloc(nc*sizeof(double));
    //Nonecesario
    double * yFut= malloc(nc*sizeof(double));
    double *yT0= malloc(nc*sizeof(double));
    double *yT1= malloc(nc*sizeof(double));
    double *yT2= malloc(nc*sizeof(double));
    double *yT3= malloc(nc*sizeof(double));
    
    
    double r=0.5, dx,dt;
    T=2.0*L/c;
    t1=T/8;
    t2=T/4;
    t3=T/2;
    
    FILE *iniciales = fopen("cond_ini_cuerda.dat", "r");
    FILE *cuerda = fopen("datosCuerda.txt","w"); //Archivo donde guardaré los datos
    FILE *sound = fopen("sonido.txt","w");
    // Guardo en 2 arreglos cond. iniciales
    for(i=0; i<nc; i++){
        fscanf(iniciales,"%lf %lf \n", &x[i], &yOld[i] );
        yT0[i]=yOld[i];
    }
    fclose(iniciales);
  
    dx=x[1]-x[0];
    dt= r*dx/c;
    nt= tf/dt;
    double *sonido = malloc(nt*sizeof(double));
    
    
    //Este print me dio los valores n1,n2 y n3
    //printf("%lf %lf %lf \n ", t1/dt,t2/dt,t3/dt);
    n1=64;
    n2=128;
    n3=256;
    

    
    //Primer paso, condiciones iniciales
    yNew[0]=0.0;
    yNew[nc-1]=0.0;
    for(j=1; j<nc-1; j++){
        yNew[j]= yOld[j] + 0.5*pow(r,2)*(yOld[j+1] - 2.0*yOld[j] + yOld[j-1]);
        //printf("%lf\n",yNew[j]);
    }
    
    sonido[0]=yOld[(nc+1)/2];
    sonido[1]=yNew[(nc-1)/2];
    
    //Lleno yFut
    for(k=2; k<nt-1; k++){
        double s;
        
        for(l=1; l<nc; l++){
            //Cond inciciales
            yFut[0]=0.0;
            yFut[nc]=0.0;
            
            yFut[l]= 2.0*(1.0-pow(r,2))*yNew[l] - yOld[l] + pow(r,2)*(yNew[l+1]+yNew[l-1]);
            
            //Llenar el arreglo sonido[k] (desde el elemento 2 hasta nt) para el punto medio
            if (l==(nc+1)/2) {
                sonido[k] = yFut[l];
                //printf("\n%lf  %d ",sonido[k],k);
            }
            if(k==n1-2){
                yT1[0]=0;
                yT1[nc-1]=0;
                yT1[l]=2.0*(1.0-pow(r,2))*yNew[l] - yOld[l] + pow(r,2)*(yNew[l+1]+yNew[l-1]);
            }
            if(k==n2-2){
                yT2[0]=0;
                yT2[nc-1]=0;
                yT2[l]=2.0*(1.0-pow(r,2))*yNew[l] - yOld[l] + pow(r,2)*(yNew[l+1]+yNew[l-1]);
            }
            if(k==n3-2){
                yT3[0]=0;
                yT3 [nc-1]=0;
                yT3[l]=2.0*(1.0-pow(r,2))*yNew[l] - yOld[l] + pow(r,2)*(yNew[l+1]+yNew[l-1]);
            }
        }
        
        //Actualizo
        for(ii=1; ii<nc-1; ii++){
            yOld[ii]=yNew[ii];
            yNew[ii]=yFut[ii];
        }
        
        
    }
    
    //Imprimo
    printf("\n\nPRIMERA PARTE\n\n    x        y0         y1       y2       y3\n\n");
    for (mm=0; mm<nc; mm++) {
        printf("%lf  %lf  %lf  %lf  %lf \n",x[mm], yT0[mm], yT1[mm], yT2[mm], yT3[mm]);
    }
    
    //Lleno los archivo con los datos para graficar y de sonido
    for (m=0; m<nc; m++) {
        fprintf(cuerda, "%lf  %lf  %lf  %lf  %lf \n",x[m], yT0[m], yT1[m], yT2[m], yT3[m]);
    }
    for (iii=0; iii<nt-1; iii++) {
        fprintf(sound,"%lf \n",sonido[iii]);
    }
    
    
    
    
    
    
    
    // SEGUNDA PARTE
    double * yOld_2= malloc(nc*sizeof(double));
    double * yNew_2= malloc(nc*sizeof(double));
    double * yFut_2= malloc(nc*sizeof(double));
    double *yT0_2= malloc(nc*sizeof(double));
    double *yT1_2= malloc(nc*sizeof(double));
    double *yT2_2= malloc(nc*sizeof(double));
    double *yT3_2= malloc(nc*sizeof(double));
    int a,b,aa,bb,cc,d,dd;
    FILE *cuerdaP = fopen("datosCuerdaPerturbada.txt","w");
    
    //Condicion inicial
    for (aa=0; aa<nc; aa++) {
        yOld_2[aa]=0.0;
    }
    
    //Primer paso
    yNew_2[0]=0.0;
    yNew_2[nc-1]=ci2(dt);
    //printf("%lf \n\n",ci2(dt));
    for(bb=1; bb<nc-1; bb++){
        yNew_2[bb]= yOld_2[bb] + 0.5*pow(r,2)*(yOld_2[bb+1] - 2.0*yOld_2[bb] + yOld_2[bb-1]);
        //printf("%lf\n",yNew_2[bb]);
    }
    for (a=2; a<nt-1; a++) {
        for (b=1; b<nc-1; b++) {
            yFut_2[0]=0.0;
            yFut_2[nc-1]=ci2(a*dt);
            
            yFut_2[b]= 2.0*(1.0-pow(r,2))*yNew_2[b] - yOld_2[b] + pow(r,2)*(yNew_2[b+1] + yNew_2[b-1]);
            
            if(a==n1-2){
                yT1_2[0]=0;
                yT1_2[nc-1]=ci2(a*n1*dt);
                yT1_2[b]= 2.0*(1.0-pow(r,2))*yNew_2[b] - yOld_2[b] + pow(r,2)*(yNew_2[b+1] + yNew_2[b-1]);
            }
            if(a==n2-2){
                yT2_2[0]=0;
                yT2_2[nc-1]=ci2(a*n2*dt);
                yT2_2[b]= 2.0*(1.0-pow(r,2))*yNew_2[b] - yOld_2[b] + pow(r,2)*(yNew_2[b+1] + yNew_2[b-1]);
            }
            if(a==n3-2){
                yT3_2[0]=0;
                yT3_2[nc-1]=ci2(a*n3*dt);
                yT3_2[b]= 2.0*(1.0-pow(r,2))*yNew_2[b] - yOld_2[b] + pow(r,2)*(yNew_2[b+1] + yNew_2[b-1]);
            }
        }
        for (cc=0; cc<nc; cc++) {
            yOld_2[cc]=yNew_2[cc];
            yNew_2[cc]=yFut_2[cc];
        }
    }
    
    printf("\n\nSEGUNDA PARTE\n\n    x        y0         y1       y2       y3\n\n");
    for (d=0; d<nc; d++) {
        printf("%lf  %lf  %lf  %lf  %lf \n",x[d], yT0_2[d], yT1_2[d], yT2_2[d], yT3_2[d]);
    }
    for (dd=0; dd<nc; dd++) {
        fprintf(cuerdaP, "%lf  %lf  %lf  %lf  %lf \n",x[dd], yT0_2[dd], yT1_2[dd], yT2_2[dd], yT3_2[d]);
    }
    
    
    
    
    
    
    
    // TERCERA PARTE
    
    
    return 0;
}

void copiar(double arr1[], double arr2[]){ //Copia el contenido de arr1 en arr2
    int i;
    double temp[nc];
    for (i=1; i<nc-1; i++) {
        temp[i]=arr1[i];
        arr2[i]=temp[i];
    }
    
}
double ci2(double tiempo){
    return sin(3.14159265*2*c*tiempo/L);
}
