#include <stdio.h>  /* directives au préprocesseur */
//#include <iostream>
#include <math.h>
#include <stdlib.h>
//#include <fstream>
//#include <string>
//#include <Rinternals.h>
//#include <gsl/gsl_sf_log.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_blas.h>

#include <R.h>
#include "Convert.h"


/******************************************************************************/
/*            Calcul de la somme intervenant dans la formule du Tau           */
/******************************************************************************/  

double CalculSommeTau(double **Pi,double **Tau,double **G,int debut, int fin,int j,int t)  
{
   //Déclaration et initialisation des variables locales;      
   int x;
   double res;
   res=0;

   //Boucle permettant de calculer la somme produit entre Pi, Tau et 1/G;
   for(x=debut; x<fin; x++){ 
             res=res+(Pi[j][x]*Tau[t][x])/G[t][x];
   }
   
   //Retourne la somme des éléments;
   return res;
}    


/******************************************************************************/
/*             Calcul de la somme du produit entre Pi et F,                   */
/*            identique à sumFPi sauf pour les indices du Pi                  */
/******************************************************************************/
    
  
  double sumFPi2(double **Pi,double **F, int debut, int fin,int j,int t)  
{
   //Déclaration et initialisation des variables locales;      
   int x;
   double res;
   res=0;

   //Boucle permettant de calculer la somme produit entre Pi et F;
   for(x=debut; x<fin; x++){ 
             res=res+(Pi[x][j])*(F[t][x]);
   }

   //Retourne la somme des éléments;
   return res;
}    

/******************************************************************************/
/*                               Calcul du G                                  */
/******************************************************************************/  

void CalculG(double **Pi,double **F, int nbClass,int nbInd,double **G)  
{
     //Déclaration et initialisation des variables locales;   
     int x,y;
     
     //Double boucle permettant de calculer G à tout les instants et pour toutes les classes;
     for(x=0; x<(nbInd-1); x++){
                     for(y=0; y<nbClass; y++){                                      
                             G[x+1][y] =sumFPi2(Pi,F,0,nbClass,y,x);   
                             G[0][y] = 0;                                
                    }         
     }   
}


/******************************************************************************/
/*                              Calcul des Tau                                */
/******************************************************************************/  

void CalculTau(double **Pi,double **F, int nbClass,int nbInd,double **G,double **Tau)  
{
     //Déclaration et initialisation des variables locales;   
     int x,y;
     double NormF = 0;
     
     //Calcul de la matrice G nécessaire au calcul des Tau;
     CalculG(Pi,F,nbClass,nbInd,G);
     
     //Double boucle permettant de calculer Tau à tout les instants t(sauf le dernier) et pour toutes les classes;     
     for(x=(nbInd-2); x>=0; x--){
                     for(y=0; y<nbClass; y++){                                   
                             Tau[x][y] =  F[x][y]*CalculSommeTau(Pi,Tau,G,0,nbClass,y,x+1);
                              
                    }         
     }
}






/******************************************************************************/
/*         Calcul du premier instant de l'étape Backward  (Tau[N][.])         */
/******************************************************************************/

void CalculTauN(double **F,int nbInd, int nbClass, double **Tau)  
{
    int y;
    for(y=0; y<nbClass; y++){ 
    Tau[nbInd-1][y]=F[nbInd-1][y];
    }
}   


/******************************************************************************/
/*       Etape Backward de l'algorithme Forward-Backward pour les HMM         */
/******************************************************************************/
    
void Backward(double *FVect, double *PiVect,int *nbInd, int *nbClass, double *TauVect, double *GVect)  
{
     int i;
    
    double **F = (double**)malloc(*nbInd*sizeof(double*));
            for(i=0;i<*nbInd;i++){ 
                   F[i] = (double*)malloc(*nbClass*sizeof(double));
            }
            
      double **Pi = (double**)malloc(*nbClass*sizeof(double*));
             for(i=0;i<*nbClass;i++){ 
                   Pi[i] = (double*)malloc(*nbClass*sizeof(double));
             }    
             
        double **Tau = (double**)malloc(*nbInd*sizeof(double*));
            for(i=0;i<*nbInd;i++){ 
                   Tau[i] = (double*)malloc(*nbClass*sizeof(double));
  	    }
  
  double **G = (double**)malloc(*nbInd*sizeof(double*));
            for(i=0;i<*nbInd;i++){ 
                   G[i] = (double*)malloc(*nbClass*sizeof(double));
            }
  
      VectToMat(PiVect,Pi, *nbClass,*nbClass);
      VectToMat(FVect,F, *nbClass,*nbInd);    
      VectToMat(TauVect,Tau, *nbClass,*nbInd); 
      VectToMat(GVect,G, *nbClass,*nbInd); 
     
    //Premier instant de l'étape Backward (instant t); 
    CalculTauN(F,*nbInd,*nbClass,Tau); 
       
    //Calcul des autres termes(de l'instant 1 à t-1);           
    CalculTau(Pi,F,*nbClass,*nbInd,G,Tau);
   
    MatToVect(F,FVect,*nbClass,*nbInd);
    MatToVect(Pi,PiVect,*nbClass,*nbClass);
    MatToVect(Tau,TauVect,*nbClass,*nbInd);
    MatToVect(G,GVect,*nbClass,*nbInd);


for(i=0;i<*nbInd;i++){
                  free(F[i]) ;
           }
   free(F);
for(i=0;i<*nbInd;i++){
                  free(Tau[i]) ;
           }
   free(Tau);
for(i=0;i<*nbInd;i++){
                  free(G[i]) ;
           }
   free(G);
for(i=0;i<*nbClass;i++){
                  free(Pi[i]) ;
           }
   free(Pi); 
}
