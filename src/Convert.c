/******************************************************************************/
/*       Convertir un vecteur en matrice (remplie par colonne)                */
/******************************************************************************/


double VectToMat(double *Vect,double **Mat, int nbcol,int nbrow)  
{
   //Déclaration et initialisation des variables locales;      
   int x,y;
  

   
   for(x=0; x<nbrow; x++){ 
            for(y=0; y<nbcol; y++){
                     
             Mat[x][y] = Vect[x+y*nbrow];

   }
}

} 



/******************************************************************************/
/*       Convertir une matrice en vecteur (par colonne de la matrice)         */
/******************************************************************************/

double MatToVect(double **Mat,double *Vect, int nbcol,int nbrow)  
{
   //Déclaration et initialisation des variables locales;      
   int x,y;
   

   
   for(x=0; x<nbrow; x++){ 
            for(y=0; y<nbcol; y++){
                     
             Vect[x+y*nbrow]=Mat[x][y];
   }
}

} 
