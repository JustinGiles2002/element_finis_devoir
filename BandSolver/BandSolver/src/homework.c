
#include"fem.h"


#ifndef NORENUMBER 
double *X_or_Y ; 
int compare_node(const void *a,const void *b ){
    int node_a = *((int*)a) ; 
    int node_b = *((int*)b) ;
    return X_or_Y[node_b] - X_or_Y[node_a] ;  

}
void femMeshRenumber(femMesh *theMesh, femRenumType renumType)
{
    int node[theMesh->nNode] ; 
    for (size_t i = 0; i < theMesh->nNode; i++)
        node[i] = i ; 
    switch (renumType) {
        case FEM_NO :
            break;
        case FEM_XNUM : 
            X_or_Y = theMesh ->X ;  
            qsort(node,theMesh->nNode,sizeof(int),compare_node) ; 
            break;
        case FEM_YNUM : 
            X_or_Y = theMesh ->Y ;  
            qsort(node,theMesh->nNode,sizeof(int),compare_node) ;  // regarde ou sont les plus grand_noeud pour après pouvoir les remplacé 
            break;            

        default : Error("Unexpected renumbering option"); 
        }
    for (size_t i = 0; i < theMesh->nNode; i++)
        theMesh->number[node[i]] = i ;  // remplace par les premier noeud 
}

#endif
#ifndef NOBAND 

int femMeshComputeBand(femMesh *theMesh)
{
    int max,min,band = 0 ; 
    int myBand = theMesh->nNode;
    int nLcolal = theMesh->nLocalNode ; 
    int list_element[nLcolal] ; 
    for (size_t i = 0; i < theMesh->nElem; i++)
    {
        int *elem = &(theMesh->elem[nLcolal*i]) ;
        for(size_t j = 0 ; j < nLcolal ; j++)
            list_element[j] = theMesh ->number[elem[j]] ;  // obliger de passer par la car, les élément sont retriée. 
        
        max = list_element[0] ; 
        min =  list_element[0]  ; 
        for(size_t j = 1; j < nLcolal ; j++){
            max = (int) fmax((double) max, (double) list_element[j]) ; 
            min = (int) fmin((double) min, (double) list_element[j]) ; 
        }
        band = (band < max - min) ?  max - min  : band   ; 
    }
    
    return(++myBand);
}


#endif
#ifndef NOBANDASSEMBLE


void femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc)  // asseemblage traingulaire supérieur 
{
    for (size_t i = 0; i < nLoc; i++){
        int row = map[i] ;
        myBandSystem->B[row] += Bloc[i] ; 
        for(size_t j = 0 ; j<nLoc ; j++) {
            int col = map[j] ; 
            if(col>= row) 
                myBandSystem->A[row][col] += Aloc[i*nLoc + j] ; 
        }
    }
}


#endif
#ifndef NOBANDELIMINATE


double  *femBandSystemEliminate(femBandSystem *myBand) // par implimentation des matrices denses 
{
    double  **A, *B, factor;
    int i, j, k, jend, size, band;
    A = myBand->A;
    B = myBand->B; 
    size = myBand->size;
    band = myBand->band;

    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-8 ) {
            printf("Pivot index %d  ",k);
            printf("Pivot value %e  ",A[k][k]);
            Error("Cannot eliminate with such a pivot"); 
        }
        // ici matrice band 
        jend = (k + band < size) ? k + band : size ; 
        for (i = k+1 ; i <  jend; i++) {
            factor = A[k][i] / A[k][k];     // changement car  a[k][i] == a[i][k] et ce a[i][k] n'est pas décalré caron n'a pris car la triangulaire supérieur 
            for (j = k+1 ; j < jend; j++) 
                A[i][j] -= A[k][j] * factor;
            B[i] -= B[k] * factor; }}
    /* Back-substitution */
    
    for (i = size-1; i >= 0 ; i--) {
        factor = 0;
        jend = (k+band<size) ? i + band : size ; 
        
        for (j = i+1 ; j < jend; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }
    return(myBand->B);
}


#endif

