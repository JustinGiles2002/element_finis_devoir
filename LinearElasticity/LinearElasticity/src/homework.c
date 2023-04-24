#include "fem.h"




void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    
    int ierr;
    double r = w/4;
    int idRect = gmshModelOccAddRectangle(0.0,0.0,0.0,w,h,-1,0.0,&ierr); 
    int idDisk = gmshModelOccAddDisk(w/2.0,h/2.0,0.0,r,r,-1,NULL,0,NULL,0,&ierr); 
    int idSlit = gmshModelOccAddRectangle(w/2.0,h/2.0-r,0.0,w,2.0*r,-1,0.0,&ierr); 
    int rect[] = {2,idRect};
    int disk[] = {2,idDisk};
    int slit[] = {2,idSlit};

    gmshModelOccCut(rect,2,disk,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    gmshModelOccCut(rect,2,slit,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    gmshModelOccSynchronize(&ierr); 

    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",8,&ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  }
  
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }
 
    return;
}


#include "fem.h"

double *femElasticitySolve(femProblem *theProblem)
{
    // Récupération des données du problème
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femMesh        *theMesh = theProblem->geometry->theElements;
    femNodes       *theNodes = theProblem->geometry->theNodes;

    // Constantes physiques
    double a = theProblem->A;
    double b = theProblem->B;
    double c = theProblem->C;
    double rho = theProblem->rho;
    double g = theProblem->g;
    int nLocal = theMesh->nLocalNode;

    // Variables locales
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4];
    // Boucle sur les éléments
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j = 0; j < nLocal; ++j) {
            map[j] = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j]+1;
            x[j] = theNodes->X[map[j]];
            y[j] = theNodes->Y[map[j]]; 
        }

    for (iInteg = 0; iInteg < theRule->n; iInteg++) {
        double xsi = theRule->xsi[iInteg];
        double eta = theRule->eta[iInteg];
        double weight = theRule->weight[iInteg];

        femDiscretePhi2(theSpace, xsi, eta, phi);
        femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

        double dxdxsi = 0, dxdeta = 0, dydxsi = 0, dydeta = 0;

        for (i = 0; i < theSpace->n; i++) {
            dxdxsi += x[i] * dphidxsi[i];
            dxdeta += x[i] * dphideta[i];
            dydxsi += y[i] * dphidxsi[i];
            dydeta += y[i] * dphideta[i];
        }

        double jac = dxdxsi * dydeta - dxdeta * dydxsi;

        for (i = 0; i < theSpace->n; i++) {
            dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
            dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
        }

        for (i = 0; i < theSpace->n; i++) {
            for (j = 0; j < theSpace->n; j++) {
                theSystem->A[mapX[i]][mapX[j]] += (dphidx[i] * dphidx[j] * a
                                         + dphidy[i] * dphidy[j] * c) * jac * weight;
                theSystem->A[mapY[i]][mapX[j]] += (dphidy[i] * dphidx[j] * b
                                         + dphidx[i] * dphidy[j] * c) * jac * weight;
                theSystem->A[mapX[i]][mapY[j]] += (dphidx[i] * dphidy[j] * b
                                         + dphidy[i] * dphidx[j] * c) * jac * weight;
                theSystem->A[mapY[i]][mapY[j]] += (dphidy[i] * dphidy[j] * a
                                         + dphidx[i] * dphidx[j] * c) * jac * weight;
            }
        }

        for (i = 0; i < theSpace->n; i++) {
            theSystem->B[mapX[i]] += 0;
            theSystem->B[mapY[i]] += phi[i] * -1 * g * rho * jac * weight;
        }
    }
}

int *theConstrainedNodes = theProblem->constrainedNodes;

  
    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem,i,value); }}
                            
    return femFullSystemEliminate(theSystem);
}