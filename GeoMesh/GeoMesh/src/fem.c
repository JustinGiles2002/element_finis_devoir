/*
 *  fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"

femGeo *theGeometry = NULL;

static int cmp_size_t(const void *a, const void *b){
    return (int) (*(const size_t*) a - * (const size_t*) b);
}

static int searchNode(size_t node, size_t *nodeTags, size_t nNodes){
    size_t *index = (size_t*) bsearch(&node, nodeTags, nNodes, sizeof(size_t), cmp_size_t);
    if (index == NULL){
        Error("Node not found\n");
    }
    return (int) (index - nodeTags);
}


static double gmshSizeCallback(int dim, int tag, double x, double y, double z, double lc, void *data){
    return geoSize(x,y,(double*)data);    
}

double geoSetSizeCallback(double *data){
    int ierr;
    gmshModelMeshSetSizeCallback(gmshSizeCallback, (void*)data, &ierr);
}

static inline int get_element_type_from_nLocalNode(const int nLocalNode){
    if (nLocalNode == 2) return 1; // edge
    if (nLocalNode == 3) return 2; // triangle
    if (nLocalNode == 4) return 3; // quadrangle
    if (nLocalNode == 1) return 15;// node
}

femGeo* geoCreate(double width, double height, double meshSize){
    int ierr;
    if (theGeometry != NULL) geoFree(theGeometry);
    theGeometry = malloc(sizeof(femGeo));
    gmshInitialize(0, NULL, 1, 0, &ierr);
    gmshModelAdd("MyGeometry", &ierr); chk(ierr);
    theGeometry->width = width;
    theGeometry->height = height;
    theGeometry->meshSize = meshSize;
    return theGeometry;
}

void geoFree(femGeo *geometry){
    int ierr;
    free(geometry);
    gmshFinalize(&ierr); chk(ierr);
}


femNodes* femNodesCreate(femGeo *geometry){
    femNodes *theNodes = malloc(sizeof(femNodes));
    int ierr;
    size_t nNodeTags, nCoord, nParametricCoord;
    size_t *nodeTags;
    double *coord, *parametricCoord;
    gmshModelMeshRenumberNodes(&ierr);
    gmshModelMeshGetNodes(&nodeTags, &nNodeTags,
                          &coord, &nCoord,
                          &parametricCoord, &nParametricCoord,
                          -1, -1, 0, 0, &ierr);
    theNodes->nNodes = (int) nNodeTags;
    theNodes->X = malloc(sizeof(double)*(theNodes->nNodes));
    theNodes->Y = malloc(sizeof(double)*(theNodes->nNodes));
    for(int i = 0; i < theNodes->nNodes; ++i){
        int tag = nodeTags[i]-1;
        theNodes->X[i] = coord[3*tag+0];
        theNodes->Y[i] = coord[3*tag+1];
    }
    gmshFree(nodeTags);
    gmshFree(coord);
    gmshFree(parametricCoord);
    return theNodes;
}

void femNodesFree(femNodes *theNodes){
    free(theNodes->X);
    free(theNodes->Y);
    free(theNodes);
}


femMesh* femMeshCreate(femGeo *geometry, femNodes *theNodes, int elementType) {
    femMesh *theMesh = malloc(sizeof(femMesh));  
    theMesh->nodes = theNodes;
    if (elementType == 1)       {theMesh->nLocalNode = 2; strncpy(theMesh->name, "mesh_edge",       256);}
    else if(elementType == 2)   {theMesh->nLocalNode = 3; strncpy(theMesh->name, "mesh_triangle",   256);}
    else if(elementType == 3)   {theMesh->nLocalNode = 4; strncpy(theMesh->name, "mesh_quad",       256);}
    else if(elementType == 15)  {theMesh->nLocalNode = 1; strncpy(theMesh->name, "mesh_node",       256);}
    else                        {Error("Unsupported elementType");}
    int ierr;
    size_t nNodeTags, nElement;
    size_t *elementTags, *nodeTags;
    gmshModelMeshGetElementsByType(elementType,
                                   &elementTags, &nElement,
                                   &nodeTags, &nNodeTags,
                                   -1, 0, 1, &ierr);
    theMesh->nElem = (int) nElement;
    theMesh->elem = malloc(sizeof(int)*theMesh->nElem*theMesh->nLocalNode);
    for(int i = 0; i < theMesh->nElem; ++i){
        for(int j = 0; j < theMesh->nLocalNode; j++){
            int index = theMesh->nLocalNode*i+j;
            theMesh->elem[index] = nodeTags[index]-1; 
        }
    }
    gmshFree(nodeTags);
    gmshFree(elementTags);
    return theMesh;
}

void femMeshFree(femMesh *theMesh) {
    free(theMesh->elem);
    free(theMesh);
}

static void femMeshExport(const femMesh *theMesh, FILE *file){
    fprintf(file, "Mesh name : %20s \n", theMesh->name);
    fprintf(file, "Number of nodes %6d \n", theMesh->nodes->nNodes);
    for (int i = 0; i < theMesh->nodes->nNodes; ++i)
        fprintf(file,"%6d : %14.7e %14.7e \n",i,theMesh->nodes->X[i],theMesh->nodes->Y[i]);
    if (theMesh->nLocalNode == 4) fprintf(file, "Number of quads %6d \n", theMesh->nElem);  
    if (theMesh->nLocalNode == 3) fprintf(file, "Number of triangles %6d \n", theMesh->nElem);  
    if (theMesh->nLocalNode == 2) fprintf(file, "Number of edges %6d \n", theMesh->nElem);  
    if (theMesh->nLocalNode == 1) fprintf(file, "Number of nodes %6d \n", theMesh->nElem);  
    for(int i=0; i < theMesh->nElem; ++i){
        int *elem = &(theMesh->elem[i*theMesh->nLocalNode]);
        fprintf(file, "%6d : ", i);
        for(int j=0; j < theMesh->nLocalNode; ++j)
            fprintf(file, "%6d ", elem[j]);
        fprintf(file, "\n");
    }
}

static void femMeshImport(femMesh *theMesh, FILE *file){
    char str[256];
    int trash;
    femNodes *theNodes = malloc(sizeof(femNodes));
    theMesh->nodes = theNodes;
    ErrorScan(fscanf(file, "Mesh name : %20s \n", str));
    strncpy(theMesh->name, str, 256);
    ErrorScan(fscanf(file, "Number of nodes %6d \n", &theMesh->nodes->nNodes));
    theMesh->nodes->X = malloc(sizeof(double)*theMesh->nodes->nNodes);
    theMesh->nodes->Y = malloc(sizeof(double)*theMesh->nodes->nNodes);
    for (int i = 0; i < theMesh->nodes->nNodes; ++i) {
        ErrorScan(fscanf(file,"%d : %le %le \n",&trash,&theMesh->nodes->X[i],&theMesh->nodes->Y[i])); 
    }
    if (fgets(str, sizeof(str), file) == NULL) Error("Corrupted mesh file !");
    if (!strncmp(str,"Number of quads",15))     {ErrorScan(sscanf(str,"Number of quads %6d \n", &theMesh->nElem)); theMesh->nLocalNode = 4;} 
    if (!strncmp(str,"Number of triangles",19)) {ErrorScan(sscanf(str,"Number of triangles %6d \n", &theMesh->nElem)); theMesh->nLocalNode = 3;} 
    if (!strncmp(str,"Number of edges",15))     {ErrorScan(sscanf(str,"Number of edges %6d \n", &theMesh->nElem)); theMesh->nLocalNode = 2;}
    if (!strncmp(str,"Number of nodes",15))     {ErrorScan(sscanf(str,"Number of nodes %6d \n", &theMesh->nElem)); theMesh->nLocalNode = 1;} 
    theMesh->elem = malloc(sizeof(int)*theMesh->nLocalNode*theMesh->nElem);
    for(int i=0; i < theMesh->nElem; ++i){
        int *elem = &(theMesh->elem[i*theMesh->nLocalNode]);
        ErrorScan(fscanf(file, "%6d : ", &trash));
        for(int j=0; j < theMesh->nLocalNode; ++j)
            ErrorScan(fscanf(file, "%6d ", &elem[j]));
        ErrorScan(fscanf(file, "\n"));
    }
}

void femMeshWrite(const femMesh *theMesh, const char *filename) {
    int i,*elem;
    FILE* file = fopen(filename,"w");
    femMeshExport(theMesh, file);
    fclose(file);
}

void femMeshAppend(const femMesh *theMesh, const char *filename) {
    int i,*elem;
    FILE* file = fopen(filename,"a");
    femMeshExport(theMesh, file);
    fclose(file);
}

femMesh *femMeshRead(const char *filename) {
    femMesh *theMesh = malloc(sizeof(femMesh));
    FILE* file = fopen(filename,"r");
    femMeshImport(theMesh, file);
    fclose(file);
    return theMesh;
}


femMeshes *femMeshesCreate(femGeo *geometry, femNodes *theNodes){
    femMeshes *theMeshes = malloc(sizeof(femMeshes));
    int ierr;
    int *elementTypes;
    size_t nElementTypes;
    gmshModelMeshGetElementTypes(&elementTypes, &nElementTypes, -1,-1, &ierr); chk(ierr);
    theMeshes->nMesh = nElementTypes;
    theMeshes->mesh = malloc(sizeof(femMesh*)*theMeshes->nMesh);
    for(size_t i=0; i < nElementTypes; ++i){
        theMeshes->mesh[i] = femMeshCreate(geometry, theNodes, elementTypes[i]);
    }
    gmshFree(elementTypes);
    return theMeshes;
}

void femMeshesFree(femMeshes *theMeshes){
    for(int i=0; i < theMeshes->nMesh; ++i){
        femMeshFree(theMeshes->mesh[i]);
    }
    free(theMeshes->mesh);
    free(theMeshes);
}

void femMeshesWrite(const femMeshes *theMeshes, const char *filename) {
    FILE* file = fopen(filename, "w");
    fprintf(file, "Number of meshes : %6d \n", theMeshes->nMesh);
    fclose(file);
    for(int i=0; i < theMeshes->nMesh; ++i){
        femMeshAppend(theMeshes->mesh[i], filename);
    }
}

femMeshes *femMeshesRead(const char *filename) {
    FILE* file = fopen(filename,"r");
    if (file == NULL) {
        printf( "Cannot open file %s\n", filename);
        exit(0);
    }
    femMeshes *theMeshes = malloc(sizeof(femMeshes));
    ErrorScan(fscanf(file, "Number of meshes : %6d \n", &theMeshes->nMesh));
    theMeshes->mesh = malloc(sizeof(femMesh*)*theMeshes->nMesh);
    int trash;
    char str[256]; 
    for(int imesh=0; imesh < theMeshes->nMesh; ++imesh){
        theMeshes->mesh[imesh] = malloc(sizeof(femMesh));
        femMeshImport(theMeshes->mesh[imesh], file);
    }
    return theMeshes;
}


femDomain* femDomainCreate(femGeo *geometry, femMesh *mesh, const char *name) {
    femDomain *theDomain = malloc(sizeof(femDomain)); 
    int ierr;
    size_t nDimTags;
    int *dimTags;
    gmshModelGetPhysicalGroups(&dimTags, &nDimTags, -1, &ierr);
    for(size_t i=0; i < nDimTags/2; ++i){
        int dim = dimTags[2*i+0];
        int tag = dimTags[2*i+1];
        char *physical_name;
        gmshModelGetPhysicalName(dim, tag, &physical_name, &ierr);
        if (strncasecmp(name, physical_name, 256)==0){
            size_t nMeshNodeTags, nMeshElement;
            size_t *meshElementTags, *meshNodeTags;
            gmshModelMeshGetElementsByType(get_element_type_from_nLocalNode(mesh->nLocalNode),
                                &meshElementTags, &nMeshElement,
                                &meshNodeTags, &nMeshNodeTags,
                                -1, 0, 1, &ierr);
            int *entityTags;
            size_t nEntityTags;
            gmshModelGetEntitiesForPhysicalGroup(dim, tag, &entityTags, &nEntityTags, &ierr);
            for(int i=0; i < nEntityTags; ++i){
                int *elementType;
                size_t nElementType, **elementTags, *nElementTags, nnElementTags, **nodesTags, *nNodesTags, nnNodesTags; 
                gmshModelMeshGetElements(&elementType, &nElementType, &elementTags, &nElementTags, &nnElementTags, &nodesTags, &nNodesTags, &nnNodesTags, dim, entityTags[i], &ierr);
                if(nElementType > 1){
                    Error("Mixed Element not supported\n");
                }
                for(int ei=0; ei < nElementType; ++ei){
                    int nLocalNode, edim, eorder, n_primal;
                    double *xi;
                    size_t n_xi;
                    char *ename;
                    gmshModelMeshGetElementProperties(elementType[ei], &ename, &edim, &eorder, &nLocalNode, &xi, &n_xi, &n_primal, &ierr); chk(ierr);
                    if(nLocalNode != mesh->nLocalNode){
                        Error("Invalid Element Type\n");
                    }
                    gmshFree(xi);
                    gmshFree(ename);

                }
                strncpy(theDomain->name, name, 256);
                theDomain->mesh = mesh;
                theDomain->nElem = (int) nElementTags[0];
                theDomain->elem = malloc(sizeof(int)*theDomain->nElem);
                for(int ie=0; ie < theDomain->nElem; ++ie){
                    size_t *index = (size_t*) bsearch(&elementTags[0][ie], meshElementTags, nMeshElement, sizeof(size_t), cmp_size_t);
                    if(index == NULL){
                        Error("Index not found\n");
                    }
                    theDomain->elem[ie] = (int)((size_t)(index - meshElementTags));
                }
                for(int ie=0; ie < nnElementTags; ++ie){
                    gmshFree(elementTags[ie]);
                }
                for(int in=0; in < nnNodesTags; ++in){
                    gmshFree(nodesTags[in]);
                }
                gmshFree(meshElementTags);
                gmshFree(meshNodeTags);
                gmshFree(nElementTags);
                gmshFree(nNodesTags);
                gmshFree(elementTags);
                gmshFree(nodesTags);
                gmshFree(elementType);
            }
            gmshFree(entityTags);
        }
        gmshFree(physical_name);
    }
    gmshFree(dimTags);
    return theDomain;
}

void femDomainFree(femDomain *theDomain){
    free(theDomain->elem);
    free(theDomain);
}

void femDomainPrint(femDomain *theDomain){
    printf("Name : %s\n", theDomain->name);
    printf("mesh : %s\n", theDomain->mesh->name);
    printf("Number of elements : %d \n", theDomain->nElem);
    for(int i=0; i < theDomain->nElem; ++i){
        int iel = theDomain->elem[i];
        printf("%6d : %6d : ", i, iel);
        for(int j=0; j < theDomain->mesh->nLocalNode; ++j){
            printf(" %6d", theDomain->mesh->elem[iel]+j);
        }
        printf("\n");
    }
}

void femDomainWrite(femDomain *theDomain, const char *filename){
    int i,*elem;
    FILE* file = fopen(filename,"w");
    fprintf(file, "Domain Name : %20s\n", theDomain->name);
    fprintf(file, "Number of elements : %6d\n", theDomain->nElem);
    for(int i=0; i < theDomain->nElem; ++i){
        int iel = theDomain->elem[i];
        fprintf(file, "%6d : %6d : ", i, iel);
        for(int j=0; j < theDomain->mesh->nLocalNode; ++j){
            fprintf(file," %6d", theDomain->mesh->elem[iel]+j);
        }
        fprintf(file, "\n");
    }
    fclose(file);
    femMeshAppend(theDomain->mesh, filename);
    fclose(file);
}

femDomain* femDomainRead(const char *filename){
    // TO DO
    return NULL;
}



femDomains* femDomainsCreate(femGeo *geometry, femMesh *mesh){
    femDomains *theDomains = malloc(sizeof(femDomains));
    int ierr;
    size_t nDimTags;
    int *dimTags;
    gmshModelGetPhysicalGroups(&dimTags, &nDimTags, -1, &ierr);
    int nDomain = 0;
    for(int i=0; i < nDimTags/2; ++i){
        int dim = dimTags[2*i+0];
        if((mesh->nLocalNode == 1 && dim==0)
        || (mesh->nLocalNode == 2 && dim==1)
        || (mesh->nLocalNode >= 3 && dim==2)){
            nDomain++;
        }
    }
    theDomains->nDomain = nDomain;
    theDomains->domain = malloc(sizeof(femDomain*)*nDomain);
    for(size_t i=0; i < nDimTags/2; ++i){
        int dim = dimTags[2*i+0];
        int tag = dimTags[2*i+1];
        if((mesh->nLocalNode == 2 && dim==1)
        || (mesh->nLocalNode >= 3 && dim==2)){
            char *physical_name;
            gmshModelGetPhysicalName(dim, tag, &physical_name, &ierr);
            theDomains->domain[i] = femDomainCreate(geometry, mesh, physical_name);
            gmshFree(physical_name);
        }
    }
    gmshFree(dimTags);
    return theDomains;
}

void femDomainsFree(femDomains *theDomains){
    for(int i=0; i<theDomains->nDomain; ++i){
        femDomainFree(theDomains->domain[i]);
    }
    free(theDomains->domain);
    free(theDomains);
}


double femMin(double *x, int n) 
{
    double myMin = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMin = fmin(myMin,x[i]);
    return myMin;
}

double femMax(double *x, int n) 
{
    double myMax = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMax = fmax(myMax,x[i]);
    return myMax;
}

void femError(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);                                                 
}

void femErrorScan(int test, int line, char *file)                                  
{ 
    if (test >= 0)  return;
    
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in fscanf or fgets in %s at line %d : \n", file, line);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");   
    exit(69);                                       
}

void femWarning(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Warning in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");                                              
}