
/*
 *  fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#ifndef _FEM_H_
#define _FEM_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "gmshc.h"


#define ErrorScan(a)   femErrorScan(a,__LINE__,__FILE__)
#define Error(a)       femError(a,__LINE__,__FILE__)
#define Warning(a)     femWarning(a,  __LINE__, __FILE__)
#define FALSE 0 
#define TRUE  1
#define chk(ierr)                                               \
  if(ierr != 0){                                                \
    fprintf(stderr, "Error on line %i in function '%s': "       \
            "gmsh function returned non-zero error code: %i\n", \
            __LINE__, __FUNCTION__, ierr);                      \
    gmshFinalize(NULL);                                         \
    exit(0 /* ierr  --- for ctest */);                          \
  }

typedef enum {FEM_TRI,FEM_QUAD} femElementType;

typedef struct {
    int nNodes;
    double *X;
    double *Y;
} femNodes;

typedef struct {
    int nLocalNode;
    int nElem;
    int *elem;
    femNodes *nodes;
    char name[256];
} femMesh;

typedef struct {
    int nMesh;
    femMesh **mesh;
} femMeshes;

typedef struct {
    femMesh *mesh;
    int nElem;
    int *elem;
    char name[256];
} femDomain;

typedef struct {
    int nDomain;
    femDomain **domain;
} femDomains;

typedef struct {
    double width;
    double height;
    double meshSize;
    double x0;
} femGeo;

extern femGeo *theGeometry;

femGeo*             geoCreate(double width, double height, double meshSize);
double              geoSize(double x, double y, double *data);
double              geoSetSizeCallback(double *data);
void                geoFree(femGeo *geometry);

femNodes*           femNodesCreate(femGeo *geometry);
void                femNodesFree(femNodes *theNodes);

femMesh*            femMeshCreate(femGeo *geometry, femNodes *theNodes, int element_type);
void                femMeshFree(femMesh *theMesh);
femMesh*            femMeshRead(const char *filename);
void                femMeshWrite(const femMesh *theMesh, const char *filename);
void                femMeshAppend(const femMesh *theMesh, const char *filename);


femMeshes*          femMeshesCreate(femGeo *geometry, femNodes *theNodes);
void                femMeshesFree(femMeshes *theMeshes);
void                femMeshesWrite(const femMeshes *theMeshes, const char *filename);
femMeshes*          femMeshesRead(const char *filename);


femDomain*          femDomainCreate(femGeo *geometry, femMesh *mesh, const char *name);
void                femDomainFree(femDomain *theDomain);
void                femDomainPrint(femDomain *theDomain);
void                femDomainWrite(femDomain *theDomain, const char *filename);
femDomain*          femDomainRead(const char *filename);

femDomains*         femDomainsCreate(femGeo *geometry, femMesh *mesh);
void                femDomainsFree(femDomains *theDomains);
void                femDomainsPrint(femDomains *theDomains);


double              femMin(double *x, int n);
double              femMax(double *x, int n);
void                femError(char *text, int line, char *file);
void                femErrorScan(int test, int line, char *file);
void                femWarning(char *text, int line, char *file);


#endif
