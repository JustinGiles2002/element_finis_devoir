/*
 *  main.c
 *  Library for EPL1110 : Finite Elements for dummies
 *  Integration of Baltic sea :-)
 *
 *  Copyright (C) 2022 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include <stdio.h>
#include <math.h>

#include "fem.h"
#include "glfem.h"
#define GRAPHICS 1

void geoGenerate(femGeo *geometry, 
                double x0, double y0, double r0, double lc0,
                double x1, double y1, double r1, double lc1);

int main(void)
{   
    double w = 1.0;
    double h = 2.0;
    double meshSize = 0.1;
    double x0 = -w/2, y0 = -h/2, r0 = w/2, lc0 = 0.05;
    double x1 =  w/4, y1 =  h/4, r1 = w/8, lc1 = 0.20;
    
    theGeometry = geoCreate(w, h, meshSize);
    geoGenerate(theGeometry, x0, y0, r0, lc0, x1, y1, r1, lc1);
    femNodes *theNodes = femNodesCreate(theGeometry);
    femMeshes* theMeshes = femMeshesCreate(theGeometry, theNodes);
    char filename[] = "../data/my_meshes.txt";
    femMeshesWrite(theMeshes, filename);
    geoFree(theGeometry);
    double data[] = {x0, y0, r0, lc0, x1, y1, r1, lc1, meshSize};
    double *meshSizeField = malloc(theNodes->nNodes*sizeof(double));
    for(int i=0; i < theNodes->nNodes; ++i){
        meshSizeField[i] = geoSize(theNodes->X[i], theNodes->Y[i], data);
    }
    
//
//  Partie graphique :-)
//  Peut etre commentee si vous n'arrivez pas a compiler la partie graphique
//  Il faut aussi retirer l'include de glfem.h et retirer la partie graphique de
//  de la compilation.  Il suffit alors de compiler main.c fem.c homework.c
//
//  So easy !
// 

#if GRAPHICS
    int mode = 1; // Change mode by pressing "j", "k", "l"
    int show_mesh=0; // Change show_mesh by pressing "i", "o"
    glfemWindowCreate("EPL1110 : GeoMesh",480,480);
    do {
        glfemChangeState(&mode, theMeshes->nMesh);
        glfemChangeBool(&show_mesh);
        int w,h;
        glfemReshape(theNodes); glfemSetLineWidth(0.002);
        glfemSetColor(GLFEM_BLACK);  glfemPlotField(theMeshes->mesh[mode], meshSizeField);
        if (show_mesh){
            glfemSetLineWidth(0.002); glfemSetColor(GLFEM_BLACK);   glfemPlotMesh(theMeshes->mesh[mode]);
        }
        glfemWindowUpdate();
    } while(!glfemWindowShouldClose());
    glfemWindowFree();
#endif
    free(meshSizeField);
    femMeshesFree(theMeshes);
    femNodesFree(theNodes);
    exit(EXIT_SUCCESS);
    return 0;  
}

 
