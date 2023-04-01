#include "fem.h"
#define ___ 0



double geoSize(double x, double y, double *data){
    double x0 = data[0], y0 = data[1], r0 = data[2], lc0=data[3]; // data from disk0
    double x1 = data[4], y1 = data[5], r1 = data[6], lc1=data[7]; // data from disk1
    double lc_global = data[8];
    // Your contribution starts here





    return lc_global;
    // Your contribution ends here
}

void geoGenerate(femGeo *geometry, double x0, double y0, double r0, double lc0, double x1, double y1, double r1, double lc1){
    int ierr;
    double w = geometry->width;
    double h = geometry->height;
    double meshSize = geometry->meshSize;
    // Use OPENCASCADE to describe the geometrical entities
    int rect_id = gmshModelOccAddRectangle(___, ___, ___, ___, ___, ___, ___, &ierr); chk(ierr);
    int disk_id[2];
    disk_id[0] = gmshModelOccAddDisk(___, ___, ___, ___, ___, ___, NULL, 0, NULL, 0, &ierr); chk(ierr);
    disk_id[1] = gmshModelOccAddDisk(___, ___, ___, ___, ___, ___, NULL, 0, NULL, 0, &ierr); chk(ierr);
    int object[] = {___, ___};
    int tool0[] = {___, ___};
    gmshModelOccCut(___,___,___,___,
                   NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr); chk(ierr);
    int tool1[] = {___, ___};
    gmshModelOccCut(___,___,___,___,
                   NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr); chk(ierr);
    // pass the OPENCASCADE data to GMSH
    gmshModelOccSynchronize(&ierr); chk(ierr);
    // Create physical group (used for boundary conditions)
    int dim = 2;
    int ntags = 1;
    int tags[] = {1};
    gmshModelAddPhysicalGroup(dim, tags, ntags, 1, "surface", &ierr); chk(ierr);
    int bnd_dim = dim-1;
    tags[0] = 1; gmshModelAddPhysicalGroup(bnd_dim, tags, ntags, 1, "outer_disk", &ierr); chk(ierr);
    tags[0] = 2; gmshModelAddPhysicalGroup(bnd_dim, tags, ntags, 2, "bottom", &ierr); chk(ierr);
    tags[0] = 3; gmshModelAddPhysicalGroup(bnd_dim, tags, ntags, 3, "left", &ierr); chk(ierr);
    tags[0] = 4; gmshModelAddPhysicalGroup(bnd_dim, tags, ntags, 4, "right", &ierr); chk(ierr);
    tags[0] = 5; gmshModelAddPhysicalGroup(bnd_dim, tags, ntags, 5, "top", &ierr); chk(ierr);
    tags[0] = 6; gmshModelAddPhysicalGroup(bnd_dim, tags, ntags, 6, "inner_disk", &ierr); chk(ierr);
    double data[] = {x0, y0, r0, lc0, x1, y1, r1, lc1, meshSize};
    geoSetSizeCallback(data);
    gmshModelMeshGenerate(2, &ierr);  
}