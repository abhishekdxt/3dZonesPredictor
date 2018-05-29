/**
 ========================================================================================
  Predict 3d zones centered around a residue in protein of a given structure
  Created 6 Nov 2010
  Abhishek Dixit
 ========================================================================================
 **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct {
	char atom[4], resid[4], chain;
	int  atomSeq, residSeq;
	double coordX, coordY, coordZ;
} ATOM;

int main (int argc, char *argv[]) {
	if( argc != 5 ) {  // Check argument length
		fprintf(stderr, "\n==== Predicts 3d zones centered around a residue in protein of a given structure ====\nUsage: ./3dZonesPredictor <PDB_file_path> <chain> <residue_seq> <threshold>\n\n");
		exit(1);
	}
	else  {
		char chain = *argv[2];
		int resid = atoi(argv[3]);
		double threshold = atof(argv[4]);
		FILE *inFile = NULL;
		inFile = fopen(argv[1], "r");

		if( inFile == NULL ) { // PDB file not found
			fprintf(stderr, "%s file not found\n", argv[1]);
			exit(1);
		}

		char line[80];
		ATOM *atom = NULL;
		int atomCount = 0, residCount = 0, atomsOut = 0, headerOut = 0, *residIndx = NULL;

		while( fgets(line, sizeof(line), inFile) != NULL ) {
			if( line[0] == 'A' && line[1] == 'T' && line[2] == 'O' && line[3] == 'M' && line[4] == ' ' ) {
				atom = (ATOM *) realloc (atom, (atomCount + 1) * sizeof(ATOM));
				sscanf(line, "ATOM   %d  %s %s %c %d\t%lf %lf %lf", &atom[atomCount].atomSeq, atom[atomCount].atom, atom[atomCount].resid, &atom[atomCount].chain, &atom[atomCount].residSeq, &atom[atomCount].coordX, &atom[atomCount].coordY, &atom[atomCount].coordZ);
				if( atom[atomCount].chain == chain && atom[atomCount].residSeq == resid ) {
					residIndx = (int *) realloc (residIndx, (residCount + 1) * sizeof(int));
					residIndx[residCount++] = atomCount;
				}
				atomCount++;
			}
		}

		fclose(inFile);

		for( int i = 0; i < residCount; i++ ) {
			for( int j = 0; j < atomCount; j++)       {
				if( residIndx[0] > j || residIndx[residCount-1] < j ) {
					double dist = sqrt(pow(atom[residIndx[i]].coordX - atom[j].coordX, 2) + pow(atom[residIndx[i]].coordY - atom[j].coordY, 2) + pow(atom[residIndx[i]].coordZ - atom[j].coordZ, 2));
					if( dist <= threshold ) {
						if ( headerOut == 0 ) {
							printf("%s", "SrcAtomSeq\tSrcAtom\tSrcResidue\tSrcChain\tSrcResidueSeq\tAtomSeq\tAtom\tResid\tChain\tResidueSeq\tDistance\n" );
							headerOut = 1;
						}
						printf("%d\t%s\t%s\t%c\t%d\t%d\t%4s\t%4s\t%c\t%d\t%.3lf\n", atom[residIndx[i]].atomSeq, atom[residIndx[i]].atom, atom[residIndx[i]].resid, atom[residIndx[i]].chain, atom[residIndx[i]].residSeq, atom[j].atomSeq, atom[j].atom, atom[j].resid, atom[j].chain, atom[j].residSeq, dist);
						atomsOut++;
					}
				}
			}
		}

		free(atom);
		free(residIndx);

		if ( atomsOut == 0 ) {
			fprintf(stderr, "No atoms found within specified radius\n");
		}

		return 0;
	}
}
