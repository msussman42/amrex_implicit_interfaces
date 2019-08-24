#include <stdlib.h>
#include <stdio.h>
#include <string>

#include "ArrayT.h"
#include "WriteAmiraMesh.h"

using namespace std;

int main()
{
	string file = "visual";
	printf("%s\n", file.c_str());



	char ext[5];
	string filename;
	char line[256];
	int nx = 256;
	int ny = 256;
	int nz = 256;

	ArrayT<float> v1(nx, ny, nz);
	ArrayT<float> v2(nx, ny, nz);
	WriteAmiraMesh wm;

	int nmax = 2320;
	int n1 = 1800;
	int n2 = 2090;
	nmax = 1950;
	for (int i=n1; i <= n2; i+=10) {
		sprintf(ext, "%.5d", i);
		filename = file + ext;
		printf("filename= %s\n", filename.c_str());
		FILE* fd = fopen(filename.c_str(), "r");
		if (!fd) continue;
		fgets(line, 255, fd);
		fgets(line, 255, fd);
		fgets(line, 255, fd);

		int ii, jj, kk;
		int ix;
		float f, g;
		float* v1f = v1.getDataPtr();
		float* v2f = v2.getDataPtr();
		int nbPts = v1.getSize();

		/** */
		for (int ix=0; ix < nbPts; ix++) {
			if (ix % (256*256) == 0) printf("plane %d\n", ix/(256*256));
			fgets(line, 255, fd);
			sscanf(line, "%d %d %d %f %f", &ii, &jj, &kk, &f, &g);
			v1f[ix] = f;
			v2f[ix] = g;
			//printf("f,g= %f, %f\n", f, g);
		}

		// Write AmiraMesh format
		/*** ***/

		string level = "level";
		string VOF = "VOF";
		string am_filename = level + ext;
		am_filename += ".am";
		printf("filename= %s\n", am_filename.c_str());

		wm.init(nx, ny, nz, 1, am_filename.c_str());
		wm.writeData(v1f, "level");

		am_filename = VOF + ext + ".am";
		printf("filename= %s\n", am_filename.c_str());

		wm.init(nx, ny, nz, 1, am_filename.c_str());
		wm.writeData(v2f, "VOF");
	}

	exit(0);
}
//----------------------------------------------------------------------
