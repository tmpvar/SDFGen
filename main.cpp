//SDFGen - A simple grid-based signed distance field (level set) generator for triangle meshes.
//Written by Christopher Batty (christopherbatty@yahoo.com, www.cs.columbia.edu/~batty)
//...primarily using code from Robert Bridson's website (www.cs.ubc.ca/~rbridson)
//This code is public domain. Feel free to mess with it, let me know if you like it.

#include "makelevelset3.h"
//#include "config.h"

#ifdef HAVE_VTK
  #include <vtkImageData.h>
  #include <vtkFloatArray.h>
  #include <vtkXMLImageDataWriter.h>
  #include <vtkPointData.h>
  #include <vtkSmartPointer.h>
#endif


#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>

unsigned int align_up(unsigned int v, unsigned int align) {
    return ((v + align - 1) / align) * align;
}

int main(int argc, char* argv[]) {
  
  if(argc != 7) {
    std::cout << "SDFGen - A utility for converting closed oriented triangle meshes into grid-based signed distance fields.\n";
    std::cout << "\nThe output file format is:";
    std::cout << "<ni> <nj> <nk>\n";
    std::cout << "<origin_x> <origin_y> <origin_z>\n";
    std::cout << "<dx>\n";
    std::cout << "<value_1> <value_2> <value_3> [...]\n\n";
    
    std::cout << "(ni,nj,nk) are the integer dimensions of the resulting distance field.\n";
    std::cout << "(origin_x,origin_y,origin_z) is the 3D position of the grid origin.\n";
    std::cout << "<dx> is the grid spacing.\n\n";
    std::cout << "<value_n> are the signed distance data values, in ascending order of i, then j, then k.\n";

    std::cout << "The output filename will match that of the input, with the OBJ suffix replaced with SDF.\n\n";

    std::cout << "Usage: SDFGen <filename> <dx> <padmin> <padmax> <align> <padmax_postalign>\n\n";
    std::cout << "Where:\n";
    std::cout << "\t<filename> specifies a Wavefront OBJ (text) file representing a *triangle* mesh (no quad or poly meshes allowed). File must use the suffix \".obj\".\n";
    std::cout << "\t<dx> specifies the length of grid cell in the resulting distance field.\n";
    std::cout << "\t<padmin> specifies the number of cells worth of padding between the object bound box min and the boundary of the distance field grid. Minimum is 1.\n\n";
	std::cout << "\t<padmax> specifies the number of cells worth of padding between the object bound box max and the boundary of the distance field grid.\n\n";
	std::cout << "\t<align> specifies the alignment of the bounding box size (round up).\n\n";
	std::cout << "\t<padmax_postalign> specifies the number of cells worth of padding added to box max after align.\n\n";
    
    exit(-1);
  }

  std::string filename(argv[1]);
  if(filename.size() < 5 || filename.substr(filename.size()-4) != std::string(".obj")) {
    std::cerr << "Error: Expected OBJ file with filename of the form <name>.obj.\n";
    exit(-1);
  }

  std::stringstream arg2(argv[2]);
  float dx;
  arg2 >> dx;
  
  std::stringstream arg3(argv[3]);
  int padmin;
  arg3 >> padmin;
  if(padmin < 1) padmin = 1;

  std::stringstream arg4(argv[4]);
  int padmax;
  arg4 >> padmax;

  std::stringstream arg5(argv[5]);
  int align;
  arg5 >> align;

  std::stringstream arg6(argv[6]);
  int padmax_postalign;
  arg6 >> padmax_postalign;

  //start with a massive inside out bound box.
  Vec3f min_box(std::numeric_limits<float>::max(),std::numeric_limits<float>::max(),std::numeric_limits<float>::max()), 
    max_box(-std::numeric_limits<float>::max(),-std::numeric_limits<float>::max(),-std::numeric_limits<float>::max());
  
  std::cout << "Reading data.\n";

  std::ifstream infile(argv[1]);
  if(!infile) {
    std::cerr << "Failed to open. Terminating.\n";
    exit(-1);
  }

  int ignored_lines = 0;
  std::string line;
  std::vector<Vec3f> vertList;
  std::vector<Vec3ui> faceList;
  while(!infile.eof()) {
    std::getline(infile, line);

    //.obj files sometimes contain vertex normals indicated by "vn"
    if(line.substr(0,1) == std::string("v") && line.substr(0,2) != std::string("vn")){
      std::stringstream data(line);
      char c;
      Vec3f point;
      data >> c >> point[0] >> point[1] >> point[2];
      vertList.push_back(point);
      update_minmax(point, min_box, max_box);
    }
    else if(line.substr(0,1) == std::string("f")) {
      std::stringstream data(line);
      char c;
      int v0,v1,v2;
      data >> c >> v0 >> v1 >> v2;
      faceList.push_back(Vec3ui(v0-1,v1-1,v2-1));
    }
    else if( line.substr(0,2) == std::string("vn") ){
      std::cerr << "Obj-loader is not able to parse vertex normals, please strip them from the input file. \n";
      exit(-2); 
    }
    else {
      ++ignored_lines; 
    }
  }
  infile.close();
  
  if(ignored_lines > 0)
    std::cout << "Warning: " << ignored_lines << " lines were ignored since they did not contain faces or vertices.\n";

  std::cout << "Read in " << vertList.size() << " vertices and " << faceList.size() << " faces." << std::endl;

  //Add padding around the box.
  Vec3f unit(1,1,1);
  min_box -= padmin*dx*unit;
  max_box += padmax*dx*unit;
  Vec3ui sizes = Vec3ui((max_box - min_box)/dx);

  sizes[0] = align_up(sizes[0], align);
  sizes[1] = align_up(sizes[1], align);
  sizes[2] = align_up(sizes[2], align);

  sizes += Vec3ui(padmax_postalign);
  
  std::cout << "Bound box size: (" << min_box << ") to (" << max_box << ") with dimensions " << sizes << "." << std::endl;

  std::cout << "Computing signed distance field.\n";
  Array3f phi_grid;
  make_level_set3(faceList, vertList, min_box, dx, sizes[0], sizes[1], sizes[2], phi_grid);

  std::string outname;

  #ifdef HAVE_VTK
    // If compiled with VTK, we can directly output a volumetric image format instead
    //Very hackily strip off file suffix.
    outname = filename.substr(0, filename.size()-4) + std::string(".vti");
    std::cout << "Writing results to: " << outname << "\n";
    vtkSmartPointer<vtkImageData> output_volume = vtkSmartPointer<vtkImageData>::New();

    output_volume->SetDimensions(phi_grid.ni ,phi_grid.nj ,phi_grid.nk);
    output_volume->SetOrigin( phi_grid.ni*dx/2, phi_grid.nj*dx/2,phi_grid.nk*dx/2);
    output_volume->SetSpacing(dx,dx,dx);

    vtkSmartPointer<vtkFloatArray> distance = vtkSmartPointer<vtkFloatArray>::New();
    
    distance->SetNumberOfTuples(phi_grid.a.size());
    
    output_volume->GetPointData()->AddArray(distance);
    distance->SetName("Distance");

    for(unsigned int i = 0; i < phi_grid.a.size(); ++i) {
      distance->SetValue(i, phi_grid.a[i]);
    }

    vtkSmartPointer<vtkXMLImageDataWriter> writer =
    vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(outname.c_str());

    #if VTK_MAJOR_VERSION <= 5
      writer->SetInput(output_volume);
    #else
      writer->SetInputData(output_volume);
    #endif
    writer->Write();

  #else
    // if VTK support is missing, default back to the original ascii file-dump.
    //Very hackily strip off file suffix.
    outname = filename.substr(0, filename.size()-4) + std::string(".sdf");
    std::cout << "Writing results to: " << outname << "\n";

    // 16 bit binary output implemented by sebbbi
    Vec3f diagonal = Vec3f(
        dx * float(phi_grid.ni),
        dx * float(phi_grid.nj),
        dx * float(phi_grid.nk));

    float diagonal_length = mag(diagonal);

    std::vector<uint16_t> voxels_u16(phi_grid.a.size());
    for (unsigned long i = 0; i < phi_grid.a.size(); ++i)
    {
        float v = phi_grid.a[i];
        float v_01 = (v / diagonal_length) * 0.5f + 0.5f;
        voxels_u16[i] = uint16_t(v_01 * 65535.0f);
    }

    std::ofstream outfile( outname.c_str(), std::ofstream::binary);
    outfile.write((char*)&phi_grid.ni, sizeof(phi_grid.ni));
    outfile.write((char*)&phi_grid.nj, sizeof(phi_grid.nj));
    outfile.write((char*)&phi_grid.nk, sizeof(phi_grid.nk));
    outfile.write((char*)&min_box, sizeof(min_box));
    outfile.write((char*)&dx, sizeof(dx));
    outfile.write((char*)voxels_u16.data(), voxels_u16.size() * sizeof(uint16_t));
    outfile.close();
  #endif

  std::cout << "Processing complete.\n";

return 0;
}
