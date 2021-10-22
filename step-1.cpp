// Translate this file with
//
// g++ -O3 assignment-code.cpp -o assignment-code
//
// Run it with
//
// ./demo-code
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2018-2020 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <iomanip>

#include <cmath>

double t          = 0;
double tFinal     = 0;
double tPlot      = 0;
double tPlotDelta = 0;

int NumberOfBodies = 0;

double colDist;

/**
 * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
 * each pointer represents one molecule/particle/body.
 */
double** x;

/**
 * Equivalent to x storing the velocities.
 */
double** v;

/**
 * One mass entry per molecule/particle.
 */
double*  mass;

/**
 * Global time step size used.
 */
double   timeStepSize = 0.0;

/**
 * Maximum velocity of all particles.
 */
double   maxV;

/**
 * Minimum distance between two elements.
 */
double   minDx;

// Array for distances

double** distCubMat;

double* force0; // int *array{ new int[length]{} };
double* force1; // Structure of arrays vs array of structures
double* force2;

/**
 * Set up scenario from the command line.
 *
 * If you need additional helper data structures, you can
 * initialise them here. Alternatively, you can introduce a
 * totally new function to initialise additional data fields and
 * call this new function from main after setUp(). Either way is
 * fine.
 *
 * This operation's semantics is not to be changed in the assignment.
 */
void setUp(int argc, char** argv) {
  NumberOfBodies = (argc-4) / 7;

  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];
  mass = new double [NumberOfBodies];
  force0 = new double [NumberOfBodies]; // int *array{ new int[length]{} };
  force1 = new double [NumberOfBodies];
  force2 = new double [NumberOfBodies];
  
  
  colDist = 0.01/NumberOfBodies; // Initialize the collision distance
  
  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

  for (int i=0; i<NumberOfBodies; i++) {
    x[i] = new double[3];
    v[i] = new double[3];

    x[i][0] = std::stof(argv[readArgument]); readArgument++;
    x[i][1] = std::stof(argv[readArgument]); readArgument++;
    x[i][2] = std::stof(argv[readArgument]); readArgument++;

    v[i][0] = std::stof(argv[readArgument]); readArgument++;
    v[i][1] = std::stof(argv[readArgument]); readArgument++;
    v[i][2] = std::stof(argv[readArgument]); readArgument++;

    mass[i] = std::stof(argv[readArgument]); readArgument++;

    if (mass[i]<=0.0 ) {
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
  }
  
  // Precalculate distances:
  distCubMat = new double*[NumberOfBodies-1];
  
  for (int i =0; i<NumberOfBodies - 1; i++) {
	  distCubMat[i] = new double[NumberOfBodies];
  	for (int j = i + 1; j<NumberOfBodies; j++) {
		  double xDiff = x[j][0]-x[i][0];
    	double yDiff = x[j][1]-x[i][1];
    	double zDiff = x[j][2]-x[i][2];
    	
    	// Calculate distance
    	const double distance = sqrt( 
		    (xDiff * xDiff) +
		    (yDiff * yDiff) +
		    (zDiff * zDiff)
		  );
		
		// dist cubed
		distCubMat[i][j] = distance * distance * distance; 
  	}
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;
  
  if (tPlotDelta<=0.0) {
    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
    std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
    tPlot = 0.0;
  }
}


std::ofstream videoFile;


/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
  videoFile.open( "result.pvd" );
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}


/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
}


/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
  static int counter = -1;
  counter++;
  std::stringstream filename;
  filename << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
//      << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

  for (int i=0; i<NumberOfBodies; i++) {
    out << x[i][0]
        << " "
        << x[i][1]
        << " "
        << x[i][2]
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}


/**
 * This is the main operation you should change in the assignment. You might
 * want to add a few more variables or helper functions, but this is where the
 * magic happens.
 */
void updateBody() {
  minDx  = std::numeric_limits<double>::max();

  // force0 = force along x direction
  // force1 = force along y direction
  // force2 = force along z direction

  for (int i = 0; i < NumberOfBodies; i++) {
	  force0[i] = 0;
	  force1[i] = 0;
	  force2[i] = 0;
  }
  
  // Main force loop
  for (int i=0; i<NumberOfBodies; i++) { // Iterate through each body, our 'current' body is i and we calculate all forces from objects j = i + 1.
  	//double massDistCubed = mass[i]*mass[j] / (distance * distance * distance);
    for (int j = i + 1; j < NumberOfBodies; j++) { // We have worked out all forces on i from objects < i (used newtons 3rd law)
    	// Calculate difference along each component
    	double xDiff = x[j][0]-x[i][0];
    	double yDiff = x[j][1]-x[i][1];
    	double zDiff = x[j][2]-x[i][2];
    	
    	// Fx = ((m1 * m2) / r^2) * (m1x - m2x)/r
		
		  // Calculate mass * dist cubed
		  double massDistCubed = mass[j] / distCubMat[i][j]; // Took out multiplication by mass[i] to save time and increase accuracy
      
		  double forceX = xDiff * massDistCubed;
		  double forceY = yDiff * massDistCubed;
		  double forceZ = zDiff * massDistCubed;
		
		  // Apply forces from j on i
		  force0[i] += forceX;
		  force1[i] += forceY;
		  force2[i] += forceZ;
		
		  // Apply force from i on j
		  force0[j] -= forceX;
		  force1[j] -= forceY;
		  force2[j] -= forceZ;
   	}
   	
   	// After calculating all forces, use equations of motion
   	// Position update
   	x[i][0] = x[i][0] + timeStepSize * v[i][0];
  	x[i][1] = x[i][1] + timeStepSize * v[i][1];
  	x[i][2] = x[i][2] + timeStepSize * v[i][2];
  	
  	// Velocity update, F = MA -> A = F/M, F = M*F1 + M*F2.. -> A = F
   	v[i][0] = v[i][0] + (timeStepSize * force0[i]);
  	v[i][1] = v[i][1] + (timeStepSize * force1[i]);
  	v[i][2] = v[i][2] + (timeStepSize * force2[i]);
  	
  	double speedSqr = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
  	
  	if (speedSqr > maxV) {
  	  maxV = speedSqr;
  	}
  }
  
  maxV = sqrt(maxV); // Save squaring till after
  
  // Check for collisions & calculate distances
  // If an object i collides with another j, it is possible that the new object at i causes additional collisions with objects that have already checked for collisions with i.
  // This causes an issue as either we have to check everywhere for collisions again (slow) or find a way to make the previous objects only check for collisisions with i (difficult).
  // As i use a triangular matrix this is even more difficult, and the fact that checking whether an object collides with iteself will cause it to be overwritten completely by the end
  // object complicates things further, so this method comprimises between speed and complexity (otherwise we would have to use a stack to simulate true recursion)
  

  for (int i = 0; i<NumberOfBodies - 1; i++) {
  	for (int j = i + 1; j<NumberOfBodies; j++) {
		  if (i == j) {
	      continue;
	    }
	    
		  double xDiff = x[j][0]-x[i][0];
      double yDiff = x[j][1]-x[i][1];
      double zDiff = x[j][2]-x[i][2];
    	
    	// Calculate distance
    	double distance = sqrt( 
		    (xDiff * xDiff) +
		    (yDiff * yDiff) +
		    (zDiff * zDiff)
		  );

		
		  // Check if collided, if so then need to move things around
		  // common subexpression elimination
		  double combMass = mass[i] + mass[j];
		  if (distance < minDx) {
		    minDx = distance;
		  }
		  //std::cout << colDist * (combMass) << "\n";
		  if (distance <= colDist * (combMass)) {
		    //std::cout << "COLLISION!!!!!!!!!!!!!!!!!!!!!";
        // Change i to be collided object, swap j with item at end of list and then reduce number of bodies by 1
        v[i][0] = ((mass[i] * v[i][0]) + (mass[j] * v[j][0])) / combMass; // Save division till end to reduce error, rather than doing it twice & adding
        v[i][1] = ((mass[i] * v[i][1]) + (mass[j] * v[j][1])) / combMass;
        v[i][2] = ((mass[i] * v[i][2]) + (mass[j] * v[j][2])) / combMass;
        
        // Set weighted pos of i 
        x[i][0] = ((mass[i] * x[i][0]) + (mass[j] * x[j][0])) / combMass;
        x[i][1] = ((mass[i] * x[i][1]) + (mass[j] * x[j][1])) / combMass;
        x[i][2] = ((mass[i] * x[i][2]) + (mass[j] * x[j][2])) / combMass;
        
        mass[i] = combMass; // Do this last!
        
        NumberOfBodies -= 1;
        
        for (int m = 0; m < i; m++) { // Move distances calculated to j
          distCubMat[m][j] = distCubMat[m][NumberOfBodies];
        } 
        
        //Duplicate item and end of list into j's position
        x[j][0] = x[NumberOfBodies][0];
        x[j][1] = x[NumberOfBodies][1];
        x[j][2] = x[NumberOfBodies][2];
        
        v[j][0] = v[NumberOfBodies][0];
        v[j][1] = v[NumberOfBodies][1];
        v[j][2] = v[NumberOfBodies][2];
        
        mass[j] = mass[NumberOfBodies];
        
        if (j < i) { // If j is less than i then we know we have collided with a particle that has already had its values calcated and needs to be removed, so we start again from there (this is fairly unlikely to happen)
          i = j;
          j = j + 1;
          continue;
        }

        j = 0; // Recursion step. We effectively 'pivot' around i, checking if when when set i to be its new mass whether it collides with anything.
		  }
		
		  // dist cubed, alternatively could use a max function but i'm not sure that's faster
		  if (i < j) { // This is slow can i get rid of this somehow?
		    distCubMat[i][j] = distance * distance * distance;
		  } else {
		    distCubMat[j][i] = distance * distance * distance;
		  }
  	}
  }


  //maxV = std::sqrt( v[0][0]*v[0][0] + v[0][1]*v[0][1] + v[0][2]*v[0][2] );

  t += timeStepSize;


}


/**
 * Main routine.
 *
 * No major changes in assignment. You can add a few initialisation
 * or stuff if you feel the need to do so. But keep in mind that you
 * may not alter what the program plots to the terminal.
 */
int main(int argc, char** argv) {
  if (argc==1) {
    std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time dt objects" << std::endl
              << "  snapshot        interval after how many time units to plot. Use 0 to switch off plotting" << std::endl
              << "  final-time      simulated time (greater 0)" << std::endl
              << "  dt              time step size (greater 0)" << std::endl
              << std::endl
              << "Examples:" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0  \t One spiralling around the other one" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0 \t Three body setup from first lecture" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0 \t Five body setup" << std::endl
              << std::endl
              << "In this naive code, only the first body moves" << std::endl;

    return -1;
  }
  else if ( (argc-4)%7!=0 ) {
    std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
    std::cerr << "got " << argc << " arguments (three of them are reserved)" << std::endl;
    std::cerr << "run without arguments for usage instruction" << std::endl;
    return -2;
  }

  std::cout << std::setprecision(15);

  setUp(argc,argv);

  openParaviewVideoFile();

  int snapshotCounter = 0;
  if (t > tPlot) {
    printParaviewSnapshot();
    std::cout << "plotted initial setup" << std::endl;
    tPlot = tPlotDelta;
  }

  int timeStepCounter = 0;
  while (t<=tFinal) {
    updateBody();
    timeStepCounter++;
    if (t >= tPlot) {
      printParaviewSnapshot();
      std::cout << "plot next snapshot"
    		    << ",\t time step=" << timeStepCounter
    		    << ",\t t="         << t
				<< ",\t dt="        << timeStepSize
				<< ",\t v_max="     << maxV
				<< ",\t dx_min="    << minDx
				<< std::endl;

      tPlot += tPlotDelta;
    }
  }
  
  delete[] force0;
  delete[] force1;
  delete[] force2;

  std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
  std::cout << "Position of first remaining object: " << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;

  closeParaviewVideoFile();

  return 0;
}