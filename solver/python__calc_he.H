#include "python_rho_cp_he_psi_alpha_mu__get_input_data.H"
#include <iostream>
#include <vector>
#include "Tensor.H"
	  

// Assuming pValueHe is a PyArrayObject* containing the data you want to write into the scalar field

// Get the size of the PyArrayObject
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*npy_intp size = PyArray_SIZE(pValueHe);

// Assuming scalar_field is a 1D array (you might need to adjust this based on your scalar field)
// and assuming pValueHe contains double values
double* scalar_field = new double[size];
//double* scalar_field[471];

// Copy values from PyArrayObject to scalar_field
for (npy_intp i = 0; i < size; ++i) {
    // Accessing the elements of pValueHe
    double value = *reinterpret_cast<double*>(PyArray_GETPTR1(pValueHe, i));
    // Assigning the value to the corresponding position in scalar_field
    scalar_field[i] = value;   //NS: This part doesn't run in parallel as well :(    

}*/
//delete[] scalar_field[i]; //think about it NS

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*pValueHe = reinterpret_cast<PyArrayObject*>(
                // PyObject_CallObject(pFunc_rho, pArgs_rho)
                PyObject_CallObject(pFunc_he, my_func_args)
                );*/
//#pragma omp parallel for


/*PyArrayObject* pArray = reinterpret_cast<PyArrayObject*>(pValueHe);
    npy_intp* shape = PyArray_DIMS( pArray);
    int num_elements = shape[0];
    double* data = static_cast<double*>(PyArray_DATA( pArray));
    Py_DECREF(pValueHe);
    Py_DECREF(data);
    Py_DECREF(pArray);
    // Finalize Python
    Py_Finalize();

    std::vector<double> values(data, data + num_elements);
    
    
    // Make a deep copy of the data
    //std::vector<double> localData(data, data + num_elements);

    // Construct the values vector using the local copy
    //std::vector<double> values(std::move(localData));
    // Write to OpenFOAM field file
    writeScalarField(values, "/home/nasrin/from_docker_test/parallel_run/from_0.001_PyFOAM/0/he");*/ // Replace "0/T" with the appropriate file path and name
  //  Py_DECREF(pValueHe);
    // Finalize Python
    //Py_Finalize();





/*
    // Define ScalarField structure
    struct ScalarField {
       std::vector<double> values;
    };
 // Assume pFunc_he and my_func_args are already defined and initialized

    // Call the Python function
    PyObject* pResult = PyObject_CallObject(pFunc_he, my_func_args);

    // Cast the result to PyArrayObject
    pValueHe = reinterpret_cast<PyArrayObject*>(pResult);

    // Get the number of elements
    npy_intp num_elements = PyArray_SIZE(pValueHe);

    // Get a pointer to the data
    double* data = static_cast<double*>(PyArray_DATA(pValueHe));

    // Create and resize the ScalarField to hold all values
    ScalarField he2;
    he2.values.resize(num_elements);

    // Populate the ScalarField in parallel using OpenMP
   //#pragma omp parallel for
    for (npy_intp i = 0; i < num_elements; ++i) {
        he2.values[i] = data[i];
    }

    #pragma omp parallel
    {
        int num_threads = omp_get_num_threads();
        int thread_id = omp_get_thread_num();
        npy_intp chunk_size = num_elements / num_threads;
        npy_intp start = thread_id * chunk_size;
        npy_intp end = (thread_id == num_threads - 1) ? num_elements : start + chunk_size;

        for (npy_intp i = start; i < end; ++i) {
            he2.values[i] = data[i];
        }
    }

    // Clean up
    Py_DECREF(pValueHe);
    Py_Finalize();

    // Print the ScalarField values (for verification)
    for (const auto& value : he2.values) {
        std::cout << value << "\n";
    }

    std::cout << "All values have been populated into the ScalarField\n";*/
    







pValueHe = reinterpret_cast<PyArrayObject*>(
                // PyObject_CallObject(pFunc_rho, pArgs_rho)
                PyObject_CallObject(pFunc_he, my_func_args)
                );
/*Py_DECREF(my_func_args);
Py_DECREF(pFunc_he);
//Py_DECREF(pModule);
npy_intp num_elements = PyArray_SIZE(pValueHe);
std::vector<double> he2;
he2.reserve(num_elements);
double arr[471];*/


//he1= {1};
 

/*for (npy_intp i = 0; i < num_elements; ++i) {
        he = *((double*)PyArray_GETPTR1(pValueHe, i));
        //he2.push_back(he[i]);
    }*/
//Py_DECREF(pValueHe);
//Py_Finalize();

forAll(he.internalField(), id) // for boundary field use u_.boundaryField()
                 {
                    
		   //double value = *reinterpret_cast<double*>(PyArray_GETPTR1(pValueHe, id));	 
		  //double_array[id] = PyFloat_AsDouble(pValueHe[id]);	 
//	         double he1= *((double*)PyArray_GETPTR1(pValueHe,id));
//                 he2.push_back(he1);i
                    he[id]= *((double*)PyArray_GETPTR1(pValueHe,id));
		  //he2.push_back(he1);
	 //          he[id]=he1[id];    
		  
		  //he[id]=arr[id];
		  //he[id]=double_array[id];
                   //he[id]=doubleValue[id];
                 };

/*forAll(he.internalField(), id) // for boundary field use u_.boundaryField()
                 {
		    he = *((double*)PyArray_GETPTR1(pValueHe,id)) ;            //Extract value from numpy array
    		    he[id] = arr[id];  // Set value to the scalar field
                 }*/

//he[470]=arr[470];

	/*for (int i = 0; i < 471; ++i) {

                 he[i]=arr[i];
                  }*/
/*forAll(he,  cellI) /i/ for boundary field use u_.boundaryField()
                 {
                    he[cellI]= *((double*)PyArray_GETPTR1(pValueHe,cellI));
                  //  he=he.push_back(he);
                 } */// for boundary field use u_.boundaryField()
// Info<<"cellI is:"<<cellI<<endl;
/*const int size = 471;
double arr[size]; // Create an array of size elements

    // Append values to the array in a loop
for (int i = 0; i < size; ++i) {
        arr[i] = *((double*)PyArray_GETPTR1(pValueHe,i)); // Assign the value of i to the ith element of the array
    }

forAll(he, cellI)
{
    // Assuming arr is a 1D array and you want to assign its values to each cell
    he[cellI] = arr[cellI];
}
he.write();*/
      /*        irank_val = PyLong_FromLong(rank);
	      PyTuple_SetItem(my_func_args, 0, array_2d);
              PyTuple_SetItem(my_func_args, 0, rank_val);*/

	    // pValueRho = (PyArrayObject*) PyObject_CallObject (pFunc_rho, pArgs_rho);
     /*       pValueHe = reinterpret_cast<PyArrayObject*>( 
                // PyObject_CallObject(pFunc_rho, pArgs_rho)
		PyObject_CallObject(pFunc_he, my_func_args)
		);*/
              
                 //alpha.setSize(174400,);



  // ////////////////////////          if (pValueHe != NULL && PyArray_NDIM(pValueHe) == 1) {
              // Get the size of the array
              /*int size = PyArray_DIM(pValueHe, 0);

              // Allocate memory for the double array
              double* doubleArray = new double[size];

              // Get a pointer to the data
              double* pData = static_cast<double*>(PyArray_DATA(pValueHe));
              //Info<<"pData is:"<<pData<<endl;
             
	      // Copy data from PyArrayObject to double array
              
	      int rank, size;
              MPI_Init(NULL, NULL);
              MPI_Comm_rank(MPI_COMM_WORLD, &rank);
              MPI_Comm_size(MPI_COMM_WORLD, &size);

             int chunk_size = size / num_processes;
             int remainder = size % num_processes;
             int start = rank * chunk_size + (rank < remainder ? rank : remainder);
             int end = start + chunk_size + (rank < remainder ? 1 : 0);

             for (int i = start; i < end; ++i) {
               doubleArray[i] = pData[i];
                 }*/
	         ///////////////////////////////////////////////////////////////////////////	     
	      

             /*void distributeData(const std::vector<double>& pValueHe, std::vector<double>& he, int rank, int num_ranks) {
                 //  int num_ranks=3;
		   int size=PyArray_DIM(pValueHe, 0);
                   int num_elements_per_rank = size / num_ranks;
                   int start_index = rank * num_elements_per_rank;
                   int end_index = (rank + 1) * num_elements_per_rank;
                   if (rank == num_ranks - 1) {
                   end_index = total_elements; // Last rank may have fewer elements
                     }

    // Copy data from pValueHe to he for the assigned range
                   for (int i = start_index; i < end_index; ++i) {
                      he[i - start_index] = pValueHe[i];
                          }
                    }*/
	   /*  int num_ranks=3;
	     //int total_elements = pValueHe.size();
	     int size = PyArray_DIM(pValueHe, 0);
             int num_elements_per_rank = size / num_ranks;
             int start_index = rank * num_elements_per_rank;
             int end_index = (rank + 1) * num_elements_per_rank;
             if (rank == num_ranks - 1) {
                 end_index = size; // Last rank may have fewer elements
                      }*/
             
             
		 //{   
               /*forAll(he.internalField(), id) // for boundary field use u_.boundaryField()
                 {
                     //  std::array<int,1> alpha;
                       
		       //double* pointerValue = (double*)PyArray_GETPTR1(pValueHe, id);
                       //double scalarValue = *pointerValue;
                       //he[id] = scalarValue;
                    
		     //he[id] = *((double*)PyArray_GETPTR1(pValueHe, id));        
                    //he[id]=1;
                      
		  //  Info<<*((double*)PyArray_GETPTR1(pValueHe,0))<<endl;
		    he[id] = *((double*)PyArray_GETPTR1(pValueHe,id));
	//	    he[id] = std::vector<double> pValueHe[id];
		      
		 }*/
//              Py_Finalize(); 
	       //char var_name[20];
               //sprintf(var_name, "he%d", mode);
               //he.rename(&var_name[0]);
               //he.write();
                
                

