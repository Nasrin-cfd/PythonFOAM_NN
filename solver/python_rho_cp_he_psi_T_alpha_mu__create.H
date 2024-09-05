//.....................................................pyhtonCreate...........................................................
            // Initialize python, load module and create some scratch space etc.
            

            #include <Python.h>
            #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
            #include <numpy/arrayobject.h>
            Info<< "Initialize python" << nl << endl;
            Py_Initialize();
            PyRun_SimpleString("import sys");
            PyRun_SimpleString("sys.path.append(\".\")");
            PyRun_SimpleString("print('Hello world from Python interpreter!')");
            Info<< "Python initialize successful" << nl << endl;
           /* void init_numpy() {
            import_array1();
            }*/  
            //initialize numpy array library
            init_numpy();
            import_array1(-1);
            Info<< "numpy initialize successful" << nl << endl;


            // Load Python modules from python scripts in the working directory (case folder)
            Info<< "    load python_modules" << nl << endl;
            PyObject *pRho = PyUnicode_DecodeFSDefault("Python_Module_NN_rho"); // Python filename for rho
            PyObject *pModuleRho = PyImport_Import(pRho);
            
            /*PyObject *pCp = PyUnicode_DecodeFSDefault("Python_Module_NN_cp");   // Python filename for Cp
            PyObject *pModuleCp = PyImport_Import(pCp);*/

            PyObject *pHe = PyUnicode_DecodeFSDefault("Python_Module_NN_he");   // Python filename for he
            PyObject *pModuleHe = PyImport_Import(pHe);

            PyObject *pPsi = PyUnicode_DecodeFSDefault("Python_Module_NN_psi"); // Python filename for psi
            PyObject *pModulePsi = PyImport_Import(pPsi);

            PyObject *pT = PyUnicode_DecodeFSDefault("Python_Module_NN_T"); // Python filename for T
            PyObject *pModuleT = PyImport_Import(pT);

	    PyObject *pAlpha = PyUnicode_DecodeFSDefault("Python_Module_NN_alpha"); // Python filename for alpha
            PyObject *pModuleAlpha = PyImport_Import(pAlpha);

	    PyObject *pMu = PyUnicode_DecodeFSDefault("Python_Module_NN_mu"); // Python filename for mu
            PyObject *pModuleMu = PyImport_Import(pMu);



           if ( (!pModuleRho) || (!pModuleHe) || (!pModulePsi) || (!pModuleT) || (!pModuleAlpha) || (!pModuleMu) )
              {
                FatalErrorInFunction
                << "Errors loading python modules (missing imports?)" << nl
                << exit(FatalError);
              }
            
            //Inform OpenFOAM about what function names are and how many arguments they take
            PyObject *pFunc_rho = PyObject_GetAttrString(pModuleRho, "my_func"); // get the name
            //PyObject *pArgs_rho = PyTuple_New(1); // Tuple for arguments to my_func + rank 
            
           // PyObject *pFunc_cp = PyObject_GetAttrString(pModuleCp, "my_func"); // get the name
            // PyObject *pArgs_cp = PyTuple_New(1); // Tuple for arguments to my_func + rank 
            
            PyObject *pFunc_he = PyObject_GetAttrString(pModuleHe, "my_func"); // get the name
            // PyObject *pArgs_he = PyTuple_New(1); // Tuple for arguments to my_func + rank 
            
            PyObject *pFunc_psi = PyObject_GetAttrString(pModulePsi, "my_func"); // get the name
            // PyObject *pArgs_psi = PyTuple_New(1); // Tuple for arguments to my_func + rank 
            
            PyObject *pFunc_T = PyObject_GetAttrString(pModuleT, "my_func"); // get the name
            // PyObject *pArgs_T = PyTuple_New(1); // Tuple for arguments to my_func + rank   
            
	    PyObject *pFunc_alpha = PyObject_GetAttrString(pModuleAlpha, "my_func"); // get the name
	    // PyObject *pArgs_alpha = PyTuple_New(1); // Tuple for arguments to my_func + rank 
	   
	    PyObject *pFunc_mu = PyObject_GetAttrString(pModuleMu, "my_func"); // get the name
	    // PyObject *pArgs_mu = PyTuple_New(1); // Tuple for arguments to my_func + rank 
	    
	    // Make array first argument
	    PyObject* my_func_args = PyTuple_New(1); // Tuple for arguments to my_func  .we assign 3 places to the arg for 3 inputs..replace all args with my_func_args...ADDED NS
            PyObject* my_func_args_T = PyTuple_New(1);

            // Numpy datastructure to obtain return value from python method
            PyArrayObject *pValueRho;
           // PyArrayObject *pValueCp;
            PyArrayObject *pValueHe;
            PyArrayObject *pValuePsi;
            PyArrayObject *pValueT;
	    PyArrayObject *pValueAlpha;
	    PyArrayObject *pValueMu;

            // Release memory
            Py_DECREF(pRho);
           // Py_DECREF(pCp);
            Py_DECREF(pHe);
            Py_DECREF(pPsi);
            Py_DECREF(pT);
	    Py_DECREF(pAlpha);
            Py_DECREF(pMu);

            
            Py_DECREF(pModuleRho);
          //  Py_DECREF(pModuleCp);
            Py_DECREF(pModuleHe);
            Py_DECREF(pModulePsi);
            Py_DECREF(pModuleT);
	    Py_DECREF(pModuleAlpha);
	    Py_DECREF(pModuleMu);


            // Placeholder to grab data before sending to Python
            int num_cells = mesh.cells().size();
            double input_vals[num_cells][3];  // [num_cells][Y1 T p] for all variables, except T 
            double input_vals_T[num_cells][3];  // [num_cells][Y1 p h] for T..........3 is the input dimension, we have 3 inputs
            npy_intp dim[] = {num_cells, 3};  // num_cells is per rank   .....3 to 2 because array is 2D   NS

            // retreive ranks
            int rank = Pstream::myProcNo();
            PyObject *rank_val = PyLong_FromLong(rank);
	    //rank_val = PyLong_FromLong(rank);
            //Info<< "retreived ranks (rank_val)" << endl;

            // Preparing input data for sending to python
            PyObject *array_2d(nullptr);        // will store Y1, T, p  --> used to calc. rho, he, cp, psi, alpha, mu
	    PyObject *array_2d_for_T(nullptr);  // will store Y1, p, h  --> used to calc. T ...........................COMMENTED NS
