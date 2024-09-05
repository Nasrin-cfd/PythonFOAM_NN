#include "python_rho_cp_he_psi_alpha_mu__get_input_data.H"

//.....................................................pyhtonReturnOutoutValues......................................................
            // Call function and get a pointer to PyArrayObject
          
	  
/*	    array_2d = PyArray_SimpleNewFromData(2, dim, NPY_DOUBLE, &input_vals[0]);//?????it 2 means that array is 2D ............................ADDED NS
            Info<< "Created numpy array (array_2d)" << endl;

	    // pValue stores the return value of the function
            PyTuple_SetItem(my_func_args, 0, array_2d);//?????............................................ADDED NS
            Info<< "Set array_2d as first argument of my_func_args tuple" << endl;*/
 
	    // pValueRho = (PyArrayObject*) PyObject_CallObject (pFunc_rho, pArgs_rho);
            pValuePsi = reinterpret_cast<PyArrayObject*>( 
		// PyObject_CallObject(pFunc_he, pArgs_psi)
                PyObject_CallObject(pFunc_psi, my_func_args) //..............................NS REPLECE my_func_args
                );
            // assign return values to the fields
	          /*forAll(psi.internalField(), id)
	           {
               double psi = *((double*)PyArray_GETPTR1(pValuePsi,id));
             }      */
           //for (int mode = 0; mode < 1; ++mode)
             //{
             // Overwrite data
               forAll(psi.internalField(), id) // for boundary field use u_.boundaryField()
                 {
                       //  std::array<int,1> psi;

		   	 psi[id] = *((double*)PyArray_GETPTR1(pValuePsi, id));

                 }

                 // Write out
              /* char var_name[20];
               sprintf(var_name, "psi%d", mode);
               //psi.rename(&var_name[0]);
               psi.write();*/

           // }	    
            	   
                  
