#include "python_rho_cp_he_psi_alpha_mu__get_input_data.H"
// Info<< "min/max(rho_after_get) = "
  //          << min(rho).value() << ", " << max(rho).value() << endl;

//.....................................................pyhtonReturnOutoutValues......................................................
            // Call function and get a pointer to PyArrayObject
            // pValue stores the return value of the function
           
/*            array_2d = PyArray_SimpleNewFromData(2, dim, NPY_DOUBLE, &input_vals[0]);//?????it 2 means that array is 2D...................ADDED NS
            Info<< "Created numpy array (array_2d)" << endl;


	    PyTuple_SetItem(my_func_args, 0, array_2d);//?????.......................................................ADDED NS
            Info<< "Set array_2d as first argument of my_func_args tuple" << endl;*/

	   
	    // pValueRho = (PyArrayObject*) PyObject_CallObject (pFunc_rho, pArgs_rho);
            pValueAlpha = reinterpret_cast<PyArrayObject*>( 
                // PyObject_CallObject(pFunc_rho, pArgs_rho)
		PyObject_CallObject(pFunc_alpha, my_func_args)
		);
            // assign return values to the fields
	       /*   forAll(rho.internalField(), id)
	           {
               double rho = *((double*)PyArray_GETPTR1(pValueRho,id));
             }*/
	    //for (int mode = 0; mode < 1; ++mode)
             //{
             // Overwrite data
	     
	      // alpha->delete();
	      // std::array<int,1> alpha;
              
                 //alpha.setSize(174400,);
                 forAll(alpha.internalField(), id) // for boundary field use u_.boundaryField()
                 {
                     //  std::array<int,1> alpha;

                       alpha[id] = *((double*)PyArray_GETPTR1(pValueAlpha, id));
                       
                 }
             /*  forAll(alpha.internalField(), id) // for boundary field use u_.boundaryField()
                 {
                       alpha[id] = alpha1[id];
                       
                 }*/

               

                 // Write out
              /* char var_name[20];
               sprintf(var_name, "alpha%d", mode);
              // alpha.rename(&var_name[0]);
               alpha.write();*/

            //}
