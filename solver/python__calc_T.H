#include "python_T__get_input_data.H"
//#include"python_rho_cp_he_psi__get_input_data.H"
//.....................................................pyhtonReturnOutoutValues......................................................
            // Call function and get a pointer to PyArrayObject
            // pValue stores the return value of the function
            
/*            array_2d_for_T = PyArray_SimpleNewFromData(2, dim, NPY_DOUBLE, &input_vals_T[0]);//?????it 2 means that array is 2D
            Info<< "Created numpy array (array_2d_for_T)" << endl;


            PyTuple_SetItem(my_func_args, 0, array_2d_for_T);//?????
            Info<< "Set array_2d as first argument of my_func_args tuple" << endl;*/

	    // pValueRho = (PyArrayObject*) PyObject_CallObject (pFunc_rho, pArgs_rho);
            pValueT = reinterpret_cast<PyArrayObject*>( 
                 // PyObject_CallObject(pFunc_T, pArgs_T)
	            PyObject_CallObject(pFunc_T, my_func_args_T)
                );
            // assign return values to the fields
            	    
           
	    //double d= *((double*)PyArray_GETPTR1(pValueT,:));
	  /* forAll(T.internalField(), id)
                 {
                   double T = *((double*)PyArray_GETPTR1(pValueT, id));

                 }*/
      //      T.write();
	   //for (int mode = 0; mode < 1; ++mode)
            //{
             // Overwrite data
               forAll(T.internalField(), id) // for boundary field use u_.boundaryField()
                 {
                        // std::array<int,1> T;

		    	 T[id] = *((double*)PyArray_GETPTR1(pValueT,id));
            
                 }

                 // Write out
              /* char var_name[20];
               sprintf(var_name, "T%d", mode);
            //   T.rename(&var_name[0]);
               T.write();*/

            //}
	    

            



	    
            
	    /*forAll(T.internalField(), id)
                 {

                    double T = *((double*)PyArray_GETPTR1(pValueT,id));
		    }*/

//           printf("T  is....................................................................... %f \n",T);
	   // printf("pValueT is....................................................................... %f \n",pValueT);

//            Info<< "this is from calc__T.........................................................."<< endl;

