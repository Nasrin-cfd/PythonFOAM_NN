//.....................................................pyhtonGetDataCommon..............................................
            // INPUT VARIABLES: Y1, T, p
            //Compile Openfoam field data in a double array
	    forAll(Y[1].internalField(), id)
             {
              input_vals[id][0] = Y[1][id];
	     }
	    forAll(T.internalField(), id)
             {
              input_vals[id][1] = T[id];
             }
            forAll(p.internalField(), id)
             {
              input_vals[id][2] = p[id];
             }
             //npy_intp dim[] = {num_cells, 3};
             array_2d = PyArray_SimpleNewFromData(2, dim, NPY_DOUBLE, &input_vals[0]);//it 2 means that array is 2D 
             //Info<< "Created numpy array (array_2d)" << endl;

	     
	   
             PyTuple_SetItem(my_func_args, 0, array_2d);
	    //PyTuple_SetItem(my_func_args, 0, array_2d_for_T);
            //Info<< "Set array_2d as first argument of my_func_args tuple" << endl;  //......................................ADDED NS


            // Pass the rank to the function
             rank_val = PyLong_FromLong(rank);
             //PyTuple_SetItem (my_func_args, 0, rank_val);
	     PyTuple_SetItem (my_func_args, 1, rank_val);
	     //PyTuple_SetItem (my_func_args, 2, rank_val);
	     //PyTuple_SetItem (my_func_args, 3, rank_val);
	     //PyTuple_SetItem (my_func_args, 4, rank_val); 
//            PyTuple_SetItem (pArgs_rho, 1, rank_val);
//            PyTuple_SetItem (pArgs_cp, 1, rank_val); 
//            PyTuple_SetItem (pArgs_he, 1, rank_val); 
//            PyTuple_SetItem (pArgs_psi, 1, rank_val);
//            Info<< "Set rank_val as second argument of func_args tuple for each outpt" << endl;


