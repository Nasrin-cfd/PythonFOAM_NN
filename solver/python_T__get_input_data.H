//.....................................................pyhtonGetDataCommon....................................................

            // INPUT VARIABLES: Y1, p, h
            //Compile Openfoam field data in a double array
	    forAll(Y[1].internalField(), id)
             {
              input_vals[id][0] = Y[1][id];
	         }
            forAll(p.internalField(), id)
             {
              input_vals[id][1] = p[id];
             }
            forAll(he.internalField(), id)
             {
              input_vals[id][2] = he[id];
             }
            
            // Cast to numpy before sharing reference with Python
            // Tell numpy about array dimensions
            // create a new array
            array_2d_for_T= PyArray_SimpleNewFromData(2, dim, NPY_DOUBLE, &input_vals[0]);//?????it 2 means that array is 2D
            Info<< "Created numpy array (array_2d_for_T)" << endl;

            PyTuple_SetItem(my_func_args_T, 0, array_2d_for_T);//?????
            Info<< "Set array_2d as first argument of my_func_args tuple" << endl;

            // Make array first argument
            // PyTuple_SetItem(pArgs_T, 0, array_2d_for_T);
            // Info<< "Set array_2d as first argument of func_args tuple for each outpt" << endl;

            // Pass the rank to the function
           // rank_val = PyLong_FromLong(rank);
	    // PyTuple_SetItem (pArgs_T, 1, rank_val);
            PyTuple_SetItem (my_func_args_T, 1, rank_val); 
           // Info<< "Set rank_val as second argument of func_args tuple for each outpt" << endl;
