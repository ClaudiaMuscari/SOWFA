      Pstream::gather(minDisGlobalEVM,minOp<List<scalar> >());
      Pstream::scatter(minDisGlobalEVM);
     // Compare the global to local lists.  Where the lists agree, this processor controls
    // the actuator line point.

        
        forAll(EVMlinesControlled, t)
     {
        
        int f = EVMlinesControlled[t]; 
        int iterEVM=0;
        for(int j=0; j <3; j++)
        {
        for(int k=0; k < numBladePoints[0]; k++)
        {
		for(int m=0; m < numEVMPoints; m++)
       	        {       
 
                	if(minDisGlobalEVM[iterEVM] != minDisLocalEVM[iterEVM])
                        {       
                        	EVMMinDisCellID[f][j][k][m] = -1;
                        }
                        iterEVM++;
                 }
         }
         }
     }
