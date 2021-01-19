// Compute meanRelVelVersor
                      
                        meanRelVel+=radius* rotOmega * alignedvector1;
                        if (mag(meanRelVel)>0)
                        {	
                             meanRelVelVersor=meanRelVel/mag(meanRelVel);
                        //Compute EVMsampDir versor orthogonal to both the EVMplaneNormal and the EVMupStreamDir versors
                        // (like a T, base of T is iEffForcePoint, T leg is stream dir and T top is EVM sampling direction).

                        	EVMsampDir=(alignedvector2^meanRelVelVersor);
                         
                         } 
                        else 
                        {   
                            EVMsampDir=-alignedvector0;
                        }               
	                    EVMsampDir/=mag(EVMsampDir); 
                    
                        // Compute the points of the EVM sampling line _CLAU
                        scalar EVMint=EVMlength/(numEVMPoints-1);
                       
                        int  i=0;  //valid if just one turbine
                        	for (int m = 0; m < numEVMPoints; m++)
                        	{  
                            		EVMPoint=point-EVMdist*meshDim*meanRelVelVersor+(EVMlength/2-EVMint*m)*meshDim*EVMsampDir;  
                            		EVMSamplePoints[i][j][k][m]=EVMPoint;
                        	}	 
