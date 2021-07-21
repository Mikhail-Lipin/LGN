#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// To get an image Cube_Im_1i, execute Batch(0,1,1). "1i" at the end stands for "1 iteration". 
// Set the first argument to either 0 (shorter time, preferred) or 1 (longer time).
// The second argument should be 1 for shapes with vertices that can be represented by 3 crossing planes (e.g., cube and triangular pyramid,
// and 4 to 20 for shapes with vertices that cannot be represented by 3 crossing planes (e.g., sphere).
// The third argument is the maximal number of iterations. To iterate the edge-enhancement 6 times, execute Batch(0,1,6).


// Making 3D- image out of 1D- array:
Function Ar2Im(W)
	Wave/D W
	Variable Ni=22, Nj=22, Nk=22
	String NW=NameOfWave(W)+"_Im"
	Variable i, j, k
	Make/O/D/N=(Ni, Nj, Nk) RLGN
	for (i=0; i<Ni; i+=1)
		for (j=0; j<Nj; j+=1)
			for (k=0; k<Nk; k+=1)
				RLGN[i][j][k]=W[k+Nk*j+Nk*Nj*i]
			endfor
		endfor
	endfor
	Duplicate/O/D RLGN, $NW
	KillWaves/Z RLGN
End

// Initiating Edge-Enhancemet of the images listed in the table ImName
Function Batch(SL,Nmeans,Niterations)
	Variable SL, Nmeans, Niterations
	Variable Junk
	Wave/T ImName
	Junk=Nmeans
	if (Junk<=1)
		Nmeans=1
	else
		Nmeans=4*ceil(Junk/4)
	endif
	Junk=Niterations
	if (Junk<=1)
		Niterations=1
	else
		Niterations=round(Junk)
	endif
	Variable bi
	Variable NIm=dimsize(ImName,0)
	if (NIm == 0)
		Abort "ATTENTION:"+" The table ImName is empty"
	endif
	Wave/D Image
	String CurName
	for (bi=0; bi < NIm; bi+=1) // <11
		CurName=ImName[bi]
		if (WaveExists($CurName) == 0)
			Abort "ATTENTION:"+" The image "+CurName+" on the line #"+num2str(bi)+" of the table ImName does not exist"
		endif
		Reconstruction(CurName, SL, Nmeans, Niterations)
	endfor
End

// The edge-enhancing function:
Function Reconstruction(CurName, SL, Nmeans, Niterations) // The input is the 22x22x22 image and the desired number of edges to sum up in the processed image
	String CurName// name of the image to be processed
	Variable SL // Short or Long
	Variable Nmeans // the number of component edges 
	Variable Niterations // the number of iterations
	Variable Resampling=1 // Resampling == 1 results in averaging of border values and interlacing.
	Wave/D Image // image to be processed
	Duplicate/O/D $CurName, Image
	Variable NmeansRec=1/Nmeans
	if (SL == 0)
		Wave/D E4DNd=root:MultiplaneUnitsShort:E3P4DNd
		Wave/D E4D27=root:MultiplaneUnitsShort:E3P4D27
		Wave/D E4D27ByVariance=root:MultiplaneUnitsShort:E3P4D27ByVariance
		Wave/D E2D27Stat=root:MultiplaneUnitsShort:E3P2D27Stat
	elseif (SL == 1)
		Wave/D E4DNd=root:MultiplaneUnitsLong:E3P4DNd
		Wave/D E4D27=root:MultiplaneUnitsLong:E3P4D27
		Wave/D E4D27ByVariance=root:MultiplaneUnitsLong:E3P4D27ByVariance
		Wave/D E2D27Stat=root:MultiplaneUnitsLong:E3P2D27Stat
	else
		Abort "The value of the first argument of Batch function must be either 0 (short) or 1 (long)"
	endif
	Variable t1=ticks/3600
	Variable Nx=dimsize(Image,0), Ny=dimsize(Image,1), Nz=dimsize(Image,2)
	Variable Nd=2
	Variable Tol=1e-16
	Variable DS=dimsize(E2D27Stat,0) // the number of edges the image is to be fit with.
	Make/O/D/N=(DS) SSRs, EdgeIDs // the columns with SSR values and EdgeID values
	Variable SSRblank, SNRblank, EdgeIDsTemp
	Variable SumOfCovariance, base, baseWeighted, baseWeightedDiv, SumOfSNRs
	Variable i, j, k, ii, jj, kk, iii, jjj, kkk, i3, j3, k3
	Variable coef, coef_Fit, base_Fit, baseOr, SSR_Fit, SSR_Or, SSRJ, itot
	Variable SumTemp, SSRmax, Nit, NitB, NitBx, NitBy, NitBz, MaxTemp=(DS-1)
	Make/O/D/N=(Nmeans) SNRs, Weights
	Make/O/D/N=(3,3,3) ImZM, ImZMSq // local 3x3x3 image with zero mean.
	Make/O/D/N=(3,3,3) E3D27Covariance, ResidualSq
// Increasing size of the Image by 2 voxels in each dimension by duplicating border voxels:
	Duplicate/O/D/R=[0][][] Image, ImageF
	Duplicate/O/D/R=[Nx-1][][] Image, ImageL
	Concatenate/O/NP=0/D {ImageF, Image, ImageL}, ImageX
	Duplicate/O/D/R=[][0][] ImageX, ImageF
	Duplicate/O/D/R=[][Ny-1][] ImageX, ImageL
	Concatenate/O/NP=1/D {ImageF, ImageX, ImageL}, ImageY	
	Duplicate/O/D/R=[][][0] ImageY, ImageF
	Duplicate/O/D/R=[][][Nz-1] ImageY, ImageL
	Concatenate/O/NP=2/D {ImageF, ImageY, ImageL}, ImageFiltered0
	Make/O/D/N=(Nx+3, Ny+3, Nz+3) ImageFiltered1=0 
	Make/O/D/N=((Nx+2)*Nd, (Ny+2)*Nd, (Nz+2)*Nd) ImageUp0=0
	Make/O/D/N=((Nx+3)*Nd, (Ny+3)*Nd, (Nz+3)*Nd) ImageUp1=0
	Make/O/D/N=((Nx+2)*Nd, (Ny+2)*Nd, (Nz+2)*Nd) ImageN0=0
	Make/O/D/N=((Nx+3)*Nd, (Ny+3)*Nd, (Nz+3)*Nd) ImageN1=0
	String ImUpName
// Sort SSRs, SSRs, EdgeIDs   will arrange SSRs in ascending order saving the correspondence SSR-EdgeID
	for (Nit=0; Nit<floor((Niterations+1.2)/2); Nit+=1) // Nit = the number of iterations
			NitB=0
			NitBx=Nx+NitB+2 // NitBx by NitBy by NitBz is the size of ImageFiltered
			NitBy=Ny+NitB+2
			NitBz=Nz+NitB+2
			ImageUp0=0
			ImageN0=0
			ImageUp1=0
			ImageN1=0
			for (i=1; i<NitBx-1; i+=1)
				for (j=1; j<NitBy-1; j+=1)
					for (k=1; k<NitBz-1; k+=1)
						Duplicate/O/D/R=[i-1,i+1][j-1,j+1][k-1,k+1] ImageFiltered0, Im
						base=mean(Im)
						ImZM=Im-base
						ImZMSq=ImZM^2
						EdgeIDs=x
						for (ii=0; ii<DS; ii+=1) // over all the edges
							Duplicate/O/D/R=[][][][ii] E4D27, JunkEdge // making 3x3x3 copy of the essential edge to be compared to.
							Redimension/D/N=(3,3,3,0) JunkEdge // This is "Edge"
							E3D27Covariance=JunkEdge*ImZM
							SumOfCovariance=sum(E3D27Covariance) // SumOfCovariance is zero for dissimilar edges.
							Duplicate/O/D/R=[][][][ii] E4D27ByVariance, JunkByVariance
							Redimension/D/N=(3,3,3,0) JunkByVariance
							ResidualSq = (ImZM - JunkByVariance*SumOfCovariance)^2
							SSRs[ii]=sum(ResidualSq)
						endfor // ii<DS
						Sort SSRs, SSRs, EdgeIDs // Topping EdgeIDs with the least SSR.
						Duplicate/O/D/R=[0,Nmeans-1] SSRs, SSRsNmeans
						SSRmax=(SSRs[Nmeans-1]-SSRs[0]+Tol)
						Duplicate/O/D/R=[0,Nmeans-1] SSRs, SSRsNmeans
						SNRs=exp(-20*(SSRsNmeans-SSRs[0])/SSRmax)
						SumOfSNRs=sum(SNRs)
						Weights=SNRs/SumOfSNRs // Nmeans weights.
						baseWeightedDiv=0
						for (iii=0; iii<Nmeans; iii+=1) // Cycling over weights
							EdgeIDsTemp=EdgeIDs[iii]
							Duplicate/O/D/R=[][][][EdgeIDsTemp] E4D27, JunkEdge // making 3x3x3 copy of the essential edge to be compared to.
							Redimension/D/N=(3,3,3,0) JunkEdge // This is "Edge"
							E3D27Covariance=JunkEdge*ImZM
							SumOfCovariance=sum(E3D27Covariance) // SumOfCovariance is zero for dissimilar edges.
							for (i3=0; i3<3; i3+=1)
								for (j3=0; j3<3; j3+=1)
									for (k3=0; k3<3; k3+=1)
										for (ii=0; ii<Nd; ii+=1) // Cycling over subvoxels 
											for (jj=0; jj<Nd; jj+=1)
												for (kk=0; kk<Nd; kk+=1)
													if ( ((i+i3) == 1) || ((j+j3) == 1) || ((k+k3) == 1) || ((i+i3) == NitBx) || ((j+j3) == NitBy) || ((k+k3) == NitBz) ) // Marginal voxels
														ImageUp0[(i+i3-1)*Nd+ii][(j+j3-1)*Nd+jj][(k+k3-1)*Nd+kk]+=(baseWeightedDiv+Weights[iii]*(base+(E4DNd[i3*Nd+ii][j3*Nd+jj][k3*Nd+kk][EdgeIDsTemp]-E2D27Stat[EdgeIDsTemp][0])*SumOfCovariance/E2D27Stat[EdgeIDsTemp][1]))
														ImageN0[(i+i3-1)*Nd+ii][(j+j3-1)*Nd+jj][(k+k3-1)*Nd+kk]+=NmeansRec
													elseif ( (i3 == 1) && (j3 == 1) && (k3 == 1) ) // Non-marginal voxels
														ImageUp0[i*Nd+ii][j*Nd+jj][k*Nd+kk]+=(baseWeightedDiv+Weights[iii]*(base+(E4DNd[Nd+ii][Nd+jj][Nd+kk][EdgeIDsTemp]-E2D27Stat[EdgeIDsTemp][0])*SumOfCovariance/E2D27Stat[EdgeIDsTemp][1]))
														ImageN0[i*Nd+ii][j*Nd+jj][k*Nd+kk]+=NmeansRec
													endif
												endfor // kk
											endfor // jj
										endfor // ii
									endfor // k3
								endfor // j3
							endfor // i3
						endfor // iii
					endfor // k
				endfor // j
			endfor // i
			// Third, divide the sum of image intensities by the number of edges in the sum:
			ImageUp0=ImageUp0/ImageN0
			ImUpName=CurName+"_"+num2str(2*(Nit+1)-1)+"i"
			Duplicate/O/D ImageUp0, $ImUpName
			if ((2*Nit+1) == Niterations)
				break
			endif
			if (Resampling == 1)
				NitB=1
				NitBx=Nx+NitB+2 // NitBx by NitBy by NitBz    is the size of ImageFiltered
				NitBy=Ny+NitB+2
				NitBz=Nz+NitB+2
				Duplicate/O/D/R=[0][][] ImageUp0, ImageF
				Duplicate/O/D/R=[NitBx*Nd-1][][] ImageUp0, ImageL
				Concatenate/O/NP=0/D {ImageF, ImageUp0, ImageL}, ImageX
				Duplicate/O/D/R=[][0][] ImageX, ImageF
				Duplicate/O/D/R=[][NitBy*Nd-1][] ImageX, ImageL
				Concatenate/O/NP=1/D {ImageF, ImageX, ImageL}, ImageY	
				Duplicate/O/D/R=[][][0] ImageY, ImageF
				Duplicate/O/D/R=[][][NitBz*Nd-1] ImageY, ImageL
				Concatenate/O/NP=2/D {ImageF, ImageY, ImageL}, ImageUp1
				for (i=0; i<NitBx; i+=1) // NitBx=24. Making ImageFiltered1 by averaging values in ImageUp1:
					for (j=0; j<NitBy; j+=1)
						for (k=0; k<NitBz; k+=1)
							Duplicate/O/D/R=[i*Nd,i*Nd+Nd-1][j*Nd,j*Nd+Nd-1][k*Nd,k*Nd+Nd-1] ImageUp1, Junk
							ImageFiltered1[i][j][k]=mean(Junk)
						endfor // k
					endfor // j
				endfor // i				
				ImageUp1=0
				for (i=1; i<NitBx-1; i+=1)
					for (j=1; j<NitBy-1; j+=1)
						for (k=1; k<NitBz-1; k+=1)
							Duplicate/O/D/R=[i-1,i+1][j-1,j+1][k-1,k+1] ImageFiltered1, Im
							base=mean(Im)
							ImZM=Im-base
							ImZMSq=ImZM^2
							EdgeIDs=x
							for (ii=0; ii<DS; ii+=1) // over all the edges
								Duplicate/O/D/R=[][][][ii] E4D27, JunkEdge // making 3x3x3 copy of the essential edge to be compared to.
								Redimension/D/N=(3,3,3,0) JunkEdge // This is "Edge"
								E3D27Covariance=JunkEdge*ImZM
								SumOfCovariance=sum(E3D27Covariance) // SumOfCovariance is zero for dissimilar edges.
								Duplicate/O/D/R=[][][][ii] E4D27ByVariance, JunkByVariance
								Redimension/D/N=(3,3,3,0) JunkByVariance
								ResidualSq = (ImZM - JunkByVariance*SumOfCovariance)^2
								SSRs[ii]=sum(ResidualSq)
							endfor // ii<DS
							Sort SSRs, SSRs, EdgeIDs // Topping EdgeIDs with the least SSR.
							Duplicate/O/D/R=[0,Nmeans-1] SSRs, SSRsNmeans
							SSRmax=(SSRs[Nmeans-1]-SSRs[0]+Tol)
							Duplicate/O/D/R=[0,Nmeans-1] SSRs, SSRsNmeans
							SNRs=exp(-20*(SSRsNmeans-SSRs[0])/SSRmax)
							SumOfSNRs=sum(SNRs)
							Weights=SNRs/SumOfSNRs // Nmeans weights.
							baseWeightedDiv=0						
							for (iii=0; iii<Nmeans; iii+=1) // Cycling over weights
								EdgeIDsTemp=EdgeIDs[iii]
								Duplicate/O/D/R=[][][][EdgeIDsTemp] E4D27, JunkEdge // making 3x3x3 copy of the essential edge to be compared to.
								Redimension/D/N=(3,3,3,0) JunkEdge // This is "Edge"
								E3D27Covariance=JunkEdge*ImZM
								SumOfCovariance=sum(E3D27Covariance) // SumOfCovariance is zero for dissimilar edges.
								for (i3=0; i3<3; i3+=1)
									for (j3=0; j3<3; j3+=1)
										for (k3=0; k3<3; k3+=1)
											for (ii=0; ii<Nd; ii+=1) // Cycling over subvoxels 
												for (jj=0; jj<Nd; jj+=1)
													for (kk=0; kk<Nd; kk+=1)
														if ( ((i+i3) == 1) || ((j+j3) == 1) || ((k+k3) == 1) || ((i+i3) == NitBx) || ((j+j3) == NitBy) || ((k+k3) == NitBz) ) // Marginal voxels
															ImageUp1[(i+i3-1)*Nd+ii][(j+j3-1)*Nd+jj][(k+k3-1)*Nd+kk]+=(baseWeightedDiv+Weights[iii]*(base+(E4DNd[i3*Nd+ii][j3*Nd+jj][k3*Nd+kk][EdgeIDsTemp]-E2D27Stat[EdgeIDsTemp][0])*SumOfCovariance/E2D27Stat[EdgeIDsTemp][1]))
															ImageN1[(i+i3-1)*Nd+ii][(j+j3-1)*Nd+jj][(k+k3-1)*Nd+kk]+=NmeansRec
														elseif ( (i3 == 1) && (j3 == 1) && (k3 == 1) ) // Non-marginal voxels
															ImageUp1[i*Nd+ii][j*Nd+jj][k*Nd+kk]+=(baseWeightedDiv+Weights[iii]*(base+(E4DNd[Nd+ii][Nd+jj][Nd+kk][EdgeIDsTemp]-E2D27Stat[EdgeIDsTemp][0])*SumOfCovariance/E2D27Stat[EdgeIDsTemp][1]))
															ImageN1[i*Nd+ii][j*Nd+jj][k*Nd+kk]+=NmeansRec
														endif
													endfor // kk
												endfor // jj
											endfor // ii
										endfor // k3
									endfor // j3
								endfor // i3
							endfor // iii
						endfor // k
					endfor // j
				endfor // i
				// Third, divide the sum of image intensities by the number of edges in the sum:
				ImageUp1=ImageUp1/ImageN1
				ImUpName=CurName+"_"+num2str(2*(Nit+1))+"i"
				Duplicate/O/D ImageUp1, $ImUpName
				for (i=0; i<(NitBx-1); i+=1) // Making ImageFiltered0 by averaging values in ImageUp1:
					for (j=0; j<(NitBy-1); j+=1)
						for (k=0; k<(NitBz-1); k+=1)
							Duplicate/O/D/R=[i*Nd+1,i*Nd+Nd][j*Nd+1,j*Nd+Nd][k*Nd+1,k*Nd+Nd] ImageUp1, Junk
							ImageFiltered0[i][j][k]=mean(Junk)
						endfor // k
					endfor // j
				endfor // i
		else // NitBx=24
			for (i=0; i<(NitBx); i+=1) // Making ImageFiltered0 by averaging values in ImageUp0:
				for (j=0; j<(NitBy); j+=1)
					for (k=0; k<(NitBz); k+=1)
						Duplicate/O/D/R=[i*Nd,i*Nd+Nd-1][j*Nd,j*Nd+Nd-1][k*Nd,k*Nd+Nd-1] ImageUp0, Junk
						ImageFiltered0[i][j][k]=mean(Junk)
					endfor // k
				endfor // j
			endfor // i
		endif
	endfor //Nit
	KillWaves/Z Junk, Image, E3D27Covariance, EdgeIDs, Im, ImageF, ImageFiltered0, ImageFiltered1, ImageL, ImageN0, ImageN1
	KillWaves/Z ImageUp0, ImageUp1, ImageX, ImageY, ImName, ImZM, ImZMSq, JunkByVariance, JunkEdge
	KillWaves/Z ResidualSq, SNRs, SSRs, SSRsNmeans, Weights, V_Flag 
	print "The elapsed time is "+num2str(ticks/3600-t1)+" minutes"
	print "The image "+CurName+" has been processed"
	Beep
	Sleep/T 30
	Beep
End


// Making 1D- array from 3D- image:
Function Im2Ar(W)
	Wave W
	Variable Nd=2 // a voxel of the original image is split into Nd*Nd*Nd subvoxels
	Variable Ni, Nj, Nk, Nstart
	print "dimsize=",dimsize(W,0)
	if (dimsize(W,0) == 48)
		Ni=dimsize(W,0)-2*Nd
		Nj=dimsize(W,1)-2*Nd
		Nk=dimsize(W,2)-2*Nd
		Nstart=Nd
	elseif (dimsize(W,1) == 50)
		Ni=dimsize(W,0)-3*Nd
		Nj=dimsize(W,1)-3*Nd
		Nk=dimsize(W,2)-3*Nd
		Nstart=Nd+1
	elseif (dimsize(W,1) == 22)
		Ni=dimsize(W,0)
		Nj=dimsize(W,1)
		Nk=dimsize(W,2)
		Nstart=0
	else
		Abort "ATTENTION: the size of the input wave should be either 50x50x50 or 48x48x48 voxels"
		Beep
	endif
		String NW=NameOfWave(W), NA=NW+"", NAtxt=NA+".txt"
		Variable i, j, k, Ntot
		Make/O/T/N=(Ni*Nj) Arr
		Make/O/T/N=0 ArrConc
		for (k=0; k<Nk; k+=1)
			for (j=0; j<Nj; j+=1)
				for (i=0; i<Ni; i+=1)
					Ntot=i+Ni*j
					Arr[Ntot]=num2str(i)+"x"+num2str(j)+"x"+num2str(k)+"vox "+num2str(W[i+Nstart][j+Nstart][k+Nstart])
				endfor
			endfor
		Duplicate/O/T ArrConc, ArrConcJunk
		Concatenate/O/T/NP {ArrConcJunk, Arr}, ArrConc
		if ( mod(k,50) == 0)
		endif
		endfor
		Save/J/M="\n" ArrConc as NAtxt
		KillWaves Arr, ArrConc, ArrConcJunk
End

