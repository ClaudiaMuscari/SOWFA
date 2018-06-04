{
    // The inverse of the momentum equation "A" operator (the matrix diagonal) and
    // its projection from cell centers to cell faces.  Because it is a diagonal matrix,
    // The inverse is just the diagonal matrix of the reciprocals, hence the "r" in the
    // name "rAU".
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("(1|A(U))", fvc::interpolate(rAU));

    // The momentum equation "H" operator is the off-diagonal times the last known U. 
    // In the equations it is always multiplied with inv(A).  This is HbyA.
    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();

    // This is the Boussinesq buoyancy term multipled by the r
    surfaceScalarField phig(-rAUf*ghf*fvc::snGrad(rhok)*mesh.magSf());

    // Project HbyA to cell faces and apply a correction for time stepping scheme.
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
       + rAUf*fvc::ddtCorr(U, phi)
    );

    // This part can be confusing.  
  //surfaceScalarField phiFixedFlux = fvc::interpolate(U) & mesh.Sf();
    surfaceScalarField phiFixedFlux = phiHbyA;

   
  //adjustPhi(phiHbyA, U, p_rgh);
    adjustPhi(phiFixedFlux, U, p_rgh);


    phiHbyA += phig;

/*
    forAll(p_rgh.boundaryField(),i)
    {
        if ( !mesh.boundary()[i].coupled() )
        {
            const fvPatchScalarField& pp = p_rgh.boundaryField()[i];
            bool isAffp = isA<fixedFluxPressureFvPatchScalarField>(pp);
            word boundaryName = mesh.boundary()[i].name();
            if (isAffp)
            {
               Info << "Patch " << i << ", " << boundaryName << ", is a fixedFluxPressure boundary for p_rgh" << endl;
            }
            else
            {
               Info << "Patch " << i << ", " << boundaryName << ", is not a fixedFluxPressure boundary for p_rgh" << endl;
            }
        }
    }
*/

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryField(),
        (
            phiHbyA.boundaryField()
          - phiFixedFlux.boundaryField()
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );




    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi = phiHbyA - p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
        }
    }

    p = p_rgh + rhok*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rhok*gh;
    }
}
