/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    poissonCGAdjointFoam

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // Disable solvers performance output
    lduMatrix::debug = 0;
    solverPerformance::debug = 0;

    // Cost function value
    scalar Jy = 0;
    scalar Ju = 0;
    scalar J = Jy + Ju;

    // Compute cost function value
    #include "costFunctionValue.H"

    // Compute b: right-hand side of Ax = b
    solve(fvm::laplacian(k, yf) + f);
    solve(fvm::laplacian(k, lambda) + yd - yf);
    volScalarField b = lambda;

    // Compute A*f0
    // Primal equation on f to solve for T + Adjoint equation on T to solve for lambda
    solve(fvm::laplacian(k, yu) + u);
    solve(fvm::laplacian(k, lambda) + yu);

    y = yf + yu;
    yua = yu;

    // Compute g0: initial gradient
    volScalarField g = lambda + beta*u - b;
    scalar gL2 = gSum( volField * g.internalField() * g.internalField() );
    scalar gL2a = 0;

    // Compute initial residual and its norm
    volScalarField r = -g;
    scalar rL2 = gL2;

    volScalarField w = g;

    scalar alpha = 0;
    scalar gamma = 0;

    std::ofstream file("results.csv");
    file << 0 << "," << Jy << "," << Ju << "," << J << "," << 0 << nl;

    while (runTime.loop() && ::sqrt( rL2 ) > tol)
    {
	r.dimensions().reset( u.dimensions() );

	// Compute w = A*r
	solve(fvm::laplacian(k, yu) + r);
	solve(fvm::laplacian(k, lambda) + yu);
	w = lambda + beta*r;

	// Update alpha
	alpha = gL2 / gSum( volField * r.internalField() * w.internalField() );

	// Update control
	u = u + alpha*r;

	// Update cost gradient and its norm
	g = g + alpha*w;
	gL2a = gL2;
	gL2 = gSum( volField * g.internalField() * g.internalField() );
	gamma = gL2/gL2a;

	// Update residual and its norm
	r.dimensions().reset( g.dimensions() );
	r = -g + gamma*r;
	rL2 = gSum( volField * r.internalField() * r.internalField() );

	yua = yua + alpha * yu;
	y = yf + yua;

	//solve(fvm::laplacian(k, y) + f + u);

	#include "costFunctionValue.H"

	Info << "Iteration " << runTime.timeName() << " - " \
		<< "Jy " << Jy << " - " \
		<< "Ju " << Ju << " - " \
		<< "J " << J << " - " \
		<< "r_L2 " << ::sqrt( rL2 ) << endl;

	file << runTime.value() << "," << Jy << "," << Ju << "," << J << "," << ::sqrt(rL2) << nl;

	runTime.write();
    }

    runTime++;
    y.write();
    lambda.write();
    u.write();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << nl << endl;
    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
	<< nl << endl;

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
