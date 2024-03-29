    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info<<"Reading tubulentSchmidt number\n"<<endl;
    dimensionedScalar Sct
    (
        transportProperties.lookup("Sct")
    );

    Info<<"Reading molecular mass diffusivity \n"<<endl;
    dimensionedScalar D
    (
        transportProperties.lookup("D")
    );

    Info<< "Reading field C\n" << endl;

    volScalarField C
    (
        IOobject
        (
            "C",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field k\n" << endl;

    volScalarField k
    (
        IOobject
        (
            "k",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Reading field epsilon\n" << endl;

    volScalarField epsilon
    (
        IOobject
        (
            "epsilon",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Reading field R\n" << endl;

    volSymmTensorField R
    (
        IOobject
        (
            "R",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );


Info<< "Reading the turbulent diffusivity"<<endl;
    volSymmTensorField Dt
    (
        IOobject
        (
            "Dt",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
	dimensionedSymmTensor( "Dt",dimensionSet(0,2,-1,0,0,0,0), symmTensor::zero), 
	"zeroGradient"
    );

    Dt = ((k/(Sct*epsilon)))*R;


    Info<< "Reading field U\n" << endl;

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

#   include "createPhi.H"

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );
