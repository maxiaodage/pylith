// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/bc/DirichletBIE.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/materials/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream>

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::DirichletBIE::DirichletBIE(void) :
    _boundaryMesh(NULL),
    _auxMaterialFactory(new pylith::materials::AuxiliaryFactory)
{ // constructor
    _description.label = "unknown";
    _description.vectorFieldType = pylith::topology::FieldBase::OTHER;
    _description.numComponents = 0;
    _description.scale = 1.0;
    _description.validator = NULL;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::DirichletBIE::~DirichletBIE(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::DirichletBIE::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    ConstraintPointwise::deallocate();

    delete _boundaryMesh; _boundaryMesh = NULL;
    delete _auxMaterialFactory; _auxMaterialFactory = NULL;

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::DirichletBIE::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.label()<<")");

    if (!solution.hasSubfield(_field.c_str())) {
        std::ostringstream msg;
        msg << "Cannot constrain field '"<< _field
            << "' in component '" << PyreComponent::identifier() << "'"
            << "; field is not in solution.";
        throw std::runtime_error(msg.str());
    } // if

    const topology::Field::SubfieldInfo& info = solution.subfieldInfo(_field.c_str());
    const int numComponents = info.description.numComponents;
    const int numConstrained = _constrainedDOF.size();
    for (int iConstrained = 0; iConstrained < numConstrained; ++iConstrained) {
        if (_constrainedDOF[iConstrained] >= numComponents) {
            std::ostringstream msg;
            msg << "Cannot constrain degree of freedom '" << _constrainedDOF[iConstrained] << "'"
                << " in component '" << PyreComponent::identifier() << "'"
                << "; solution field '" << _field << "' contains only " << numComponents << " components.";
            throw std::runtime_error(msg.str());
        } // if
    } // for

    PYLITH_METHOD_END;
} // verifyConfiguration


// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::bc::DirichletBIE::initialize(const pylith::topology::Field& solution)
{ // initialize
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("initialize(solution="<<solution.label()<<")");

    const topology::Field::SubfieldInfo& info = solution.subfieldInfo(_field.c_str());
    _description = info.description;

    delete _boundaryMesh; _boundaryMesh = new pylith::topology::Mesh(solution.mesh(), _label.c_str()); assert(_boundaryMesh);
    PetscDM dmBoundary = _boundaryMesh->dmMesh(); assert(dmBoundary);
    pylith::topology::CoordsVisitor::optimizeClosure(dmBoundary);

    delete _auxField; _auxField = new pylith::topology::Field(*_boundaryMesh); assert(_auxField);
    _auxField->label("DirichletBIE auxiliary");
    _auxFieldSetup(solution);
    _auxField->subfieldsSetup();
    _auxField->allocate();
    _auxField->zeroLocal();

    assert(_normalizer);
    pylith::feassemble::AuxiliaryFactory* factory = _auxFactory(); assert(factory);
    factory->initializeSubfields();

    //_auxField->view("AUXILIARY FIELD"); // :DEBUG: TEMPORARY
    writeInfo();

    const PetscDM dmSoln = solution.dmMesh(); assert(dmSoln);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dmSoln, &prob); PYLITH_CHECK_ERROR(err);

    void* context = NULL;
    const int labelId = 1;
    const PylithInt numConstrained = _constrainedDOF.size();
    err = PetscDSAddBoundary(prob, DM_BC_ESSENTIAL, label(), label(), info.index, numConstrained, &_constrainedDOF[0],
                             NULL, 1, &labelId, context); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Set constrained values in solution field.
void
pylith::bc::DirichletBIE::setSolution(pylith::topology::Field* solution,
                                   const double t)
{ // setSolution
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setSolution(solution="<<solution->label()<<", t="<<t<<")");

    assert(solution);
    assert(_auxField);

    PetscErrorCode err;
    PetscDM dmSoln = solution->dmMesh();
    PetscDM dmAux = _auxField->dmMesh();

    // Get label for constraint.
    PetscDMLabel dmLabel;
    err = DMGetLabel(dmSoln, _label.c_str(), &dmLabel); PYLITH_CHECK_ERROR(err);

    // Set auxiliary data
    err = PetscObjectCompose((PetscObject) dmSoln, "dmAux", (PetscObject) dmAux); PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmSoln, "A", (PetscObject) _auxField->localVector()); PYLITH_CHECK_ERROR(err);

    IS points;
    PetscInt num_points;
    const PetscInt *pts;
    PetscSection section;
    PetscScalar *array;
    const int labelId = 1;
    const PylithInt numConstrained = _constrainedDOF.size();
    const int fieldIndex = solution->subfieldInfo(_field.c_str()).index;

    err = DMLabelGetStratumIS(dmLabel,labelId,&points);PYLITH_CHECK_ERROR(err);
    err = ISGetLocalSize(points,&num_points);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(points,&pts);PYLITH_CHECK_ERROR(err);
    err = DMGetDefaultSection(dmSoln,&section);PYLITH_CHECK_ERROR(err);
    err = VecGetArray(solution->localVector(),&array);PYLITH_CHECK_ERROR(err);
    for (PetscInt p=0;p<num_points; ++p)
    {
        const PetscInt point = pts[p];
        PetscInt dof, off;
        err = PetscSectionGetFieldDof(section,point,fieldIndex,&dof);PYLITH_CHECK_ERROR(err);
        err = PetscSectionGetFieldOffset(section,point,fieldIndex,&off);PYLITH_CHECK_ERROR(err);
        if (!dof) continue;
        assert(numConstrained<=dof);
        for (PetscInt d=0; d<numConstrained;++d) array[off+_constrainedDOF[d]]=1.0;

    }
    err = VecRestoreArray(solution->localVector(),&array);PYLITH_CHECK_ERROR(err);
    err = ISRestoreIndices(points,&pts);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&points);PYLITH_CHECK_ERROR(err);



#if 0
    void* context = NULL;
    const int labelId = 1;
    const int fieldIndex = solution->subfieldInfo(_field.c_str()).index;
    const PylithInt numConstrained = _constrainedDOF.size();
    assert(solution->localVector());
    err = DMPlexLabelAddCells(dmSoln, dmLabel); PYLITH_CHECK_ERROR(err);
    err = DMPlexInsertBoundaryValuesEssentialField(dmSoln, t,
                                                   solution->localVector(), fieldIndex, numConstrained, &_constrainedDOF[0], dmLabel, 1, &labelId, _bcKernel,
                                                   context, solution->localVector());
    PYLITH_CHECK_ERROR(err);
    err = DMPlexLabelClearCells(dmSoln, dmLabel); PYLITH_CHECK_ERROR(err);
#endif

    //solution->view("SOLUTION at end of setSolution()"); // :DEBUG: TEMPORARY

    PYLITH_METHOD_END;
} // setSolution

// ----------------------------------------------------------------------
// Setup auxiliary subfields (discretization and query fns).
void
pylith::bc::DirichletBIE::_auxFieldSetup(const pylith::topology::Field& solution)
{ // _auxFieldsSetup
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_auxFieldsSetup(solution="<<solution.label()<<")");

    assert(_auxMaterialFactory);
    assert(_normalizer);
    _auxMaterialFactory->initialize(_auxField, *_normalizer, solution.spaceDim(),
                                         &solution.subfieldInfo(_field.c_str()).description);

    // :ATTENTION: The order of the factory methods must match the order of the auxiliary subfields in the FE kernels.
    _auxMaterialFactory->density(); // 0
    _auxMaterialFactory->shearModulus(); // 1
    _auxMaterialFactory->bulkModulus(); // 2

    PYLITH_METHOD_END;
}     // _auxFieldSetup

// ----------------------------------------------------------------------
// Get factory for setting up auxliary fields.
pylith::feassemble::AuxiliaryFactory*
pylith::bc::DirichletBIE::_auxFactory(void)
{ // _auxFactory
    return _auxMaterialFactory;
} // auxFactory


// End of file
