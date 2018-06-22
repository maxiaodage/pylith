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
    _boundaryMesh(NULL)
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
    _auxField->subfieldsSetup();
    _auxField->allocate();
    _auxField->zeroLocal();

    assert(_normalizer);
    // pylith::feassemble::AuxiliaryFactory* factory = _auxFactory(); assert(factory);
    // factory->initializeSubfields();

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

/* Calculate stress for 2-D plane strain isotropic linear
 * elasticity WITHOUT a reference stress and strain.
 *
 * Used to output the stress field.
 *
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
 */
void
pylith::bc::DirichletBIE::_stress(const PylithInt dim,
                                  const PylithInt numS,
                                  const PylithInt numA,
                                  const PylithInt sOff[],
                                  const PylithInt sOff_x[],
                                  const PylithScalar s[],
                                  const PylithScalar s_t[],
                                  const PylithScalar s_x[],
                                  const PylithInt aOff[],
                                  const PylithInt aOff_x[],
                                  const PylithScalar a[],
                                  const PylithScalar a_t[],
                                  const PylithScalar a_x[],
                                  const PylithReal t,
                                  const PylithScalar x[],
                                  const PylithInt numConstants,
                                  const PylithScalar constants[],
                                  PylithScalar stress[])
{
    PetscInt     paramOff[3] = {0, 1, 2};
    PylithScalar param[3]    = {0.0, _shearModulus, _bulkModulus};
    pylith::fekernels::IsotropicLinearElasticityPlaneStrain::stress(dim, numS, 3,
        sOff, sOff_x, s, s_t, s_x, paramOff, NULL, param, NULL, NULL, t, x, numConstants, constants, stress);
} // stress

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

    PetscVec stressLocal;
    PetscIS points;
    PetscInt num_points;
    const PetscInt *pts;
    PetscSection section;
    PetscScalar *array;
    const int labelId = 1;
    const PylithInt numConstrained = _constrainedDOF.size();
    const int fieldIndex = solution->subfieldInfo(_field.c_str()).index;
    PetscPointFunc *stressKernels;

    // Calculate stress on the Boundary
    err = DMGetLocalVector(dmSoln, &stressLocal); PYLITH_CHECK_ERROR(err);
    err = PetscCalloc1(solution->subfieldNames().size(), &stressKernels); PYLITH_CHECK_ERROR(err);
    stressKernels[fieldIndex] = _stressKernel;
    err = DMProjectFieldLabelLocal(dmSoln, t, dmLabel, 1, &labelId, numConstrained, &_constrainedDOF[0], solution->localVector(), stressKernels, INSERT_VALUES, stressLocal); PYLITH_CHECK_ERROR(err);
    err = PetscFree(stressKernels); PYLITH_CHECK_ERROR(err);

    // BIE setSolution
    err = DMLabelGetStratumIS(dmLabel,labelId,&points);PYLITH_CHECK_ERROR(err);
    err = ISGetLocalSize(points,&num_points);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(points,&pts);PYLITH_CHECK_ERROR(err);
    err = DMGetDefaultSection(dmSoln,&section);PYLITH_CHECK_ERROR(err);

    // // Get coordiantes
    // PetscDM cdm;
    // PetscVec coordinates;
    // PetscSection csection;
    // PetscInt vStart, vEnd;
    //
    // err = DMGetDepthStratum(dmSoln, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
    // err = DMGetCoordinateDM(dmSoln, &cdm);PYLITH_CHECK_ERROR(err);
    // err = DMGetDefaultSection(cdm, csection);PYLITH_CHECK_ERROR(err);
    // err = DMGetCoordinatesLocal(dmSoln, &coordinates);PYLITH_CHECK_ERROR(err);
    // err = VecGetArray(coordinates,&array);PYLITH_CHECK_ERROR(err);
    // for (PetscInt p=0;p<num_points; ++p)
    // {
    //     const PetscInt point = pts[p];
    //     PetscReal coord[3];
    //     PetscInt dof, off;
    //
    //     if (point < vStart || point >= vEnd) continue;
    //     err = PetscSectionGetDof(section,point,&dof);PYLITH_CHECK_ERROR(err);
    //     err = PetscSectionGetOffset(section,point,&off);PYLITH_CHECK_ERROR(err);
    //     for (PetscInt d=0; d<dof;++d) coord[d] = array[off+d];
    //
    // }
    // err = VecRestoreArray(coordinates,&array);PYLITH_CHECK_ERROR(err);

    //SBIE apply value
    err = VecGetArray(solution->localVector(),&array);PYLITH_CHECK_ERROR(err);
    for (PetscInt p=0;p<num_points; ++p)
    {
        const PetscInt point = pts[p];
        PetscInt dof, off;
        err = PetscSectionGetFieldDof(section,point,fieldIndex,&dof);PYLITH_CHECK_ERROR(err);
        err = PetscSectionGetFieldOffset(section,point,fieldIndex,&off);PYLITH_CHECK_ERROR(err);
        if (!dof) continue;
        assert(numConstrained<=dof);
        for (PetscInt d=0; d<numConstrained;++d) array[off+_constrainedDOF[d]]=100.0;

    }
    err = VecRestoreArray(solution->localVector(),&array);PYLITH_CHECK_ERROR(err);
    err = ISRestoreIndices(points,&pts);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&points);PYLITH_CHECK_ERROR(err);

    err = DMRestoreLocalVector(dmSoln, &stressLocal); PYLITH_CHECK_ERROR(err);

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
// Compute stress values in solution field.
PetscVec
 _computeStress(pylith::topology::Field* solution,
                 const double t)
 {

 }


// ----------------------------------------------------------------------
// Get factory for setting up auxliary fields.
pylith::feassemble::AuxiliaryFactory*
pylith::bc::DirichletBIE::_auxFactory(void)
{ // _auxFactory
    return NULL;
} // auxFactory


// End of file
