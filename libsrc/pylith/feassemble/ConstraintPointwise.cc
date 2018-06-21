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

#include "pylith/feassemble/ConstraintPointwise.hh" // implementation of object methods

#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/meshio/OutputManager.hh" // USES OutputManager
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()
#include <stdexcept> \
    // USES std::runtime_error

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::ConstraintPointwise::ConstraintPointwise(void) :
    _normalizer(new spatialdata::units::Nondimensional),
    _auxField(NULL),
    _output(NULL),
    _logger(NULL)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::feassemble::ConstraintPointwise::~ConstraintPointwise(void)
{ // destructor
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::ConstraintPointwise::deallocate(void)
{ // deallocate
    PYLITH_METHOD_BEGIN;

    delete _normalizer; _normalizer = NULL;
    delete _logger; _logger = NULL;
    delete _auxField; _auxField = NULL;

    _output = NULL; // :KLUDGE: Memory managed by Python object. :TODO: Use shared pointer.

    PYLITH_METHOD_END;
} // deallocate


// ----------------------------------------------------------------------
// Set indices of constrained degrees of freedom at each location.
void
pylith::feassemble::ConstraintPointwise::constrainedDOF(const int* flags,
                                                        const int size)
{ // constrainedDOF
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("constrainedDOF(flags="<<flags<<", size="<<size<<")");

    assert((size > 0 && flags) || (!size && !flags));

    _constrainedDOF.resize(size);
    for (int i = 0; i < size; ++i) {
        if (flags[i] < 0) {
            std::ostringstream msg;
            msg << "Constrained DOF '" << flags[i] << "' must be nonnegative in component '" << PyreComponent::identifier() << "'.";
            throw std::runtime_error(msg.str());
        } // if
        _constrainedDOF[i] = flags[i];
    } // for

    PYLITH_METHOD_END;
} // constrainedDOF


// ----------------------------------------------------------------------
// Get indices of constrained degrees of freedom.
const pylith::int_array&
pylith::feassemble::ConstraintPointwise::constrainedDOF(void) const
{ // constrainedDOF
    return _constrainedDOF;
} // constrainedDOF

// ----------------------------------------------------------------------
// Return auxiliary subfields for this problem.
const pylith::topology::Field&
pylith::feassemble::ConstraintPointwise::auxField(void) const
{ // auxField
    PYLITH_METHOD_BEGIN;

    assert(_auxField);

    PYLITH_METHOD_RETURN(*_auxField);
} // auxField

// ----------------------------------------------------------------------
// Set database for filling auxiliary subfields.
void
pylith::feassemble::ConstraintPointwise::auxFieldDB(spatialdata::spatialdb::SpatialDB* value)
{ // auxFieldDB
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("auxFieldDB(value="<<value<<")");


    pylith::feassemble::AuxiliaryFactory* factory = _auxFactory();
    if (!factory) {
        PYLITH_COMPONENT_ERROR("Constraint does not contain an auxiliary factory.");
        throw std::logic_error("Constraint does not contain an auxiliary factory.");
    } else {
        assert(factory);
        factory->queryDB(value);
    } // if/else
    PYLITH_METHOD_END;
} // auxFieldDB

// ----------------------------------------------------------------------
// Set discretization information for auxiliary subfield.
void
pylith::feassemble::ConstraintPointwise::auxSubfieldDiscretization(const char* name,
                                                                   const int basisOrder,
                                                                   const int quadOrder,
                                                                   const bool isBasisContinuous,
                                                                   const pylith::topology::FieldBase::SpaceEnum feSpace)
{ // auxSubfieldDiscretization
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("auxSubfieldDiscretization(name="<<name<<", basisOrder="<<basisOrder<<", quadOrder="<<quadOrder<<", isBasisContinuous="<<isBasisContinuous<<")");

    pylith::feassemble::AuxiliaryFactory* factory = _auxFactory();
    if (!factory) {
        PYLITH_COMPONENT_ERROR("Constraint does not contain an auxiliary factory.");
        throw std::logic_error("Constraint does not contain an auxiliary factory.");
    } else {
        assert(factory);
        factory->subfieldDiscretization(name, basisOrder, quadOrder, isBasisContinuous, feSpace);
    } // if/else

    PYLITH_METHOD_END;
} // auxSubfieldDiscretization


// ----------------------------------------------------------------------
// Set manager of scales used to nondimensionalize problem.
void
pylith::feassemble::ConstraintPointwise::normalizer(const spatialdata::units::Nondimensional& dim)
{ // normalizer
    PYLITH_COMPONENT_DEBUG("normalizer(dim="<<typeid(dim).name()<<")");

    if (!_normalizer) {
        _normalizer = new spatialdata::units::Nondimensional(dim);
    } else {
        *_normalizer = dim;
    } // if/else
} // normalizer


// ----------------------------------------------------------------------
// Set output manager.
void
pylith::feassemble::ConstraintPointwise::output(pylith::meshio::OutputManager* manager) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("output(manager="<<manager<<")");

    _output = manager;

    PYLITH_METHOD_END;
} // output


// ----------------------------------------------------------------------
// Write information (auxiliary field) output.
void
pylith::feassemble::ConstraintPointwise::writeInfo(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("writeInfo(void)");

    if (_output) {
        assert(_auxField);
        _output->writeInfo(*_auxField);
    } // if

    PYLITH_METHOD_END;
} // writeInfo


// ----------------------------------------------------------------------
// Update auxiliary fields at beginning of time step.
void
pylith::feassemble::ConstraintPointwise::prestep(const double t,
                                                 const double dt)
{ // prestep
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("prestep(t="<<t<<", dt="<<dt<<") empty method");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // prestep


// End of file
