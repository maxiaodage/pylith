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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/bc/TestDirichletBIE.hh
 *
 * @brief C++ TestDirichletBIE object.
 *
 * C++ unit testing for DirichletBIEBC.
 */

#if !defined(pylith_bc_testdirichletbie_hh)
#define pylith_bc_testdirichletbie_hh

#include <cppunit/extensions/HelperMacros.h>
#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/bc/bcfwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // forward declarations

#include "pylith/topology/Field.hh" // HOLDSA Field

#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA UserFunctionDB
#include "spatialdata/geocoords/geocoordsfwd.hh" // HOLDSA CoordSys
#include "spatialdata/units/unitsfwd.hh" // HOLDSA Nondimensional


/// Namespace for pylith package
namespace pylith {
    namespace bc {
        class TestDirichletBIE;
        class TestDirichletBIE_Data;
    } // bc
} // pylith

/// C++ unit testing for DirichletBC.
class pylith::bc::TestDirichletBIE : public CppUnit::TestFixture, public pylith::utils::GenericComponent {

    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestDirichletBIE);

    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testAccessors);
    //CPPUNIT_TEST(testAuxFieldDiscretization);
    //CPPUNIT_TEST(testAuxFieldDB);
    CPPUNIT_TEST(testNormalizer);
    CPPUNIT_TEST(testVerifyConfiguration);
    CPPUNIT_TEST(testInitialize);
    CPPUNIT_TEST(testPrestep);
//    CPPUNIT_TEST(testSetSolution);
    // Adding a test for computeStress()
    CPPUNIT_TEST(testcomputeStress);
    //CPPUNIT_TEST(testAuxFieldSetup);

    CPPUNIT_TEST_SUITE_END_ABSTRACT();

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Tear down testing data.
    void tearDown(void);

    /// Test constructor.
    void testConstructor(void);

    /// Test accessors (field, constrainedDOF, dbTimeHistory, useInitial, useRate, useTimeHistory).
    void testAccessors(void);

    /// Test auxFieldDiscretization().
    //void testAuxFieldDiscretization(void);

    /// Test auxFieldDB().
    //void testAuxFieldDB(void);

    /// Test normalizer().
    void testNormalizer(void);

    /// Test verifyConfiguration().
    void testVerifyConfiguration(void);

    /// Test initialize().
    void testInitialize(void);

    /// Test prestep().
    void testPrestep(void);

    /// Test setSolution().
    void testSetSolution(void);

    /// Test _computeStress().
    void testcomputeStress(void);

    /// Test _auxFieldsSetup().
    //void testAuxFieldSetup(void);

    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    TestDirichletBIE_Data* _data; ///< Data for testing

    pylith::bc::DirichletBIE* _bc; /// Test subject.
    pylith::topology::Mesh* _mesh; /// Mesh used in testing.
    pylith::topology::Field* _solution; ///< Solution field used in testing.

    static const double FILL_VALUE; ///< Fill value for unconstrained values.

    // PRIVATE METHODS ////////////////////////////////////////////////////
private:

    /// Initializer boundary condition for testing.
    void _initialize(void);

    /// Setup solution field.
    void _setupSolutionField(void);

}; // class TestDirichletTimeDependent


// ======================================================================
class pylith::bc::TestDirichletBIE_Data {

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    TestDirichletBIE_Data(void);

    /// Destructor
    ~TestDirichletBIE_Data(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    const char* meshFilename; ///< Name of file with ASCII mesh.
    const char* bcLabel; ///< Label defining cells associated with material.

    spatialdata::geocoords::CoordSys* cs; ///< Coordinate system.
    spatialdata::units::Nondimensional* normalizer; ///< Scales for nondimensionalization.

    const char* field; ///< Name of solution field constrained.
    pylith::topology::FieldBase::VectorFieldEnum vectorFieldType; ///< Vector field type for constrained field.
    PylithReal scale; ///< Scale of constrained field.
    int numConstrainedDOF; ///< Number of constrained DOF;
    const int* constrainedDOF; ///< Array of constrained DOF.

    int numAuxSubfields; ///< Number of auxiliary subfields.
    const char** auxSubfields; ///< Names of auxiliary subfields.
    pylith::topology::Field::Discretization* auxDiscretizations; ///< Discretizations for auxiliary fields.
    spatialdata::spatialdb::UserFunctionDB* auxDB; ///< Spatial database with auxiliary field.

    spatialdata::spatialdb::UserFunctionDB* tractionDB; ///< Spatial database with traction field.


    PylithReal t; ///< Time associated with setting solution.
    PylithReal dt; ///< Time step associated with setting solution.
    int solnNumSubfields; ///< Number of solution fields.
    pylith::topology::FieldBase::Discretization* solnDiscretizations; ///< Discretizations for solution fields.
    spatialdata::spatialdb::UserFunctionDB* solnDB; ///< Spatial database with solution.

}; // class TestDirichletBIE_Data


#endif // pylith_bc_dirichletbie_hh


// End of file
