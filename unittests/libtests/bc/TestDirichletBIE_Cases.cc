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

#include <portinfo>

#include "TestDirichletBIE.hh" // Implementation of cases

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional


// ----------------------------------------------------------------------
namespace pylith {
    namespace bc {

        // --------------------------------------------------------------
        class TestDirichletBIE_InitialDisp2D : public TestDirichletBIE {

protected:

            static const char* disp_units(void) {
                return "m";
            }
            static const char* vel_units(void) {
                return "m/s";
            }
            static const char* pressure_units(void) {
                return "Pa";
            }

            // Velocity and fluid pressure solution fields at time t.

            static double vel_x(const double x,
                                const double y) {
                return FILL_VALUE;
            } // vel_x
            static double vel_y(const double x,
                                const double y) {
                return FILL_VALUE;
            } // vel_y
            static double fluid_press(const double x,
                                      const double y) {
                return FILL_VALUE;
            } // fluid_press

            void setUp(void) {
                TestDirichletBIE::setUp();
                _data = new TestDirichletBIE_Data();CPPUNIT_ASSERT(_data);
                _data->cs = new spatialdata::geocoords::CSCart();CPPUNIT_ASSERT(_data->cs);
                _data->cs->setSpaceDim(2);
                _data->cs->initialize();

                CPPUNIT_ASSERT(_data->normalizer);
                _data->normalizer->lengthScale(1000.0);
                _data->normalizer->timeScale(10.0);
                _data->normalizer->pressureScale(0.1);
                _data->normalizer->densityScale(2.0);

                _data->field = "displacement";
                _data->vectorFieldType = pylith::topology::Field::VECTOR;
                _data->scale = _data->normalizer->lengthScale();
                _data->numConstrainedDOF = 2;
                static const int constrainedDOF[2] = { 0, 1 };
                _data->constrainedDOF = const_cast<int*>(constrainedDOF);

            } // setUp
        }; // class TestDirichletTimeDependent_TriP1


        // --------------------------------------------------------------
        class TestDirichletBIE_TriP1 : public TestDirichletBIE_InitialDisp2D {

            // Spatial database user functions for auxiliary subfields.



            // Displacement solution field at time t.

            static double disp_x(const double x,
                                 const double y) {
                return 100.0e+3;
            } // disp_x
            static double disp_y(const double x,
                                 const double y) {
                return 100.0e+3;
            } // disp_y

protected:

            CPPUNIT_TEST_SUB_SUITE(TestDirichletBIE_TriP1, TestDirichletBIE_InitialDisp2D);
            CPPUNIT_TEST_SUITE_END();
            void setUp(void) {
                TestDirichletBIE_InitialDisp2D::setUp();
                _data->meshFilename = "data/tri_small.mesh";
                _data->bcLabel = "boundary_bottom";

                _data->t = 1.23;
                _data->dt = 0.1;
                _data->solnNumSubfields = 3;
                static const pylith::topology::Field::Discretization solnDiscretizations[3] = {
                    {1, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // displacement
                    {1, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // velocity
                    {1, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // fluid_pressure
                };
                _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(solnDiscretizations);

                CPPUNIT_ASSERT(_data->solnDB);
                _data->solnDB->coordsys(*_data->cs);
                _data->solnDB->addValue("displacement_x", disp_x, disp_units());
                _data->solnDB->addValue("displacement_y", disp_y, disp_units());
                _data->solnDB->addValue("velocity_x", vel_x, vel_units());
                _data->solnDB->addValue("velocity_y", vel_y, vel_units());
                _data->solnDB->addValue("fluid_pressure", fluid_press, pressure_units());

            } // setUp
        }; // class TestDirichletTimeDependent_TriP1
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDirichletBIE_TriP1);

        // --------------------------------------------------------------
        class TestDirichletBIE_TriP2 : public TestDirichletBIE_InitialDisp2D {

            // Spatial database user functions for auxiliary subfields.



            // Displacement solution field at time t.

            static double disp_x(const double x,
                                 const double y) {
                return 100.0e+3;
            } // disp_x
            static double disp_y(const double x,
                                 const double y) {
                return 100.0e+3;
            } // disp_y

protected:

            CPPUNIT_TEST_SUB_SUITE(TestDirichletBIE_TriP2, TestDirichletBIE_InitialDisp2D);
            CPPUNIT_TEST_SUITE_END();
            void setUp(void) {
                TestDirichletBIE_InitialDisp2D::setUp();
                _data->meshFilename = "data/tri_small.mesh";
                _data->bcLabel = "boundary_bottom";

                static const pylith::topology::Field::Discretization auxDiscretizations[1] = {
                    {2, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // initial_amplitude
                };
                _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(auxDiscretizations);

                _data->t = 1.23;
                _data->dt = 0.1;
                _data->solnNumSubfields = 3;
                static const pylith::topology::Field::Discretization solnDiscretizations[3] = {
                    {2, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // displacement
                    {2, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // velocity
                    {2, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // fluid_pressure
                };
                _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(solnDiscretizations);

                CPPUNIT_ASSERT(_data->solnDB);
                _data->solnDB->coordsys(*_data->cs);
                _data->solnDB->addValue("displacement_x", disp_x, disp_units());
                _data->solnDB->addValue("displacement_y", disp_y, disp_units());
                _data->solnDB->addValue("velocity_x", vel_x, vel_units());
                _data->solnDB->addValue("velocity_y", vel_y, vel_units());
                _data->solnDB->addValue("fluid_pressure", fluid_press, pressure_units());

            } // setUp
        }; // class TestDirichletTimeDependent_TriP2
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDirichletBIE_TriP2);

#if 0
        // --------------------------------------------------------------
        class TestDirichletTimeDependent_QuadP1 : public TestDirichletTimeDependent {
            CPPUNIT_TEST_SUB_SUITE(TestDirichletTimeDependent_QuadP1, TestDirichletTimeDependent);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestDirichletTimeDependent::setUp();
                _data = new TestDirichletTimeDependent_Data();CPPUNIT_ASSERT(_data);
                _data->meshFilename = "data/quad_small.mesh";
                _data->bcLabel = "boundary_top";
            } // setUp
        }; // class TestDirichletTimeDependent_QuadP1
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDirichletTimeDependent_QuadP1);


        // --------------------------------------------------------------
        class TestDirichletTimeDependent_TetP1 : public TestDirichletTimeDependent {
            CPPUNIT_TEST_SUB_SUITE(TestDirichletTimeDependent_TetP1, TestDirichletTimeDependent);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestDirichletTimeDependent::setUp();
                _data = new TestDirichletTimeDependent_Data();CPPUNIT_ASSERT(_data);
            } // setUp
        }; // class TestDirichletTimeDependent_TetP1
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDirichletTimeDependent_TetP1);


        // --------------------------------------------------------------
        class TestDirichletTimeDependent_Hex : public TestDirichletTimeDependent {
            CPPUNIT_TEST_SUB_SUITE(TestDirichletTimeDependent_Hex, TestDirichletTimeDependent);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestDirichletTimeDependent::setUp();
                _data = new TestDirichletTimeDependent_Data();CPPUNIT_ASSERT(_data);
            } // setUp
        }; // class TestDirichletTimeDependent_Hex
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDirichletTimeDependent_Hex);
#endif

    } // namespace bc
} // namespace pylith


// End of file
