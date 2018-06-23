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
            static const char* stress_units(void) {
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
        };
        // --------------------------------------------------------------
        class TestDirichletBIE_QuadP1 : public TestDirichletBIE_InitialDisp2D {

            // Spatial database user functions for auxiliary subfields.
            // Displacement solution field at time t.

            static double disp_x(const double x,
                                 const double y) {
                return x;
            } // disp_x
            static double disp_y(const double x,
                                 const double y) {
                return y;
            } // disp_y
            static double stress_xx(const double x,
                                      const double y) {
                return 4.0;

            } // stress_xx
            static double stress_yy(const double x,
                                      const double y) {
                return 4.0;
            } // stress_yy
            static double stress_xy(const double x,
                                      const double y) {
                return 0.0;
            } // stress_xy
            static double stress_zz(const double x,
                                      const double y) {
                return 2;
            } // stress_zz



protected:

            CPPUNIT_TEST_SUB_SUITE(TestDirichletBIE_QuadP1, TestDirichletBIE_InitialDisp2D);
            CPPUNIT_TEST_SUITE_END();
            void setUp(void) {
                TestDirichletBIE_InitialDisp2D::setUp();
                _data->meshFilename = "data/quad_small.mesh";
                _data->bcLabel = "boundary_top";
                _data->t = 1.23;
                _data->dt = 0.1;
                _data->solnNumSubfields = 2;
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

                CPPUNIT_ASSERT(_data->tractionDB);
                _data->tractionDB->coordsys(*_data->cs);
                _data->tractionDB->addValue("stress_xx", stress_xx, stress_units());
                 _data->tractionDB->addValue("stress_yy", stress_yy, stress_units());
                 _data->tractionDB->addValue("stress_xy", stress_xy, stress_units());
                 _data->tractionDB->addValue("stress_zz", stress_zz, stress_units());

            } // setUp
        }; // class TestDirichletBIE_QuadP1
        CPPUNIT_TEST_SUITE_REGISTRATION(TestDirichletBIE_QuadP1);





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
