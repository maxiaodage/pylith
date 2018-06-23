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

/** @file libsrc/bc/DirichletBIE.hh
 *
 * @brief C++ implementation of Dirichlet using Spectral Boundary integral
 * implementation defined in [Xiao et al, 2018, IJNAMG].
 */

#if !defined(pylith_bc_dirichletbie_hh)
#define pylith_bc_dirichletbie_hh

// Include directives ---------------------------------------------------
#include "BoundaryCondition.hh" // ISA BoundaryCondition
#include "pylith/feassemble/ConstraintPointwise.hh" // ISA ConstraintPointwise
#include "pylith/materials/materialsfwd.hh" // HODSA AuxillaryFactory
#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/petscfwd.h" // USES Vec
#include "pylith/fekernels/IsotropicLinearElasticityPlaneStrain.hh" // USES stress
// Dirichlet ----------------------------------------------------
/// @brief Dirichlet (prescribed values at degrees of freedom) boundary
/// conditions with points on a boundary.
class pylith::bc::DirichletBIE :
    public pylith::bc::BoundaryCondition,
    public pylith::feassemble::ConstraintPointwise {

    friend class DirichletAuxiliaryFactory; // factory for auxiliary fields
    friend class TestDirichletBIE;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    DirichletBIE(void);

    /// Destructor.
    ~DirichletBIE(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    /** Initialize boundary condition.
     *
     * @param[in] solution Solution field.
     */
    void initialize(const pylith::topology::Field& solution);

    /** Set constrained values in solution field.
     *
     * @param[out] solution Solution field.
     * @param[in] t Current time.
     */
    void setSolution(pylith::topology::Field* solution,
                     const double t);



    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    pylith::topology::Mesh* _boundaryMesh;   ///< Boundary mesh.
    pylith::topology::FieldBase::Description _description; ///< Description for constrained field.
    pylith::feassemble::AuxiliaryFactory* _auxFactory(void);

    // PRIVATE METHODS //////////////////////////////////////////////////
private:
    /** Compute stress values on the boundary
     *
     * @param[out] solution Solution field.
     * @param[in] t Current time.
     */
    void  _computeStress(pylith::topology::Field* stress,const pylith::topology::Field& solution,
                                             const double t);

   /** Compute SBIE solution on the boundary
    *
    * @param[out] solution Solution field.
    * @param[in] t Current time.
    */
     void _computeSBIEsolution(PetscVec stressLocal,
                      PetscScalar *array);


     /** Set SBIE solution in solution field.
      *
      * @param[out] solution Solution field.
      * @param[in] t Current time.
      */
     void _setsolutionfromSBIEsolution(pylith::topology::Field* solution,
                      const double t,
                      PetscScalar *array);
    // PRIVATE MEMBERS //////////////////////////////////////////////////
private:

    PylithReal _bulkModulus; ///< Bulk modulus of the far field.
    PylithReal _shearModulus; ///< Shear modulus of the far field.
    PetscPointFunc _stressKernel; ///< Kernel for stress calculation.
    static const char* _pyreComponent; ///< Name of Pyre component.


    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    DirichletBIE(const DirichletBIE&); ///< Not implemented.
    const DirichletBIE& operator=(const DirichletBIE&); ///< Not implemented.
    void
    _stress(const PylithInt dim,
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
            PylithScalar stress[]);

}; // class DirichletBIE

#endif // pylith_bc_dirichletbie_hh


// End of file
