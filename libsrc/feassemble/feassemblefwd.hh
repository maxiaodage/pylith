// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

/** @file libsrc/feassemble/feassemblefwd.hh
 *
 * @brief Forward declarations for PyLith feassemble objects.
 *
 * Including this header file eliminates the need to use separate
 * forward declarations.
 */

#if !defined(pylith_feassemble_feassemblefwd_hh)
#define pylith_feassemble_feassemblefwd_hh

namespace pylith {
  namespace feassemble {

    class CellGeometry;
    class GeometryPoint1D;
    class GeometryPoint2D;
    class GeometryPoint3D;
    class GeometryLine1D;
    class GeometryLine2D;
    class GeometryLine3D;
    class GeometryTri2D;
    class GeometryTri3D;
    class GeometryQuad2D;
    class GeometryQuad3D;
    class GeometryTet3D;
    class GeometryHex3D;

    class QuadratureBase;
    template<typename mesh_type> class Quadrature;
    template<typename mesh_type> class Quadrature0D;
    template<typename mesh_type> class Quadrature1D;
    template<typename mesh_type> class Quadrature1Din2D;
    template<typename mesh_type> class Quadrature1Din3D;
    template<typename mesh_type> class Quadrature2D;
    template<typename mesh_type> class Quadrature2Din3D;
    template<typename mesh_type> class Quadrature3D;

    class Constraint;
    template<typename quadrature_type> class Integrator;

    class IntegratorElasticity;
    class ElasticityImplicit;
    class ElasticityExplicit;

  } // feassemble
} // pylith


#endif // pylith_feassemble_feassemblefwd_hh


// End of file 
