// -*- C++ -*-
// $Id: simple_baryon_seqsrc_w.h,v 3.3 2006-11-28 20:00:49 edwards Exp $
/*! \file
 *  \brief Construct baryon sequential sources.
 */

#ifndef __simple_baryon_seqsrc_w_h_v2__
#define __simple_baryon_seqsrc_w_h_v2__

#include "meas/hadron/simple_baryon_seqsrc_w.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup hadron */
  namespace SimpleBaryonSeqSourceEnv
  {
    bool registerAll_extra();

    class WeakCurrentgAgV : public BaryonSeqSourceBase
    {
    public:
      //! Full constructor
      WeakCurrentgAgV(const Params& p, const SpinMatrix& spinT, const SpinMatrix& spinCg5) :
        params(p), T(spinT), Cg5(spinCg5) {}

      //! Default destructor
      ~WeakCurrentgAgV() {}

      //! Construct the source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
                                   const multi1d<ForwardProp_t>& forward_headers,
                                   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
                        const multi1d<ForwardProp_t>& forward_headers,
                        const multi1d<LatticePropagator>& forward_props,
                        int gamma_insertion);

    protected:
      //! Set bc
      multi1d<int>& getBC() {return bc;}

      //! Get bc
      const multi1d<int>& getBC() const {return bc;}

      //! Set t_srce
      multi1d<int>& getTSrce() {return t_srce;}

      //! Get t_srce
      const multi1d<int>& getTSrce() const {return t_srce;}

      //! Get t_sink
      int getTSink() const {return params.t_sink;}

      //! Get sink_mom
      const multi1d<int>& getSinkMom() const {return params.sink_mom;}

      //! Get decay_dir
      const int getDecayDir() const {return params.j_decay;}

    private:
      //! Hide partial constructor
      WeakCurrentgAgV() {}

    private:
      multi1d<int>  t_srce;   /*<! Must come from propagator headers */
      multi1d<int>  bc;       /*<! Must come from propagator headers */
      Params  params;   /*!< Seqsource params */
      SpinMatrix T;     /*!< The spin projector matrix */
      SpinMatrix Cg5;   /*!< The Cg5 at the source and sink */
    };


  }  // end namespace

}  // end namespace Chroma


#endif
