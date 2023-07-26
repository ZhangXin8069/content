// $Id: simple_baryon_seqsrc_w.cc,v 3.7 2009-03-19 17:17:20 mcneile Exp $
/*! \file
 *  \brief Construct baryon sequential sources.
 */

#include "qdp_config.h"

#include "meas/hadron/simple_baryon_seqsrc_w.h"
#include "simple_baryon_seqsrc_w.h"
#include "meas/hadron/seqsource_factory_w.h"
#include "meas/hadron/barspinmat_w.h"
#include "meas/hadron/barhqlq_w.h"
#include "meas/hadron/baryon_spinmat_funcmap_w.h"
#include "util/ft/sftmom.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{

  //! Baryon sequential sources
  /*! \ingroup hadron */
  namespace SimpleBaryonSeqSourceEnv
  { 

    //! Anonymous namespace
    namespace
    {

       void checkArgs(const char* name, const multi1d<LatticePropagator>& quark_propagators)
       {
         if (quark_propagators.size() != 1)
         {
//         QDPIO::cerr << __func__ << ": expect only 1 prop" << std::endl;
//         QDP_abort(1);

           std::ostringstream s;
           s << name << ": expecting 1 prop, instead passed = " << quark_propagators.size() << std::endl;
           throw s.str();
         }
       }

       // For baryon case: 2 spectator
       void checkArgs2(const char* name, const multi1d<LatticePropagator>& quark_propagators)
       {
         if (quark_propagators.size() != 2)
         {
           std::ostringstream s;
           s << name << ": expecting 2 prop, instead passed = " << quark_propagators.size() << std::endl;
           throw s.str();
         }
       }
    }

//! Nucleon-Nucleon ISO piece with general projector and Cg5
    LatticePropagator
    WeakCurrentgAgV::operator()(const multi1d<LatticeColorMatrix>& u,
                             const multi1d<ForwardProp_t>& forward_headers,
                             const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();
      if ( Nc != 3 ){    /* Code is specific to Ns=4 and Nc=3. */
        QDPIO::cerr<<" code only works for Nc=3 and Ns=4\n";
        QDP_abort(111) ;
      }
#if QDP_NC == 3

      checkArgs2("WeakCurrentgAgV", quark_propagators);

      LatticePropagator src_prop_tmp;

      src_prop_tmp += T * transposeSpin( quarkContract13(quark_propagators[0] * Cg5, Cg5 * quark_propagators[1]) );
      src_prop_tmp += T * quarkContract24(quark_propagators[0] * Cg5, Cg5 * quark_propagators[1]);
      src_prop_tmp += quarkContract13(quark_propagators[0] * Cg5, Cg5 * quark_propagators[1]) * T;
      src_prop_tmp += quarkContract14(quark_propagators[0] * Cg5, Cg5 * quark_propagators[1] * T);


      END_CODE();

      return projectBaryon(src_prop_tmp,
                           forward_headers);
#else
      LatticePropagator q1_tmp;

      q1_tmp = zero ;
      return q1_tmp ;
#endif
    }

    Complex
    WeakCurrentgAgV::twoPtSink(const multi1d<LatticeColorMatrix>& u,
                            const multi1d<ForwardProp_t>& forward_headers,
                            const multi1d<LatticePropagator>& quark_propagators,
                            int gamma_insertion)
    {
      //checkArgs("BarNucl_ISO_TCg5", quark_propagators);
      //setTSrce(forward_headers);
      //setBC(forward_headers);

      // Constructor the 2pt
      //LatticeComplex b_prop = Baryon2PtContractions::sigma2pt(quark_propagators[0],
                                                              //quark_propagators[0],
                                                              //T, Cg5);

      // Extract the sink at the appropriate momenta
      //SftMom sft(0, getTSrce(), getSinkMom(), false, getDecayDir());
      //multi2d<DComplex> hsum;
      //hsum = sft.sft(b_prop);

      // D-quark rescaling factor is 1 for this type of seqsource
      //return Real(1) * hsum[0][getTSink()];
      return 0;
      
    }



    //! Anonymous namespace
    namespace
    {

      HadronSeqSource<LatticePropagator>* WeakCurrentgAgV_cur(XMLReader& xml_in,
							  const std::string& path)
      {
        SpinMatrix g_one = 1.0;        
        return new WeakCurrentgAgV(SimpleBaryonSeqSourceEnv::Params(xml_in, path), 
                   SpinMatrix(0.5*(g_one + Gamma(8)*g_one)*(g_one + Gamma(11)*g_one)),
                   BaryonSpinMats::Cg5());//BaryonSpinMats::Cg5()
      }

      
      bool registered_extra = false;

    }// end anonymous namespace

    //! Register all the factories
    bool registerAll_extra() 
    {
      bool success = true; 
      if (! registered_extra)
      {
	//! Register needed stuff

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(std::string("WEAK_CURRENT_gAgV_CUR"),
										      WeakCurrentgAgV_cur);

	registered_extra = true;
      }
      return success;
    }

  } // namespace BaryonSeqSourceCallMapEnv


}  // end namespace Chroma

