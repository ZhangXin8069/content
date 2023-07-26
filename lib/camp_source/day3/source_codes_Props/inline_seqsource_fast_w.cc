// -*- C++ -*-
/*! \file
 * \brief Inline construction of sequential sources
 *
 * Sequential source construction
 */

#include "handle.h"
#include "inline_seqsource_fast_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/hadron/seqsource_factory_w.h"
#include "meas/hadron/seqsource_aggregate_w.h"
#include "meas/sinks/sink_smearing_factory.h"
#include "meas/sinks/sink_smearing_aggregate.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "util/info/unique_id.h"
#include "util/ferm/transf.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  //! Propagator input
  void read(XMLReader& xml, const std::string& path, InlineSeqSourceFastEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "prop_ids", input.prop_ids);
    read(inputtop, "seqsource_id", input.seqsource_id);

    // single element case without <elem>
    if(input.seqsource_id.size() == 0)
    {
      std::string m_seqsrc_id;
      read(inputtop, "seqsource_id", m_seqsrc_id);
      input.seqsource_id.resize(1);
      input.seqsource_id[0] = m_seqsrc_id;
    }
  }

  //! Propagator output
  void write(XMLWriter& xml, const std::string& path, const InlineSeqSourceFastEnv::Params::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "prop_ids", input.prop_ids);
    write(xml, "seqsource_id", input.seqsource_id);

    pop(xml);
  }


  namespace InlineSeqSourceFastEnv 
  { 
    // Anonymous namespace
    namespace TSF
    {
      // Standard Time Slicery
      class TimeSliceFunc : public SetFunc
      {
      public:
        TimeSliceFunc(int dir): dir_decay(dir) {}
        int operator() (const multi1d<int>& coordinate) const {return coordinate[dir_decay];}
        int numSubsets() const {return Layout::lattSize()[dir_decay];}
        int dir_decay;
      private:
        TimeSliceFunc() {}  // hide default constructor
      };
    }

    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
          const std::string& path) 
      {
        return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "SEQSOURCE_FAST";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {	
        success &= QuarkSinkSmearingEnv::registerAll();
        success &= HadronSeqSourceEnv::registerAll();
        success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
        registered = true;
      }
      return success;
    }


    // Param stuff
    Params::Params() { frequency = 0; }

    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      try 
      {
        XMLReader paramtop(xml_in, path);

        if (paramtop.count("Frequency") == 1)
          read(paramtop, "Frequency", frequency);
        else
          frequency = 1;

        // smeared_props is true if sink-smeared props are given; default=false
        if (paramtop.count("SmearedProps") == 1)
          read(paramtop, "SmearedProps", smeared_props);
        else
          smeared_props = false;

        // The parameters holds the version number
        read(paramtop, "Param", param);

        // If multi_t_sinks are given, use those t_sinks instead of
        // the t_sink in <SeqSource>
        if (paramtop.count("multi_tSinks") == 1)
          read(paramtop, "multi_tSinks", multi_tSinks);
        else
        {
          multi_tSinks.resize(1);
          multi_tSinks[0] = param.t_sink; 
        }

        // The parameters for smearing the sink
        read(paramtop, "PropSink", sink_header);

        // Read in the forward_prop/seqsource info
        read(paramtop, "NamedObject", named_obj);

        // Sanity check
        if(multi_tSinks.size() != named_obj.seqsource_id.size())
        {
          std::ostringstream s;
          s << "Number of t_sinks in multi_tSinks should be the same as the number of seqsource_id elements but, "
            << multi_tSinks.size() << " != " << named_obj.seqsource_id.size() << std::endl;
          throw std::string(s.str());
        }

        // Sanity check 2
        if(multi_tSinks.size() > 1)
        {
          std::vector<int> m_tsinks;
          for(int i=0; i<multi_tSinks.size(); ++i)
            m_tsinks.push_back(multi_tSinks[i]);

          std::sort(m_tsinks.begin(), m_tsinks.end());

          for(int i=0; i<multi_tSinks.size()-1; ++i)
            if(m_tsinks[i+1] - m_tsinks[i] < Ns*Nc)
            {
              std::ostringstream s;
              s << "Separation of multi_tSinks should be larger than Ns*Nc = " << Ns*Nc;
              throw std::string(s.str());
            }
        }
      }
      catch(const std::string& e) 
      {
        QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
        QDP_abort(1);
      }
    }


    void
      Params::writeXML(XMLWriter& xml_out, const std::string& path) 
      {
        push(xml_out, path);

        write(xml_out, "Param", param);
        write(xml_out, "PropSink", sink_header);
        write(xml_out, "NamedObject", named_obj);

        pop(xml_out);
      }


    // Function call
    void 
      InlineMeas::operator()(unsigned long update_no,
          XMLWriter& xml_out) 
      {
        START_CODE();

        StopWatch snoop;
        snoop.reset();
        snoop.start();

        // Test and grab a reference to the gauge field
        XMLBufferWriter gauge_xml;
        try
        {
          TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
          TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
        }
        catch( std::bad_cast ) 
        {
          QDPIO::cerr << name << ": caught dynamic cast error" 
            << std::endl;
          QDP_abort(1);
        }
        catch (const std::string& e) 
        {
          QDPIO::cerr << name << ": map call failed: " << e 
            << std::endl;
          QDP_abort(1);
        }
        const multi1d<LatticeColorMatrix>& u = 
          TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

        push(xml_out, "seqsource");
        write(xml_out, "update_no", update_no);

        QDPIO::cout << name << ": propagator sequential source constructor" << std::endl;
        StopWatch swatch;

        proginfo(xml_out);    // Print out basic program info

        // Write out the input
        params.writeXML(xml_out, "Input");

        // Write out the config header
        write(xml_out, "Config_info", gauge_xml);

        push(xml_out, "Output_version");
        write(xml_out, "out_version", 1);
        pop(xml_out);

        // Calculate some gauge invariant observables just for info.
        MesPlq(xml_out, "Observables", u);

        // Sanity check
        if (params.named_obj.prop_ids.size() == 0)
        {
          QDPIO::cerr << name << ": sanity error: " << std::endl;
          QDP_abort(1);
        }

        //
        // Read the quark propagator and extract headers
        //
        PropSinkSmear_t prop_sink_header;  //if smeared_props
        multi1d<LatticePropagator> forward_props(params.named_obj.prop_ids.size());
        multi1d<ForwardProp_t> forward_headers(params.named_obj.prop_ids.size());
        push(xml_out, "Forward_prop_infos");
        for(int loop=0; loop < params.named_obj.prop_ids.size(); ++loop)
        {
          push(xml_out, "elem");
          try
          {
            // Snarf the data into a copy
            forward_props[loop] =
              TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_ids[loop]);

            // Snarf the source info. This is will throw if the source_id is not there
            XMLReader prop_file_xml, prop_record_xml;
            TheNamedObjMap::Instance().get(params.named_obj.prop_ids[loop]).getFileXML(prop_file_xml);
            TheNamedObjMap::Instance().get(params.named_obj.prop_ids[loop]).getRecordXML(prop_record_xml);

            // Try to invert this record XML into a ChromaProp struct
            // If sink-smeared propagators are input
            if(params.smeared_props)
            {
              ForwardProp_t header;
              read(prop_record_xml, "/SinkSmear", header);

              prop_sink_header                    = header.sink_header;
              forward_headers[loop].prop_header   = header.prop_header;
              forward_headers[loop].source_header = header.source_header;
              forward_headers[loop].gauge_header  = header.gauge_header;
            }
            else
            {
              Propagator_t  header;
              read(prop_record_xml, "/Propagator", header);

              forward_headers[loop].prop_header   = header.prop_header;
              forward_headers[loop].source_header = header.source_header;
              forward_headers[loop].gauge_header  = header.gauge_header;
            }

            // Save prop input
            write(xml_out, "Propagator_info", prop_record_xml);
          }
          catch( std::bad_cast ) 
          {
            QDPIO::cerr << name << ": caught dynamic cast error" 
              << std::endl;
            QDP_abort(1);
          }
          catch (const std::string& e) 
          {
            QDPIO::cerr << name << ": map call failed: " << e 
              << std::endl;
            QDP_abort(1);
          }
          pop(xml_out);
        }
        pop(xml_out);

        QDPIO::cout << "Forward propagator successfully read and parsed" << std::endl;

        // Derived from input prop
        int j_decay  = forward_headers[0].source_header.j_decay;

        // Initialize the slow Fourier transform phases
        SftMom phases(0, true, j_decay);

        // Sanity check - write out the norm2 of the forward prop in the j_decay direction
        // Use this for any possible verification
        push(xml_out, "Forward_prop_correlators");
        for(int loop=0; loop < params.named_obj.prop_ids.size(); ++loop)
        {
          multi1d<Double> forward_prop_corr = sumMulti(localNorm2(forward_props[loop]),
              phases.getSet());

          push(xml_out, "elem");
          write(xml_out, "forward_prop_corr", forward_prop_corr);
          pop(xml_out);
        }
        pop(xml_out);

        // A sanity check
        if (params.param.t_sink < 0 || params.param.t_sink >= QDP::Layout::lattSize()[j_decay]) 
        {
          QDPIO::cerr << "Sink time coordinate incorrect." << std::endl;
          QDPIO::cerr << "t_sink = " << params.param.t_sink << std::endl;
          QDP_abort(1);
        }


        //------------------ Start main body of calculations -----------------------------

        multi1d<LatticePropagator> quark_prop_src(params.multi_tSinks.size());

        try
        {
          // Sink smear the forward propagators
          // NOTE: The smearing construction is pulled outside the loop
          // for efficiency. However, I'm anticipating that we will 
          // have different smearings at the sink of the forward props.
          // In that case, the loop needs to be in inverted.
          std::istringstream  xml_s(params.sink_header.sink.xml);
          XMLReader  sinktop(xml_s);
          QDPIO::cout << "Sink = " << params.sink_header.sink.id << std::endl;

          Handle< QuarkSourceSink<LatticePropagator> >
            sinkSmearing(ThePropSinkSmearingFactory::Instance().createObject(params.sink_header.sink.id,
                  sinktop,
                  params.sink_header.sink.path,
                  u));

          // Do the sink smearing BEFORE the interpolating operator
          for(int loop=0; loop < params.named_obj.prop_ids.size(); ++loop)
          {
            // If propagators are not sink-smeared
            if(!params.smeared_props) {
              forward_headers[loop].sink_header = params.sink_header;
              (*sinkSmearing)(forward_props[loop]);
            }
            else {
             forward_headers[loop].sink_header = prop_sink_header;
            }
          }

          for(int iter_tsink=0; iter_tsink<params.multi_tSinks.size(); ++iter_tsink)
          {
            // Make seqsrc xml for the correct tSink
            params.param.t_sink = params.multi_tSinks[iter_tsink];

            std::ostringstream s;
            s << "<t_sink>" << params.multi_tSinks[iter_tsink];
            std::string new_str(s.str());

            size_t sbeg = params.param.seqsrc.xml.find("<t_sink>");
            size_t send = params.param.seqsrc.xml.find("</t_sink>");
            params.param.seqsrc.xml.replace(sbeg, send-sbeg, new_str); 

            //
            // Construct the sequential source
            //
            QDPIO::cout << "Sequential source = " << params.param.seqsrc.xml << std::endl;

            std::istringstream  xml_seq(params.param.seqsrc.xml);
            XMLReader  seqsrctop(xml_seq);
            QDPIO::cout << "SeqSource = " << params.param.seqsrc.id << std::endl;

            Handle< HadronSeqSource<LatticePropagator> >
              hadSeqSource(TheWilsonHadronSeqSourceFactory::Instance().createObject(params.param.seqsrc.id,
                    seqsrctop,
                    params.param.seqsrc.path));

            swatch.reset();
            swatch.start();
            quark_prop_src[iter_tsink] = (*hadSeqSource)(u, forward_headers, forward_props);

            swatch.stop();

            QDPIO::cout << "Hadron sequential source computed: time= " 
              << swatch.getTimeInSeconds() 
              << " secs" << std::endl;
          } // End of loop over iter_tsink

            
          //-------------------------------------------------------
          // Do the sink smearing AFTER the interpolating operator
          //-------------------------------------------------------
          // Sink smearing needs to be done only for the timeslice t_sink
          // Hence here we reduce the problem size from LatticePropagator
          // to LatticeFermion. The 12 spin x color components of the
          // LatticePropagator are redistributed into 12 timeslices of
          // the LatticeFermion. Gauge links on t_sink are also replicated 
          // to corresponding 12 timeslices.
         
          int j_dec  = params.param.j_decay;

          // Check j_decay
          if(j_decay != 3 && Nd != 4)
          {
            QDPIO::cerr << __func__ << ": j_decay " << j_decay
              << ", Nd " << Nd
              << " not tested." << std::endl;
            QDP_abort(1);
          }

          Set TS;
          TS.make(TSF::TimeSliceFunc(Nd-1));

          multi1d<LatticeColorMatrix> u_smr_d(Nd);
          for(int mu=0; mu<Nd; ++mu)  u_smr_d[mu] = zero;

          // Gauge links only for the timeslice of t_sink
          multi1d<LatticeColorMatrix> u_smr_1(Nd);
          for(int mu=0; mu<Nd; ++mu)  u_smr_1[mu] = zero;

          // block: Make gauge field that has duplicated timeslices
          {
            for(int iter_tsink=0; iter_tsink<params.multi_tSinks.size(); ++iter_tsink)
            {
              int t_sink = params.multi_tSinks[iter_tsink];

              // Loop over \mu in gauge link
              for(int mu=0; mu<Nd; ++mu)
              {
                // Do not need mu in time direction for the klein_gord
                if(mu == j_decay) continue;

                // Replicate one timeslice to u_smr_1
                LatticeColorMatrix u_tmp =  zero;
                u_tmp[TS[t_sink]]        =  u[mu];
                u_smr_1[mu]              += u_tmp;
              }
            }

            // Gauge links whose 12 timeslices are the t_sink.
            // They are copied to t=t_sink, (t_sink+1)%LT, ..., (t_sink+11)%LT
            for(int mu=0; mu<Nd; ++mu)
            {
              // Do not need mu in time direction for the klein_gord
              if(mu == j_decay) continue;
              for(int virtual_t_src=0; virtual_t_src<Ns*Nc; ++virtual_t_src)
              {
                u_smr_d[mu]+= u_smr_1[mu];
                if(virtual_t_src != Ns*Nc-1)
                  u_smr_1[mu] = shift(u_smr_1[mu],BACKWARD,j_decay);
              }
            }
          } // End of block: Make gauge field that has duplicated timeslices

          // Reduce quark_prop_src to LatticeFermion
          LatticeFermion    quark_ferm_src       = zero;

          for(int iter_tsink=0; iter_tsink<params.multi_tSinks.size(); ++iter_tsink)
          {
            int t_sink = params.multi_tSinks[iter_tsink];
          
            LatticePropagator quark_prop_src_tsink = zero;

            quark_prop_src_tsink[TS[t_sink]] = quark_prop_src[iter_tsink];

            for(int color_source = 0; color_source < Nc; ++color_source)
            for(int spin_source = 0; spin_source < Ns; ++spin_source)
            {
              LatticeFermion chi = zero;
              PropToFerm(quark_prop_src_tsink, chi, color_source, spin_source); 
              quark_ferm_src += chi;

              // First source is at t=(t_sink+11)%LT, second source is at t=(t_sink+10)%LT,
              // ... last source is at t=t_sink
              if(color_source!=Nc-1 || spin_source!= Ns-1)
                quark_ferm_src = shift(quark_ferm_src,BACKWARD,j_decay);
            } // loop over color_source & spin_source

            // shift it back to 0 for the next source insertion
            if(iter_tsink != params.multi_tSinks.size()-1)
              for(int i=0; i<(Ns*Nc)-1; ++i)
                quark_ferm_src = shift(quark_ferm_src,FORWARD,j_decay);
          } // End of loop over iter_tsink

          Handle< QuarkSourceSink<LatticeFermion> >
            sinkSmearingFerm(TheFermSinkSmearingFactory::Instance().createObject(params.sink_header.sink.id,
                  sinktop,
                  params.sink_header.sink.path,
                  u_smr_d));

          // Do the sink smearing AFTER the interpolating operator
          (*sinkSmearingFerm)(quark_ferm_src);

          for(int iter_tsink=0; iter_tsink<params.multi_tSinks.size(); ++iter_tsink)
          {
            int t_sink = params.multi_tSinks[iter_tsink];

            // Reconstruct smeared source in LatticePropagator form of quark_prop_src
            // Loop inverse order because smeared quark_ferm_src has the source position inverse order
            quark_prop_src[iter_tsink] = zero;
            for(int color_source = Nc-1; color_source >=0; --color_source)
            for(int spin_source = Ns-1; spin_source >=0; --spin_source)
            {
              LatticeFermion quark_ferm_src_t = zero;
              LatticePropagator eta = zero;
              quark_ferm_src_t[TS[t_sink]] = quark_ferm_src;
              FermToProp(quark_ferm_src_t, eta, color_source, spin_source);
              quark_prop_src[iter_tsink] += eta;
            
              // shift back so that next spin/color component is in t_sink
              if(color_source!=0 || spin_source!=0)
                quark_ferm_src = shift(quark_ferm_src,FORWARD,j_decay);
            }

            // shift it back to 0 for the next t_sink position 
            if(iter_tsink != params.multi_tSinks.size()-1)
              for(int i=0; i<(Ns*Nc)-1; ++i)
                quark_ferm_src = shift(quark_ferm_src,BACKWARD,j_decay);
          } // End of loop over iter_tsink
        }
        catch(const std::string& e) 
        {
          QDPIO::cerr << name << ": Caught Exception in sink: " << e << std::endl;
          QDP_abort(1);
        }


        // Sanity check - write out the norm2 of the propagator source in the j_decay direction
        // Use this for any possible verification
        for(int iter_tsink=0; iter_tsink<params.multi_tSinks.size(); ++iter_tsink)
        {
          multi1d<Double> seqsource_corr = sumMulti(localNorm2(quark_prop_src[iter_tsink]), 
              phases.getSet());

          push(xml_out, "SeqSource_correlator");
          write(xml_out, "seqsource_corr", seqsource_corr);
          pop(xml_out);
        }


        /*
         *  Write the sequential source out to a named buffer
         */
        try
        {
          for(int iter_tsink=0; iter_tsink<params.multi_tSinks.size(); ++iter_tsink)
          {
            // Make seqsrc xml for the correct tSink
            params.param.t_sink = params.multi_tSinks[iter_tsink];

            std::ostringstream s;
            s << "<t_sink>" << params.multi_tSinks[iter_tsink];
            std::string new_str(s.str());

            size_t sbeg = params.param.seqsrc.xml.find("<t_sink>");
            size_t send = params.param.seqsrc.xml.find("</t_sink>");
            params.param.seqsrc.xml.replace(sbeg, send-sbeg, new_str); 

            QDPIO::cout << "Attempt to store sequential source" << std::endl;

            XMLBufferWriter file_xml;
            push(file_xml, "seqsource");
            write(file_xml, "id", uniqueId());  // NOTE: new ID form
            pop(file_xml);
              
            // Sequential source header
            // Header composed of all forward prop headers
            SequentialSource_t new_header;
            new_header.sink_header      = params.sink_header;
            new_header.seqsource_header = params.param;
            new_header.forward_props    = forward_headers;
            new_header.gauge_header     = gauge_xml.printCurrentContext();

            XMLBufferWriter record_xml;
            write(record_xml, "SequentialSource", new_header);

            // Store the seqsource
            TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.seqsource_id[iter_tsink]);
            TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.seqsource_id[iter_tsink]) = quark_prop_src[iter_tsink];
            TheNamedObjMap::Instance().get(params.named_obj.seqsource_id[iter_tsink]).setFileXML(file_xml);
            TheNamedObjMap::Instance().get(params.named_obj.seqsource_id[iter_tsink]).setRecordXML(record_xml);

            QDPIO::cout << "Sequential source successfully stored"  << std::endl;
          }
        }
        catch (std::bad_cast)
        {
          QDPIO::cerr << name << ": dynamic cast error" 
            << std::endl;
          QDP_abort(1);
        }
        catch (const std::string& e) 
        {
          QDPIO::cerr << name << ": error storing seqsource: " << e << std::endl;
          QDP_abort(1);
        }

        pop(xml_out);    // seqsource

        snoop.stop();
        QDPIO::cout << name << ": total time = "
          << snoop.getTimeInSeconds() 
          << " secs" << std::endl;

        QDPIO::cout << name << ": ran successfully" << std::endl;

        END_CODE();
      } 

  }

}
