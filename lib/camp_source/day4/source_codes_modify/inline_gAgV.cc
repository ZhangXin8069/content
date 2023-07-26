#include "inline_gAgV.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/hadron/BuildingBlocks_w.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "actions/ferm/fermstates/ferm_createstate_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_aggregate_w.h"
#include "io_general_class.h"
#include "util/ferm/transf.h"

namespace Chroma
{

  // Environment in which the measurement lives (holds params and such)
  namespace InlineMeasgAgVEnv
  {  
    namespace
    {
      // Function to register with a factory
      AbsInlineMeasurement* createMeasurement( XMLReader& xml_in,
                                       const std::string& path )
      {
        return new InlineMeasgAgV(
                        InlineMeasgAgVParams(xml_in, path) );
      }

      //! Local registration flag
      bool registered = false;
    }

    // The name of the measurement for the XML file
    const std::string name = "Measure_gAgV";

    //! Register all the factories
    bool registerAll()
    {
      bool success = true;
      if (!registered)
      {
        success &= CreateFermStateEnv::registerAll(); 
        success &= TheInlineMeasurementFactory::Instance().
          registerObject(name, createMeasurement);
        registered = true;
      }
      return success;
    }
  } // end namespace InlineMeasgAgVEnv

  //------------------------------------------------------------------------------
  // Parameter reading, writing, handling

  //! Reader for parameters
  void read(                         XMLReader& xml,
                             const std::string& path,
           InlineMeasgAgVParams::Param_t& param )
  {
    XMLReader paramtop(xml, path);
    read(paramtop, "cfg_serial", param.cfg_serial);
    read(paramtop, "gAgV_curr", param.gAgV_curr);
    read(paramtop, "l_prop", param.l_prop);
    read(paramtop, "seq_prop", param.seq_prop);
    read(paramtop, "file_name", param.file_name); 
  }

  //! Writer for parameters
  void write(                              XMLWriter& xml,
                                   const std::string& path,
           const InlineMeasgAgVParams::Param_t& param )
  {
    push(xml, path);

    write(xml, "cfg_serial", param.cfg_serial);
    write(xml, "gAgV_curr", param.gAgV_curr);
    write(xml, "l_prop", param.l_prop);
    write(xml, "seq_prop", param.seq_prop);
    write(xml, "file_name", param.file_name);
  
    pop(xml);
  }


  // Construct params from XML
  InlineMeasgAgVParams::InlineMeasgAgVParams(
            XMLReader& xml_in,
    const std::string& path )
  {
    try
    {
      XMLReader paramtop(xml_in, path);
      frequency = 1;
      read(paramtop, "Param", param);  // Read in the parameters
    }
    catch(const std::string& e)
    {
      QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
      QDP_abort(1);
    }
  }


  // Write out the parameters we constructed
  void InlineMeasgAgVParams::write( XMLWriter& xml_out,
                                  const std::string& path )
  {
    push(xml_out, path);
    
    Chroma::write(xml_out, "Param", param);
    QDP::write(xml_out, "xml_file", xml_file);

    pop(xml_out);
  }

  // Set up the XML and invoke func, which does the acual work
  void InlineMeasgAgV::operator()( unsigned long update_no, XMLWriter& xml_out )
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      std::string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "BuildingBlocks");
      write(xml_out, "update_no", update_no);
      write(xml_out, "xml_file", xml_file);
      pop(xml_out);

      XMLFileWriter xml(xml_file);
      func(update_no, xml);
    }
    else
    {
      func(update_no, xml_out);
    }
  }

  //------------------------------------------------------------------------------
  // Real work done here
  void InlineMeasgAgV::func( unsigned long update_no, XMLWriter& xml_out)
  {
    START_CODE();

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    //--------------------------------------------------------------------------
    // Start building the output XML
    push(xml_out, "BuildingBlocks");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << " BuildingBlocks" << std::endl;
    QDPIO::cout << "     volume: " << QDP::Layout::lattSize()[0];
    for (int i=1; i<Nd; ++i) {
      QDPIO::cout << " x " << QDP::Layout::lattSize()[i];
    }
    QDPIO::cout << std::endl;

    proginfo(xml_out); // Print out basic program info
    push(xml_out, "Output_version");
    write(xml_out, "out_version", 2);
    pop(xml_out);

    //--------------------------------------------------------------------------
    // Grab propagator
    LatticePropagator L_prop, Seq_prop;
    PropSourceConst_t source_header;
    QDPIO::cout << "Attempt to parse forward propagator" << std::endl;

    try
    {
      // Grab the forward propagator
      L_prop = TheNamedObjMap::Instance().getData<LatticePropagator>(params.param.l_prop);
      Seq_prop = TheNamedObjMap::Instance().getData<LatticePropagator>(params.param.seq_prop);
    }
    catch (const std::string& e)
    {
      QDPIO::cerr << InlineMeasgAgVEnv::name
                  << ": propagators: error message: " << e << std::endl;
      QDP_abort(1);
    }
    QDPIO::cout << "All propagators successfully parsed" << std::endl;

    //--------------------------------------------------------------------------
    // Grab the current type
    multi1d<std::string> current_list = params.param.gAgV_curr;
    for(int i = 0; i < current_list.size(); i++)
    {
       QDPIO::cout << "Calculate currents include " << current_list[i] << std::endl;
    }

    QDPIO::cout << "Total current number: " << current_list.size() << std::endl;
   
    int operator_no = current_list.size();
    int j_decay = 3;

    LatticeReal  Phases = 1.;

    // Keep a copy of the phases with no momenta
    SftMom phases_nomom(0, true, j_decay);
    Set timeslice= phases_nomom.getSet();
    int tlen = Layout::lattSize()[j_decay];

    general_data_base res(params.param.file_name.c_str());
    res.add_dimension(dim_conf, 1, &params.param.cfg_serial);
    res.add_dimension(dim_operator, operator_no);
    res.add_dimension(dim_t,tlen);
    res.add_dimension(dim_complex, 2);
    if(Layout::primaryNode()) res.initialize();


    int offset = 0;
    for ( int i = 0; i < operator_no; i++ )
    {
      LatticeComplex corr = zero;

      if (current_list[i] == "gV")
      {
        corr = trace(adj(Gamma(15) * Seq_prop * Gamma(15)) * Gamma(8) * L_prop);
      }
      else if (current_list[i] == "gA")
      {
        corr = trace(adj(Gamma(15) * Seq_prop * Gamma(15)) * Gamma(11) * L_prop);
      }
      else
      {
        QDPIO::cerr << "Unknown current name: " << current_list[i] << std::endl;
        QDP_abort(1);
      }

      multi1d<DComplex> hsum = sumMulti( Phases * corr, timeslice );
      if(Layout::primaryNode())
      for (int t=0; t < tlen; ++t) 
      {
        res.data[offset*tlen*2 + 2*t] = hsum[t].elem().elem().elem().real();
        res.data[offset*tlen*2 + 2*t + 1] = hsum[t].elem().elem().elem().imag();
      }
      offset++;
    }

    if(Layout::primaryNode()) res.save();

    snoop.stop();
    QDPIO::cout << InlineMeasgAgVEnv::name << ": total time = "
                << snoop.getTimeInSeconds()
                << " secs" << std::endl;

    QDPIO::cout << InlineMeasgAgVEnv::name << ": ran successfully"
                << std::endl;

    END_CODE();
  } // end of InlineMeasgAgV::func

}; // end of namespace Chroma
