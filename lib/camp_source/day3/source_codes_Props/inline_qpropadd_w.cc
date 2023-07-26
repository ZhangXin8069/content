/*! \file
 * \brief Inline measurement of qpropadd
 *
 * Addition of props
 */

#include "inline_qpropadd_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  //! Propagator parameters
  void read(XMLReader& xml, const std::string& path, InlineQpropAddCohenEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "factorA", input.factorA);
    read(inputtop, "propA", input.propA);
    read(inputtop, "factorB", input.factorB);
    read(inputtop, "propB", input.propB);
    read(inputtop, "propApB", input.propApB);


    // Read timeslices to add of each propagator
    if (inputtop.count("j_decay") == 1)
      read(inputtop, "j_decay", input.j_decay);
    else
      input.j_decay = Nd-1;

    int N_T = Layout::lattSize()[input.j_decay];

    // Default value is add all timeslices
    multi1d<int> tDefault(2);
    tDefault[0] = 0;
    tDefault[1] = N_T-1;

    if (inputtop.count("tA") == 1)
      read(inputtop, "tA", input.tA);
    else
      input.tA = tDefault;

    if (inputtop.count("tB") == 1)
      read(inputtop, "tB", input.tB);
    else
      input.tB = tDefault;

    if ( input.tA[0] < 0  || input.tA[1] < 0  || input.tB[0] < 0  || input.tB[1] < 0
      || input.tA[0]>N_T-1 || input.tA[1]>N_T-1 || input.tB[0]>N_T-1 || input.tB[1]>N_T-1)
    {
      QDPIO::cerr << __func__ << ": elements of tA or tB should be in [0,NT-1] but they are"
        << "[" << input.tA[0] << ", " << input.tA[1] << "],  "
        << "[" << input.tB[0] << ", " << input.tB[1] << "],  "
        << std::endl;
      QDP_abort(1);
    }

  }

  //! Propagator parameters
  void write(XMLWriter& xml, const std::string& path, const InlineQpropAddCohenEnv::Params::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "j_decay", input.j_decay);
    write(xml, "tA", input.tA);
    write(xml, "factorA", input.factorA);
    write(xml, "propA", input.propA);
    write(xml, "tB", input.tB);
    write(xml, "factorB", input.factorB);
    write(xml, "propB", input.propB);
    write(xml, "propApB", input.propApB);

    pop(xml);
  }


  namespace InlineQpropAddCohenEnv 
  { 
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

    const std::string name = "QPROPADD_cohen";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
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

	// Parameters for source construction
	// Read in the output propagator/source configuration info
	read(paramtop, "NamedObject", named_obj);
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
    
      // Write out the output propagator/source configuration info
      write(xml_out, "NamedObject", named_obj);

      pop(xml_out);
    }

    //--------------------------------------------------------------


    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();

      push(xml_out, "qpropadd");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << "QPROPADD: propagator transformation utility" << std::endl;

      // Write out the input
      params.writeXML(xml_out, "Input");

      //
      // Read in the source along with relevant information.
      // 
      XMLReader propA_file_xml, propA_record_xml;
    
      LatticePropagator propA ;
      LatticePropagator propB ;
      LatticePropagator propApB ;
      QDPIO::cout << "Snarf the props from a named buffer" << std::endl;
      try
      {
	// Try the cast to see if this is a valid source
	propA = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.propA);

	TheNamedObjMap::Instance().get(params.named_obj.propA).getFileXML(propA_file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.propA).getRecordXML(propA_record_xml);

	propB = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.propB);
      }    
      catch (std::bad_cast)
      {
	QDPIO::cerr << name << ": caught dynamic cast error" 
		    << std::endl;
	QDP_abort(1);
      }
      catch (const std::string& e) 
      {
	QDPIO::cerr << name << ": error extracting source_header: " << e << std::endl;
	QDP_abort(1);
      }


        //----------------------------
        // Add the props with weights
        //----------------------------
        int j_decay = params.named_obj.j_decay;
        int N_T = Layout::lattSize()[j_decay];
        if(  params.named_obj.tA[0] == 0 && params.named_obj.tA[1] == N_T
          && params.named_obj.tB[0] == 0 && params.named_obj.tB[1] == N_T )
        {
          propApB = params.named_obj.factorA*propA  +  params.named_obj.factorB*propB;
        }
        else
        {
          // Machinery to do timeslice replication with
          propApB = zero;
          LatticeReal facA=where(
               (((Layout::latticeCoordinate(3))-params.named_obj.tA[0])>=0)&&
               (((Layout::latticeCoordinate(3))-params.named_obj.tA[1])<=0),
               LatticeReal(params.named_obj.factorA),LatticeReal(0.0));
          if(params.named_obj.tA[0] > params.named_obj.tA[1])
          {
            facA=where(
               (((Layout::latticeCoordinate(3))-params.named_obj.tA[0])>=0)||
               (((Layout::latticeCoordinate(3))-params.named_obj.tA[1])<=0),
               LatticeReal(params.named_obj.factorA),LatticeReal(0.0));
          }
          LatticeReal facB=where(
               (((Layout::latticeCoordinate(3))-params.named_obj.tB[0])>=0)&&
               (((Layout::latticeCoordinate(3))-params.named_obj.tB[1])<=0),
               LatticeReal(params.named_obj.factorB),LatticeReal(0.0));
          if(params.named_obj.tB[0] > params.named_obj.tB[1])
          {
            facB=where(
               (((Layout::latticeCoordinate(3))-params.named_obj.tB[0])>=0)||
               (((Layout::latticeCoordinate(3))-params.named_obj.tB[1])<=0),
               LatticeReal(params.named_obj.factorB),LatticeReal(0.0));
          }
          propApB +=facA*propA;
          propApB +=facB*propB;
        }


      /*
       *  Write the a source out to a named buffer
       */
      try
      {
	QDPIO::cout << "Attempt to store sequential source" << std::endl;


	// Store the seqsource
	TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.propApB);
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.propApB) = propApB ;
	TheNamedObjMap::Instance().get(params.named_obj.propApB).setFileXML(propA_file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.propApB).setRecordXML(propA_record_xml);

	QDPIO::cout << "Propagator sum successfully stored"  << std::endl;
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


      pop(xml_out);   // qpropadd
        
      QDPIO::cout << "QpropAddCohen ran successfully" << std::endl;

      END_CODE();
    }

  }

}  
